library(tidyverse)
library(lubridate)
library(streamMetabolizer)
library(RMariaDB)
library(glue)

setwd('~/git/streampulse/other_projects/nwqp_ingest/')

# 0 setup ####

cur_time = Sys.time()
attr(cur_time, 'tzone') = 'UTC'

conf = read_lines('../../server_copy/sp/config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}

#connect to mysql via mariadb
pw = extract_from_config('MYSQL_PW')
con = dbConnect(RMariaDB::MariaDB(),
                dbname = 'sp',
                username = 'root',
                password = pw)

# 1 wrangle site data ####

site_info = read_csv('Metabolism.site.info.csv') %>%
    mutate(STATE = c('IL', 'IN', 'IN', 'IA', 'CT', 'CT', 'CT', 'NY', 'NY', 'OR',
                     'OR', 'WA', 'WA', 'WA', 'GA', 'GA', 'GA', 'NC', 'NC', 'NC'),
           SHORT_NAME = ifelse(SHORT_NAME == 'NFDeep', 'NFDeepCreek', SHORT_NAME),
           addDate = cur_time,
           embargo = 0,
           by = -904,
           datum = 'WGS84',
           contact = 'https://doi.org/10.5066/P9YHB00S',
           contactEmail = NA_character_) %>%
    rename(region = STATE,
           site = SHORT_NAME,
           name = STATION_NM,
           latitude = DEC_LAT,
           longitude = DEC_LONG,
           usgs = SITE_NO) %>%
    select(region, site, name, latitude, longitude, datum, usgs, addDate, embargo,
           by, contact, contactEmail)

# 2 harmonize input datasets (and prepare sections of RDS files in section 4) ####

#group data files by variant
site_files = dir('site_files/')

sites_type1 = site_files %>%
    str_match('(.+)?[_\\.]inputs?\\.csv') %>%
    .[, 2] %>%
    na.omit()

sites_type2 = site_files %>%
    str_match('metab_input_(.+)?\\.csv') %>%
    .[, 2] %>%
    na.omit()

#process type 1
d1 = tibble()
d1_raw_list = list()
for(i in seq_along(sites_type1)){

    site = sites_type1[i]
    this_site_info = site_info[site_info$site == site, ]

    input_file_ind = which(grepl(site, site_files) & grepl('input', site_files))
    d1_raw_list[[i]] = read.csv(file.path('site_files',
                                          site_files[input_file_ind]),
                                stringsAsFactors = FALSE) %>%
        mutate(utc.time = streamMetabolizer::convert_solartime_to_UTC(
            any.solar.time = as.POSIXct(solar.time,
                                        tz = 'UTC'),
            longitude = site_info$longitude[site_info$site == !!site]),
            date = as.Date(utc.time)) %>%
        select(solar.time, DO.obs, DO.sat, depth, temp.water, light, date) %>%
        tidyr::drop_na(everything())
    names(d1_raw_list)[i] = site

    d1 = as_tibble(d1_raw_list[[i]]) %>%
        mutate(DateTime_UTC = streamMetabolizer::convert_solartime_to_UTC(
            any.solar.time = as.POSIXct(solar.time,
                                        tz = 'UTC'),
            longitude = site_info$longitude[site_info$site == !!site])) %>%
        select(-solar.time, -any_of(c('X', 'date'))) %>%
        rename(DO_mgL = DO.obs,
               satDO_mgL = DO.sat,
               Depth_m = depth,
               WaterTemp_C = temp.water,
               Light_PAR = light) %>%
        pivot_longer(cols = c(DO_mgL, satDO_mgL, Depth_m, WaterTemp_C, Light_PAR),
                     names_to = 'variable',
                     values_to = 'value',
                     values_drop_na = TRUE) %>%
        mutate(region = this_site_info$region,
               site = this_site_info$site,
               flag = NA_integer_,
               upload_id = -904) %>%
        select(region, site, DateTime_UTC, variable, value, flag, upload_id) %>%
        bind_rows(d1)
}

#process type 2
d2 = tibble()
d2_raw_list = list()
for(i in seq_along(sites_type2)){

    site = sites_type2[i]
    this_site_info = site_info[site_info$site == site, ]

    input_file_ind = which(grepl(site, site_files) & grepl('input', site_files))
    d2_raw_list[[i]] = read.csv(file.path('site_files',
                                          site_files[input_file_ind]),
                                stringsAsFactors = FALSE) %>%
        mutate(utc.time = lubridate::parse_date_time(utc.time,
                orders = c('%m/%d/%Y %H:%M', '%Y-%m-%d %H:%M:%S')),
            solar.time = streamMetabolizer::convert_UTC_to_solartime(
                date.time = utc.time,
                longitude = this_site_info$longitude),
            date = as.Date(utc.time),
            DO.sat = LakeMetabolizer::o2.at.sat.base(
                temp = temp.water,
                baro = pressure.air,
                salinity = 0,
                model = 'garcia-benson')) %>%
        select(solar.time, DO.obs, DO.sat, depth, temp.water, light, date,
               utc.time, pressure.air) %>%
        tidyr::drop_na(everything())
    names(d2_raw_list)[i] = site

    d2 = as_tibble(d2_raw_list[[i]]) %>%
        select(-any_of(c('X', 'date', 'solar.time'))) %>%
        rename(DateTime_UTC = utc.time,
               DO_mgL = DO.obs,
               Depth_m = depth,
               WaterTemp_C = temp.water,
               AirPres_kPa = pressure.air,
               Light_PAR = light) %>%
        pivot_longer(cols = c(DO_mgL, Depth_m, WaterTemp_C, AirPres_kPa, Light_PAR),
                     names_to = 'variable',
                     values_to = 'value',
                     values_drop_na = TRUE) %>%
        mutate(region = this_site_info$region,
               site = this_site_info$site,
               flag = NA_integer_,
               upload_id = -904) %>%
        select(region, site, DateTime_UTC, variable, value, flag, upload_id) %>%
        bind_rows(d2)

    d2_raw_list[[i]] = select(d2_raw_list[[i]], -utc.time, -pressure.air)
}

d = bind_rows(d1, d2)
d_list = c(d1_raw_list, d2_raw_list)

# 3 insert datasets into db ####

dbWriteTable(con, 'nwqp', d, append = TRUE)
dbWriteTable(con, 'site', site_info, append = TRUE)

# 4 wrangle results data and write to shiny app folder ####

shiny_lists_folder = '../../server_copy/sp/shiny/model_viz/nwqp_data/shiny_lists'
dir.create(shiny_lists_folder, recursive = TRUE, showWarnings = FALSE)

results_type1 = site_files %>%
    str_match('(.+)?[_\\.]metab.+?\\.csv') %>%
    .[, 2] %>%
    na.omit()

results_type2 = site_files %>%
    str_match('metab_predict_lightobs_meank_(.+)?\\.csv') %>%
    .[, 2] %>%
    na.omit()

for(i in seq_along(results_type1)){

    sm_out = list()
    site = results_type1[i]
    this_site_info = site_info[site_info$site == site, ]

    results_file_ind = which(grepl(site, site_files) & grepl('metab.+?\\.csv', site_files))
    sm_out$predictions = read.csv(file.path('site_files',
                                            site_files[results_file_ind]),
                  stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(site_name = site,
               date = as.Date(date),
               # resolution = NA_integer_,
               # GPP.n_eff = NA_real_,
               # GPP.Rhat = NA_real_,
               # ER.n_eff = NA_real_,
               # ER.Rhat = NA_real_,
               K600.lower = NA_real_,
               K600.upper = NA_real_) %>%
               # K600.n_eff = NA_real_,
               # K600.Rhat = NA_real_,
               # DO.obs = NA_real_,
               # DO.sat = NA_real_,
               # DO.amp = NA_real_,
               # DO.psat = NA_real_,
               # depth = NA_real_,
               # temp.water = NA_real_,
               # day.length = NA_real_,
               # discharge = NA_real_,
               # shortwave = NA_real_,
               # velocity = NA_real_,
               # DO.tdist80 = NA_real_) %>%
        # select(site_name, resolution, date, GPP, GPP.lower, GPP.upper,
        #        GPP.n_eff, GPP.Rhat, ER, ER.lower, ER.upper, ER.n_eff, ER.Rhat,
        #        K600 = K600.fixed,
        #        K600.lower, K600.upper, K600.n_eff, K600.Rhat, DO.obs,
        #        DO.sat, DO.amp, DO.psat, depth, temp.water, day.length, discharge,
        #        shortwave, velocity, DO.tdist80)
        rename_all(dplyr::recode, K600.fixed = 'K600') %>%
        select(site_name, date, GPP, GPP.lower, GPP.upper, ER, ER.lower, ER.upper,
               any_of('K600'),
               K600.lower, K600.upper)

    if(! 'K600' %in% colnames(sm_out$predictions)){
        sm_out$predictions$K600 = NA
        sm_out$predictions = select(sm_out$predictions,
               site_name, date, GPP, GPP.lower, GPP.upper, ER, ER.lower, ER.upper,
               K600, K600.lower, K600.upper)
    }

    # zz = readRDS('../../../streampulse/server_copy/sp/shiny/model_viz/powell_data/shiny_lists/AK_15298040_2010.rds')
    # print(as_tibble(zz$data), n=2)
    sm_out$data = d_list[[site]]
    # zz = readRDS('../../../streampulse/server_copy/sp/shiny/model_viz/powell_data/shiny_lists/AR_07075250_2013.rds')
    # print(as_tibble(zz$data_daily), n=2)
    sm_out$data_daily = sm_out$data %>%
        group_by(date) %>%
        mutate(discharge = NA_real_) %>%
        ungroup() %>%
        select(date, discharge) %>%
        as.data.frame()
    # print(as_tibble(zz$fit$daily), n=2)
    sm_out$fit = list(daily = data.frame(date = sm_out$predictions$date,
                                         K600_daily_mean = sm_out$predictions$K600,
                                         K600_daily_97.5pct = NA_real_,
                                         K600_daily_2.5pct = NA_real_,
                                         ER_daily_mean = sm_out$predictions$ER,
                                         GPP_daily_mean = sm_out$predictions$GPP))
                                         # discharge.daily = NA_real_))

    yr = unique(lubridate::year(d_list[[site]]$date))
    if(length(yr) > 1) stop('oi')

    saveRDS(sm_out,
            glue('{slf}/{r}_{s}_{y}.rds',
                 slf = shiny_lists_folder,
                 r = this_site_info$region,
                 s = this_site_info$site,
                 y = yr))
}

for(i in seq_along(results_type2)){

    sm_out = list()
    site = results_type2[i]
    this_site_info = site_info[site_info$site == site, ]

    results_file_ind = which(grepl(site, site_files) &
                                 grepl('metab.+?\\.csv', site_files) &
                                 ! grepl('input', site_files))
    sm_out$predictions = read.csv(file.path('site_files',
                                            site_files[results_file_ind]),
                  stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(site_name = site,
               date = as.Date(date)) %>%
        select(site_name, date, GPP, GPP.lower, GPP.upper, ER, ER.lower, ER.upper,
               K600, K600.lower, K600.upper)

    sm_out$data = d_list[[site]]

    sm_out$data_daily = sm_out$data %>%
        group_by(date) %>%
        mutate(discharge = NA_real_) %>%
        ungroup() %>%
        select(date, discharge) %>%
        as.data.frame()

    sm_out$fit = list(daily = data.frame(date = sm_out$predictions$date,
                                         K600_daily_mean = sm_out$predictions$K600,
                                         K600_daily_97.5pct = sm_out$predictions$K600.upper,
                                         K600_daily_2.5pct = sm_out$predictions$K600.lower,
                                         ER_daily_mean = sm_out$predictions$ER,
                                         GPP_daily_mean = sm_out$predictions$GPP))
                                         # discharge.daily = NA_real_))

    yr = unique(lubridate::year(d_list[[site]]$date))
    if(length(yr) > 1) stop('oi')

    saveRDS(sm_out,
            glue('{slf}/{r}_{s}_{y}.rds',
                 slf = shiny_lists_folder,
                 r = this_site_info$region,
                 s = this_site_info$site,
                 y = yr))
}

