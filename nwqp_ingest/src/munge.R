library(tidyverse)
library(lubridate)
library(streamMetabolizer)
library(RMariaDB)

setwd('~/git/streampulse/other_projects/nwqp_ingest/')

#setup ####

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

# wrangle site data ####

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

# harmonize input datasets ####

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
for(i in seq_along(sites_type1)){

    site = sites_type1[i]
    this_site_info = site_info[site_info$site == site, ]

    input_file_ind = which(grepl(site, site_files) & grepl('input', site_files))
    d1 = read.csv(file.path('site_files', site_files[input_file_ind]),
                 stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(DateTime_UTC = streamMetabolizer::convert_solartime_to_UTC(
            any.solar.time = as.POSIXct(solar.time,
                                        tz = 'UTC'),
            longitude = site_info$longitude[site_info$site == !!site])) %>%
        select(-solar.time) %>%
        rename(DO_mgL = DO.obs,
               satDO_mgL = DO.sat,
               Depth_m = depth,
               WaterTemp_C = temp.water,
               Light_PAR = light) %>%
        pivot_longer(cols = c(DO_mgL, satDO_mgL, Depth_m, WaterTemp_C, Light_PAR),
                     names_to = 'variable',
                     values_to = 'value',
                     values_drop_na = TRUE) %>%
        mutate(region = this_site_info$STATE,
               site = this_site_info$SHORT_NAME,
               flag = NA_integer_,
               upload_id = -904) %>%
        select(region, site, DateTime_UTC, variable, value, flag, upload_id) %>%
        bind_rows(d1)
}

#process type 2
d2 = tibble()
for(i in seq_along(sites_type2)){

    site = sites_type2[i]
    this_site_info = site_info[site_info$site == site, ]

    input_file_ind = which(grepl(site, site_files) & grepl('input', site_files))
    d2 = read.csv(file.path('site_files', site_files[input_file_ind]),
                 stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(DateTime_UTC = lubridate::parse_date_time(
            utc.time,
            orders = c('%m/%d/%Y %H:%M', '%Y-%m-%d %H:%M:%S'))) %>%
        select(-utc.time) %>%
        rename(DO_mgL = DO.obs,
               Depth_m = depth,
               WaterTemp_C = temp.water,
               AirPres_kPa = pressure.air,
               Light_PAR = light) %>%
        pivot_longer(cols = c(DO_mgL, Depth_m, WaterTemp_C, AirPres_kPa, Light_PAR),
                     names_to = 'variable',
                     values_to = 'value',
                     values_drop_na = TRUE) %>%
        mutate(region = this_site_info$STATE,
               site = this_site_info$SHORT_NAME,
               flag = NA_integer_,
               upload_id = -904) %>%
        select(region, site, DateTime_UTC, variable, value, flag, upload_id) %>%
        bind_rows(d2)
}

d = bind_rows(d1, d2)


# insert datasets into db ####

dbWriteTable(con, 'nwqp', d, append = TRUE)
dbWriteTable(con, 'site', site_info, append = TRUE)

# wrangle results data and write to shiny app folder ####

results_type1 = site_files %>%
    str_match('(.+)?[_\\.]metab.+?\\.csv') %>%
    .[, 2] %>%
    na.omit()

results_type2 = site_files %>%
    str_match('metab_predict_lightobs_meank_(.+)?\\.csv') %>%
    .[, 2] %>%
    na.omit()

for(i in seq_along(results_type1)){

    site = results_type1[i]
    this_site_info = site_info[site_info$site == site, ]

    results_file_ind = which(grepl(site, site_files) & grepl('metab', site_files))
    #HERE: D1 SHOULD MATCH zz$predictions (or some other element)
    d1 = read.csv(file.path('site_files', site_files[results_file_ind]),
                  stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(DateTime_UTC = streamMetabolizer::convert_solartime_to_UTC(
            any.solar.time = as.POSIXct(solar.time,
                                        tz = 'UTC'),
            longitude = site_info$longitude[site_info$site == !!site])) %>%
        select(-solar.time) %>%
        rename(DO_mgL = DO.obs,
               satDO_mgL = DO.sat,
               Depth_m = depth,
               WaterTemp_C = temp.water,
               Light_PAR = light) %>%
        pivot_longer(cols = c(DO_mgL, satDO_mgL, Depth_m, WaterTemp_C, Light_PAR),
                     names_to = 'variable',
                     values_to = 'value',
                     values_drop_na = TRUE) %>%
        mutate(region = this_site_info$STATE,
               site = this_site_info$SHORT_NAME,
               flag = NA_integer_,
               upload_id = -904) %>%
        select(region, site, DateTime_UTC, variable, value, flag, upload_id) %>%
        bind_rows(d1)
}

for(i in 1:length(powIn)){

    shinylist = list(predictions=Es,
        # fit=list(daily=PO, KQ_binned=PE$KQ_binned),
        fit=list(daily=PO), #cant do KQ nodes because data_index field is all NA
        data=PI, #missing DO.mod (not published)
        data_daily=Es[,c('date','discharge')])

    zz = shinylist
    tibble::as_tibble(zz$predictions)
    tibble::as_tibble(zz$data)
    tibble::as_tibble(zz$data_daily)
    tibble::as_tibble(zz$fit$daily)

    saveRDS(shinylist, paste0('powell_data/shiny_lists/', Pc[1], '_', Pc[2],
        '_', Pc[3], '.rds'))
}
