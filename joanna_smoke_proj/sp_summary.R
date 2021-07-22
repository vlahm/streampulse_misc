library(tidyverse)
library(glue)
library(lubridate)

setwd('~/git/streampulse/other_projects/joanna_smoke_proj/')

site_data = read_csv('all_basic_site_data.csv') %>%
    select(region = regionID,
           site = siteID,
           site_longname = siteName,
           latitude,
           longitude,
           usgs_gauge_id = USGSgageID,
           first_record = firstRecord,
           last_record = lastRecord,
           variable_list = variableList)

ests = read_csv('all_daily_model_results.csv')

out = ests %>%
    filter(! region %in% c('AU', 'AT', 'GL', 'PR', 'SE')) %>%
    mutate(region = ifelse(region == 'co', 'CO', region)) %>%
    group_by(region, site, year = year(date)) %>%
    summarize(mean_GPP = mean(GPP, na.rm = TRUE),
              coverage_pct = n() / 365,
              .groups = 'drop') %>%
    group_by(region, site) %>%
    summarize(mean_ann_gpp = round(mean(mean_GPP, na.rm = TRUE), 2),
              mean_ann_coverage_pct = round(mean(coverage_pct, na.rm = TRUE), 2),
              nyears = n(),
              .groups = 'drop') %>%
    left_join(site_data, by = c('region', 'site')) %>%
    select(region, site, site_longname, mean_ann_gpp, mean_ann_coverage_pct,
           nyears, first_record, last_record, latitude, longitude, usgs_gauge_id,
           variable_list)

write_csv(out, 'ann_gpp_summary_streampulse_sites.csv')
