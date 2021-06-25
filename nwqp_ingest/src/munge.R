library(tidyverse)
setwd('~/git/streampulse/other_projects/nwqp_ingest/')
sites = dir('site_files/') %>%
    str_match('metab_input_(.+)?.csv') %>%
    .[, 2]
na.omit(str_match(dir('site_files'), 'metab_input_(.+)?.csv')[, 2])

#why are some sites not included?
