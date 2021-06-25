site_deets = read.csv('~/git/streampulse/model/site_deets.csv',
    stringsAsFactors=FALSE)
sd2 = read.csv('~/git/streampulse/other_projects/big_stats_paper/sp_siteyears_full.csv',
    stringsAsFactors=FALSE)

zq = paste(site_deets$site_code, site_deets$start_date, sep='_')
zq2 = paste(sd2$sitecode, paste0(sd2$year, '-01-01'), sep='_')
setdiff(zq, zq2)
setdiff(zq2, zq)

sd2$start_date = paste0(sd2$year, '-01-01')
sd2$end_date = paste0(sd2$year, '-12-31')
sd2$year = NULL
sd2 = site_deets %>%
    select(site_code, int) %>%
    right_join(sd2, c('site_code'='sitecode')) %>%
    mutate(int=ifelse(is.na(int), '15 min', int))

site_deets = site_deets %>%
    full_join(sd2, by=c('site_code', 'start_date', 'end_date', 'int')) %>%
    distinct(site_code, start_date, end_date, skip, run_status,
        .keep_all=TRUE) %>%
    arrange(site_code, start_date) %>%
    select(-tok, everything()) %>%
    mutate(skip=ifelse(is.na(skip), '', skip),
        run_status=ifelse(is.na(run_status), '', run_status))


write.csv(site_deets, '~/git/streampulse/model/site_deets_postAZ.csv',
    row.names=FALSE)
