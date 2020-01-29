library(StreamPULSE)
library(tidyverse)
rm(list=ls()); cat('/014')

setwd('~/git/streampulse/other_projects/big_stats_paper/')

sites = as_tibble(readRDS('sites_COMID.rds'))
nwis_ind = substr(sites$Site_ID, 1, 4) == 'nwis'
sites$Site_ID[nwis_ind] = gsub('_', '-', sites$Site_ID[nwis_ind])
sites$site = unname(sapply(sites$Site_ID, function(x) strsplit(x, '_')[[1]][2]))
sites$site[is.na(sites$site)] = sites$Site_ID[is.na(sites$site)]
sites$Site_ID = NULL

mods = query_available_results('all')[[1]] %>%
    as_tibble() %>%
    mutate_all(as.character) %>%
    # left_join(select(sites, Site_ID), c('site'='Site_ID')) %>%
    left_join(sites, 'site') %>%
    mutate(Site_ID=paste(region, site, sep='_')) %>%
    filter(! is.na(Name)) %>%
    distinct() %>%
    arrange(Site_ID, year)
# nwis_ind = substr(mods$Site_ID, 4, 7) == 'nwis'
# nwis_split = strsplit(mods$Site_ID[nwis_ind], '_')
# mods$Site_ID[nwis_ind] = sapply(nwis_split, function(x) x[2])


#gpp-er biplots ####
cnt = 0
hold_cnt = FALSE
errtypes = errmods = errs = list()

pdf(file='streampulse_gppXer_all.pdf', onefile=TRUE)
par(mfrow=c(3, 3), mar=c(3, 3, 1, 1), oma=c(1, 1, 1, 1))

for(i in 1:nrow(mods)){

    m = mods[i,]
    res = try( request_results(sitecode=m$Site_ID,
                year=as.numeric(m$year)) )

    if(class(res) == 'try-error'){
        errtypes = append(errtypes, 'request_err')
        errmods = append(errmods, paste(m$Site_ID, i))
        errs = append(errs, res)
        hold_cnt = TRUE
        next
    }

    if(nrow(res$predictions) == 0){
        errtypes = append(errtypes, 'no_predictions')
        errmods = append(errmods, paste(m$Site_ID, i))
        errs = append(errs, 'NA')
    }

    tryCatch({

        mout = lm(res$predictions$ER ~ res$predictions$GPP)
        mres = summary(mout)
        coeffs = round(unname(mres$coefficients[,'Estimate']), 2)
        icept = sprintf('%+.2f', coeffs[1])
        icept = paste(substr(icept, 0, 1), substr(icept, 2, nchar(icept)))
        coeff_lab = paste0('y = ', coeffs[2], 'x ', icept, ' (Adj. R^2: ',
            round(mres$adj.r.squared, 2), ')')

        # plot(res$predictions$GPP, res$predictions$ER, main=m$Site_ID,
        plot(res$predictions$GPP, res$predictions$ER, main='', xlab='GPP',
            ylab='ER', bty='l', col=alpha('gray30', 0.7), pch=20)
        abline(mout, lty=2, col='red', lwd=2)
        mtext(paste0(m$Site_ID, ' (', m$year, ')\n', coeff_lab),
            side=3, line=0, cex=0.7)


    }, error=function(e){
            errtypes = append(errtypes, 'plot/stat_err')
            errmods = append(errmods, paste(m$Site_ID, i))
            errs = append(errs, e)
        }
    )

    cnt = cnt + 1
    gc()
}

dev.off()


# time series plots without K ####
cnt = 0
hold_cnt = FALSE
errtypes2 = errmods2 = errs2 = list()

pdf(file='streampulse_metab_ts_all_noK.pdf', onefile=TRUE)
par(mfrow=c(3, 3), mar=c(3, 3, 1, 1), oma=c(1, 1, 1, 1))

# for(i in 1:nrow(mods)){
for(i in 1:10){

    m = mods[i,]
    res = try( request_results(sitecode=m$Site_ID,
        year=as.numeric(m$year)) )

    if(class(res) == 'try-error'){
        errtypes2 = append(errtypes2, 'request_err')
        errmods2 = append(errmods2, paste(m$Site_ID, i))
        errs2 = append(errs2, res)
        hold_cnt = TRUE
        next
    }

    if(nrow(res$predictions) == 0){
        errtypes2 = append(errtypes2, 'no_predictions')
        errmods2 = append(errmods2, paste(m$Site_ID, i))
        errs2 = append(errs2, 'NA')
    }

    tryCatch({

        o = res$model_results$fit$daily

        llim = min(c(o$GPP_daily_2.5pct, o$ER_daily_2.5pct), na.rm=TRUE)
        ulim = max(c(o$GPP_daily_97.5pct, o$ER_daily_97.5pct), na.rm=TRUE)
        # maxmin_day = range(doy, na.rm=TRUE)
        plot(o$date, o$GPP_daily_mean, type='l', col='red', xlab='', las=0,
            ylab='', xaxs='i', yaxs='i', ylim=c(llim, ulim), bty='l', yaxt='n')
        axis(2, tck=-.02, labels=FALSE)
        axis(2, tcl=0, col='transparent', line=-0.5)
        mtext(paste0(m$Site_ID, ' (', m$year, ')'), side=3, line=0, cex=0.7)
        lines(o$date, o$ER_daily_mean, col='blue')
        mtext(expression(paste("g"~O[2]~"m"^"-2"~" d"^"-1")), side=2,
            line=1.5, font=2, cex=0.7)


        rl = rle(is.na(o$GPP_daily_2.5pct))
        vv = !rl$values
        chunkfac = rep(cumsum(vv), rl$lengths)
        chunkfac[chunkfac == 0] = 1
        chunks = split(o, chunkfac)
        noNAchunks = lapply(chunks, function(x) x[!is.na(x$GPP_daily_2.5pct),] )

        for(i in 1:length(noNAchunks)){
            polygon(x=c(noNAchunks[[i]]$date, rev(noNAchunks[[i]]$date)),
                y=c(noNAchunks[[i]]$GPP_daily_2.5pct,
                    rev(noNAchunks[[i]]$GPP_daily_97.5pct)),
                col=adjustcolor('red', alpha.f=0.3), border=NA)
            polygon(x=c(noNAchunks[[i]]$date, rev(noNAchunks[[i]]$date)),
                y=c(noNAchunks[[i]]$ER_daily_2.5pct,
                    rev(noNAchunks[[i]]$ER_daily_97.5pct)),
                col=adjustcolor('blue', alpha.f=0.3), border=NA)
        }

        abline(h=0, lty=3, col='gray50')

    }, error=function(e){
            errtypes2 = append(errtypes2, 'plot/stat_err')
            errmods2 = append(errmods2, paste(m$Site_ID, i))
            errs2 = append(errs2, e)
    }, warning=function(w) NULL)

    cnt = cnt + 1
    gc()
}

dev.off()


# time series plots with K ####
cnt = 0
hold_cnt = FALSE
errtypes_3 = errmods_3 = errs_3 = list()

pdf(file='streampulse_metab_ts_all.pdf', onefile=TRUE)
par(mfrow=c(3, 3), mar=c(3, 3, 1, 3), oma=c(1, 1, 1, 1))

# for(i in 1:nrow(mods)){
for(i in 1:10){

    m = mods[i,]
    res = try( request_results(sitecode=m$Site_ID,
        year=as.numeric(m$year)) )

    if(class(res) == 'try-error'){
        errtypes_3 = append(errtypes_3, 'request_err')
        errmods_3 = append(errmods_3, paste(m$Site_ID, i))
        errs_3 = append(errs_3, res)
        hold_cnt = TRUE
        next
    }

    if(nrow(res$predictions) == 0){
        errtypes_3 = append(errtypes_3, 'no_predictions')
        errmods_3 = append(errmods_3, paste(m$Site_ID, i))
        errs_3 = append(errs_3, 'NA')
    }

    tryCatch({

        o = res$model_results$fit$daily

        llim = min(c(o$GPP_daily_2.5pct, o$ER_daily_2.5pct), na.rm=TRUE)
        ulim = max(c(o$GPP_daily_97.5pct, o$ER_daily_97.5pct), na.rm=TRUE)
        # maxmin_day = range(doy, na.rm=TRUE)
        plot(o$date, o$GPP_daily_mean, type='l', col='red', xlab='', las=0,
            ylab='', xaxs='i', yaxs='i', ylim=c(llim, ulim), bty='l', yaxt='n')
        axis(2, tck=-.02, labels=FALSE)
        axis(2, tcl=0, col='transparent', line=-0.5)
        mtext(paste0(m$Site_ID, ' (', m$year, ')'), side=3, line=0, cex=0.7)
        lines(o$date, o$ER_daily_mean, col='blue')
        mtext(expression(paste("g"~O[2]~"m"^"-2"~" d"^"-1")), side=2,
            line=1.5, font=2, cex=0.7)


        rl = rle(is.na(o$GPP_daily_2.5pct))
        vv = !rl$values
        chunkfac = rep(cumsum(vv), rl$lengths)
        chunkfac[chunkfac == 0] = 1
        chunks = split(o, chunkfac)
        noNAchunks = lapply(chunks, function(x) x[!is.na(x$GPP_daily_2.5pct),] )

        for(i in 1:length(noNAchunks)){
            polygon(x=c(noNAchunks[[i]]$date, rev(noNAchunks[[i]]$date)),
                y=c(noNAchunks[[i]]$GPP_daily_2.5pct,
                    rev(noNAchunks[[i]]$GPP_daily_97.5pct)),
                col=adjustcolor('red', alpha.f=0.3), border=NA)
            polygon(x=c(noNAchunks[[i]]$date, rev(noNAchunks[[i]]$date)),
                y=c(noNAchunks[[i]]$ER_daily_2.5pct,
                    rev(noNAchunks[[i]]$ER_daily_97.5pct)),
                col=adjustcolor('blue', alpha.f=0.3), border=NA)
        }

        abline(h=0, lty=3, col='gray50')

        par(new=TRUE)
        llim2 = min(o$K600_daily_2.5pct, na.rm=TRUE)
        ulim2 = max(o$K600_daily_97.5pct, na.rm=TRUE)
        plot(o$date, o$K600_daily_mean, col='orange',
            type='l', xlab='', las=0, ylab='', xaxs='i', yaxs='i',
            xaxt='n', bty='u', yaxt='n', ylim=c(llim2, ulim2))
        axis(4, tck=-.02, labels=FALSE)
        axis(4, tcl=0, col='transparent', line=-0.5, las=1)
        mtext(expression(paste('K600 (d'^'-1' * ')')), side=4, line=2,
            font=2, cex=0.7)
        # for(i in 1:length(noNAchunks)){
        #     polygon(x=c(noNAchunks[[i]]$date, rev(noNAchunks[[i]]$date)),
        #         y=c(noNAchunks[[i]]$K600_daily_2.5pct,
        #             rev(noNAchunks[[i]]$K600_daily_97.5pct)),
        #         col=adjustcolor('orange', alpha.f=0.3), border=NA)
        # }

    }, error=function(e){
        errtypes_3 = append(errtypes_3, 'plot/stat_err')
        errmods_3 = append(errmods_3, paste(m$Site_ID, i))
        errs_3 = append(errs_3, e)
    }, warning=function(w) NULL)

    cnt = cnt + 1
    gc()
}

dev.off()
