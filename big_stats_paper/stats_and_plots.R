library(plyr)
library(StreamPULSE)
library(tidyverse)
rm(list=ls()); cat('/014')

#setup ####

#set mode: 'run' will rerun data generators; 'retrieve' will collect outputs
mode = 'retrieve'
# mode = 'run'

setwd('~/git/streampulse/other_projects/big_stats_paper/')
source('helpers.R')
phil_srcs = list.files('phil_stuff/metab_synthesis/R/functions/',
    full.names=TRUE)
for(x in phil_srcs) source(x)

#model results from streampulse portal
mods = query_available_results('all')[[1]] %>%
    as_tibble() %>%
    mutate_all(as.character)

#datasets from Phil
metab_d = readRDS('phil_stuff/output/synthesis_standardized.rds')
diag = as_tibble(readRDS('phil_stuff/metab_synthesis/output/yearly_diagnostics.rds'))
metr = as_tibble(readRDS('phil_stuff/output/site_metrics.rds')) %>%
    arrange(desc(ann_GPP_C))
ws_metr = as_tibble(readRDS('phil_stuff/output/metrics_bound.rds'))
ws_metr = phil_to_mike_format(ws_metr, mods)
ws_metr = select(ws_metr, -Source, -Lat, -Lon, -VPU, -site, -region, -COMID, -Name)
# filled = readRDS('phil_stuff/output/synthesis_gap_filled.rds')

#subset of streampulse + powell sites used in this analysis
sites = as_tibble(readRDS('sites_COMID.rds'))
sites = phil_to_mike_format(sites, mods)

WGS84 = 4326 #EPSG code for coordinate reference system

#bind NHDPlusV2 data ####

if(mode == 'run'){

    errflag = FALSE
    sites$reach_proportion = NA
    ni = nrow(sites)

    for(i in 1:ni){

        print(paste(i, ni, sep='/'))
        tryCatch({
            rp = calc_reach_prop(sites$VPU[i], sites$COMID[i], sites$Lat[i],
                sites$Lon[i], WGS84, quiet=TRUE)
        }, error=function(e) {print('error'); errflag <<- TRUE} )

        if(errflag) {
            errflag = FALSE
            sites$reach_proportion[i] = NA
            next
        }

        sites$reach_proportion[i] = rp
    }

    saveRDS(sites$reach_proportion, 'output/reach_prop_col.rds')

} else {
    sites$reach_proportion = readRDS('output/reach_prop_col.rds')
}

#separate NHD-able sites from non; will recombine later
no_comid_site_inds = which(is.na(sites$COMID))
no_comid_sites = sites[no_comid_site_inds, ]
sites = sites[-no_comid_site_inds, ]

if(mode == 'run'){

    #construct list of DSN=component pairs to acquire. see NHDPlus docs for more.
    setlist = list('NHDPlusAttributes'='PlusFlowlineVAA',
        'NHDPlusAttributes'='ElevSlope')

    #retrieve NHDPlusV2 data
    nhdplusv2_data = nhdplusv2_bulk(sites, setlist, quiet=TRUE)

    #nhd variable names do not have consistent naming conventions. sometimes they're
    #all caps; other times camel case. here's a crude way to deal with that.
    colnames(nhdplusv2_data) = toupper(colnames(nhdplusv2_data))
    nhdplusv2_data = nhdplusv2_data[, ! duplicated(colnames(nhdplusv2_data))]

    #choose variables to join
    nhdplusv2_data = select(nhdplusv2_data, COMID, STREAMORDE, FROMMEAS, TOMEAS, SLOPE,
        REACHCODE, AREASQKM, TOTDASQKM, MAXELEVSMO, MINELEVSMO)

    saveRDS(nhdplusv2_data, 'output/nhdplusv2_data.rds')

} else {
    nhdplusv2_data = readRDS('output/nhdplusv2_data.rds')
}


sites = left_join(sites, nhdplusv2_data, by='COMID')
sites = sites[! duplicated(sites$Name),]

#correct catchment area (AREASQKM) based on where each site falls within its reach.
#use this to correct watershed area (TOTDASQKM) and to determine an areal
#correction factor that can be multiplied with any areal summary data.
sites$AREASQKM_corr = round(sites$AREASQKM * sites$reach_proportion, 5)
sites$TOTDASQKM_corr = sites$TOTDASQKM - (sites$AREASQKM - sites$AREASQKM_corr)
sites$areal_corr_factor = sites$TOTDASQKM_corr / sites$TOTDASQKM

#bind StreamCat data (superfluous since Phil added ws_metr) ####

if(mode == 'run'){

    #construct vector of streamcat datasets to acquire (check variable list for deets)
    setlist2 = c('Elevation', 'PRISM_1981_2010', 'NLCD2011', 'Runoff',
        'STATSGO_Set2', 'NADP', 'GeoChemPhys1', 'GeoChemPhys2', 'BFI')

    streamcat_data1 = streamcat_bulk(sites[1:100,], setlist2)
    saveRDS(streamcat_data1, 'output/sc1.rds')
    streamcat_data2 = streamcat_bulk(sites[101:200,], setlist2)
    saveRDS(streamcat_data2, 'output/sc2.rds')
    streamcat_data3 = streamcat_bulk(sites[201:300,], setlist2)
    saveRDS(streamcat_data3, 'output/sc3.rds')
    streamcat_data4 = streamcat_bulk(sites[301:404,], setlist2)
    saveRDS(streamcat_data4, 'output/sc4.rds')

    #select and join variables
    streamcat_data = select(streamcat_data, COMID, ElevWs, Precip8110Ws, Tmin8110Ws,
        Tmax8110Ws, Tmean8110Ws, RunoffWs, matches('^Pct[a-zA-z]+2011Ws$'),
        PermWs, RckDepWs, OmWs, WtDepWs, matches('^[a-zA-z0-9]_2008Ws$'), BFIWs,
        NWs, Al2O3Ws, CaOWs, Fe2O3Ws, K2OWs, MgOWs, Na2OWs, P2O5Ws, SWs, SiO2Ws) %>%
        mutate(precip_runoff_ratio=Precip8110Ws / RunoffWs)

    saveRDS(streamcat_data, 'output/streamcat_data.rds')

} else {
    streamcat_data = readRDS('output/streamcat_data.rds')
}

sites = left_join(sites, streamcat_data, by='COMID')
sites = sites[! duplicated(sites$sitecode),]
sites = arrange(sites, region, sitecode)

#bind Phil's StreamCat dataset ####
sites = left_join(sites, ws_metr, by='sitecode')
sites = sites[! duplicated(sites$sitecode),]

#generate/retrieve watershed boundaries ####

if(mode == 'run'){

    library(rgdal)
    library(streamstats)
    library(htmlwidgets)

    unlink('spatial/indiv_sp_ws/*', recursive=TRUE)
    unlink('spatial/interactive_maps/*', recursive=TRUE)

    sites_sp = filter(sites, Source == 'StreamPULSE')

    for(i in 1:nrow(sites_sp)){

        current_sitecode = sites_sp$sitecode[i]

        ws_bound = try( streamstats::delineateWatershed(sites_sp$Lon[i],
            sites_sp$Lat[i], crs=WGS84) )
        if('try-error' %in% class(ws_bound)) next

        #save watershed boundary as shapefile
        tpf = tempfile(fileext='.geojson')
        streamstats::writeGeoJSON(ws_bound, file=tpf, what='boundary')
        spatialdf = rgdal::readOGR(tpf)
        unlink(tpf)
        rgdal::writeOGR(spatialdf, dsn='spatial/indiv_sp_ws',
            layer=current_sitecode, driver='ESRI Shapefile')

        #thin to just 30 points and save copy for earthexplorer
        crds = spatialdf@polygons[[1]]@Polygons[[1]]@coords
        filt_ind = round(seq(1, nrow(crds), length.out=30))
        spatialdf@polygons[[1]]@Polygons[[1]]@coords = crds[filt_ind, ]

        rgdal::writeOGR(spatialdf, dsn='spatial/indiv_sp_ws_30pts',
            layer=current_sitecode, driver='ESRI Shapefile')

        #save interactive webpage for viewing boundary on basemap
        webfile = try( streamstats::leafletWatershed(ws_bound) )
        if('try-error' %in% class(webfile)) next

        setwd('spatial/interactive_maps/')
        htmlwidgets::saveWidget(webfile, selfcontained=FALSE,
            file=paste0(current_sitecode, '.html'))
        setwd('../..')
    }

}

#filter and zip boundaries ####

shedfiles = list.files('spatial/indiv_sp_ws/')
sheds = unique(unname(sapply(shedfiles, function(x){
        strsplit(x, '\\.')[[1]][1]
    })))

for(s in sheds){
    setwd('spatial/indiv_sp_ws/')
    zip(paste0(s, '.zip'), list.files(pattern=s))
    file.rename(paste0(s, '.zip'), paste0('../indiv_sp_ws_zips/', s, '.zip'))
    setwd('../..')
    # wshed = rgdal::readOGR('spatial/indiv_sp_ws/', s)
}

shedfiles_30pt = list.files('spatial/indiv_sp_ws_30pts/')
sheds_30pt = unique(unname(sapply(shedfiles_30pt, function(x){
        strsplit(x, '\\.')[[1]][1]
    })))

for(s in sheds_30pt){
    setwd('spatial/indiv_sp_ws_30pts/')
    zip(paste0(s, '.zip'), list.files(pattern=s))
    file.rename(paste0(s, '.zip'), paste0('../indiv_sp_ws_30pts_zips/', s, '.zip'))
    setwd('../..')
    # wshed = rgdal::readOGR('spatial/indiv_sp_ws/', s)
}

# HERE: wget and unzip powell datasets

#bind IGBP landcover classifications (MCD12Q1v006 LC_Type1) ####

library(gdalUtils)
library(raster)

setwd('spatial/igbp_classes/hdfs')
hdfs = list.files(pattern='.hdf')

#extract igbp classification layer from hdf4s and convert to geotiffs
for(f in hdfs){
    sitecode = strsplit(f, '\\.')[[1]][1]
    sds = gdalUtils::get_subdatasets(f)
    igbp_ind = grep('Type1', sds)
    gtfile = paste0('../geotiffs/', sitecode, '.tif')
    gdalUtils::gdal_translate(sds[igbp_ind], dst_dataset=gtfile)
}

#HERE: CUT GEOTIFF RASTER OF ECOREGIONS BY WS BOUNDS
    rast = raster::raster(gtfile)


#summarize raster by wshed boundary
pgon = sp::SpatialPolygons(wsboundary@polygons,
    # proj4string=sp::CRS(PROJ4)) %>%
    #only one CRS is accepted? should probs reproject wsboundary tp WGS84
    proj4string=sp::CRS('+proj=longlat +datum=WGS84')) %>%
    simplegeom()
spset_job = geoknife(stencil=pgon, fabric=spset, wait=TRUE)
spset_cutout = result(spset_job, with.units=TRUE)
check(spset_job)
head(spset_cutout)

#summarize by point
station = as.data.frame(t(site[,c('longitude', 'latitude')]))
colnames(station) = paste(site$region, site$sitecode, sep='_')
station = simplegeom(station)

spset_job = geoknife(stencil=station, fabric=spset, wait=TRUE)
spset_summ = result(spset_job, with.units=TRUE)
check(spset_job)
head(spset_summ)


#bind other MODIS data? ####
#bind Fluxnet data? ####

#bind site data to model output data; separate sp/powell ####

mods = mods %>%
    left_join(sites, 'site'='sitecode') %>%
    filter(! is.na(Name)) %>%
    distinct() %>%
    arrange(sitecode, year)

spmods = filter(mods, Source == 'StreamPULSE')
powmods = filter(mods, Source == 'USGS (Powell Center)')

#plot functions ####

gpp_er_biplot = function(mods, outfile){

    cnt = 0
    hold_cnt = FALSE
    errtypes = errmods = errs = list()

    pdf(file=outfile, onefile=TRUE)
    par(mfrow=c(3, 3), mar=c(3, 3, 1, 1), oma=c(1, 2, 1, 1))

    for(i in 1:nrow(mods)){

        m = mods[i,]
        res = try( request_results(sitecode=m$sitecode,
                    year=as.numeric(m$year)) )

        if(class(res) == 'try-error'){
            errtypes = append(errtypes, 'request_err')
            errmods = append(errmods, paste(m$sitecode, i))
            errs = append(errs, res)
            hold_cnt = TRUE
            next
        }

        if(nrow(res$predictions) == 0){
            errtypes = append(errtypes, 'no_predictions')
            errmods = append(errmods, paste(m$sitecode, i))
            errs = append(errs, 'NA')
        }

        er_vals = res$predictions$ER
        er_vals[er_vals > 0] = NA
        er_vals = gO2_to_gC(er_vals, PQ=1.25)
        er_vals = temp20_std(er_vals)
        res$predictions$ER = er_vals * -1
        gpp_vals = res$predictions$GPP
        gpp_vals[gpp_vals < 0] = NA
        res$predictions$GPP = gO2_to_gC(gpp_vals, PQ=1.25)
        n_days = sum(complete.cases(res$predictions[, c('GPP', 'ER')]))

        tryCatch({

            mout = lm(res$predictions$ER ~ res$predictions$GPP)
            mres = summary(mout)
            coeffs = round(unname(mres$coefficients[,'Estimate']), 2)
            icept = sprintf('%+.2f', coeffs[1])
            icept = paste(substr(icept, 0, 1), substr(icept, 2, nchar(icept)))
            coeff_lab = paste0('y = ', coeffs[2], 'x ', icept, ' (Adj. R^2: ',
                round(mres$adj.r.squared, 2), '; n days: ', n_days, ')')

            # plot(res$predictions$GPP, res$predictions$ER, main=m$sitecode,
            plot(res$predictions$GPP, res$predictions$ER, main='', xlab='GPP',
                ylab='ER', bty='l', col=alpha('gray30', 0.7), pch=20,
                xaxt='n', yaxt='n')
            axis(1, tck=-.02, labels=FALSE)
            axis(1, tcl=0, col='transparent', line=-0.5)
            axis(2, tck=-.02, labels=FALSE)
            axis(2, tcl=0, col='transparent', line=-0.5)
            abline(mout, lty=2, col='red', lwd=2)
            mtext(paste0(m$sitecode, ' (', m$year, ')\n', coeff_lab),
                side=3, line=0, cex=0.7)
            if(cnt %% 9 == 0){
                mtext(expression(paste(GPP[20]~"(gC"~"m"^"-2"~" d"^"-1"*')')),
                    side=1, outer=TRUE)
                mtext(expression(paste(ER[20]~"(gC"~"m"^"-2"~" d"^"-1"*')')),
                    side=2, outer=TRUE)
            }

        }, error=function(e){
                errtypes = append(errtypes, 'plot/stat_err')
                errmods = append(errmods, paste(m$sitecode, i))
                errs = append(errs, e)
            }
        )

        cnt = cnt + 1
        gc()
    }

    dev.off()

    return(list(errtypes=errtypes, errmosd=errmods, errs=errs))
}

ts_plot = function(mods, outfile, with_K){

    cnt = 0
    hold_cnt = FALSE
    errtypes2 = errmods2 = errs2 = list()

    pdf(file=outfile, onefile=TRUE)

    if(with_K){
        par(mfrow=c(3, 3), mar=c(3, 3, 1, 3), oma=c(1, 1, 1, 1))
    } else {
        par(mfrow=c(3, 3), mar=c(3, 3, 1, 1), oma=c(1, 1, 1, 1))
    }

    for(i in 1:nrow(mods)){

        m = mods[i,]
        res = try( request_results(sitecode=m$sitecode,
            year=as.numeric(m$year)) )

        if(class(res) == 'try-error'){
            errtypes2 = append(errtypes2, 'request_err')
            errmods2 = append(errmods2, paste(m$sitecode, i))
            errs2 = append(errs2, res)
            hold_cnt = TRUE
            next
        }

        if(nrow(res$predictions) == 0){
            errtypes2 = append(errtypes2, 'no_predictions')
            errmods2 = append(errmods2, paste(m$sitecode, i))
            errs2 = append(errs2, 'NA')
        }

        tryCatch({

            o = res$model_results$fit$daily

            llim = min(c(o$GPP_daily_2.5pct, o$ER_daily_2.5pct), na.rm=TRUE)
            ulim = max(c(o$GPP_daily_97.5pct, o$ER_daily_97.5pct), na.rm=TRUE)

            fulldateseq = seq(as.Date(paste0(m$year, '-01-01')),
                as.Date(paste0(m$year, '-12-31')), by='day')
            o = left_join(tibble(date=fulldateseq), o)
            plot(o$date, o$GPP_daily_mean, type='l', col='red', xlab='', las=0,
                ylab='', xaxs='i', yaxs='i', ylim=c(llim, ulim), bty='l', yaxt='n',
                xaxt='n')
            month_starts = o$date[substr(o$date, 9, 10) == '01']
            # all_month_starts_abb = paste0(sprintf('%02d', 1:12), '-01')
            # which_starts = all_month_starts_abb %in% substr(month_starts, 6, 10)
            axis(1, at=month_starts, labels=substr(month.abb, 1, 1),
                cex.axis=0.6)
            axis(2, tck=-.02, labels=FALSE)
            axis(2, tcl=0, col='transparent', line=-0.5)
            mtext(paste0(m$sitecode, ' (', m$year, ')'), side=3, line=0, cex=0.7)
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

            if(with_K){

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
            }


        }, error=function(e){
                errtypes2 = append(errtypes2, 'plot/stat_err')
                errmods2 = append(errmods2, paste(m$sitecode, i))
                errs2 = append(errs2, e)
        }, warning=function(w) NULL)

        cnt = cnt + 1
        gc()
    }

    dev.off()

    return(list(errtypes=errtypes2, errmosd=errmods2, errs=errs2))
}

#plots ####
gpp_er_biplot(spmods, 'output/gppXer_sp.pdf')
gpp_er_biplot(powmods, 'output/_gppXer_powell.pdf')
ts_plot(spmods, 'output/metab_ts_sp_noK.pdf', FALSE)
ts_plot(powmods, 'output/metab_ts_powell_noK.pdf', FALSE)
ts_plot(spmods, 'output/metab_ts_sp.pdf', TRUE)
ts_plot(powmods, 'output/metab_ts_powell.pdf', TRUE)

# barplot(metr$ann_GPP_C, ylab='', yaxt='n', yaxs='i',
#     width=0.03, space=0.5, xlim=c(0, 10))
# axis(2, las=2, line=-1.8, at=seq500)
# axis(2, las=2, line=-1.8, at=seq500,
#     tcl=0, col='white', lwd=2)
# mtext('Mean annual GPP (C)', side=2, line=2)

pdf(file='output/metab_dist.pdf', width=12, height=7)

par(mfrow=c(2, 1), mar=c(0, 3, 1, 1), oma=c(0, 1, 0, 0))
seq500 = seq(500, 3000, 500)

plot(metr$ann_GPP_C, ylab='', yaxt='n', yaxs='i', type='n', bty='n', xaxt='n')
segments(x0=1:nrow(metr), y0=rep(0, nrow(metr)), y1=metr$ann_GPP_C, lwd=3,
    lend=2)
axis(2, las=2, line=-1.8, at=seq500, xpd=NA)
axis(2, las=2, line=-1.8, tcl=0, col='white', lwd=2, at=seq500)
mtext(expression(paste("Mean annual GPP (g"~O[2]~"m"^"-2"~" d"^"-1"*')')),
    side=2, line=2)

par(mar=c(3, 3, 0, 1))

plot(metr$ann_ER_C, ylab='', yaxt='n', yaxs='i', type='n', bty='n', xaxt='n')
segments(x0=1:nrow(metr), y0=rep(0, nrow(metr)), y1=metr$ann_ER_C, lwd=3,
    lend=2, col='gray50')
axis(2, las=2, line=-1.8, at=seq500 * -1, xpd=NA, labels=rep('', 6))
axis(2, las=2, line=-1.8, at=seq500 * -1, labels=seq500 * -1,
    tcl=0, col='white', lwd=2)
mtext(expression(paste("Mean annual ER (g"~O[2]~"m"^"-2"~" d"^"-1"*')')),
    side=2, line=2)

dev.off()


#apply Philters and gapPhills; generate sub-datasets ####

filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.6')
saveRDS(filt, 'output/filtered_dsets/ERKunder60.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.4')
saveRDS(filt, 'output/filtered_dsets/ERKunder40.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 165')
saveRDS(filt, 'output/filtered_dsets/daysOver165.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 250')
saveRDS(filt, 'output/filtered_dsets/daysOver250.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 165', 'ER_K < 0.6')
saveRDS(filt, 'output/filtered_dsets/daysOver165_ERKunder60.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 165', 'ER_K < 0.4')
saveRDS(filt, 'output/filtered_dsets/daysOver165_ERKunder40.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 250', 'ER_K < 0.6')
saveRDS(filt, 'output/filtered_dsets/daysOver250_ERKunder60.rds')
filt = filter_and_impute(diag, models=metab_d, 'num_days > 250', 'ER_K < 0.4')
saveRDS(filt, 'output/filtered_dsets/daysOver250_ERKunder40.rds')

filt = readRDS('output/filtered_dsets/daysOver165_ERKunder40.rds')

smry = consolidate_list(filt) %>%
    as_tibble() %>%
    group_by(DOY) %>%
    summarize_all(list(median=~median(., na.rm=TRUE),
        quant25=~quantile(.)[2], quant75=~quantile(.)[4]))

#lips plots ####

pdf(file='output/lips_daysOver165_ERKunder40.pdf', width=12, height=7)

par(mfrow=c(2, 1), mar=c(0, 3, 1, 2), oma=c(0, 1, 0, 0))

gpplim = c(0, max(smry$GPP_C_filled_quant75, na.rm=TRUE))
plot(smry$DOY, smry$GPP_C_filled_median, ylab='', yaxs='i', type='l',
    bty='n', lwd=2, xlab='', ylim=gpplim, xaxs='i', xaxt='n', yaxt='n')
polygon(x=c(smry$DOY, rev(smry$DOY)),
    y=c(smry$GPP_C_filled_quant25, rev(smry$GPP_C_filled_quant75)),
    border=NA, col=alpha('red', alpha=0.6))
axis(2, las=2, line=0, xpd=NA, at=round(seq(0, gpplim[2], length.out=5), 1))
mtext(expression(paste("GPP (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
abline(h=0, lty=1, lwd=2, col='gray60')
legend('topright', legend='> 165 days; ER * K < 0.4', bty='n')

par(mar=c(4, 3, 0, 1))

erlim = c(min(smry$ER_C_filled_quant25, na.rm=TRUE), 0)
plot(smry$DOY, smry$ER_C_filled_median, ylab='', yaxs='i', type='l',
    bty='n', lwd=2, xlab='', ylim=erlim, xaxs='i', xaxt='n', yaxt='n')
polygon(x=c(smry$DOY, rev(smry$DOY)),
    y=c(smry$ER_C_filled_quant25, rev(smry$ER_C_filled_quant75)),
    border=NA, col=alpha('blue', alpha=0.6))
axis(2, las=2, line=0, at=round(seq(0, erlim[1], length.out=5), 1))
axis(1, line=0, at=seq(0, max(smry$DOY), 30))
mtext(expression(paste("ER (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
mtext('DOY', side=1, line=2.5)

dev.off()
