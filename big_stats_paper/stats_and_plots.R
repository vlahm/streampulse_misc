library(plyr)
library(StreamPULSE)
library(tidyverse)
library(RColorBrewer)
library(plotrix)
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

# model results from streampulse portal
mods = query_available_results('all')[[1]] %>%
    as_tibble() %>%
    mutate_all(as.character)

#datasets from Phil
metab_d = readRDS('phil_stuff/output/synthesis_standardized.rds')
# metab_d = metab_d[! grepl('nwis', names(metab_d))]
names(metab_d) = phil_to_mike_format(tibble(Site_ID=names(metab_d)), mods,
        arrange=FALSE) %>%
    pull(sitecode)
erXk_p_vals = sapply(metab_d, function(x){
    summary(lm(x$ER~x$K600))$coefficients[2,4]
})
erXk_p_vals = tibble(sitecode=names(erXk_p_vals), er_k_pval=erXk_p_vals)
diag = as_tibble(readRDS('phil_stuff/metab_synthesis/output/yearly_diagnostics.rds'))
diag = filter(diag, ! is.na(ER_K))
diag$ER_K = abs(diag$ER_K)
diag = phil_to_mike_format(diag, mods) %>%
    select(-site, -region)
# diag = filter(diag, ! grepl('nwis', Site_ID)) %>%
#     rename(sitecode=Site_ID)
metr = as_tibble(readRDS('phil_stuff/output/site_metrics.rds'))
metr = readRDS('phil_stuff/output/metrics_bound.rds') %>%
    select(-one_of('Name', 'Source', 'Lat', 'Lon', 'COMID', 'VPU')) %>%
    as_tibble() %>%
    left_join(metr, by='Site_ID')
metr = phil_to_mike_format(metr, mods) %>%
    select(-site, -region)
# metr = filter(metr, ! grepl('nwis', Site_ID)) %>%
#     rename(sitecode=Site_ID)

metr = diag %>%
    group_by(sitecode) %>%
    summarize(ER_K=mean(ER_K, na.rm=TRUE), GPP_neg=mean(GPP_neg, na.rm=TRUE),
        ER_pos=mean(ER_pos, na.rm=TRUE), num_days=sum(num_days, na.rm=TRUE)) %>%
    ungroup() %>%
    right_join(metr, by='sitecode') %>%
    left_join(erXk_p_vals, by='sitecode') %>%
    select(-ER_K, -GPP_neg, -ER_pos, -num_days, -er_k_pval, everything())

mods = mutate(mods, sitecode=paste(region, site, sep='_'))

#subset of streampulse + powell sites used in this analysis
sites = as_tibble(readRDS('sites_COMID.rds'))
sites = phil_to_mike_format(sites, mods)
    # rename(sitecode=Site_ID)

WGS84 = 4326 #EPSG code for coordinate reference system

# 0 bind NHDPlusV2 data ####

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

#bind StreamCat data (superfluous since Phil added his own streamcat pull) ####

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

# 0 bind Phil's StreamCat dataset####
sites = left_join(sites, metr, by='sitecode')
sites = sites[! duplicated(sites$sitecode),]

#write variable key table ####

varnames = colnames(sites)
varnames = varnames[! varnames %in%
    c('Name', 'Source', 'Lat', 'Lon', 'COMID', 'VPU', 'site', 'region', 'sitecode')]

vartypes = rep('watershed', length(varnames))
vartypes[46:65] = 'metabolism'
vartypes[66:93] = 'stream phys'
vartypes[7:45] = 'watershed'
vartypes[grep('Cat', varnames)] = 'catchment'
vartypes[grep('Rp100', varnames)] = 'riparian100'
vartypes[c(1:6, 14, 29, 30)] = 'stream channel'
vartypes[94:98] = 'diagnostic'

varsources = rep('StreamCat', length(varnames))
varsources[48:49] = 'MODIS'
varsources[2:12] = 'NHDPlusV2'
varsources[c(1, 13, 14)] = 'NHD derived'
varsources[c(46:65, 94:98)] = 'model derived'
varsources[66:93] = 'other'

vardesc = rep(paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
    'NHDPlusLandscapeAttributes/StreamCat/Documentation/VariableList-',
    'QuickReference.html'), length(varnames))
vardesc[varsources == 'NHDPlusV2'] = paste0('user guide pdf link here: ',
    'http://www.horizon-systems.com/NHDPlus/NHDPlusV2_documentation.php',
    '#NHDPlusV2%20User%20Guide')
vardesc[1] = 'location of site along reach (in downstream direction) as proportion of reach length'
vardesc[13] = 'shrinking factor applied to areal metrics (based on reach_proportion)'
vardesc[c(14, 48:98)] = ''
vardesc[46] = 'annual gpp (gC m-2 d-1), converted from gO2 via PQ=1.25'
vardesc[47] = 'annual er (gC m-2 d-1), converted from gO2 via PQ=1.25'

variable_key = tibble(name=varnames, type=vartypes, source=varsources,
        description=vardesc) %>%
    arrange(vartypes, varsources)

write.csv(variable_key, 'output/variable_key.csv', row.names=FALSE)

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

# 0 bind site data to model output data; separate sp/powell ####

mods = mods %>%
    left_join(sites, by='sitecode') %>%
    filter(! is.na(Name)) %>%
    distinct() %>%
    arrange(sitecode, year)

spmods = filter(mods, Source == 'StreamPULSE')
powmods = filter(mods, Source == 'USGS (Powell Center)')

#apply Philters and gapPhills; generate sub-datasets ####

filt = filter_and_impute(diag, models=metab_d, 'ER_K <= 1')
# sum(sapply(filt, function(x) ! is.null(x)))
saveRDS(filt, 'output/filtered_dsets/no_filter.rds')

# filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.1')
# saveRDS(filt, 'output/filtered_dsets/ERKunder10.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.2')
saveRDS(filt, 'output/filtered_dsets/ERKunder20.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.3')
saveRDS(filt, 'output/filtered_dsets/ERKunder30.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.4')
saveRDS(filt, 'output/filtered_dsets/ERKunder40.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.5')
saveRDS(filt, 'output/filtered_dsets/ERKunder50.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.6')
saveRDS(filt, 'output/filtered_dsets/ERKunder60.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.7')
saveRDS(filt, 'output/filtered_dsets/ERKunder70.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.8')
saveRDS(filt, 'output/filtered_dsets/ERKunder80.rds')
filt = filter_and_impute(diag, models=metab_d, 'ER_K < 0.9')
saveRDS(filt, 'output/filtered_dsets/ERKunder90.rds')

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

# 0 plot functions ####

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
        res$predictions$ER = er_vals * -1
        # daily_mean_temp = res$model_results$data %>%
        #     group_by(date) %>%
        #     summarize_at(vars(one_of('temp.water')), list(~mean(., na.rm=TRUE))) %>%
        #     ungroup() %>%
        #     right_join(select(res$predictions, date), by='date') %>%
        #     pull(temp.water)
        # res$predictions$ER = temp20_std(er_vals, daily_mean_temp)

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
                round(mres$adj.r.squared, 2), ')')

            # plot(res$predictions$GPP, res$predictions$ER, main=m$sitecode,
            plot(res$predictions$GPP, res$predictions$ER, main='', xlab='GPP',
                ylab='ER', bty='l', col=alpha('gray30', 0.7), pch=20,
                xaxt='n', yaxt='n')
            axis(1, tck=-.02, labels=FALSE)
            axis(1, tcl=0, col='transparent', line=-0.5)
            axis(2, tck=-.02, labels=FALSE)
            axis(2, tcl=0, col='transparent', line=-0.5)
            abline(mout, lty=2, col='red', lwd=2)
            mtext(paste0(m$sitecode, ' (', m$year, '; n days: ', n_days,
                ')\n', coeff_lab),
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

        res$predictions$ER[res$predictions$ER > 0] = NA
        res$predictions$GPP[res$predictions$GPP < 0] = NA
        n_days = sum(complete.cases(res$predictions[, c('GPP', 'ER')]))

        tryCatch({

            o = select(res$model_results$fit$daily, date, ER_daily_mean,
                ER_daily_2.5pct, ER_daily_97.5pct, GPP_daily_mean,
                GPP_daily_2.5pct, GPP_daily_97.5pct, K600_daily_mean,
                K600_daily_2.5pct, K600_daily_97.5pct)

            bogus_bool = o$ER_daily_mean > 0 | o$GPP_daily_mean < 0
            bogus_bool[is.na(bogus_bool)] = TRUE
            o[bogus_bool, -1] = rep(NA, ncol(o) - 1)
            o$ER_daily_mean = gO2_to_gC(o$ER_daily_mean, PQ=1.25)
            o$ER_daily_2.5pct = gO2_to_gC(o$ER_daily_2.5pct, PQ=1.25)
            o$ER_daily_97.5pct = gO2_to_gC(o$ER_daily_97.5pct, PQ=1.25)
            # daily_mean_temp = res$model_results$data %>%
            #     group_by(date) %>%
            #     summarize_at(vars(one_of('temp.water')), list(~mean(., na.rm=TRUE))) %>%
            #     ungroup() %>%
            #     right_join(select(res$predictions, date), by='date') %>%
            #     pull(temp.water)
            # res$predictions$ER = temp20_std(er_vals, daily_mean_temp)

            o$GPP_daily_mean = gO2_to_gC(o$GPP_daily_mean, PQ=1.25)
            o$GPP_daily_2.5pct = gO2_to_gC(o$GPP_daily_2.5pct, PQ=1.25)
            o$GPP_daily_97.5pct = gO2_to_gC(o$GPP_daily_97.5pct, PQ=1.25)

            llim = min(c(o$GPP_daily_2.5pct, o$ER_daily_2.5pct), na.rm=TRUE)
            ulim = max(c(o$GPP_daily_97.5pct, o$ER_daily_97.5pct), na.rm=TRUE)

            fulldateseq = seq(as.Date(paste0(m$year, '-01-01')),
                as.Date(paste0(m$year, '-12-31')), by='day')
            o = left_join(tibble(date=fulldateseq), o)
            plot(o$date, o$GPP_daily_mean, type='l', col='forestgreen', xlab='',
                las=0, ylab='', xaxs='i', yaxs='i', ylim=c(llim, ulim), bty='l',
                yaxt='n', xaxt='n')
            month_starts = o$date[substr(o$date, 9, 10) == '01']
            # all_month_starts_abb = paste0(sprintf('%02d', 1:12), '-01')
            # which_starts = all_month_starts_abb %in% substr(month_starts, 6, 10)
            axis(1, at=month_starts, labels=substr(month.abb, 1, 1),
                cex.axis=0.6)
            axis(2, tck=-.02, labels=FALSE)
            axis(2, tcl=0, col='transparent', line=-0.5)
            mtext(paste0(m$sitecode, ' (', m$year, '; n days: ', n_days, ')'),
                side=3, line=0, cex=0.7)
            lines(o$date, o$ER_daily_mean, col='sienna')
            mtext(expression(paste("gC"~"m"^"-2"~" d"^"-1")), side=2,
                line=1.5, font=2, cex=0.7)

            rl = rle(is.na(o$GPP_daily_2.5pct))
            vv = !rl$values
            chunkfac = rep(cumsum(vv), rl$lengths)
            chunkfac[chunkfac == 0] = 1
            chunks = split(o, chunkfac)
            noNAchunks = lapply(chunks, function(x) x[!is.na(x$GPP_daily_2.5pct),] )

            for(j in 1:length(noNAchunks)){
                polygon(x=c(noNAchunks[[j]]$date, rev(noNAchunks[[j]]$date)),
                    y=c(noNAchunks[[j]]$GPP_daily_2.5pct,
                        rev(noNAchunks[[j]]$GPP_daily_97.5pct)),
                    col=adjustcolor('forestgreen', alpha.f=0.3), border=NA)
                polygon(x=c(noNAchunks[[j]]$date, rev(noNAchunks[[j]]$date)),
                    y=c(noNAchunks[[j]]$ER_daily_2.5pct,
                        rev(noNAchunks[[j]]$ER_daily_97.5pct)),
                    col=adjustcolor('sienna', alpha.f=0.3), border=NA)
            }

            abline(h=0, lty=3, col='gray50')

            if(with_K){

                par(new=TRUE)
                # llim2 = min(o$K600_daily_2.5pct, na.rm=TRUE)
                # ulim2 = max(o$K600_daily_97.5pct, na.rm=TRUE)
                ylims = range(o$K600_daily_mean, na.rm=TRUE)
                ylim_expansion = ylims[2] - mean(ylims)
                ylims = c(ylims[1] - ylim_expansion, ylims[2] + ylim_expansion)
                plot(o$date, o$K600_daily_mean, col='black',
                    type='l', xlab='', las=0, ylab='', xaxs='i', yaxs='i',
                    xaxt='n', bty='u', yaxt='n', ylim=ylims)
                    # xaxt='n', bty='u', yaxt='n', ylim=c(llim2, ulim2))
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

dist_plot = function(m, outfile, plot_comm_resp=FALSE){

    m = arrange(m, desc(ann_GPP_C))

    pdf(file=outfile, width=12, height=7)

    if(plot_comm_resp){
        par(mfrow=c(3, 1), mar=c(0, 3, 1, 1), oma=c(0, 1, 0, 0))
    } else {
        par(mfrow=c(2, 1), mar=c(0, 3, 1, 1), oma=c(0, 1, 0, 0))
    }

    seq500 = seq(500, 3000, 500)

    plot(m$ann_GPP_C, ylab='', yaxt='n', yaxs='i', type='n', bty='n', xaxt='n')
    segments(x0=1:nrow(m), y0=rep(0, nrow(m)), y1=m$ann_GPP_C, lwd=3,
        lend=2, col='forestgreen')
    axis(2, las=2, line=-1.8, at=seq500, xpd=NA)
    axis(2, las=2, line=-1.8, tcl=0, col='white', lwd=2, at=seq500)
    mtext(expression(paste("Mean annual GPP (gC"~"m"^"-2"~" d"^"-1"*')')),
        side=2, line=2)

    par(mar=c(3, 3, 0, 1))

    plot(m$ann_ER_C, ylab='', yaxt='n', yaxs='i', type='n', bty='n', xaxt='n')
    segments(x0=1:nrow(m), y0=rep(0, nrow(m)), y1=m$ann_ER_C, lwd=3,
        lend=2, col='sienna')
    axis(2, las=2, line=-1.8, at=seq500 * -1, xpd=NA, labels=rep('', 6))
    axis(2, las=2, line=-1.8, at=seq500 * -1, labels=seq500 * -1,
        tcl=0, col='white', lwd=2)
    mtext(expression(paste("Mean annual ER (gC"~"m"^"-2"~" d"^"-1"*')')),
        side=2, line=2)

    if(plot_comm_resp){
        comm_resp = m$ann_ER_C - m$ann_GPP_C
        plot(rev(sort(comm_resp, na.last=FALSE)), ylab='', yaxt='n', yaxs='i',
            type='n', bty='n', xaxt='n')
        segments(x0=1:nrow(m), y0=rep(0, nrow(m)), y1=comm_resp, lwd=3,
            lend=2, col='black')
        axis(2, las=2, line=-1.8)
        axis(2, las=2, line=-1.8,
        # axis(2, las=2, line=-1.8, at=seq500 * -1, xpd=NA, labels=rep('', 6))
        # axis(2, las=2, line=-1.8, at=seq500 * -1, labels=seq500 * -1,
            tcl=0, col='white', lwd=2)
        mtext(expression(paste("Mean annual CR (gC"~"m"^"-2"~" d"^"-1"*')')),
            side=2, line=2)
    }

    dev.off()
}

er_k_filter_plot = function(diagnostics, outfile){

    pdf(file=outfile, width=7, height=6)

    r2_seq = seq(1, 0.05, -0.05)
    unfilt_sites = rep(NA, length(r2_seq))

    for(i in 1:length(r2_seq)){
        unfilt_sites[i] = nrow(filter(diagnostics, ER_K < !!(r2_seq[i])))
    }

    plot(rev(r2_seq), unfilt_sites, ylab='Modeled Siteyears',
        xlab=expression(paste('Maximum |R'^2 * '|')), xaxt='n', yaxt='n')
    xseq = seq(0.1, 1, 0.1)
    axis(1, at=xseq - 0.05, labels=rev(xseq))
    axis(2, cex.axis=0.8)

    dev.off()
}

lips_plot = function(readfile, datlist, diagnostics, sitedata, quant_filt=NULL,
    standalone, outfile, ...){

    #quant_filt example: 'width_calc > 0.25'
    if(standalone){
        filt = readRDS(readfile)
    } else {
        filt = datlist
    }

    nsites_included = sum(sapply(filt, nrow) != 0)

    var_quant_filt = NULL
    if(! is.null(quant_filt)){
        quant_comp = strsplit(quant_filt, ' ')[[1]]
        qf = quantile(sitedata[, quant_comp[1], drop=TRUE], na.rm=TRUE,
            probs=as.numeric(quant_comp[3]))
        filt_sites = sitedata %>%
            filter_(paste(quant_comp[1], quant_comp[2], qf)) %>%
            pull(sitecode)
        filt = filt[names(filt) %in% filt_sites]

        var_quant_filt = paste0(quant_comp[1], ' ', quant_comp[2], ' ',
            as.numeric(quant_comp[3]) * 100, '%')
    }

    smry = consolidate_list(filt) %>%
        as_tibble() %>%
        group_by(DOY) %>%
        summarize_all(list(median=~median(., na.rm=TRUE),
            quant25=~quantile(., na.rm=TRUE)[2],
            quant75=~quantile(., na.rm=TRUE)[4]))

    if(standalone){
        pdf(file=outfile, width=12, height=7)
        par(mfrow=c(2, 1), oma=c(0, 1, 0, 0))
    }
    par(mar=c(0, 3, 1, 1), lend=2)

    gpplim = c(0, max(smry$GPP_C_filled_quant75, na.rm=TRUE))
    erlim = c(min(smry$ER_C_filled_quant25, na.rm=TRUE), 0)
    gpplim = c(0, max(gpplim[2], abs(erlim[1])))
    erlim = c(min(abs(gpplim[2]) * -1, erlim[1]), 0)

    plot(smry$DOY, smry$GPP_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=2, xlab='', ylim=gpplim, xaxs='i', xaxt='n', yaxt='n')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$GPP_C_filled_quant25, rev(smry$GPP_C_filled_quant75)),
        border=NA, col=alpha('forestgreen', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, gpplim[2], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, gpplim[2], length.out=5), 1))
    abline(h=0, lty=1, lwd=2, col='gray60')
    medsums = round(colSums(select(smry, contains('median'))), 1)

    if(standalone){
        mtext(expression(paste("GPP (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
        legend('topleft', legend=c('Median', '', '25-75%', '', 'NEP Median'),
            col=c('darkgreen', 'sienna4', alpha('forestgreen', alpha=0.6),
                alpha('sienna', alpha=0.6), 'black'),
            bty='n', lty=1, lwd=c(2, 2, 10, 10, 2))
        legend('topright', title='Filters', bty='n', title.col='gray30',
            lty=1, seg.len=0.2, lwd=2, legend=c(..., var_quant_filt))
    }

    legend('right', title='Cumulative\nMedian Sums', bty='n',
        legend=c(paste('GPP:', medsums[1]), paste('ER:', medsums[2]),
            paste('NEP:', medsums[3])), title.col='gray30')
    legend('left', paste('Sites included:', nsites_included), bty='n')

    par(mar=c(3, 3, 0, 1))

    plot(smry$DOY, smry$ER_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=2, xlab='', ylim=erlim, xaxs='i', xaxt='n', yaxt='n')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$ER_C_filled_quant25, rev(smry$ER_C_filled_quant75)),
        border=NA, col=alpha('sienna', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(1, line=0, tck=-.02, labels=FALSE, at=seq(0, max(smry$DOY), 30))
    axis(1, line=-0.5, tcl=0, col='transparent', at=seq(0, max(smry$DOY), 30))
    lines(smry$DOY, smry$NEP_C_filled_median, col='black', lwd=2, xpd=NA)

    if(standalone){
        mtext(expression(paste("ER (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
        mtext('DOY', side=1, line=2)
        dev.off()
    }
}

lips_facet = function(readfile, diagnostics, sitedata, split_vars,
    outfile, ...){

    pdf(file=outfile, width=12, height=8)

    layoutmat = matrix(c(
        17,17,18,18,0,5,5,5,5,
        1,1,1,1,13,6,6,6,6,
        2,2,2,2,14,7,7,7,7,
        0,19,19,0,21,8,8,8,8,
        0,0,0,0,0,9,9,9,9,
        3,3,3,3,15,10,10,10,10,
        4,4,4,4,16,11,11,11,11,
        0,20,20,0,22,12,12,12,12),
        nrow=8, ncol=9, byrow=TRUE)
    layout(layoutmat)

    filt = readRDS(readfile)
    fpath = strsplit(readfile, '/')[[1]]
    filters = strsplit(fpath[length(fpath)], '\\.')[[1]][1]
    filt_vec = strsplit(filters, '_')[[1]]
    filt_disp = c()

    dayfilt = grep('days', filt_vec)
    if(length(dayfilt)){
        dispy = paste('>', c(str_match(filt_vec[dayfilt], '[0-9]+')), 'days')
        filt_disp = append(filt_disp, dispy)
    }

    erkfilt = grep('ERK', filt_vec)
    if(length(erkfilt)){
        dispy = paste('ER * K <',
            as.numeric(c(str_match(filt_vec[erkfilt], '[0-9]+'))) / 100)
        filt_disp = append(filt_disp, dispy)
    }

    # quant_comp = strsplit(quant_filt, ' ')[[1]]
    q75 = quantile(sitedata[, split_vars[1], drop=TRUE], na.rm=TRUE, probs=0.75)
    upper75sites = sitedata %>%
        filter_(paste(split_vars[1], '>', q75)) %>%
        pull(sitecode)
    upper75 = filt[names(filt) %in% upper75sites]
    upper75 = upper75[sapply(upper75, length) != 0]
    upper75sitedata = filter(sitedata, sitecode %in% upper75sites)

    q25 = quantile(sitedata[, split_vars[1], drop=TRUE], na.rm=TRUE, probs=0.25)
    lower25sites = sitedata %>%
        filter_(paste(split_vars[1], '>', q25)) %>%
        pull(sitecode)
    lower25 = filt[names(filt) %in% lower25sites]
    lower25 = lower25[sapply(lower25, length) != 0]
    lower25sitedata = filter(sitedata, sitecode %in% lower25sites)

    lips_plot(readfile='placeholder', datlist=filt, diagnostics=diagnostics,
        sitedata=sitedata, quant_filt=paste(split_vars[1], '> 0.75'),
        standalone=FALSE, outfile='placeholder')
    lips_plot(readfile='placeholder', datlist=filt, diagnostics=diagnostics,
        sitedata=sitedata, quant_filt=paste(split_vars[1], '< 0.25'),
        standalone=FALSE, outfile='placeholder')
    #split 2a
    lips_plot(readfile='placeholder', datlist=upper75, diagnostics=diagnostics,
        sitedata=upper75sitedata, quant_filt=paste(split_vars[2], '> 0.75'),
        standalone=FALSE, outfile='placeholder')
    lips_plot(readfile='placeholder', datlist=upper75, diagnostics=diagnostics,
        sitedata=upper75sitedata, quant_filt=paste(split_vars[2], '< 0.25'),
        standalone=FALSE, outfile='placeholder')
    #split 2b
    lips_plot(readfile='placeholder', datlist=lower25, diagnostics=diagnostics,
        sitedata=lower25sitedata, quant_filt=paste(split_vars[2], '> 0.75'),
        standalone=FALSE, outfile='placeholder')
    lips_plot(readfile='placeholder', datlist=lower25, diagnostics=diagnostics,
        sitedata=lower25sitedata, quant_filt=paste(split_vars[2], '< 0.25'),
        standalone=FALSE, outfile='placeholder')

    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    segments(0, 0.25, 1, 1)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    segments(0, 0.75, 1, 0)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    segments(0, 0.25, 1, 1)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    segments(0, 0.75, 1, 0)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    if(length(filt_disp) == 1){
        text('Prefilters:', x=0.2, y=0.5, pos=4, font=2)
        text(filt_disp, x=0.2, y=0.3, pos=4)
    } else if(length(filt_disp) == 2){
        text('Prefilters:', x=0.2, y=0.6, pos=4, font=2)
        text(filt_disp[1], x=0.2, y=0.4, pos=4)
        text(filt_disp[2], x=0.2, y=0.2, pos=4)
    }
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    legend('left', legend=c('Median', '', '25-75%', '', 'NEP Median'),
        col=c('darkgreen', 'sienna4', alpha('forestgreen', alpha=0.6),
            alpha('sienna', alpha=0.6), 'black'),
        bty='n', lty=1, lwd=c(2, 2, 7, 7, 2))
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    text(paste('Upper quartile of', split_vars[1]), x=0.5, y=0.8, cex=1.2, font=2)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    text(paste('Lower quartile of', split_vars[1]), x=0.5, y=0.8, cex=1.2, font=2)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    text(paste('Upper and\nlower quartiles of\n', split_vars[2]),
        x=0.5, y=0.5, cex=0.8, font=2)
    plot(1, 1, ann=FALSE, axes=FALSE, type='n', xlim=c(0,1), ylim=c(0,1))
    text(paste('Upper and\nlower quartiles of\n', split_vars[2]),
        x=0.5, y=0.5, cex=0.8, font=2)

    dev.off()
}

light_gpp_biplot = function(sitedata, outfile, colvar, quantvar=NULL,
        xl, xr, yb, yt, txt_x, txt_y){

    pdf(file=outfile, width=8, height=8)

    if(! is.null(quantvar)){

        s = select(sitedata, sitecode, Stream_PAR_sum, ann_GPP_C,
                one_of(quantvar, colvar)) %>%
            filter_at(vars(sitecode, Stream_PAR_sum, ann_GPP_C,
                one_of(quantvar, colvar)),
                all_vars(! is.na(.)))

        colvar_rng = round(range(s[[colvar]], na.rm=TRUE), 2)
        cols = brewer.pal(9, 'Blues')[2:9]
        colfac = cut(s[[colvar]], 8)
        levels(colfac) = cols

        qfs = quantile(sitedata[[quantvar]], na.rm=TRUE,
            probs=c(0.25, 0.75))
        upquant_ind = s[[quantvar]] > qfs[2]
        dnquant_ind = s[[quantvar]] < qfs[1]
        s_up = s[upquant_ind, ]
        s_dn = s[dnquant_ind, ]
        cols_up = as.character(colfac[upquant_ind])
        cols_dn = as.character(colfac[dnquant_ind])

        par(mfrow=c(2, 1), mar=c(4, 1, 1, 1), oma=c(0, 4, 0, 0))

        plot(s_up$Stream_PAR_sum, s_up$ann_GPP_C, col=cols_up, xlab='', ylab='',
            main=paste('Site filter:', quantvar, '> 75th pctile'))
        color.legend(xl=xl, xr=xr, yb=yb, yt=yt, align='rb', pos=4, offset=0.2,
            # color.legend(xl=2, xr=2.5, yb=2400, yt=2900, align='rb', pos=4, offset=0.2,
            rect.col=cols, gradient='y', legend=c(colvar_rng[1], colvar_rng[2]),
            cex=0.9)
        text(colvar, x=txt_x, y=txt_y, pos=4, offset=0, cex=0.9)
        # text('AR-1 of Mean Annual', x=2, y=3170, pos=4, offset=0, cex=0.9)
        # text(expression(paste('Discharge (m'^3 * 's)')),
        #     x=2, y=3030, pos=4, offset=0, cex=0.9)
        plot(s_dn$Stream_PAR_sum, s_dn$ann_GPP_C, col=cols_dn, xlab='', ylab='',
            main=paste('Site filter:', quantvar, '< 25th pctile'))
        mtext(expression(paste('Mean Annual GPP (gC'~'m'^'-2'~' d'^'-1'*')')),
            side=2, line=2.5, outer=TRUE)
        mtext('Mean Annual Stream PAR', side=1, line=2.5)

    } else {

        par(mar=c(4, 4, 1, 1))

        s = select(sitedata, Stream_PAR_sum, ann_GPP_C, one_of(colvar)) %>%
            filter_at(vars(Stream_PAR_sum, ann_GPP_C, one_of(colvar)),
                all_vars(! is.na(.)))

        colvar_rng = round(range(s[[colvar]], na.rm=TRUE), 2)
        cols = brewer.pal(9, 'Blues')[2:9]
        colfac = cut(s[[colvar]], 8)
        levels(colfac) = cols

        plot(s$Stream_PAR_sum, s$ann_GPP_C, xlab='', ylab='',
            col=as.character(colfac))
        text(colvar, x=txt_x, y=txt_y, pos=4, offset=0, cex=0.9)
        mtext(expression(paste('Mean Annual Stream GPP (gC'~'m'^'-2'~' d'^'-1'*')')),
            side=2, line=2.5)
        mtext('Mean Annual Stream PAR', side=1, line=2.5)

        colvar_rng = round(range(s[[colvar]], na.rm=TRUE), 2)
        color.legend(xl=xl, xr=xr, yb=yb, yt=yt, align='rb', pos=4, offset=0.2,
            rect.col=cols, gradient='y', legend=c(colvar_rng[1], colvar_rng[2]),
            cex=0.9)
    }

    dev.off()
}

wtemp_er_biplot = function(sitedata, outfile, signif=TRUE){

    if(signif){
        filtstr = 'er_k_pval < 0.05'
    } else {
        filtstr = 'er_k_pval >= 0.05'
    }

    s = sitedata %>%
        filter_(filtstr) %>%
        select(Wtemp_mean, ann_ER_C) %>%
        filter_at(vars(Wtemp_mean, ann_ER_C), all_vars(! is.na(.)))

    pdf(file=outfile, width=8, height=8)
    par(mar=c(4, 4, 3, 1))

    plot(s$Wtemp_mean, s$ann_ER_C, xlab='', ylab='',
        main='Filter: K600~ER p < 0.05')
    mtext(expression(paste('Mean Annual ER (gC'~'m'^'-2'~' d'^'-1'*')')),
        side=2, line=2.5)
    mtext('Mean Annual Water Temp (C)', side=1, line=2.5)

    dev.off()
}
crude_plot = function(x, y, outfile){
    pdf(file=outfile, width=8, height=8)
    plot(x, y)
    dev.off()
}
crude_plot(metr$width_calc, metr$Disch_ar1, 'output/width_discharge.pdf')
crude_plot(metr$MOD_ann_ER, metr$er_C_mean, 'output/MODISer_streamER.pdf')
crude_plot(log(metr$MOD_ann_GPP), log(metr$gpp_C_mean), 'output/logMODISgpp_logStreamGPP.pdf')
gpp_modis_stream_biplot = function(sitedata, outfile, quantvar=NULL){

    pdf(file=outfile, width=8, height=8)
    par(mar=c(4, 4, 1, 1))

    if(! is.null(quantvar)){

        s = select(sitedata, sitecode, MOD_ann_GPP, ann_GPP_C, one_of(quantvar)) %>%
            filter_at(vars(sitecode, MOD_ann_GPP, ann_GPP_C, one_of(quantvar)),
                all_vars(! is.na(.)))

        qfs = quantile(sitedata[, quantvar, drop=TRUE], na.rm=TRUE,
            probs=c(0.25, 0.75))
        upper75 = sitedata %>%
            filter_(paste(quantvar, '>', qfs[2])) %>%
            pull(sitecode)
        lower25 = sitedata %>%
            filter_(paste(quantvar, '<', qfs[1])) %>%
            pull(sitecode)

        s_up = filter(s, sitecode %in% upper75)
        s_dn = filter(s, sitecode %in% lower25)

        par(mfrow=c(2, 1), mar=c(4, 1, 1, 1), oma=c(0, 4, 0, 0))
        plot(s_up$MOD_ann_GPP, s_up$ann_GPP_C, xlab='', ylab='',
            main=paste('Site filter:', quantvar, '> 75th pctile'))
        plot(s_dn$MOD_ann_GPP, s_dn$ann_GPP_C, xlab='', ylab='',
            main=paste('Site filter:', quantvar, '< 25th pctile'))
        mtext(expression(paste('Mean Annual Stream GPP (gC'~'m'^'-2'~' d'^'-1'*')')),
            side=2, line=2.5, outer=TRUE)
        mtext('Mean Annual MODIS GPP', side=1, line=2.5)

    } else {
        s = select(sitedata, MOD_ann_GPP, ann_GPP_C) %>%
            filter_at(vars(MOD_ann_GPP, ann_GPP_C),
                all_vars(! is.na(.)))

        plot(s$MOD_ann_GPP, s$ann_GPP_C, xlab='', ylab='')
        mtext(expression(paste('Mean Annual Stream GPP (gC'~'m'^'-2'~' d'^'-1'*')')),
            side=2, line=2.5)
        mtext('Mean Annual MODIS GPP', side=1, line=2.5)
    }

    dev.off()
}

#plots ####

gpp_er_biplot(spmods, 'output/gppXer_sp.pdf')
gpp_er_biplot(powmods, 'output/gppXer_powell.pdf')
# ts_plot(spmods, 'output/metab_ts_sp_noK.pdf', FALSE)
# ts_plot(powmods, 'output/metab_ts_powell_noK.pdf', FALSE)
ts_plot(spmods, 'output/metab_ts_sp.pdf', TRUE)
ts_plot(powmods, 'output/metab_ts_powell.pdf', TRUE)
dist_plot(metr, 'output/metab_dist.pdf')
er_k_filter_plot(diag, 'output/erXk_corr_filter.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='width_calc <= 1',
    standalone=TRUE, outfile='output/lips_overall.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Disch_ar1 > 0.75',
    standalone=TRUE, outfile='output/lips_Qar1_75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='BFIWs <= 1',
    standalone=TRUE, outfile='output/lips_BFI.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='BFIWs <= 1',
    standalone=TRUE, outfile='output/lips_BFI.pdf')
light_gpp_biplot(sites, 'output/parXgpp2.pdf', colvar='Disch_ar1',
    quantvar='Disch_ar1', xl=2, xr=2.5, yb=2400, yt=2900, txt_x=10, txt_y=3000)
light_gpp_biplot(sites, 'output/parXgpp2.pdf', colvar='Disch_ar1',
    quantvar=NULL, xl=10, xr=10.2, yb=2400, yt=2900, txt_x=10, txt_y=3000)
wtemp_er_biplot(sites, 'output/wtempXer_signif.pdf', signif=TRUE)
# wtemp_er_biplot(sites, 'output/wtempXer_nonsignif.pdf', signif=FALSE)
lips_facet('output/filtered_dsets/no_filter.rds', diag, sites,
    c('MOD_ann_GPP', 'ann_GPP_C'), 'output/lipset_terrGPP_aqGPP.pdf')
lips_facet('output/filtered_dsets/no_filter.rds', diag, sites,
    c('width_calc', 'Disch_ar1'), 'output/lipset_width_Qar1.pdf')
lips_facet('output/filtered_dsets/no_filter.rds', diag, sites,
    c('width_calc', 'BFIWs'), 'output/lipset_width_bfi.pdf')
lips_facet('output/filtered_dsets/daysOver165.rds', diag, sites,
    c('width_calc', 'BFIWs'), 'output/lipset_daysOver165_width_bfi.pdf')
lips_facet('output/filtered_dsets/ERKunder40.rds', diag, sites,
    c('width_calc', 'BFIWs'), 'output/lipset_ERKUnder40_width_bfi.pdf')
lips_facet('output/filtered_dsets/daysOver165_ERKunder40.rds', diag, sites,
    c('width_calc', 'BFIWs'), 'output/lipset_daysOver165_ERKUnder40_width_bfi.pdf')
gpp_modis_stream_biplot(sitedata=sites, outfile='output/modis_stream_biplot.pdf',
    quantvar=NULL)
gpp_modis_stream_biplot(sites, 'output/modis_stream_biplot_width.pdf', 'width_calc')


pcaset = metr %>%
    select(gpp_C_mean, er_C_mean, PAR_sum, Disch_ar1, width_calc, MOD_ann_GPP,
        MOD_ann_ER) %>%
    mutate(log_gpp_C_mean=log(gpp_C_mean), log_er_C_mean=log(gpp_C_mean)) %>%
    select(-gpp_C_mean, er_C_mean)
pcaset = pcaset[complete.cases(pcaset),]
out = prcomp(pcaset, scale=TRUE)

pdf(file='output/candidate_var_pca.pdf', width=8, height=8)
autoplot(out, loadings=TRUE, loadings.label=TRUE)
dev.off()
