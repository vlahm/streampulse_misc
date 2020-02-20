library(plyr)
library(StreamPULSE)
library(tidyverse)
library(RColorBrewer)
library(plotrix)
rm(list=ls()); cat('/014')

# setup ####

#set mode: 'run' will rerun data generators; 'retrieve' will collect outputs
mode = 'retrieve'
# mode = 'run'

terrcolor = 'sienna4'
aqcolor = 'cadetblue4'
gppcolor = 'forestgreen'
ercolor = 'sienna'

setwd('~/git/streampulse/other_projects/big_stats_paper/')
source('helpers.R')
phil_srcs = list.files('phil_stuff/metab_synthesis0/R/functions/',
    full.names=TRUE)
for(x in phil_srcs) source(x)

# model results from streampulse portal
mods = query_available_results('all')[[1]] %>%
    as_tibble() %>%
    mutate_all(as.character)

#datasets from Phil
fnet_list = readRDS('phil_stuff/FLUXNET_standardized.rds')
fnames = names(fnet_list)
for(i in 1:length(fnet_list)){
    fnet_list[[i]]$sitecode = fnames[i]
}
fnet = Reduce(bind_rows, fnet_list)
fnet = fnet %>%
    select(sitecode, GPP, ER) %>%
    group_by(sitecode) %>%
    summarize_all(mean, na.rm=TRUE) %>%
    rename(GPP_terr=GPP, ER_terr=ER)

metab_d = readRDS('phil_stuff/output/synthesis_standardized.rds')
# metab_d = metab_d[! grepl('nwis', names(metab_d))]
names(metab_d) = phil_to_mike_format(tibble(Site_ID=names(metab_d)), mods,
        arrange=FALSE) %>%
    pull(sitecode)
mn = names(metab_d)
for(i in 1:length(mn)){
    metab_d[[i]]$sitecode = mn[i]
}
magg = Reduce(bind_rows, metab_d)
magg = magg %>%
    select(sitecode, GPP, ER) %>%
    group_by(sitecode) %>%
    summarize_all(mean, na.rm=TRUE) %>%
    rename(GPP_aq=GPP, ER_aq=ER)
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

# junk below? ####

m = arrange(m, desc(ann_GPP_C))

pdf(file=outfile, width=12, height=7)

par(mfrow=c(2, 1), mar=c(0, 3, 1, 1), oma=c(0, 1, 0, 0))

# seq500 = seq(500, 3000, 500)

plot(distset$, ylab='', yaxt='n', yaxs='i', type='n', bty='n', xaxt='n')
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

# bind StreamCat data (superfluous since Phil added his own streamcat pull) ####

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

# write variable key table ####

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

# generate/retrieve watershed boundaries ####

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

# filter and zip boundaries ####

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

# bind IGBP landcover classifications (MCD12Q1v006 LC_Type1) ####

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


# bind other MODIS data? ####
# bind Fluxnet data? ####

# 0 bind site data to model output data; separate sp/powell ####

mods = mods %>%
    left_join(sites, by='sitecode') %>%
    filter(! is.na(Name)) %>%
    distinct() %>%
    arrange(sitecode, year)

spmods = filter(mods, Source == 'StreamPULSE')
powmods = filter(mods, Source == 'USGS (Powell Center)')

# apply Philters and gapPhills; generate sub-datasets ####

fyears = lapply(fnet_list, function(x){
    unique(x$Year)
})
flengths = sapply(fyears, length)
fake_diag = tibble(sitecode=rep(fnames, times=flengths),
    ER_K=rep(0, sum(flengths)),
    Year=unname(unlist(fyears)))
filt_terr = filter_and_impute(fake_diag, models=fnet_list, 'ER_K <= 1', terr=TRUE)
saveRDS(filt_terr, 'output/filtered_dsets/no_filter_terr.rds')

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

gpp_er_biplot = function(mods, outfile, generate_r2s=FALSE){

    cnt = 0
    hold_cnt = FALSE
    errtypes = errmods = errs = list()
    corr_coeffs = c()

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
            corr_coeff = mres$adj.r.squared
            corr_coeffs = append(corr_coeffs, corr_coeff)
            names(corr_coeffs)[length(corr_coeffs)] = m$sitecode
            coeff_lab = paste0('y = ', coeffs[2], 'x ', icept, ' (Adj. R^2: ',
                round(corr_coeff, 2), ')')

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

    if(generate_r2s){
        return(corr_coeffs)
    } else {
        return(list(errtypes=errtypes, errmosd=errmods, errs=errs))
    }
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
    standalone, outfile, filter_label=TRUE, ...){

    #fixed ylims to 5, -5 for standalone; didn't change standalone=F

    #quant_filt example: 'width_calc > 0.25'
    if(standalone){
        filt = readRDS(readfile)
    } else {
        filt = datlist
    }

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

    # nsites_included = sum(names(filt) %in%
    #     sitedata$sitecode[! is.na(sitedata$gpp_C_amp)])
    nsites_included = sum(sapply(filt, nrow) > 0)

    smry = consolidate_list(filt) %>%
        as_tibble() %>%
        group_by(DOY) %>%
        summarize_all(list(median=~median(., na.rm=TRUE),
            quant25=~quantile(., na.rm=TRUE)[2],
            quant75=~quantile(., na.rm=TRUE)[4]))

    if(standalone){
        pdf(file=outfile, width=10, height=10)
        par(mfrow=c(2, 1), oma=c(1, 1, 0, 0))
    }
    par(mar=c(0, 3, 1, 1), lend=2)

    # gpplim = c(0, max(smry$GPP_C_filled_quant75, na.rm=TRUE))
    # erlim = c(min(smry$ER_C_filled_quant25, na.rm=TRUE), 0)
    # gpplim = c(0, max(gpplim[2], abs(erlim[1])))
    # erlim = c(min(abs(gpplim[2]) * -1, erlim[1]), 0)
    gpplim=c(0, 5)
    erlim=c(-5, 0)

    plot(smry$DOY, smry$GPP_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=4, xlab='', ylim=gpplim, xaxs='i', xaxt='n', yaxt='n',
        col='gray30')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$GPP_C_filled_quant25, rev(smry$GPP_C_filled_quant75)),
        border=NA, col=alpha('forestgreen', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, gpplim[2], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, gpplim[2], length.out=5), 1))
    abline(h=0, lty=2, lwd=2, col='gray60')
    medsums = round(colSums(select(smry, contains('median'))), 1)

    if(standalone){
        mtext(expression(paste(bold("gC") ~ bold("m") ^ bold("-2") ~
                bold(" d") ^ bold("-1"))), side=2,
            line=-0.5, outer=TRUE)
        # mtext(expression(paste("GPP (gC"~"m"^"-2"~" d"^"-1"*')')), side=2,
        #     line=2.5)
        # legend('topleft', legend=c('Median', '', '25-75%', '', 'NEP Median'),
        #     col=c('forestgreen', 'sienna4', alpha('forestgreen', alpha=0.6),
        #         alpha('sienna', alpha=0.6), 'black'),
        #     bty='n', lty=1, lwd=c(4, 4, 10, 10, 4))
        if(filter_label){
            legend('topright', title='Filters', bty='n', title.col='gray30',
                lty=1, seg.len=0.2, lwd=2, legend=c(..., var_quant_filt))
        }
    }

    legend('right', title='Cumulative\nMedian Sums', bty='n',
        legend=c(paste('GPP:', medsums[1]), paste('ER:', medsums[2]),
            paste('NEP:', medsums[3])), title.col='gray30')
    legend('left', paste('Sites included:', nsites_included), bty='n')

    par(mar=c(3, 3, 0, 1))

    plot(smry$DOY, smry$ER_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=4, xlab='', ylim=erlim, xaxs='i', xaxt='n', yaxt='n')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$ER_C_filled_quant25, rev(smry$ER_C_filled_quant75)),
        border=NA, col=alpha('sienna', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(1, line=0, tck=-.02, labels=FALSE, at=seq(0, max(smry$DOY), 30))
    axis(1, line=-0.5, tcl=0, col='transparent', at=seq(0, max(smry$DOY), 30))
    lines(smry$DOY, smry$NEP_C_filled_median, col='black', lwd=4, xpd=NA, lend=1)

    if(standalone){
        # mtext(expression(paste("ER (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
        mtext('DOY', side=1, line=2, font=2)
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

# crude_plot(metr$width_calc, metr$Disch_ar1, 'output/width_discharge.pdf')
# crude_plot(metr$MOD_ann_ER, metr$er_C_mean, 'output/MODISer_streamER.pdf')
# crude_plot(log(metr$MOD_ann_GPP), log(metr$gpp_C_mean), 'output/logMODISgpp_logStreamGPP.pdf')
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

pdf_plot = function(outfile, var){

    pdf(outfile, height=4, width=7)

    vv = na.omit(sort(metr[[var]]))
    dens = density(vv)
    vq = quantile(vv, probs=c(0.25, 0.75))
    densdf = tibble(x=dens$x, y=dens$y)
    dens25 = dens75 = densdf
    # dens25[densdf$x > vq[1], ] = NA
    # dens75[densdf$x < vq[2], ] = NA
    dens25 = dens25[densdf$x <= vq[1], ]
    dens75 = dens75[densdf$x >= vq[2], ]
    plot(densdf$x, densdf$y, type='l', xlab='width', ylab='density', bty='l',
        col='gray50', lwd=2)
    polygon(x=c(dens25$x, rev(dens25$x)),
        y=c(dens25$y, rep(0, nrow(dens25))), col='gray50', border='gray50')
    polygon(x=c(dens75$x, rev(dens75$x)),
        y=c(dens75$y, rep(0, nrow(dens75))), col='gray50', border='gray50')

    dev.off()
}

pdf_plot_er_gpp = function(outfile, var='gpp', xlims){

    pdf(outfile, height=4, width=7)

    par(mar=c(1,1,1,1), oma=c(0,0,0,0))

    if(var == 'gpp'){
        vaq = na.omit(sort(metr$gpp_C_mean))
        vterr = na.omit(sort(fnet$GPP_terr))
    } else {
        vaq = na.omit(sort(metr$er_C_mean))
        vterr = na.omit(sort(fnet$ER_terr))
    }
    dterr = density(vterr)
    densterr = tibble(x=dterr$x, y=dterr$y)
    daq = density(vaq)
    densaq = tibble(x=daq$x, y=daq$y)
    ylims = range(c(daq$y, dterr$y), na.rm=TRUE)
    plot(densaq$x, densaq$y, type='n', bty='n', ylim=ylims, xlim=xlims,
        # col='gray50', lwd=2, xlab='', ylab='')
        col='gray50', lwd=2, xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i',
        yaxs='i')
    terrcol = alpha('sienna', alpha=0.5)
    polygon(x=c(densterr$x, rev(densterr$x)),
        y=c(densterr$y, rep(0, nrow(densterr))), col=terrcol, border='sienna')
    aqcol = alpha('skyblue4', alpha=0.5)
    polygon(x=c(densaq$x, rev(densaq$x)),
        y=c(densaq$y, rep(0, nrow(densaq))), col=aqcol, border='skyblue4')

    dev.off()
}

# 0 generate some additional metrics for four corners analyses ####

if(mode == 'run'){
    sp_corr_coeffs = gpp_er_biplot(spmods, 'output/gppXer_sp.pdf', generate_r2s=TRUE)
    powell_corr_coeff = gpp_er_biplot(powmods, 'output/gppXer_powell.pdf',
        generate_r2s=TRUE)
    corr_coeffs = c(sp_corr_coeffs, powell_corr_coeff)
    corr_coeffs = data.frame(sitecode=names(corr_coeffs), corr_coeff=corr_coeffs)
    write.csv(corr_coeffs, 'corr_coeffs.csv', row.names=FALSE)
} else {
    corr_coeffs = read.csv('corr_coeffs.csv', stringsAsFactors=FALSE)
}

interann_gppOverEr_cv = function(x){

    if(length(unique(x$Year)) >= 3){

        x = x %>%
            select(GPP, ER, Year) %>%
            group_by(Year) %>%
            summarize(GPP=mean(GPP, na.rm=TRUE), ER=mean(ER, na.rm=TRUE)) %>%
            ungroup()

        GPP_ER_quotient = x$GPP / x$ER
        GPP_ER_quotient_cv = sd(GPP_ER_quotient, na.rm=TRUE) /
            mean(GPP_ER_quotient, na.rm=TRUE) * 100

        return(GPP_ER_quotient_cv)
    }

}

cv_quotients = sapply(metab_d, interann_gppOverEr_cv)
cv_quotients = cv_quotients[sapply(cv_quotients, function(x) ! is.null(x))]
cv_quotients = data.frame(sitecode=names(cv_quotients),
    quotient_cv=unname(unlist(cv_quotients)), stringsAsFactors=FALSE)

metr_extras = full_join(corr_coeffs, cv_quotients)
metr = left_join(metr, metr_extras)

# plots ####

gpp_er_biplot(spmods, 'output/gppXer_sp.pdf')
gpp_er_biplot(powmods, 'output/gppXer_powell.pdf')
# ts_plot(spmods, 'output/metab_ts_sp_noK.pdf', FALSE)
# ts_plot(powmods, 'output/metab_ts_powell_noK.pdf', FALSE)
ts_plot(spmods, 'output/metab_ts_sp.pdf', TRUE)
ts_plot(powmods, 'output/metab_ts_powell.pdf', TRUE)
dist_plot(metr, 'output/metab_dist.pdf')
er_k_filter_plot(diag, 'output/erXk_corr_filter.pdf')
#for reals lips plots
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='width_calc <= 1',
    standalone=TRUE, outfile='output/final/lips_overall_stream.pdf',
    filter_label=FALSE)
lips_plot(readfile='output/filtered_dsets/no_filter_terr.rds',
    diagnostics=diag, sitedata=sites, quant_filt=NULL,
    standalone=TRUE, outfile='output/final/lips_overall_terr.pdf',
    filter_label=FALSE)

lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='width_calc > 0.75',
    standalone=TRUE, outfile='output/day2/lips_width75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='width_calc < 0.25',
    standalone=TRUE, outfile='output/day2/lips_width25.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Disch_ar1 > 0.75',
    standalone=TRUE, outfile='output/day2/lips_Qar1_75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Disch_ar1 < 0.25',
    standalone=TRUE, outfile='output/day2/lips_Qar1_25.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='MOD_ann_GPP > 0.75',
    standalone=TRUE, outfile='output/day2/lips_modisGPP_75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='MOD_ann_GPP < 0.25',
    standalone=TRUE, outfile='output/day2/lips_modisGPP_25.pdf')

lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Stream_PAR_sum > 0.75',
    standalone=TRUE, outfile='output/day2/lips_streamPAR_75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Stream_PAR_sum < 0.25',
    standalone=TRUE, outfile='output/day2/lips_streamPAR_25.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Disch_cv > 0.75',
    standalone=TRUE, outfile='output/day2/lips_Qcv_75.pdf')
lips_plot(readfile='output/filtered_dsets/no_filter.rds',
    diagnostics=diag, sitedata=sites, quant_filt='Disch_cv < 0.25',
    standalone=TRUE, outfile='output/day2/lips_Qcv_25.pdf')
# end for reals
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
pdf_plot('output/day2/probdens_width.pdf', 'width_calc')
pdf_plot('output/day2/probdens_Qar1.pdf', 'Disch_ar1')
pdf_plot('output/day2/probdens_modisGPP.pdf', 'MOD_ann_GPP')
pdf_plot('output/day2/probdens_streamPAR.pdf', 'Stream_PAR_sum')
pdf_plot('output/day2/probdens_Qcv.pdf', 'Disch_cv')
pdf_plot_er_gpp('output/day2/GPP_PDFs.pdf', 'gpp', gpplim)
pdf_plot_er_gpp('output/day2/ER_PDFs.pdf', 'er', erlim)

# GPP-ER biplot and dist plots (Figure 1) ####

pdf(file='output/final/gpp_er_biplot.pdf', width=8, height=8)

log_gpp_terr = log(fnet$GPP_terr)
log_er_terr = -1 * log(fnet$ER_terr * -1)
log_gpp_aq = log(metr$gpp_C_mean)
log_er_aq = -1 * log(metr$er_C_mean * -1)
plot(log_gpp_terr, log_er_terr, col=alpha(terrcolor, alpha=0.5),
    xlab='GPP (gC)',# xaxs='i', yaxs='i',
    ylab='ER (gC)', cex=2, cex.lab=1.2, cex.axis=1.2,
    pch=20, yaxt='n', xaxt='n')
points(log_gpp_aq, log_er_aq, col=alpha(aqcolor, alpha=0.5),
    cex=2, pch=20)
legend('topright', legend=c('FLUXNET', 'StreamPULSE'), pch=20, bty='n', pt.cex=2,
    col=c(alpha(terrcolor, alpha=0.5), alpha(aqcolor, alpha=0.5)))

all_gpp = c(fnet$GPP_terr, metr$gpp_C_mean)
all_gpp[all_gpp <= 0] = NA
gpprng = range(all_gpp, na.rm=TRUE)
all_er = c(fnet$ER_terr, metr$er_C_mean)
all_er[all_er >= 0] = NA
errng = range(all_er, na.rm=TRUE)

gpptck = c(0.04, 0.25, 1, 2, 4, 8, 11.8)
gpptck_log = log(gpptck)
axis(1, at=gpptck_log, labels=gpptck)
# ertck = seq(errng[1], errng[2], length.out=20)
# ertck_log = log(ertck * -1) * -1
ertck = rev(c(14.3, 8, 4, 2, 1, 0.25, 0.06))
ertck_log = log(ertck) * -1
axis(2, at=ertck_log, labels=ertck * -1)

abline(a=0, b=-1, lty=2)
dev.off()

#---

pdf(file='output/final/distplots.pdf', width=8, height=8)
par(mfrow=c(2, 1), mar=c(1,1,1,1), oma=c(0, 0, 0, 0))

dens = density(na.omit(log_gpp_terr))
gpp_dens_terr = tibble(x=dens$x, y=dens$y)
dens = density(na.omit(log_gpp_aq))
gpp_dens_aq = tibble(x=dens$x, y=dens$y)

plot(gpp_dens_terr$x, gpp_dens_terr$y, type='n', ann=FALSE, xaxt='n', yaxt='n',
    bty='n')
mtext('GPP', 1, line=0)
polygon(x=c(gpp_dens_terr$x, rev(gpp_dens_terr$x)),
    y=c(gpp_dens_terr$y, rep(0, nrow(gpp_dens_terr))),
    col=alpha(terrcolor, alpha=0.7),
    border=alpha(terrcolor, alpha=0.7))
polygon(x=c(gpp_dens_aq$x, rev(gpp_dens_aq$x)),
    y=c(gpp_dens_aq$y, rep(0, nrow(gpp_dens_aq))),
    col=alpha(aqcolor, alpha=0.7),
    border=alpha(aqcolor, alpha=0.7))

dens = density(na.omit(log_er_terr))
er_dens_terr = tibble(x=dens$x, y=dens$y)
dens = density(na.omit(log_er_aq))
er_dens_aq = tibble(x=dens$x, y=dens$y)

plot(er_dens_terr$x, er_dens_terr$y, type='n', ann=FALSE, xaxt='n', yaxt='n',
    bty='n')
mtext('ER', 1, line=0)
polygon(x=c(er_dens_terr$x, rev(er_dens_terr$x)),
    y=c(rev(er_dens_terr$y), rep(0, nrow(er_dens_terr))),
    col=alpha(terrcolor, alpha=0.7),
    border=alpha(terrcolor, alpha=0.7))
polygon(x=c(er_dens_aq$x, rev(er_dens_aq$x)),
    y=c(rev(er_dens_aq$y), rep(0, nrow(er_dens_aq))),
    col=alpha(aqcolor, alpha=0.7),
    border=alpha(aqcolor, alpha=0.7))

dev.off()

# plot(sort(na.omit(log_gpp_terr), decreasing=FALSE),
#     1:length(na.omit(log_gpp_terr)), type='l',
#     # ylab='', xlab='GPP (gC)',
#     bty='n', col='sienna', lwd=2, ann=FALSE)
#     xaxt='n', yaxt='n')#, xlim=c(0, 366))
# GPP_terr = sort(log_gpp_terr, decreasing=TRUE)
# terrcol = alpha('sienna', alpha=0.5)
# polygon(x=log(c(rep(0.0001, 166), 166:1)),
#     y=c(GPP_terr, rev(GPP_terr)), lwd=2,
#     col=terrcol, border=terrcol)
# par(new=TRUE)
# plot(sort(magg$GPP_aq, decreasing=TRUE), type='n', xlab='',
#     ylab='GPP (gC)', bty='n', col='sienna', lwd=2, ann=FALSE,
#     xaxt='n', yaxt='n', xlim=c(0, 366))
# GPP_aq = sort(magg$GPP_aq, decreasing=TRUE)
# aqcol = alpha('skyblue3', alpha=0.5)
# polygon(x=c(rep(0, 407), 407:1), y=c(GPP_aq, rev(GPP_aq)), lwd=2,
#     col=aqcol, border=aqcol)


# plot(sort(fnet$ER_terr * -1, decreasing=TRUE), type='l', xlab='',
#     ylab='ER (gC)', bty='l', col='forestgreen', lwd=2)
# lines(sort(magg$ER_aq * -1, decreasing=TRUE), lwd=2, col='sienna')
# mtext('DOY', 1, outer=TRUE)
# dev.off()

# depth acquisition ####
plot(metr$gpp_C_mean, metr$er_C_mean)
metr = metr %>%
    arrange(desc(gpp_C_mean)) %>%
    select(gpp_C_mean, everything())
head(metr)
metr$sitecode[1:5]
StreamPULSE::query_available_results('TX', site='nwis-08437710')
x = StreamPULSE::request_results('OR_nwis-13173600', 2011)
dd = data.frame(sitecode=rep(NA, 3), depth_m=rep(NA, 3))
dd[1,] = c('OR_nwis-13173600', mean(x$model_results$data$depth, na.rm=TRUE))
x = StreamPULSE::request_results('TX_nwis-08446500', 2011)
dd[2,] = c('TX_nwis-08446500', mean(x$model_results$data$depth, na.rm=TRUE))
x = StreamPULSE::request_results('TX_nwis-08437710', 2011)
dd[3,] = c('TX_nwis-08437710', mean(x$model_results$data$depth, na.rm=TRUE))
write.csv(dd, 'output/day2/extreme_gpp_site_depths.csv', row.names=FALSE)


plot(log(metr$MOD_ann_GPP), log(metr$gpp_C_mean))

# peak day and week of productivity ####

sagg = Reduce(bind_rows, metab_d)

doymaxes_s = sagg %>%
    select(sitecode, DOY, GPP) %>%
    group_by(sitecode) %>%
    filter(GPP == max(GPP, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(DOY)

doypeaks_s = doymaxes_s %>%
    group_by(DOY) %>%
    summarize(meanGPP=mean(GPP, na.rm=TRUE), nDOY=length(DOY),
        pctDOY=length(DOY) / nrow(doymaxes_s) * 100) %>%
    filter(nDOY == max(nDOY, na.rm=TRUE)) %>%
    mutate(type='stream') %>%
    select(type, everything())

woymaxes_s = sagg %>%
    mutate(WOY=floor(sagg$DOY / 7 - 0.1) + 1) %>%
    select(sitecode, WOY, GPP) %>%
    group_by(sitecode) %>%
    filter(GPP == max(GPP, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(WOY)

woypeaks_s = woymaxes_s %>%
    group_by(WOY) %>%
    summarize(meanGPP=mean(GPP, na.rm=TRUE), nWOY=length(WOY),
        pctWOY=length(WOY) / nrow(woymaxes_s) * 100) %>%
    filter(nWOY == max(nWOY, na.rm=TRUE)) %>%
    mutate(type='stream') %>%
    select(type, everything())

tagg = Reduce(bind_rows, fnet_list)

doymaxes_t = tagg %>%
    select(sitecode, DOY, GPP) %>%
    group_by(sitecode) %>%
    filter(GPP == max(GPP, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(DOY)

doypeaks_t = doymaxes_t %>%
    group_by(DOY) %>%
    summarize(meanGPP=mean(GPP, na.rm=TRUE), nDOY=length(DOY),
        pctDOY=length(DOY) / nrow(doymaxes_t) * 100) %>%
    filter(nDOY == max(nDOY, na.rm=TRUE)) %>%
    mutate(type='terr') %>%
    select(type, everything())

woymaxes_t = tagg %>%
    mutate(WOY=floor(tagg$DOY / 7 - 0.1) + 1) %>%
    select(sitecode, WOY, GPP) %>%
    group_by(sitecode) %>%
    filter(GPP == max(GPP, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(WOY)

woypeaks_t = woymaxes_t %>%
    group_by(WOY) %>%
    summarize(meanGPP=mean(GPP, na.rm=TRUE), nWOY=length(WOY),
        pctWOY=length(WOY) / nrow(woymaxes_t) * 100) %>%
    filter(nWOY == max(nWOY, na.rm=TRUE)) %>%
    mutate(type='terr') %>%
    select(type, everything())

most_productive_DOYs = bind_rows(doypeaks_s, doypeaks_t)
most_productive_WOYs = bind_rows(woypeaks_s, woypeaks_t)
write.csv(most_productive_DOYs, 'output/final/most_productive_DOYs.csv',
    row.names=FALSE)
write.csv(most_productive_WOYs, 'output/final/most_productive_WOYs.csv',
    row.names=FALSE)

#plot---

pdf('output/final/peak_productivity_dists.pdf', width=8, height=8)

tdens_doy = density(doymaxes_t$DOY)
tdens_woy = density(woymaxes_t$WOY)
tdf = tibble(x_woy=tdens_woy$x, y_woy=tdens_woy$y,
    x_doy=tdens_doy$x, y_doy=tdens_doy$y)
sdens_doy = density(doymaxes_s$DOY)
sdens_woy = density(woymaxes_s$WOY)
sdf = tibble(x_woy=sdens_woy$x, y_woy=sdens_woy$y,
    x_doy=sdens_doy$x, y_doy=sdens_doy$y)

# woypeaks = pracma::findpeaks(woymaxes_s$GPP, npeaks=3)
# abline(v=tdf$x_doy[woypeaks])

plot(tdf$x_doy, tdf$y_doy, xlab='DOY', main='', type='l',
    col=terrcolor, lwd=2, xlim=range(c(tdf$x_doy, sdf$x_doy), na.rm=TRUE),
    ylab='Kernel density')
lines(sdf$x_doy, sdf$y_doy, col=aqcolor, lwd=2)
legend('topright', legend=c('terr', 'aq'), col=c(terrcolor, aqcolor),
    lty=1, lwd=2, bty='n')

# plot(tdf$x_woy, tdf$y_woy, xlab='WOY', main='', type='l',
#     col=terrcolor, lwd=2, xlim=range(c(tdf$x_woy, sdf$x_woy), na.rm=TRUE),
#     ylab='Kernel density')
# lines(sdf$x_woy, sdf$y_woy, col=aqcolor, lwd=2)

dev.off()

write.csv(doymaxes_t, file='output/final/peak_productivity_terr.csv',
    row.names=FALSE)
write.csv(doymaxes_s, file='output/final/peak_productivity_aq.csv',
    row.names=FALSE)
# write.csv(woymaxes, file='output/final/peak_productivity_woy.csv', row.names=FALSE)

# slope distribution for terr and aq siteyears ####

get_metab_slope = function(x){

    x = as_tibble(x) %>%
        mutate(ER=ER * -1) %>%
        select(Date, GPP, ER) %>%
        group_by(substr(Date, 6, 10)) %>%
        summarize_all(mean, na.rm=TRUE) %>%
        ungroup()

    slope = summary(lm(x$ER ~ x$GPP))$coefficients[2, 1]

    return(slope)
}

fslopes = sapply(fnet_list, get_metab_slope)
sslopes = sapply(metab_d, get_metab_slope)

GPPxER_slopeDists = data.frame(terr=round(quantile(fslopes), 2),
    aq=round(quantile(sslopes), 2))
write.csv(GPPxER_slopeDists, 'output/final/GPPxER_slopeDists.csv')

fdens = density(fslopes)
fdf = tibble(x=fdens$x, y=fdens$y)
sdens = density(sslopes)
sdf = tibble(x=sdens$x, y=sdens$y)

pdf('output/final/gpp_er_slope_dists.pdf', width=8, height=8)

plot(fdf$x, fdf$y, xlab='-ER vs. GPP slope coeff', main='', type='l',
    col=terrcolor, lwd=2, xlim=range(c(fdf$x, sdf$x), na.rm=TRUE),
    ylab='Kernel density')
lines(sdf$x, sdf$y, col=aqcolor, lwd=2)
legend('topright', legend=c('terr', 'aq'), col=c(terrcolor, aqcolor),
    lty=1, lwd=2, bty='n')

dev.off()

# intra- and inter-annual cvs of productivity ####

fprodcvs_intra = sapply(fnet_list, function(x){
    cv = sd(x$GPP, na.rm=TRUE) / mean(x$GPP, na.rm=TRUE) * 100
    return(cv)
})

sprodcvs_intra = sapply(metab_d, function(x){
    cv = sd(x$GPP, na.rm=TRUE) / mean(x$GPP, na.rm=TRUE) * 100
    return(cv)
})

inter_cv = function(x){

    if(length(unique(x$Year)) >= 3){

        x = x %>%
            select(GPP, Year) %>%
            group_by(Year) %>%
            summarize(GPP=mean(GPP, na.rm=TRUE)) %>%
            ungroup()

        cv = sd(x$GPP, na.rm=TRUE) / mean(x$GPP, na.rm=TRUE) * 100

        return(cv)
    }

}

fprodcvs_inter = sapply(fnet_list, inter_cv)
fprodcvs_inter = unlist(Filter(function(x) ! is.null(x), fprodcvs_inter))
sprodcvs_inter = sapply(metab_d, inter_cv)
sprodcvs_inter = unlist(Filter(function(x) ! is.null(x), sprodcvs_inter))

pdf('output/final/gpp_CVs.pdf', width=5, height=5)

par(mfrow=c(2, 1), mar=c(4, 4, 1, 2))

plot(density(fprodcvs_intra[fprodcvs_intra > 0]), xlim=c(0, 1000),
    ylim=c(0, 0.01), col='darkgreen', lty=2, lwd=2,
    xlab='Intra-annual CV of GPP', main='', xaxs='i', yaxs='i')
lines(density(sprodcvs_intra[sprodcvs_intra > 0]), col='cadetblue4', lwd=2)
legend('topright', legend=c('stream', 'terr'), lty=c(2, 1), lwd=2,
    col=c('darkgreen', 'cadetblue4'), bty='n')

plot(density(na.omit(fprodcvs_inter[fprodcvs_inter > 0])),# xlim=c(0, 1000),
    col='darkgreen', lty=2, lwd=2,
    xlab='Interannual CV of GPP ( > 2 siteyrs)', main='', xaxs='i', yaxs='i')
lines(density(sprodcvs_inter[sprodcvs_inter > 0]), col='cadetblue4', lwd=2)
legend('topright', legend=c('stream', 'terr'), lty=c(2, 1), lwd=2,
    col=c('darkgreen', 'cadetblue4'), bty='n')

dev.off()

# NWIS (USGS) sites included in analysis ####
nn = names(metab_d)
nn = nn[grepl('nwis', nn)]
nn = substr(nn, 9, nchar(nn))
write.csv(data.frame(usgs_gage_id=nn), row.names=FALSE,
    file='output/final/usgs_gage_ids.csv')

# 0 identify the "four corners" groups ####

par_z = scale(metr$Stream_PAR_sum)
q_z = scale(metr$Disch_ar1)
par_z_inv = par_z * -1
q_z_inv = q_z * -1
zframe = data.frame(sitecode=metr$sitecode, par=par_z, qar1=q_z,
        par_inv=par_z_inv, qar1_inv=q_z_inv) %>%
    mutate(primary_axis=par_z + q_z, quad_2=par_inv + qar1,
        quad_4=qar1_inv + par)

# quartile_n = sum(complete.cases(zframe[, c('par', 'qar1')])) / 4
# bright_stable = head(order(zframe$primary_axis),
bright_stable_thresh = quantile(zframe$primary_axis, probs=0.75 + 0.125, na.rm=TRUE)
bright_stable = which(zframe$primary_axis >= bright_stable_thresh)
er_thresh = quantile(zframe$primary_axis, probs=0.25 - 0.125, na.rm=TRUE)
dark_stormy = which(zframe$primary_axis <= er_thresh)
bright_erratic_thresh = quantile(zframe$quad_4, probs=0.75 + 0.125, na.rm=TRUE)
bright_erratic = which(zframe$quad_4 >= bright_erratic_thresh)
dark_dull_thresh = quantile(zframe$quad_2, probs=0.75 + 0.125, na.rm=TRUE)
dark_dull = which(zframe$quad_2 >= dark_dull_thresh)

bs_sites = as.character(zframe[bright_stable, 1])
ds_sites = as.character(zframe[dark_stormy, 1])
be_sites = as.character(zframe[bright_erratic, 1])
dd_sites = as.character(zframe[dark_dull, 1])

# four corners plots ####

bs_metr = metr[metr$sitecode %in% bs_sites, ]
ds_metr = metr[metr$sitecode %in% ds_sites, ]
be_metr = metr[metr$sitecode %in% be_sites, ]
dd_metr = metr[metr$sitecode %in% dd_sites, ]

# gpp_quants = data.frame(bs=quantile(bs_metr$gpp_C_ar1),
#     ds=quantile(ds_metr$gpp_C_ar1),
#     be=quantile(be_metr$gpp_C_ar1),
#     dd=quantile(dd_metr$gpp_C_ar1))

# plot(1:4, unlist(gpp_quants['50%',]), ylim=ylims, pch=15,
#     ylab='GPP AR-1')
# segments(1:4, unlist(gpp_quants['0%',]), 1:4, unlist(gpp_quants['25%',]))

boxout_gppAR1 = boxplot(list(bs_metr$gpp_C_ar1, dd_metr$gpp_C_ar1,
    ds_metr$gpp_C_ar1, be_metr$gpp_C_ar1))
boxout_erAR1 = boxplot(list(bs_metr$er_C_ar1, dd_metr$er_C_ar1,
    ds_metr$er_C_ar1, be_metr$er_C_ar1))
boxout_erXgppR2 = boxplot(list(bs_metr$corr_coeff, dd_metr$corr_coeff,
    ds_metr$corr_coeff, be_metr$corr_coeff))
bs_metr_sub = bs_metr$quotient_cv[bs_metr$quotient_cv < 4000]
boxout_quotientCV = boxplot(list(bs_metr_sub, dd_metr$quotient_cv,
    ds_metr$quotient_cv, be_metr$quotient_cv))

pdf('output/final/4corners.pdf', width=10, height=10)

par(mfrow=c(2, 2))

#gpp ar1
ylims = range(c(boxout_gppAR1$stats, boxout_gppAR1$out))
plot(1:4, boxout_gppAR1$stats[3,], ylim=ylims, pch=15, ylab='GPP AR-1', col=gppcolor,
    xlab='', xaxt='n', cex=0.6)
segments(1:4, boxout_gppAR1$stats[1,], 1:4, boxout_gppAR1$stats[2,], lwd=3,
    col=gppcolor, lend=1)
segments(1:4, boxout_gppAR1$stats[4,], 1:4, boxout_gppAR1$stats[5,], lwd=3,
    col=gppcolor, lend=1)
points(boxout_gppAR1$group, boxout_gppAR1$out, pch=4, cex=1, col=gppcolor)
axis(1, at=1:4, labels=c('I', 'II', 'III', 'IV'))

#er ar1
ylims = range(c(boxout_erAR1$stats, boxout_erAR1$out))
plot(1:4, boxout_erAR1$stats[3,], ylim=ylims, pch=15, ylab='ER AR-1', col=ercolor,
    xlab='', xaxt='n', cex=0.6)
segments(1:4, boxout_erAR1$stats[1,], 1:4, boxout_erAR1$stats[2,], lwd=3,
    col=ercolor, lend=1)
segments(1:4, boxout_erAR1$stats[4,], 1:4, boxout_erAR1$stats[5,], lwd=3,
    col=ercolor, lend=1)
points(boxout_erAR1$group, boxout_erAR1$out, pch=4, cex=1, col=ercolor)
axis(1, at=1:4, labels=c('I', 'II', 'III', 'IV'))

#gpp x er corr coeff
ylims = range(c(boxout_erXgppR2$stats, boxout_erXgppR2$out))
plot(1:4, boxout_erXgppR2$stats[3,], ylim=ylims, pch=15, ylab='ER vs. GPP R^2',
    col='black', xlab='', xaxt='n', cex=0.6)
segments(1:4, boxout_erXgppR2$stats[1,], 1:4, boxout_erXgppR2$stats[2,], lwd=3,
    col='black', lend=1)
segments(1:4, boxout_erXgppR2$stats[4,], 1:4, boxout_erXgppR2$stats[5,], lwd=3,
    col='black', lend=1)
points(boxout_erXgppR2$group, boxout_erXgppR2$out, pch=4, cex=1, col='black')
axis(1, at=1:4, labels=c('I', 'II', 'III', 'IV'))

#GPP/ER CV for sites with >=3 yrs results
ylims = range(c(boxout_quotientCV$stats, boxout_quotientCV$out))
plot(1:4, boxout_quotientCV$stats[3,], ylim=ylims, pch=15, ylab='GPP/ER CV',
    col='gray50', xlab='', xaxt='n', cex=0.6)
segments(1:4, boxout_quotientCV$stats[1,], 1:4, boxout_quotientCV$stats[2,], lwd=3,
    col='gray50', lend=1)
segments(1:4, boxout_quotientCV$stats[4,], 1:4, boxout_quotientCV$stats[5,], lwd=3,
    col='gray50', lend=1)
points(boxout_quotientCV$group, boxout_quotientCV$out, pch=4, cex=1, col='gray50')
axis(1, at=1:4, labels=c('I', 'II', 'III', 'IV'))

dev.off()

# one more four corners plot for Alice ####

bs_inter = sprodcvs_inter[names(sprodcvs_inter) %in% bs_sites]
bs_intra = sprodcvs_intra[names(sprodcvs_intra) %in% bs_sites]
dd_inter = sprodcvs_inter[names(sprodcvs_inter) %in% dd_sites]
dd_intra = sprodcvs_intra[names(sprodcvs_intra) %in% dd_sites]
ds_inter = sprodcvs_inter[names(sprodcvs_inter) %in% ds_sites]
ds_intra = sprodcvs_intra[names(sprodcvs_intra) %in% ds_sites]
be_inter = sprodcvs_inter[names(sprodcvs_inter) %in% be_sites]
be_intra = sprodcvs_intra[names(sprodcvs_intra) %in% be_sites]

pdf('output/final/4corners_inter_intra_CV_dists.pdf', height=8, width=8)

boxplot(list(bs_inter, bs_intra, dd_inter, dd_intra, ds_inter, ds_intra,
    be_inter, be_intra), ylab='CV', col=c('orange', 'cadetblue3'), xaxt='n',
    boxwex=.3)
axis(1, at=seq(1.5, 7.5, 2), labels=c('I', 'II', 'III', 'IV'))
legend('topright', legend=c('interannual', 'intra-annual'),
    fill=c('orange', 'cadetblue3'), bty='n')

dev.off()

# bright and stable site data for Audrey ####

bs_fullnames = as.character(zframe[bright_stable, 1])
xx = metab_d[bs_fullnames[grep('nwis', bs_fullnames)]]
yr_list = lapply(xx, function(x) unique(x$Year))

bs_sites = na.omit(unname(sapply(bs_fullnames,
    function(x){
        strsplit(x, '-')[[1]][2]
    })))

names(yr_list) = bs_sites

saveRDS(yr_list, 'siteyears.rds')
# write.csv(data.frame(nwis_code=bs_sites), 'bright_stable_sites.csv',
#     row.names=FALSE)

bs_siteframe = tibble(sitecode=rep(names(yr_list),
    times=sapply(yr_list, length)),
    year=unname(unlist(yr_list)))

write.csv(bs_siteframe, 'nwis_siteyears_sub.csv', row.names=FALSE)

#full sitelist
fullnames = as.character(zframe[c(dark_stormy, bright_stable,
    bright_erratic, dark_dull), 1])
xx = metab_d[fullnames[grep('nwis', fullnames)]]
yr_list = lapply(xx, function(x) unique(x$Year))

sites = na.omit(unname(sapply(fullnames,
    function(x){
        strsplit(x, '-')[[1]][2]
    })))

names(yr_list) = sites

full_siteframe = tibble(sitecode=rep(names(yr_list),
    times=sapply(yr_list, length)),
    year=unname(unlist(yr_list)))

write.csv(full_siteframe, 'nwis_siteyears_full.csv', row.names=FALSE)

# PAR vs Qar1 by gpp ####

pdf('output/final/light_vs_flow_by_gpp.pdf', height=8, width=8)
par(mar=c(5, 5, 4, 6))
gpprng = range(metr$gpp_C_mean, na.rm=TRUE)
rescaled = ((metr$gpp_C_mean - gpprng[1]) / (gpprng[2] - gpprng[1])) * (4 - 1) + 1
plot(metr$Stream_PAR_sum, metr$Disch_ar1, pch=20,
    xlab='Light Availability (Mean Annual Surface PAR)',
    ylab='Predictability of Flow (Discharge AR-1 Coeff.)',
    col=alpha('darkgreen', alpha=0.5), bty='o',
    xpd=NA, main='', cex=rescaled, font.lab=2)
legend('right', legend=c(gpprng[1], '', '', gpprng[2]), pch=20, bty='n',
    pt.cex=c(1, 2, 3, 4), col=alpha('darkgreen', alpha=0.5), xpd=NA,
    inset=c(-0.15, 0), title=expression(paste(bold('Mean\nAnnual\nGPP'))))
dev.off()


# PAR vs Qar1 by gpp with corners ####

pdf('output/final/light_vs_flow_by_gpp_corners.pdf', height=8, width=8)
# lmod = lm(metr$Disch_ar1 ~ metr$Stream_PAR_sum)
# zz = scales::rescale(metr$gpp_C_mean, to=c(1, 3))
par(mar=c(5, 5, 4, 6))
gpprng = range(metr$gpp_C_mean, na.rm=TRUE)
rescaled = ((metr$gpp_C_mean - gpprng[1]) / (gpprng[2] - gpprng[1])) * (4 - 1) + 1
# plot(metr$Stream_PAR_sum, metr$Disch_ar1, pch=20, xlab='Light Availability',
plot(metr$Stream_PAR_sum, metr$Disch_ar1, pch=20,
    xlab='Light Availability (Mean Annual Surface PAR)',
    ylab='Predictability of Flow (Discharge AR-1 Coeff.)',
    col=alpha('gray60', alpha=0.5), bty='o',
    # col=alpha('darkgreen', alpha=0.5), bty='u',
    xpd=NA, main='', cex=rescaled, font.lab=2)
points(metr[bright_stable, 'Stream_PAR_sum', drop=TRUE], xpd=NA,
    metr[bright_stable, 'Disch_ar1', drop=TRUE], cex=rescaled[bright_stable],
    # pch=20, col=alpha('black', alpha=0.5))
    pch=20, col='black')
points(metr[bright_erratic, 'Stream_PAR_sum', drop=TRUE], xpd=NA,
    metr[bright_erratic, 'Disch_ar1', drop=TRUE], cex=rescaled[bright_erratic],
    pch=20, col='blue')
points(metr[er, 'Stream_PAR_sum', drop=TRUE], xpd=NA,
    metr[er, 'Disch_ar1', drop=TRUE], cex=rescaled[er],
    pch=20, col='orange')
points(metr[dark_dull, 'Stream_PAR_sum', drop=TRUE], xpd=NA,
    metr[dark_dull, 'Disch_ar1', drop=TRUE], cex=rescaled[dark_dull],
    pch=20, col='red')
# main=paste(paste('R^2:', round(summary(lmod)$adj.r.squared, 2), '\n'),
# expression(paste('size'~alpha~'GPP'))),
# main=paste0('R^2: ', round(summary(lmod)$adj.r.squared, 2), '; size=GPP'),
# cex.main=0.8)
# abline(lmod, lty=2, col='blue')
legend('right', legend=c(gpprng[1], '', '', gpprng[2]), pch=20, bty='n',
    pt.cex=c(1, 2, 3, 4), col=alpha('gray60', alpha=0.5), xpd=NA,
    # pt.cex=c(1, 2, 3, 4), col=alpha('darkgreen', alpha=0.5), xpd=NA,
    inset=c(-0.15, 0), title=expression(paste(bold('Mean\nAnnual\nGPP'))))
dev.off()


# PAR vs Qar1 by er ####

pdf('output/final/light_vs_flow_by_er.pdf', height=8, width=8)
par(mar=c(5, 5, 4, 6))
# gpprng = range(metr$gpp_C_mean, na.rm=TRUE)
# rescaled = ((-1 * metr$er_C_mean - gpprng[1]) / (gpprng[2] - gpprng[1])) * (4 - 1) + 1
errng = range(metr$er_C_mean * -1, na.rm=TRUE)
rescaled = ((-1 * metr$er_C_mean - errng[1]) / (errng[2] - errng[1])) * (4 - 1) + 1
plot(metr$Stream_PAR_sum, metr$Disch_ar1, pch=20,
    xlab='Light Availability (Mean Annual Surface PAR)',
    ylab='Predictability of Flow (Discharge AR-1 Coeff.)',
    col=alpha('sienna4', alpha=0.5), bty='o',
    xpd=NA, main='', cex=rescaled, font.lab=2)
legend('right', legend=c(errng[1], '', '', errng[2]), pch=20, bty='n',
# legend('right', legend=c(gpprng[1], '', '', gpprng[2]), pch=20, bty='n',
    pt.cex=c(1, 2, 3, 4), col=alpha('sienna4', alpha=0.5), xpd=NA,
    inset=c(-0.15, 0), title=expression(paste(bold('Mean\nAnnual\nER'))))
dev.off()


# P:R vs stream order ####
metr = left_join(metr, select(sites, sitecode, STREAMORDE))

pdf('output/final/PR_ratio_vs_strahler_order.pdf', width=8, height=8)
plot(metr$STREAMORDE, metr$gpp_C_mean / metr$er_C_mean, ylab='P:R',
    xlab='Strahler Order', pch=20, cex=1.5, col=alpha('cadetblue4', alpha=0.1))
dev.off()
