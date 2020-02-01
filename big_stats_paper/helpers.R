comid_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    COMID = nhdplusTools::discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}

vpu_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    VPU = nhdR::find_vpu(ptc)
    return(VPU)
}

#this calculates how far along a reach any given point falls. That way when we pull in
#watershed summary data for a reach, we can adjust it according to how much
#of the total upstream area actually contributes to the point in question.
# A value of 0 means upstream end; 1 means downstream end.
calc_reach_prop = function(VPU, COMID, lat, long, CRS, quiet=FALSE){

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    fl = nhdR::nhd_plus_load(vpu=VPU, component='NHDSnapshot',
        dsn='NHDFlowline', approve_all_dl=TRUE)
    fl_etc = nhdR::nhd_plus_load(vpu=VPU, component='NHDPlusAttributes',
        dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

    colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
    colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
    fl = fl[fl$COMID == COMID,]
    fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
        by=c('COMID'='ComID'))

    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    ptct = sf::st_transform(ptc, crs=4269) #CRS for NAD 83
    x = suppressWarnings(nhdplusTools::get_flowline_index(fl, points=ptct))
    out = 1 - x$REACH_meas / 100 #0=upstream end; 1=downstream end

    return(out)
}

#this acquires nhdplusv2 data for a single site by COMID.
#it's just a thin wrapper around nhdR::nhd_plus_load
nhdplusv2_from_comid = function(VPU, COMID, component, DSN, quiet=FALSE) {

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    data = nhdR::nhd_plus_load(vpu=VPU, component=component, dsn=DSN,
        approve_all_dl=TRUE)

    colnames(data)[colnames(data) == 'ComID'] = 'COMID'
    colnames(data)[colnames(data) == 'ReachCode'] = 'REACHCODE'
    data = data[data$COMID == COMID,]

    return(data)
}

#this calls nhdplusv2_from_comid repeatedly to get data for all your sites.
#the dataframe must include COMID and VPU columns
nhdplusv2_bulk = function(site_df, nhdplusv2_sets, quiet=FALSE){

    nhdplus_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(setlist)){
            print(paste(j, nhdplusv2_sets[[i]]))

            if(i == 1 || initerr){
                row_base = try(nhdplusv2_from_comid(site_df$VPU[j],
                    site_df$COMID[j], names(setlist[i]), setlist[[i]],
                    quiet=quiet))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(nhdplusv2_from_comid(site_df$VPU[j],
                    site_df$COMID[j], names(setlist[i]), setlist[[i]],
                    quiet=quiet))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        nhdplus_data = rbind.fill(nhdplus_data, row_base)
    }

    return(nhdplus_data)
}

query_streamcat_datasets = function(keyword=NULL){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')

    url_list = getURL(ftpdir, dirlistonly=TRUE)
    url_list = strsplit(url_list, split='\n')[[1]]

    if(! is.null(keyword)){
        url_list = url_list[grep(keyword, url_list, ignore.case=TRUE)]
    }

    return(url_list)
}

#this function acquires streamcat data for a single site by NHDPlusV2 COMID.
streamcat_from_comid = function(USstate, COMID, dataset){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')
    zip_name = paste0(dataset, '_', USstate, '.zip')

    csv_name = gsub('.zip', '.csv', zip_name)
    temp = tempfile()
    download.file(paste0(ftpdir, zip_name), temp)
    data = read.csv(unz(temp, csv_name), stringsAsFactors=FALSE)
    data = data[data$COMID == COMID,]

    return(data)
}

#this calls streamcat_from_comid repeatedly to get data for all your sites
#the dataframe must include COMID and region columns, where "region" refers to
#each state's 2-letter abbreviation.
streamcat_bulk = function(site_df, streamcat_sets){

    streamcat_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(streamcat_sets)){
            print(paste(j, streamcat_sets[i]))

            if(i == 1 || initerr){
                row_base = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        streamcat_data = rbind.fill(streamcat_data, row_base)
    }

    return(streamcat_data)
}