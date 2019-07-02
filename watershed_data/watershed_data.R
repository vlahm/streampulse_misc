library(sf)#
library(mapview)###
# library(mapedit)
# library(rayshader)
library(elevatr)#
library(raster)#
# devtools::install_github("giswqs/whiteboxR")
library(whitebox) #
library(stars) #
# library(rgl)
# library(rgdal)
library(RMariaDB) #
library(stringr) #
library(dplyr) #
# knitr::knit_hooks$set(webgl = hook_webgl)

# setwd('~/git/streampulse/other_projects/watershed_data/')
setwd('~/Desktop/untracked/sp_watershed_delin') #this can be any empty directory

# #read in mysql pw
# conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
# # conf = readLines('/home/aaron/sp/config.py')
# ind = which(lapply(conf, function(x) grepl('MYSQL_PW', x)) == TRUE)
# pw = str_match(conf[ind], '.*\\"(.*)\\"')[2]
#
# #read in site and results tables from mysql
# con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)
# sites = as_tibble(dbReadTable(con, "site")) %>%
#     arrange(latitude, longitude) %>%
#     mutate(regionsite=paste(region, site, sep='_'))

# saveRDS(sites, '~/Desktop/site_data.rds')
sites = readRDS('site_data.rds')

# projection_table = list('equatorial'='+proj=laea +lon_0=-125.859375',
#     'temperate'='+proj=laea +lat_0=37.71859032558816 +lon_0=-86.484375',
#     'polar'='+proj=laea +lat_0=90.0 +lon_0=-117.0703125')

regions = unique(sites$region)
for(i in 1:length(regions)){

    subset = sites[sites$region == regions[i],
        c('regionsite', 'longitude', 'latitude')]

    #choose reasonable projection
    meanlat = mean(subset$latitude, na.rm=TRUE)
    meanlon = mean(subset$longitude, na.rm=TRUE)
    if(meanlat <= 15 && meanlat >= -15){ #equatorial
        PROJ4 = paste0('+proj=laea +lon_0=', meanlon)
    # } else if(meanlat >= 75 || meanlat <= -75){ #polar
    } else { #temperate or polar
        PROJ4 = paste0('+proj=laea +lat_0=', meanlat, ' +lon_0=', meanlon)
    }

    #convert to spatial object and project from WGS 84
    subset = subset %>%
        st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
        st_transform(PROJ4)

    #get DEM, 14 is highest res, broadest area; 1 is lowest res
    dem = get_elev_raster(subset, z=4)
    mapview(dem) + mapview(subset)

    #generate a box and check topo basemap for full watershed capture
    box = st_bbox(dem) %>% st_as_sfc()

    # Double check that we captured the whole watershed
    mapview(box) +  mapview(subset)

    #Save files so that whitebox can call the data
    writeRaster(dem, filename='dem.tif', overwrite=TRUE)
    st_write(subset, 'sites.shp', delete_layer=TRUE)

    #Fill single cell pits (for hydrologic correctness)
    # fill_single_cell_pits('./dem.tif', './breach.tif')

    #uses lindasy's algorithm (better than depression filling)
    # breach_depressions('./breach.tif', './filled.tif', fill_pits=TRUE)
    breach_depressions('./dem.tif', './breached.tif', fill_pits=TRUE)

    #D8 pointer (what is a pointer? a flow direction grid? probably)
    d8_pointer('./breached.tif', './d8_pntr.tif')

    #D8 flow
    d8_flow_accumulation('./breached.tif', './d8_flow.tif',
        out_type='catchment area')

    #snap_points
    snap_pour_points('./sites.shp', './d8_flow.tif', './snapped_sites.shp', 100)

    #Watershed delineation as "whole watersheds'
    unnest_basins('./d8_pntr.tif','./snapped_sites.shp','./basins/sheds.tif')

    #Read in flow accumulation algorithm
    fac = raster('./d8_flow.tif')

    #Get a list of the watershed created by `unnest_basins`
    sheds = list.files('./basins', full.names=TRUE)

    #Create a function that uses the stars package to transform
    # the raster watershed outlines into shapefiles
    shed_stacker = function(x){
        read_stars(sheds[x]) %>%
            st_as_sf(merge=TRUE, use_integer=TRUE) %>%
            rename(id=1) %>%
            group_by(id) %>%
            summarize()
    }

    ## Use purrr::map to apply the raster-shapefile transformation to all
    ## rasters to a list of shapefiles (map_dfr doesn't play nice with sf for
    ## unknown reasons)
    s = purrr::map(1:length(sheds), shed_stacker)
}

#Use do.call to bind these sf objects into a single one
shape_sheds = do.call('rbind', s) %>% arrange(id)

#subset flow accumulation
fac_sub = crop(fac, shape_sheds)

mapview(fac_sub) +
    mapview(shape_sheds)
    # mapview(streamstats)

#crop and mask elevatr DEM
only = dem %>%
    crop(., shape_sheds) %>%
    mask(., shape_sheds)

# #Convert to matrix so rayshader is happy
# pmat <- matrix(raster::extract(only,raster::extent(only),buffer=300),
#     nrow=ncol(only),ncol=nrow(only))
#
# #Generate a hillshade
# raymat = ray_shade(pmat,sunangle=330)
#
# #use rayshader commands to generate map
# #rglwidget embeds output in html
# pmat %>%
#     sphere_shade(texture='desert') %>%
#     add_shadow(raymat) %>%
#     plot_3d(pmat,zscale=10,fov=0,theta=135,zoom=0.75,phi=45,
#         windowsize=c(750,750))
# rglwidget()
