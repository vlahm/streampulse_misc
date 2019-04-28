airpres_NCEP = function(lat, long, startyr, endyr=startyr){


   station = data.frame(lon=long, lat=lat)

   #format site data for use with geoknife package
   station = geoknife::simplegeom(as.data.frame(t(station)))

   years = as.numeric(startyr):as.numeric(endyr)
   cat('Acquiring air pressure',
      'data for', length(years),
      'year(s).\n\tEach year may take a few minutes.\n')

   #retrieve air pressure data from noaa
   pres = data.frame(datetime=.POSIXct(character()), pres=numeric())
   for(i in 1:length(years)){

      fabric = geoknife::webdata(url=paste0('https://www.esrl.noaa.gov/psd/th',
         'redds/dodsC/Datasets/ncep.reanalysis/surface/pres.sfc.',
         years[i], '.nc'), variables='pres')
      noaa_job = geoknife::geoknife(stencil=station, fabric=fabric, wait=TRUE)
      noaa_data = geoknife::result(noaa_job, with.units=TRUE)

      datcol = ifelse('1' %in% colnames(noaa_data), '1', 'V1')
      pres = rbind(pres, noaa_data[,c('DateTime', datcol)])

      cat('Year', i, 'complete.\n')
   }

   pres = data.frame(pres)

   colnames(pres)[colnames(pres) == 'X1'] = 'V1'
   df_out = pres %>% mutate(AirPres_kPa = V1 / 1000,
      DateTime_UTC=DateTime) %>% select(AirPres_kPa, DateTime_UTC)

   return(df_out)
}

airpres = airpres_NCEP(lat=35.9925, long=-79.046, startyr=2018)
