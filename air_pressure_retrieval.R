library(devtools)
install_github('streampulse/StreamPULSE', dependencies=TRUE)
library(StreamPULSE)

#establish datetime range for which to retrieve air pressure.
#can also use e.g. 'EDT' or 'EST' if the rest of the data you're dealing with
#are in local time
start = as.POSIXct('2018-09-01 00:00:00', tz='UTC')
end = as.POSIXct('2018-11-18 23:45:00', tz='UTC')

#this uses lat and long for NC_NHC.
# the `:::` operator accesses internal functions of packages (those not
#normally intended to be accessed by users), whereas `::`
#can only access the collection of functions that are exposed when you attach
#a package.
air_pres = StreamPULSE:::FindandCollect_airpres(lat=35.9925, long=-79.046,
   start_datetime=start, end_datetime=end)
