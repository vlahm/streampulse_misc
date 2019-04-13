library(tidyr)
library(dplyr)

#select site of interest
focus = 'NHC'

#read in full streampulse dataset (not including NEON or Powell)
sp = read.csv('~/Downloads/all_sp_data.csv', stringsAsFactors=FALSE)

#filter all but NC sites
sp = sp[sp$regionID == 'NC',]

#get rid of regionID since this is a huge dataset and it's just wasting memory
sp$regionID = NULL

#grab the NHC subset
subset = sp[sp$siteID == focus,]

#remove the main dataset from the global environment if youre low on RAM
rm(sp)

#remove duplicate rows
subset = subset[! duplicated(subset),]

#convert datetime from string to POSIXct
subset$dateTimeUTC = as.POSIXct(subset$dateTimeUTC, tz='UTC')

#remove bad data
subset = subset[subset$flagID != 'Bad Data',]

#get indices of all points that have been flagged with the word "storm".
#i just flagged all likely storms for Eno 2017 and NHC 2017
storm_related_rows = grep('storm', subset$flagComment, ignore.case=TRUE)

#filter all but storm-related rows (jk)
# subset = subset[storm_related_rows,]

#get set of all storm related flag comments
unique(subset$flagComment[storm_related_rows])

#get datetimes associated with those i just flagged (comment == 'storm')
#(you could include others, but look at them with the cleaning tool first)
storm_rows = grepl('storm', subset$flagComment)
storm_dt = subset$dateTimeUTC[storm_rows]

#filter all but storm rows (jk)
# subset = subset[subset$dateTimeUTC %in% storm_dt,]

#average rows with the same siteID, datetime, and variable (for shape format step below)
subset = aggregate(value ~ siteID + dateTimeUTC + variable, mean,
    data=subset, na.action=NULL)

# create boolean column for storm or not storm
subset$storm = as.numeric(subset$dateTimeUTC %in% storm_dt)

# convert from long to wide format
stormdata = tidyr::spread(subset, variable, value)

#bring in metabolism estimates; deal with \\N characters from MySQL
modres = read.csv('~/Downloads/all_daily_model_results.csv',
    stringsAsFactors=FALSE)
modres[modres == '\\N'] = NA
numind = grepl('(GPP|ER|K600)', colnames(modres))
modres[,numind] = apply(modres[,numind], 2, as.numeric)

#convert datetime to date (make new column)
modres$date = as.Date(as.POSIXct(modres$solar_date))

#clean
modres = filter(modres, site == focus & year != 1907) %>%
    select(-region, -site, -year, -solar_date)

#average storm data by day and merge with metab estimates
stormdaily = group_by(stormdata, 'date'=as.Date(substr(dateTimeUTC, 1, 10))) %>%
    summarize_if(is.numeric, mean, na.rm=TRUE) %>%
    left_join(modres, by='date') %>%
    as.data.frame()

#see what you're dealing with
head(stormdaily)

#plot
ylims = range(stormdaily$GPP, stormdaily$ER, na.rm=TRUE)
plot(stormdaily$date, stormdaily$GPP, type='l', col='red', ylim=ylims,
    xlab='', ylab='g O2 m^-2 d^-1', main='NHC storms, 2017',
    xlim=c(as.Date('2017-01-01'), as.Date('2017-12-31')))
lines(stormdaily$date, stormdaily$ER, col='blue')
abline(v=stormdaily$date[stormdaily$storm > 0], col='gray', lty=2)

#save csv
write.csv(stormdata, '~/Desktop/NHC_stormdata.csv', row.names=FALSE)
write.csv(stormdaily, '~/Desktop/NHC_stormdaily.csv', row.names=FALSE)

#task: think about how to identify and compare pre-post storm data.
#how should we define the beginning and end of a storm?
#how much time should we compare on either side?
#what summary statistics (mean, max, etc) should we use?
#looking at the data using the qaqc tool will help with thinking about this stuff.
#code it up if you feel like it! We can go through it together too.
