library(tidyr)

#read in full streampulse dataset (not including NEON or Powell)
sp = read.csv('~/Downloads/all_sp_data.csv', stringsAsFactors=FALSE)

#filter all but NC sites
sp = sp[sp$regionID == 'NC',]

#get rid of regionID since this is a huge dataset and it's just wasting memory
sp$regionID = NULL

#grab the Eno subset
subset = sp[sp$siteID == 'NHC',]

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

#see what you're dealing with
head(stormdata)

#a useless plot, but it verifies that stuff worked
plot(stormdata$dateTimeUTC, stormdata$DO_mgL, type='l')
points(stormdata$dateTimeUTC[stormdata$storm == 1],
    stormdata$DO_mgL[stormdata$storm == 1], col='red', pch=20, cex=0.1)

#save csv
write.csv(stormdata, '~/Desktop/NHC_stormdata.csv', row.names=FALSE)

#task: think about how to identify and compare pre-post storm data.
#how should we define the beginning and end of a storm?
#how much time should we compare on either side?
#what summary statistics (mean, max, etc) should we use?
#looking at the data using the qaqc tool will help with thinking about this stuff.
#code it up if you feel like it! We can go through it together too.
