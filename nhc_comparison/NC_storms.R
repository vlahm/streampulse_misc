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

stormdaily = read.csv('~/Downloads/NHC_stormdaily.csv', stringsAsFactors=FALSE,
    colClasses=c('date'='Date'))

# stormdaily$GPP = any(is.nan(stormdaily$GPP))
during = as.logical(stormdaily$storm)
day1pre = lead(during)
day2pre = lead(during, 2)
day3pre = lead(during, 3)
day1post = lag(during)
day2post = lag(during, 2)
day3post = lag(during, 3)
storm_summary = matrix(NA, 5, 7,
    dimnames=list(c('GPPmn', 'ERmn', 'GPP:ERmn', 'GPPse', 'ERse'),
    c('day3pre', 'day2pre', 'day1pre', 'during', 'day1post', 'day2post', 'day3post')))
for(i in 1:ncol(storm_summary)){
    gppmn = mean(stormdaily$GPP[get(colnames(storm_summary)[i])], na.rm=TRUE)
    storm_summary[1, i] = gppmn
    ermn = mean(stormdaily$ER[get(colnames(storm_summary)[i])], na.rm=TRUE)
    storm_summary[2, i] = ermn
    storm_summary[3, i] = gppmn / ermn
    gppsd = sd(stormdaily$GPP[get(colnames(storm_summary)[i])], na.rm=TRUE)
    storm_summary[4, i] = gppsd
    ersd = sd(stormdaily$ER[get(colnames(storm_summary)[i])], na.rm=TRUE)
    storm_summary[5, i] = ersd
}

error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}

barplot(storm_summary[1:3,], beside=TRUE, ylim=c(-7, 5))
s1 = seq(1.5, 3.5, 1)
# s2 = seq(4, 22, 3)
xv = c(s1, s1 + 4, s1 + 7, s1 + 10, s1 + 13, s1 + 16, s1 + 19)
error.bars(xv, c(storm_summary[1:3,]), upper=rep(0, length(xv)),
    lower=rep(-2, length(xv)), cap.length=0.01)
#task: think about how to identify and compare pre-post storm data.
#how should we define the beginning and end of a storm?
#how much time should we compare on either side?
#what summary statistics (mean, max, etc) should we use?
#looking at the data using the qaqc tool will help with thinking about this stuff.
#code it up if you feel like it! We can go through it together too.
