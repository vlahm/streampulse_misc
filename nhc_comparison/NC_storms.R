library(tidyr)

#read in full streampulse dataset (not including NEON or Powell)
sp = read.csv('~/Downloads/all_sp_data.csv', stringsAsFactors=FALSE)

#filter all but NC sites
sp = sp[sp$regionID == 'NC',]

#get rid of regionID since this is a huge dataset and it's just wasting memory
sp$regionID = NULL

#grab the Eno subset
eno = sp[sp$siteID == 'Eno',]

#remove duplicate rows
eno = eno[! duplicated(eno),]

#get indices of all points that have been flagged with the word "storm".
#i just flagged all likely storms for Eno 2017 and NHC 2017
storm_related_rows = grep('storm', eno$flagComment, ignore.case=TRUE)

#filter all but storm-related rows
eno = eno[storm_related_rows,]

#get set of all storm related flag comments
unique(eno$flagComment)

#get datetimes associated with those i just flagged.
#(you could include others, but look at them with the cleaning tool first)
storm_rows = grep('', eno$flagComment)
storm_dt = eno$dateTimeUTC[storm_rows]
eno = eno[eno$dateTimeUTC %in% storm_dt,]

#average rows with the same siteID, datetime, and variable in preparation for...
eno = aggregate(value ~ siteID + dateTimeUTC + variable, mean,
    data=eno, na.action=NULL)

# converting from long to wide format
eno = tidyr::spread(eno, variable, value)

#see what you're dealing with
head(eno)

#task: think about how to identify and compare pre-post storm data.
#how should we define the beginning and end of a storm?
#how much time should we compare on either side?
#what summary statistics (mean, max, etc) should we use?
#looking at the data using the qaqc tool will help with thinking about this stuff.
#code it up if you feel like it! We can go through it together too.
