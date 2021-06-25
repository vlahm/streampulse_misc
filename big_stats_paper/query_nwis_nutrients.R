siteyrs = readRDS('siteyears.rds')
sites = names(siteyrs)

#to be loopified:
i = j = 1

site = sites[i]
yr = siteyrs[[i]][j]

startdate = paste0(yr, '-01-01')
enddate = paste0(yr, '-12-31')
vcds = '00060,00065'

url = paste0('https://nwis.waterservices.usgs.gov/nwis/iv/?format=json&sites=',
    site, '&startDT=', startdate, 'T00:00Z&endDT=', enddate,
    'T23:59Z&parameterCd=', vcds, '&siteStatus=all')

resp = httr::GET(url)
json = httr::content(resp, as='text', encoding='UTF-8')
data = jsonlite::fromJSON(json)

#now we must traverse the data labyrinth. good luck! Seeing this, I think I'll
#email Jordan and Alison after all. Maybe they know some shortcuts.
str(data) #!!!
head(data$value$timeSeries$values[[1]])
View(data$value$timeSeries$variable)
data$value$timeSeries$variable$
