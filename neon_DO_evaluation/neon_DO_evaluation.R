library(lubridate)

setwd('~/git/streampulse/other_projects/neon_DO_evaluation/')

d = read.csv('NEON_DO2.csv', stringsAsFactors=FALSE)
colnames(d) = c('id','reg','site','time','var','val','flag','upload')
head(d)
str(d)
# a = d$time[640235:640240]
# a = d$time[600000:700000]
# identical(a, as.character(as.POSIXct(a, tz='UTC')))
# a; as.POSIXct(a, tz='UTC')
d$time = as.POSIXct(d$time, tz='UTC')

#get minmax datetimes and deployment time ####
mins = tapply(d$time, d$site, min)
maxes = tapply(d$time, d$site, max)
mms = cbind(mins, maxes)
mms = apply(mms, 2, function(x) as.POSIXct(x, origin=as.Date('1970-01-01')))
out = as.data.frame(mms)
out$mins = as.POSIXct(out$mins, origin=as.Date('1970-01-01'))
tz(out$mins) = 'UTC'
out$maxes = as.POSIXct(out$maxes, origin=as.Date('1970-01-01'))
tz(out$maxes) = 'UTC'

out$maxes
out$mins
out$deployTime_days = floor(as.numeric(out$maxes - out$mins))

#explore issue with most days having either 1440 or 1200 records.
#(turns out there are regular gaps in the data, maybe for maintenance or correction?)
sum(is.na(d$val))
o = tapply(d$time[d$site == 'WALK-down'],
    substr(d$time[d$site == 'WALK-down'], 1, 10), length)
table(o)

pdf(width=5, height=7, file='~/Dropbox/streampulse/figs/NEON_DO_gap_mystery2.pdf',
    compress=FALSE)
par(mfrow=c(4,1), mar=c(0,4,0,0), oma=c(6,4,4,4))
dd = d[d$site == 'WALK-down' & substr(d$time, 1, 10) == '2016-08-13',]
plot(dd$time, dd$val, type='p', pch=20, xaxt='n', xlab='', ylab='')
legend('topright', '2016-08-13', bty='n')
dd = d[d$site == 'WALK-down' & substr(d$time, 1, 10) == '2016-08-18',]
plot(dd$time, dd$val, type='p', pch=20, xaxt='n', xlab='', ylab='')
legend('topright', '2016-08-18', bty='n')
dd = d[d$site == 'WALK-down' & substr(d$time, 1, 10) == '2016-08-26',]
plot(dd$time, dd$val, type='p', pch=20, xaxt='n', xlab='', ylab='')
legend('topright', '2016-08-26', bty='n')
dd = d[d$site == 'WALK-down' & substr(d$time, 1, 10) == '2016-08-31',]
plot(dd$time, dd$val, type='p', pch=20, xlab='', ylab='')
legend('topright', '2016-08-31', bty='n')
mtext('DO (mgL)', 2, outer=TRUE)
mtext('Time', 1, outer=TRUE, line=3)
dev.off()

#get n days of O2 coverage ####
coverage = sapply(unique(d$site), FUN=function(x){
    samp_per_day = tapply(d$time[d$site == x], substr(d$time[d$site == x],1,10), length)
    nFullDays = sum(samp_per_day == 1440)
    coverage_full = nFullDays / length(samp_per_day)
    nFullAnd240SampGapDays = sum(samp_per_day == 1200 | samp_per_day == 1440)
    coverage_fullAndPartial = nFullAnd240SampGapDays / length(samp_per_day)
    return(list(nFullDays, coverage_full, nFullAnd240SampGapDays,
        coverage_fullAndPartial))
})

coverage = t(coverage)
colnames(coverage) = c('nFullDays', 'coverage_full', 'nFullAnd240SampGapDays',
    'coverage_fullAndPartial')

#final touches on table ####
out = merge(out, coverage, by='row.names')
out$mins = NULL
out$maxes = NULL
colnames(out)[1] = 'site'

# out = ooo
out = apply(out, 2, function(x) I(x))
out$coverage_full = sapply(out$coverage_full, round, 2)
out$coverage_fullAndPartial = sapply(out$coverage_fullAndPartial, round, 2)

#write and read fixes list issue
write.csv(out, 'NEON_DO_coverage_summary_utc2.csv', row.names=FALSE) #not proper df format
out = read.csv('NEON_DO_coverage_summary_utc2.csv', stringsAsFactors=FALSE) #proper format
out1 = read.csv('NEON_DO_coverage_summary_utc.csv', stringsAsFactors=FALSE) #first run

#plot coverage ####
out$site = sapply(out$site, function(x) return(x))
out$deployTime_days = sapply(out$deployTime_days, function(x) return(x))
out$nFullDays = sapply(out$nFullDays, function(x) return(x))
out$nFullAnd240SampGapDays = sapply(out$nFullAnd240SampGapDays, function(x) return(x))

zz = t(as.matrix(as.data.frame(out[c(2,3,5)])))

#first run
out1$site = sapply(out1$site, function(x) return(x))
out1$deployTime_days = sapply(out1$deployTime_days, function(x) return(x))
out1$nFullDays = sapply(out1$nFullDays, function(x) return(x))
out1$nFullAnd240SampGapDays = sapply(out1$nFullAnd240SampGapDays, function(x) return(x))

zz1 = t(as.matrix(as.data.frame(out1[c(2,3,5)])))

pdf(width=8, height=9, file='~/Dropbox/streampulse/figs/NEON_DO_coverage2.pdf',
    compress=FALSE)
hline_seq = seq(100, 500, 100)
par(mfrow=c(2,1))
barplot(zz1, beside=TRUE, names.arg=out1$site, las=2, main='Sept 2018',
    legend.text=c('deploy time (days)', '# full coverage days',
        '# full and partial coverage days (regular gap of 240 samples)'),
    args.legend=list(x='top', bty='n', cex=0.8), cex.names=0.7, ylim=c(0,650))
abline(h=hline_seq, lty=3, col=1:length(hline_seq))
barplot(zz, beside=TRUE, names.arg=out$site, las=2, main='Nov 2018',
    args.legend=list(x='top', bty='n', cex=0.8), cex.names=0.7, ylim=c(0,650))
abline(h=hline_seq, lty=3, col=1:length(hline_seq))
mtext('NEON DO Coverage', 3, outer=TRUE)
dev.off()

#get coverage across sites by day (38 sites total) ####
d$day = substr(d$time, 1, 10)
cov2 = tapply(d$val, list(d$day, d$site), function(x) length(!is.na(x)))
breakdown_by_day = data.frame(prop_max_coverage_across_sites = apply(cov2, 1, function(x){
    round(sum(x, na.rm=TRUE) / (1440 * 38), 2) #max samples per day * nsites
}))
breakdown_by_day$nsites_with_data = apply(cov2, 1, function(x){
    sum(!is.na(x))
})
breakdown_by_day$nsites_full_coverage = apply(cov2, 1, function(x){
    sum(x == 1440, na.rm=TRUE)
})
breakdown_by_day$nsites_1440_or_1200 = apply(cov2, 1, function(x){
    sum(x %in% c(1200,1440), na.rm=TRUE)
})

breakdown_by_day = breakdown_by_day[order(breakdown_by_day$prop_max_coverage,
    decreasing=TRUE),]
write.csv(breakdown_by_day, 'coverage_by_day2.csv')

#see if there are any years with decent coverage across sites ####
library(tidyr)
d2 = d[,c('site', 'time', 'val')]
d2 = d2[-c(5655161, 5723437),] #remove duplicates
dd = spread(d2, key='site', value='val')
sites = unique(substr(colnames(dd)[-1], 1, 4))

png(width=10, height=7, filename='~/Dropbox/streampulse/figs/NEON_DO_series2.png',
    units='in', type='cairo', res=300)
if(length(sites) %% 2 == 1) print('WARNING: odd number of plots - switch commenting')
par(mfcol=c(ceiling(length(sites)/2) + 1, 2), mar=rep(0,4), oma=c(2,0,0,0)) #
# par(mfcol=c(ceiling(length(sites)/2), 2), mar=rep(0,4), oma=c(2,0,0,0))
plot(1,1,axes=FALSE, xlab='', ylab='', main='', type='n')#
legend('center', legend=c('upstream', 'downstream'), lwd=c(4,1),#
    col=c('orange','purple'), bty='n')#
xmin = min(d2$time, na.rm=TRUE)
xmax = max(d2$time, na.rm=TRUE)
for(i in 1:length(sites)){
    ymin = min(d2$val[d2$site == paste(sites[i], 'up', sep='-')], na.rm=TRUE)
    ymax = max(d2$val[d2$site == paste(sites[i], 'up', sep='-')], na.rm=TRUE)
    if(i==14){ #
        plot(1,1,axes=FALSE, xlab='', ylab='', main='', type='n')#
    } #
    if(i %in% c(ceiling(length(sites) / 2),length(sites))){
        plot(d2$time[d2$site == paste(sites[i], 'up', sep='-')],
            d2$val[d2$site == paste(sites[i], 'up', sep='-')],
            type='p', col='orange', pch=20, cex=1,
            ylim=c(ymin, ymax), ylab='', xlab='', yaxt='n',
            xlim=c(xmin, xmax))
    } else {
        plot(d2$time[d2$site == paste(sites[i], 'up', sep='-')],
            d2$val[d2$site == paste(sites[i], 'up', sep='-')],
            type='p', col='orange', pch=20, cex=1,
            ylim=c(ymin, ymax), ylab='', xlab='', yaxt='n',
            xlim=c(xmin, xmax), xaxt='n')
    }
    legend('topleft', sites[i], bty='n')
    try(points(d2$time[d2$site == paste(sites[i], 'down', sep='-')],
        d2$val[d2$site == paste(sites[i], 'down', sep='-')], col='purple', pch='.'))
}
# plot(1,1,axes=FALSE, xlab='', ylab='', main='', type='n')
# legend('center', legend=c('upstream', 'downstream'), lwd=c(4,1),
#     col=c('orange','purple'), bty='n')

dev.off()

par(mfrow=c(1,1))
sitesfull = colnames(dd)[-1]

plot(d2$time[d2$site == sitesfull[1]], d2$val[d2$site == sitesfull[1]],
    type='l')


