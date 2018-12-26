library(lubridate)
library(plotrix)
library(RColorBrewer)

setwd('~/git/streampulse/other_projects/neon_coverage_analysis/')

#read in all variables required for metabolism modeling; extract to CSV like so:
    # select * from data where upload_id=-900 and variable='Discharge_m3s' and
    # DateTime_UTC like '2017%' into outfile '/var/lib/mysql-files/NEON_Q_2017.csv'
    # fields terminated by ',' enclosed by '"' lines terminated by '\n';
cols = c('id','reg','site','time','var','val','flag','upload')
DO = read.csv('NEON_DO_2017.csv', stringsAsFactors=FALSE, col.names=cols)
depth = read.csv('NEON_depth_2017.csv', stringsAsFactors=FALSE, col.names=cols)
sat = read.csv('NEON_DOsat_2017.csv', stringsAsFactors=FALSE, col.names=cols)
# Q = read.csv('NEON_Q_2017.csv', stringsAsFactors=FALSE, col.names=cols)
temp = read.csv('NEON_temp_2017.csv', stringsAsFactors=FALSE, col.names=cols)

# vars = list('DO'=DO, 'depth'=depth, 'sat'=sat, 'Q'=Q, 'temp'=temp)
vars = list('DO'=DO, 'depth'=depth, 'sat'=sat, 'temp'=temp)

#get sample count for each day of year, site, and variable (as 3d array)
allSites = sort(Reduce(union, lapply(vars, function(x) unique(x[,'site']))))
sampByDay = array(NA, dim=c(365, length(allSites), length(vars)),
    dimnames=list('DOY'=1:365, 'site'=allSites, 'var'=names(vars)))
for(i in 1:length(vars)){
    v = vars[[i]]
    sites = unique(v$site)
    for(s in sites){
        sitedf = v[v[,'site'] == s, ]
        DOYcount = data.frame('nsamp'=tapply(sitedf$val, yday(sitedf$time), length))
        z = merge(DOYcount, data.frame('NAs'=rep(NA, 365)), by=0, all=TRUE)
        DOYcountFill = z[order(as.numeric(z$Row.names)), 'nsamp']
        sampByDay[, dimnames(sampByDay)$site == s, i] = as.numeric(DOYcountFill)
    }
}

#mutiple watertemp and depth counts by 5, as these are collected on 5m intervals
sampByDay[,,c(2,4)] = sampByDay[,,c(2,4)] * 5

#assign coverage indices to each DOY-site pair (sum across variables)
getCoverageIndex = function(arr, maxDailySamp){
    m = maxDailySamp
    out = apply(arr, MARGIN=c('DOY', 'site'), sum, na.rm=TRUE)
    out[is.na(out)] = 0
    out[out >= (m * 4)] = 1
    out[out < (m * 4) & out >= (m * 3)] = .75
    out[out < (m * 3) & out >= (m * 2)] = .5
    out[out < (m * 2) & out >= (m * 1)] = .25
    out[out < (m * 1) & out >= 1.1] = 0
    # out[! is.na(out)] = 1
    if(all(out != 1)) out[which(out == 0)[1]] = 1
    if(all(out != .75)) out[which(out == 0)[1]] = .75
    if(all(out != .5)) out[which(out == 0)[1]] = .5
    if(all(out != .25)) out[which(out == 0)[1]] = .25
    out = t(out)
}
coverage_index = getCoverageIndex(sampByDay, 1440)

#plot

png(width=10, height=7, file='~/Dropbox/streampulse/figs/NEON_coverage_2017.png',
    units='in', type='cairo', res=300)
cols = brewer.pal(5, 'YlGnBu')
cols[1] = 'white'
par(mar=c(3,5,4,1))
color2D.matplot(coverage_index, extremes=cols, border=NA,
    axes=FALSE, xlab='', ylab='', yrev=FALSE)
abline(h=1:(length(allSites)-1), col='gray60', lty=2)
# color.legend(166, 20, 176, 25, rect.col=cols, align='rb', gradient='y',
#     legend=c('[0-25)','[25-50)','[50-75)','[75-100)','100%'))
color.legend(286, 28.9, 356, 29.7, rect.col=cols, xpd=NA,
    legend=c('[0-25)','','[50-75)','','100%'), align='lt', cex=0.8)
color.legend(286, 28.9, 356, 29.7, rect.col=cols, xpd=NA,
    legend=c('','[25-50)','','[75-100)',''), align='rb', cex=0.8)
axis(1, seq(1, 365, 5), rep('', length(seq(1, 365, 5))), cex.axis=0.8, tck=-.008)
axis(1, seq(1, 365, 5), seq(1, 365, 5), cex.axis=0.8, tcl=0, col='transparent', line=-.5)
axis(2, seq(0.5, nrow(coverage_index) - 0.5, 1), rownames(coverage_index),
    las=2, cex.axis=0.8)
mtext('DOY', 1, line=1.5, font=2)
mtext(paste('Daily % coverage across all variables required for metabolism ',
    'modeling (NEON aquatic stations; 2017)'), 3, line=.5, font=2, adj=0, cex=0.8)
dev.off()
