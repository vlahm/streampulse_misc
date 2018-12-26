library(lubridate)
library(plotrix)

setwd('~/git/streampulse/other_projects/neon_coverage_analysis/')

cols = c('id','reg','site','time','var','val','flag','upload')
DO = read.csv('NEON_DO_2017.csv', stringsAsFactors=FALSE, col.names=cols)
depth = read.csv('NEON_depth_2017.csv', stringsAsFactors=FALSE, col.names=cols)
sat = read.csv('NEON_DOsat_2017.csv', stringsAsFactors=FALSE, col.names=cols)
# Q = read.csv('NEON_Q_2017.csv', stringsAsFactors=FALSE, col.names=cols)
temp = read.csv('NEON_temp_2017.csv', stringsAsFactors=FALSE, col.names=cols)

# vars = list('DO'=DO, 'depth'=depth, 'sat'=sat, 'Q'=Q, 'temp'=temp)
vars = list('DO'=DO, 'depth'=depth, 'sat'=sat, 'temp'=temp)

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

allVarsPresent = apply(sampByDay, MARGIN=c('DOY', 'site'), sum)
allVarsPresent[! is.na(allVarsPresent)] = 1
allVarsPresent[is.na(allVarsPresent)] = 0
allVarsPresent = t(allVarsPresent)
par(mar=c(4,5,1,1))
color2D.matplot(allVarsPresent, extremes=c('white','black'), border=NA,
    axes=FALSE, xlab='', ylab='')
abline(h=1:length(allSites), col='gray70', lty=1)
axis(1, seq(1, 365, 5), seq(1, 365, 5))
axis(2, seq(0.5, nrow(allVarsPresent) - 0.5, 1), rownames(allVarsPresent),
    las=2, cex.axis=0.8)
mtext('DOY', 1, line=2.5, font=2)
