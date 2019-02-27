#read in all variables required for metabolism modeling; extract to CSV like so:
    # select * from data where upload_id != -900 and variable='Discharge_m3s' and
    # DateTime_UTC like '2017%' into outfile '/var/lib/mysql-files/SP_Q_2017.csv'
    # fields terminated by ',' enclosed by '"' lines terminated by '\n';

library(lubridate)
library(plotrix)
library(RColorBrewer)

# y = '2015'
y = '2016'
# y = '2017'
# y = '2018'

setwd(paste0('~/git/streampulse/other_projects/SP_coverage_analysis/', y))


cols = c('id','reg','site','time','var','val','flag','upload')
DO = read.csv(paste0('SP_DO_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)
depth = read.csv(paste0('SP_depth_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)
watPres = read.csv(paste0('SP_waterPres_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)
sat = read.csv(paste0('SP_DOsat_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)
Q = read.csv(paste0('SP_Q_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)
watTemp = read.csv(paste0('SP_waterTemp_', y, '.csv'), stringsAsFactors=FALSE,
    col.names=cols)

vars = list('DO'=DO, 'depth'=depth, 'sat'=sat, 'Q'=Q, 'temp'=watTemp,
    'watPres'=watPres)

#get sample count for each day of year, site, and variable (as 3d array)
allSites = sort(Reduce(union, lapply(vars, function(x){
    unique(paste0(x$reg, '_', x$site))
})))

#filter embargoed sites
# allSites = allSites[! allSites %in%

#count up number of samples per variable per day per site
nday = ifelse(y %in% as.character(seq(2000, 2060, 4)), 366, 365)
sampByDay = array(NA, dim=c(nday, length(allSites), length(vars) + 1),
    dimnames=list('DOY'=1:nday, 'site'=allSites,
    'var'=c(names(vars), 'derived_sat')))
for(i in 1:(length(vars) - 1)){
    v = vars[[i]]
    sites = unique(paste0(v$reg, '_', v$site))
    for(s in sites){
        sitedf = v[paste0(v$reg, '_', v$site) == s, ]
        DOYcount = data.frame('nsamp'=tapply(sitedf$val, yday(sitedf$time),
            length))
        z = merge(DOYcount, data.frame('NAs'=rep(NA, nday)), by=0, all=TRUE)
        DOYcountFill = z[order(as.numeric(z$Row.names)), 'nsamp']
        sampByDay[, dimnames(sampByDay)$site == s, i] = as.numeric(DOYcountFill)
    }
}

#collapse waterpres into depth where depth is missing (wp + airpres can be proxy)
sampByDay[,,'depth'] = apply(sampByDay[,,c('depth', 'watPres')], MAR=1:2,
    FUN=function(x){
        m = suppressWarnings(min(x, na.rm=TRUE))
        ifelse(is.infinite(m), NA, m)
    })

#count airpres + watertemp + DO as DOsat (airpres can be modeled)
# sat_tempo = sampByDay[,,'sat']
sampByDay[,,'derived_sat'] = apply(sampByDay[,,c('temp', 'DO')], MAR=1:2,
    FUN=function(x){
        m = suppressWarnings(min(x, na.rm=TRUE))
        ifelse(is.infinite(m), NA, m)
    })
sampByDay[,,'sat'] = apply(sampByDay[,,c('derived_sat', 'sat')], MAR=1:2,
    FUN=function(x){
        m = suppressWarnings(max(x, na.rm=TRUE))
        ifelse(is.infinite(m), NA, m)
    })
# sampByDay[,,'sat'] = sat_tempo[!is.na(sat_tempo)]
sampByDay = sampByDay[,,! dimnames(sampByDay)$var %in% c('derived_sat',
    'watPres')]

#mutiply watertemp and depth counts by 5, as these are collected on 5m intervals
# sampByDay[,,c(2,4)] = sampByDay[,,c(2,4)] * 5

#assign coverage indices to each DOY-site pair (sum across variables)
getCoverageIndex = function(arr, maxDailySamp){
    m = maxDailySamp
    out = apply(arr, MARGIN=c('DOY', 'site'), sum, na.rm=TRUE)
    out[is.na(out)] = 0
    out[out >= (m * 5)] = 1 #5 here is the total number of variables in the array
    out[out < (m * 5) & out >= (m * 4)] = .75
    out[out < (m * 4) & out >= (m * 3)] = .5
    out[out < (m * 3) & out >= (m * 2)] = .25
    out[out < (m * 2) & out >= 1.1] = 0

    #artificially make sure all colors are present (for correct plotting)
    if(all(out != 1)) out[which(out == 0)[1]] = 1
    if(all(out != .75)) out[which(out == 0)[1]] = .75
    if(all(out != .5)) out[which(out == 0)[1]] = .5
    if(all(out != .25)) out[which(out == 0)[1]] = .25
    out = t(out)
}
coverage_index = getCoverageIndex(sampByDay, 96)

#find least abundant variable for each site
sumsAcrossDays = apply(sampByDay, MAR=2:3, sum, na.rm=TRUE)
colnames(sumsAcrossDays) = c('O','Z','S','Q','T')
varAbundanceOrder = apply(sumsAcrossDays, 1, order)
varAbundanceRank = apply(varAbundanceOrder, 2, function(x){
    paste(colnames(sumsAcrossDays)[x], collapse='')
})

#remove corkbrook. it has none of these core variables for 2017
# coverage_index = coverage_index[rownames(coverage_index) != 'RI_CorkBrk',]
# varAbundanceRank = varAbundanceRank[names(varAbundanceRank) != 'RI_CorkBrk']
# coverage_index = coverage_index[rownames(coverage_index) != '_',]
# varAbundanceRank = varAbundanceRank[names(varAbundanceRank) != '_']

#remove embargoed sites (manually for now)
regex = '(KS_|SE_|IN_|MC[0-9]|NHC[0-9]|ArbSeep)'
coverage_index = coverage_index[! grepl(regex, rownames(coverage_index)),]
varAbundanceRank = varAbundanceRank[! grepl(regex, names(varAbundanceRank))]

#plot

png(width=3600, height=2100, units='px', type='cairo', res=300, antialias='default',
    file=paste0('~/Dropbox/streampulse/figs/sp_coverage_analyses/SP_coverage_',
    y, '.png'))
cols = brewer.pal(5, 'YlGnBu')
cols[1] = 'white'
par(mar=c(3,5,4,4))
color2D.matplot(coverage_index, extremes=cols, border=NA,
    axes=FALSE, xlab='', ylab='', yrev=FALSE)
abline(h=1:(length(allSites)-1), col='gray60', lty=2)
li = nrow(coverage_index)
color.legend(240, li + 1.3, 310, li + 2, rect.col=cols, xpd=NA,
    legend=c('[0-25)','','[50-75)','','100%'), align='lt', cex=0.8)
color.legend(240, li + 1.3, 310, li + 2, rect.col=cols, xpd=NA,
    legend=c('','[25-50)','','[75-100)',''), align='rb', cex=0.8)
axis(1, seq(1, 365, 5), rep('', length(seq(1, 365, 5))), cex.axis=0.8, tck=-.008)
axis(1, seq(1, 365, 5), seq(1, 365, 5), cex.axis=0.8, tcl=0, col='transparent', line=-.5)
axis(2, seq(0.5, nrow(coverage_index) - 0.5, 1), rownames(coverage_index),
    las=2, cex.axis=0.7)
mtext('DOY', 1, line=1.5, font=2)
mtext(paste0('Daily % coverage across all variables required for metabolism ',
    'modeling (StreamPULSE, ', y, ')'), 3, line=2.5, font=2, adj=0, cex=0.8)
mtext(paste('Variables: DO, DOsat and/or waterTemp + airPres + DO, discharge,',
    'depth and/or waterpres + airPres, waterTemp.'), 3, line=1.5, font=1,
    adj=0, cex=0.8)
mtext(paste('This analysis does not include air pressure or light,',
    'which can be modeled.'), 3, line=0.5, font=1, adj=0, cex=0.8)
mtext(paste('O = DO\nS = DO sat\nZ = depth\nQ = discharge\nT ',
    '= water temp'), 3, line=0.5, font=1, at=345, cex=0.7)
mtext('Variable\nabundance\nrank:\nL -> H', 3, line=0.5, font=2, at=377,
    cex=0.8, xpd=NA)
mtext('Site', 3, line=0.5, font=2, at=-13, cex=0.8, xpd=NA)
axis(4, seq(0.5, nrow(coverage_index) - 0.5, 1), varAbundanceRank,
    las=2, cex.axis=0.8, tck=0)
dev.off()
