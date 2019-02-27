#questions:
#combine data from historic sites? NHC UNHC %in% Concrete, Blackwood, Wood Bridge?
#combine years?

library(StreamPULSE)
library(viridis)

nhc_68_70 = read.csv('~/git/streampulse/other_projects/nhc_comparison/hall_table_15.csv')

#subset historic data by site and year
# gpp_concrete = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Concrete']
# gpp_blackwood = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Blackwood']
# gpp_wb = nhc_68_70$GPP_gO2m2d[nhc_68_70$site == 'Wood Bridge']

gpp_68_70 = nhc_68_70$GPP_gO2m2d
gpp_68 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1968']
gpp_69 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1969']
gpp_70 = nhc_68_70$GPP_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1970']

er_68_70 = nhc_68_70$ER_gO2m2d
er_68 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1968']
er_69 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1969']
er_70 = nhc_68_70$ER_gO2m2d[substr(nhc_68_70$date, 1, 4) == '1970']

#retrieve contemporary data by year
query_available_results('NC', 'NHC')
nhc_17 = request_results(sitecode='NC_NHC', '2017')
nhc_18 = request_results(sitecode='NC_NHC', '2018')

nhc_17 = as.data.frame(nhc_17$predictions[,c('date','GPP','ER')])
nhc_18 = as.data.frame(nhc_18$predictions[,c('date','GPP','ER')])

# gpp_17 = ts(nhc_17$GPP)
# gpp_18 = ts(nhc_18$GPP)
gpp_17 = nhc_17$GPP
gpp_18 = nhc_18$GPP
gpp_17_18 = c(gpp_17, gpp_18)

er_17 = nhc_17$ER
er_18 = nhc_18$ER
er_17_18 = c(er_17, er_18)

#time-series comparison of means then and now? ####
plot(gpp_17, type='l')
acf(gpp_17, na.action=na.pass)
pacf(gpp_17, na.action=na.pass)
#strong autocorrelation and partial autocorr;
#will have to model error as an autoregressive process
#not stationary; can't use pure AR
qqnorm(gpp_17); abline(1, 1, col='red', lty=2)
#normalish; no need to go bayesian

#nhc_17 is irregular, so can't use arima; only option would be GAM


#non-timeseries comparison of distributions then and now ####

# plot(density(gpp_17, na.rm=TRUE), xlim=c(-3, 10), main='', xlab='GPP',
#     ylim=c(0,1.1))
# par(new=TRUE)
# plot(density(gpp_concrete, na.rm=TRUE), xlim=c(-3, 10), main='', xlab='',
#     ylab='', ylim=c(0,1.1))

# hist(gpp_17, breaks=15, main='2017; n=264', xlab='GPP')
# hist(gpp_concrete, breaks=15, main='1968-70 (Concrete); n=60', xlab='GPP')

par(mfrow=c(2,2))

#plot GPP dists, then and now
plot(density(gpp_68_70, na.rm=TRUE), xlim=c(-3, 10), bty='l', col='sienna3',
    main='GPP 1968-70 vs. 2017-18', xlab='GPP', ylim=c(0,1.1))
lines(density(gpp_17_18, na.rm=TRUE), col='blue')
legend('topright', legend=c('68-70; n=79', '17-18; n=475'),
    col=c('sienna3','blue'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot GPP dists by year
cols = viridis(5)
plot(density(gpp_68, na.rm=TRUE), xlim=c(-3, 10), bty='l', col=cols[1],
    main='GPP by year', xlab='GPP', ylim=c(0,1.3))
lines(density(gpp_69, na.rm=TRUE), col=cols[2])
lines(density(gpp_70, na.rm=TRUE), col=cols[3])
lines(density(gpp_17, na.rm=TRUE), col=cols[4])
lines(density(gpp_18, na.rm=TRUE), col=cols[5])
legend('topright',
    legend=c('68; n=18', '69; n=49', '70; n=12', '17; n=264', '18; n=211'),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot ER dists, then and now
plot(density(er_17_18, na.rm=TRUE), xlim=c(-15, 1), bty='l', col='sienna3',
    main='ER 1968-70 vs. 2017-18', xlab='ER', ylim=c(0,0.7))
lines(density(er_68_70 * -1, na.rm=TRUE), col='blue')
legend('topleft', legend=c('68-70; n=79', '17-18; n=475'),
    col=c('sienna3','blue'), lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

#plot ER dists by year
cols = viridisyy(6)
plot(density(er_68 * -1, na.rm=TRUE), xlim=c(-18, 2), bty='l', col=cols[1],
    main='ER by year', xlab='ER', ylim=c(0,0.7))
lines(density(er_69 * -1, na.rm=TRUE), col=cols[2])
lines(density(er_70 * -1, na.rm=TRUE), col=cols[3])
lines(density(er_17, na.rm=TRUE), col=cols[4])
lines(density(er_18, na.rm=TRUE), col=cols[5])
legend('topleft',
    legend=c('68; n=18', '69; n=49', '70; n=12', '17; n=264', '18; n=211'),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)
