# library(animation)
library(stringr)
library(RMariaDB)
library(DBI)
library(ks)

#read in mysql pw
conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
# conf = readLines('/home/aaron/sp/config.py')
ind = which(lapply(conf, function(x) grepl('MYSQL_PW', x)) == TRUE)
pw = str_match(conf[ind], '.*\\"(.*)\\"')[2]

#read in site and results tables from mysql
con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)
site = dbReadTable(con, "site")
results = dbReadTable(con, "results")
doy = as.numeric(strftime(results$solar_date, format='%j'))

#get list of fitted model names available on server
modnames = dir('../model_viz/data', pattern='modOut')
sitenm_all = str_match(modnames, 'modOut_(\\w+_\\w+)_[0-9]{4}')[,2]

#isolate those that are public (no embargo or past embargo period)
public_sites = difftime(Sys.time(), site$addDate,
    units='days') > site$embargo * 365
site = site[public_sites, c('region','site')]
sitenames_public = paste(site[,1], site[,2], sep='_')

#filter available models so that only public ones can be viewed.
#legacy code, so unnecessarily convoluted
modelnames_public = intersect(sitenm_all, sitenames_public)
sitenames = sitenm_all[sitenm_all %in% modelnames_public]
sitenames = unique(sitenames)

#compute overall kernel density so that it need not always be recomputed
overall_kernel = kde(na.omit(results[, c('GPP','ER')]))

#plot stuff ####
setwd('/home/mike/Dropbox/streampulse/figs/gifs/overall_kernel_density/')

plot_gif_slide = function(tmin, tmax, overlay, mon){

    par(mar=c(4,4,2,0), oma=rep(0,4))

    #get plot limits and overlay info
    overlay = overlay[overlay != 'None']
    regionsite = sapply(overlay, strsplit, '_')
    regions = unique(sapply(regionsite, function(x) x[1]))
    sites = unique(sapply(regionsite, function(x) x[2]))
    ylims = c(-15, 10)
    xlims = c(-5, 10)

    #overall plot
    overall_kernel = kde(na.omit(results[doy > tmin & doy < tmax,
        c('GPP','ER')]))

    plot(overall_kernel, xlab='', las=1, xaxt='n', ylab='', yaxt='n',
        ylim=ylims, xlim=xlims, display='filled.contour',
        col=c(NA, "gray70", "gray50", "gray30"))

    #overlay
    results_sub = na.omit(results[results$region %in% regions &
            results$site %in% sites & doy > tmin & doy < tmax,
        c('GPP', 'ER')])

    site_kernel = kde(results_sub)

    par(new=TRUE)
    plot(site_kernel, xlab='', las=1, xaxt='n', ylab='', yaxt='n',
        ylim=ylims, xlim=xlims, display='filled.contour',
        col=c(NA, adjustcolor('red', alpha.f=0.5),
            adjustcolor('orange', alpha.f=0.5),
            adjustcolor('yellow', alpha.f=0.5)))

    #peripheral plot items
    abline(0, -1, col='black', lty=3)
    axis(1, tcl=-0.2, padj=-1)
    axis(2, tcl=-0.2, hadj=0.5, las=1)
    mtext(expression(paste("GPP (g"~O[2]~"m"^"-2"~" d"^"-1"*")")),
        1, line=1.8)
    mtext(expression(paste("ER (g"~O[2]~"m"^"-2"~" d"^"-1"*")")),
        2, line=2)
    mtext('Metabolism of Brewery Creek, WI relative to all StreamPULSE sites', 3)

    #legends
    legend("bottomleft", c("75%", "50%", "25%"), bty="n", bg='white',
        lty=c(1,1,1), lwd=4, col=c("gray70", "gray50", "gray30"),
        seg.len=1, box.col='transparent', horiz=TRUE, title='Overall density')

    legend('bottom', '1:1', bg='white', lty=3, bty='n')

    legend("bottomright", c("75%", "50%", "25%"), bty="n", bg='white',
        lty=c(1,1,1), lwd=4, col=c("red", "orange", "yellow"),
        seg.len=1, box.col='transparent', horiz=TRUE,
        title='Density of subset')

    legend('topright', mon, col='transparent', bty='n', cex=1.5, text.font=2)
}

#gif stuff ####
tminseq = seq(1, 365, 31)
tmaxseq = c(seq(31, 365, 31), 365)

for(i in 1:12){
    png(width=8, height=8, units='in', res=300,
        filename=paste0('slide', sprintf('%02d', i), '.png'), type='cairo')
    plot_gif_slide(tminseq[i], tmaxseq[i], 'WI_BRW', mon=month.abb[i])
    dev.off()
    print(paste(i, 'done'))
}

system("convert -delay 150 slide*.png overall_kernel_density.gif")

dbDisconnect(con)
