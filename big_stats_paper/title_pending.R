library(tidyverse)
rm(list=ls()); cat('/014')

setwd('~/git/streampulse/other_projects/big_stats_paper/')

metab = readRDS('phil_stuff/output/synthesis_standardized.rds')
diag = readRDS('phil_stuff/metab_synthesis/output/yearly_diagnostics.rds')
metr = as.data.frame(readRDS('phil_stuff/output/site_metrics.rds'))
# filled = readRDS('phil_stuff/output/synthesis_gap_filled.rds')

srcs = list.files('phil_stuff/metab_synthesis/R/functions/', full.names=TRUE)
for(x in srcs) source(x)

#filter metabolism data by number of days in siteyear; then gapfill
filters = paste0("num_days >= ", 365 * 0.6)
filt = filter_metab(diag, filters, metab)
imp = synthesis_gapfill(filt, PQ=1.25, block=150, pmiss=50)

str(imp[[1]])
