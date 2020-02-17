sitesub = sites %>%
    select(sitecode, Lat, Lon, TOTDASQKM, TOTDASQKM_corr)

pcaset = metr %>%
    select(sitecode, gpp_C_mean, er_C_mean, Stream_PAR_sum, PAR_sum, Disch_ar1,
        Disch_cv, width_calc, PAR_mean, MOD_ann_GPP, MOD_ann_ER, ER_K, num_days) %>%
    mutate(log_gpp_C_mean=log(gpp_C_mean), log_er_C_mean=log(er_C_mean * -1),
        log_MOD_ann_GPP=log(MOD_ann_GPP), log_MOD_ann_ER=log(MOD_ann_ER * -1)) %>%
    left_join(sitesub, by='sitecode')
    # select(-gpp_C_mean, -er_C_mean, -MOD_ann_GPP, -MOD_ann_ER) %>%
    # filter_at(vars(-sitecode, -width_calc), any_vars(! is.na(.)))
write.csv(pcaset, 'output/candidate_vars_etc.csv', row.names=FALSE)
corrmat = cor(pcaset[, -1], use='na.or.complete', method='kendall')
write.csv(round(corrmat, 2), 'output/day2/candidate_var_corrmat_kendall.csv')
# pcaset = pcaset[complete.cases(pcaset),]
out = prcomp(pcaset, scale=TRUE)

pdf(file='output/candidate_var_pca.pdf', width=8, height=8)
autoplot(out, loadings=TRUE, loadings.label=TRUE)
dev.off()
