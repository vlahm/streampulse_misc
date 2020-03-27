sitesub = sites %>%
    select(sitecode, Lat, Lon, TOTDASQKM, TOTDASQKM_corr)

# pcaset = metr %>%
#     # left_join(aq_metab_cov$
#     select(sitecode, gpp_C_mean, er_C_mean, Stream_PAR_sum, PAR_sum, Disch_ar1,
#         Disch_cv, width_calc, PAR_mean, MOD_ann_GPP, MOD_ann_ER) %>%
#     mutate(log_gpp_C_mean=log(gpp_C_mean), log_er_C_mean=log(er_C_mean * -1),
#         log_MOD_ann_GPP=log(MOD_ann_GPP), log_MOD_ann_ER=log(MOD_ann_ER * -1)) %>%
#     left_join(sitesub, by='sitecode')
ER_K, num_days

pcaset = aq_metab_cov %>%
    left_join(metr, by='sitecode') %>%
    select(sitecode, gpp_C_mean=GPP_aq_mean, er_C_mean=ER_aq_mean,
        gpp_C_sum=GPP_aq_sum, er_C_sum=ER_aq_sum,
        Stream_PAR_sum, PAR_sum, Disch_ar1,
        Disch_cv, width_calc, PAR_mean, MOD_ann_GPP, MOD_ann_ER) %>%
    left_join(sitesub, by='sitecode')
all_vars = aq_metab_cov %>%
    left_join(metr, by='sitecode') %>%
    select(-gpp_C_mean, -er_C_mean) %>%
    rename(gpp_C_mean=GPP_aq_mean, er_C_mean=ER_aq_mean,
        gpp_C_sum=GPP_aq_sum, er_C_sum=ER_aq_sum) %>%
    left_join(sitesub, by='sitecode')

    # select(-gpp_C_mean, -er_C_mean, -MOD_ann_GPP, -MOD_ann_ER) %>%
    # filter_at(vars(-sitecode, -width_calc), any_vars(! is.na(.)))
write.csv(pcaset, 'output/final/candidate_vars.csv', row.names=FALSE)
write.csv(all_vars, 'output/final/all_vars.csv', row.names=FALSE)
write.csv(diag, 'output/final/model_diagnostics.csv', row.names=FALSE)
corrmat = cor(pcaset[, -1], use='na.or.complete', method='kendall')
write.csv(round(corrmat, 2), 'output/final/candidate_var_corrmat_kendall.csv')
# pcaset = pcaset[complete.cases(pcaset),]
out = prcomp(pcaset, scale=TRUE)

pdf(file='output/candidate_var_pca.pdf', width=8, height=8)
autoplot(out, loadings=TRUE, loadings.label=TRUE)
dev.off()
