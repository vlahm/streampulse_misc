pcaset = metr %>%
    select(sitecode, gpp_C_mean, er_C_mean, Stream_PAR_sum, PAR_sum, Disch_ar1,
        Disch_cv, width_calc, PAR_mean, MOD_ann_GPP, MOD_ann_ER) %>%
    mutate(log_gpp_C_mean=log(gpp_C_mean), log_er_C_mean=log(gpp_C_mean),
        log_MOD_ann_GPP=log(MOD_ann_GPP), log_MOD_ann_ER=log(MOD_ann_ER))
    # select(-gpp_C_meaCn, er_C_mean)
write.csv(pcaset, 'output/candidate_vars.csv', row.names=FALSE)
pcaset = pcaset[complete.cases(pcaset),]
out = prcomp(pcaset, scale=TRUE)

pdf(file='output/candidate_var_pca.pdf', width=8, height=8)
autoplot(out, loadings=TRUE, loadings.label=TRUE)
dev.off()
