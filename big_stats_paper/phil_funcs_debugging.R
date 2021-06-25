# diagnostics=diag; models=metab_d
filter_and_impute = function(diagnostics, models, ...){

    filt = filter_metab(diag=diagnostics, 'ER_K < 0.1', metab_rds=models)
    filt = filter_metab(diag=diagnostics, 'ER_K <= 1', metab_rds=models)
    sum(sapply(filt, nrow) == 0)
    imp = synthesis_gapfill(filt, PQ=1.25, block=Inf, pmiss=99)

    return(imp)
}

# metab_rds=models;
sitecode=names(site_split)[1]
filter_metab <- function(diag, ..., metab_rds){
    sy_filter <- dplyr::filter_(diag, 'ER_K < 0.1')
    sy_filter <- dplyr::filter_(diag, 'ER_K <= 1')
    site_split <- split(sy_filter, sy_filter[, "sitecode"])
    # sapply(site_split, nrow)
    retrieve_metab <- function(sitecode){
        years <- site_split[[sitecode]][, "Year"]
        metab_SOI <- metab_rds[[sitecode]]
        return(metab_SOI[metab_SOI[, "Year"] %in% years, ])
    }
    filtered <- lapply(names(site_split), retrieve_metab)
    # length(filtered)
    names(filtered) <- names(site_split)
    return(filtered)
} #End filter_metab function

synthesis_filtered=filt; PQ=1.25; block=Inf; pmiss=99
SOI=synthesis_filtered[[1]]
cnt_sg = 1
synthesis_gapfill <- function(synthesis_filtered, PQ, block, pmiss){
    site_gapfill <- function(SOI, PQ){
        print(cnt_sg)
        cnt_sg <<- cnt_sg + 1
        year_split <- split(SOI, SOI[, "Year"])
        filled <- do.call(rbind, lapply(year_split, FUN = gapfill_norm, PQ = PQ,
            block = block, pmiss = pmiss))
        return(filled)
    }
    # SOI = synthesis_filtered[[14]]
    gapfilled_data <- lapply(synthesis_filtered, FUN = site_gapfill, PQ = PQ)
    return(gapfilled_data)
} #End synthesis_gapfill function

ts = year_split[[1]]
gapfill_norm <- function(ts, PQ, block, pmiss){
    #If there is negative GPP, then replace with NA
    if(length(ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"]) !=0){
        ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"] <- NA
    } #End if statement
    if(length(ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"]) !=0){
        ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"] <- NA
    } #End if statement
    ts$NEP <- ts[, "GPP"] + ts[, "ER"]

    ts$GPP_filled <- fillMiss3(ts, "GPP", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)
    ts$ER_filled <- fillMiss3(ts, "ER", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)
    ts$NEP_filled <- fillMiss3(ts, "NEP", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)
    ts$Wtemp_filled <- fillMiss3(ts, "temp_water", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)
    ts$PAR_filled <- fillMiss3(ts, "PAR_sum", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)
    ts$Disch_filled  <- fillMiss3(ts, "discharge", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

    znorm <- function(timeseries, var){
        ts.mean <- mean(timeseries[, var], na.rm = TRUE)
        ts.dev <- sd(timeseries[, var], na.rm = TRUE)
        (timeseries[, var] - ts.mean)/ts.dev
    }
    ts$GPP_norm <- znorm(ts, "GPP_filled")
    ts$ER_norm <- znorm(ts, "ER_filled")
    ts$NEP_norm <- znorm(ts, "NEP_filled")
    ts$Wtemp_norm <- znorm(ts, "Wtemp_filled")
    ts$PAR_norm <- znorm(ts, "PAR_filled")

    if(all(is.na(ts[, "MOD_GPP"])) == "FALSE" &
            length(unique(na.omit(ts[, "MOD_GPP"]))) != 1){
        ts$MOD_GPP_filled <- fillMiss3(ts, "MOD_GPP", block = 150, pmiss = 99,
            model = "trend", smooth = TRUE, plot = FALSE)
        ts$MOD_NPP_filled <- fillMiss3(ts, "MOD_NPP", block = 150, pmiss = 99,
            model = "trend", smooth = TRUE, plot = FALSE)
        ts$MOD_GPP_norm <- znorm(ts, "MOD_GPP_filled")
        ts$MOD_NPP_norm <- znorm(ts, "MOD_NPP_filled")
    } #End if statement

    if(all(is.na(ts[, "MOD_GPP"])) == "TRUE" |
            length(unique(na.omit(ts[, "MOD_GPP"]))) == 1){
        ts$MOD_GPP_filled <- NA
        ts$MOD_NPP_filled <- NA
        ts$MOD_GPP_norm <- NA
        ts$MOD_NPP_norm <- NA
    } #End if statement

    ts$GPP_C <- ts[, "GPP"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$ER_C <- ts[, "ER"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$NEP_C <- ts[, "NEP"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$GPP_C_filled <- ts[, "GPP_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$ER_C_filled <- ts[, "ER_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$NEP_C_filled <- ts[, "NEP_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107

    final <- ts[, c("Date", "U_ID", "Year", "DOY", "GPP", "ER", "NEP", "GPP_C",
        "ER_C", "NEP_C", "K600", "DO_obs", "DO_sat", "temp_water", "discharge",
        "PAR_sum", "Stream_PAR_sum", "LAI_proc", "MOD_GPP", "MOD_NPP", "GPP_filled", "ER_filled",
        "NEP_filled", "GPP_C_filled", "ER_C_filled", "NEP_C_filled", "Wtemp_filled",
        "Disch_filled", "PAR_filled", "MOD_GPP_filled", "MOD_NPP_filled", "GPP_norm",
        "ER_norm", "NEP_norm","Wtemp_norm", "PAR_norm","MOD_GPP_norm", "MOD_NPP_norm")]
} #End gapfill_norm function

# dataset=ts; var="GPP";
# model = "trend"; smooth = TRUE; log = "y"; plot = FALSE
fillMiss3 <- function(dataset, var, block = 30, pmiss = 40, model = "trend",
    smooth = TRUE, plot = FALSE, ... ){
    pck <- is.na(dataset[, var])

    #Catch if the dataset contains only NA values
    if(all(pck == TRUE) == FALSE){
        if(sum(pck) > 0){
            #Calculate the percentage of missing values in the dataset
            percent <- round((sum(pck) / length(dataset[, var])) * 100,
                digits = 2)

            #Calculate the largest consecutive set of missing data
            rles <- rle(is.na(dataset[, var]))
            max.mis <- max(rles$lengths[rles$values])

        } else {
            percent <- 0
            max.mis <- 0
        }
        #If there are not too many missing values then proceed
        if(percent < pmiss & max.mis < block){
            my.series <- window(dataset[, var])
            #First value can't be NA for StructTS, replace NA with nearest value P.S. 2019
            if(is.na(my.series)[1]){my.series[1] <- zoo::na.locf(my.series,
                option = "nocb", na.remaining = "rev")[1]}

            tryCatch({
                my.struct <- StructTS(my.series, type = model)
            }, error=function(e){
                my.struct <- StructTS(my.series, type = "level")
            }, warning=function(w){
                my.struct <- StructTS(my.series, type = "level")
            })

            # struct_try <- try(StructTS(my.series, type = model), silent = TRUE)
            # if(class(struct_try)[1] != 'try-error'){
            #     my.struct <- struct_try
            # } #End if statement
            # if(class(struct_try)[1] == 'try-error'){
            #     message(paste0("type = level was used for ", var))
            #     my.struct <- StructTS(my.series, type = "level")
            # } #End if statement
            #Perform smoothing (if chosen)
            if(smooth) fit <- tsSmooth(my.struct) else fit <- fitted(my.struct)
            if(plot == TRUE){
                plot(my.series, typ = "l", lwd = 4, xlab = "Observation",
                    ylab = "Observed and estimated times series")
                lines(fit[, 1], col = "green")
                leg.txt <- c("Observed values", "New time series")
                legend("topleft", leg.txt, col = c("black", "green"), lwd = c(4, 1),
                    bty = "n", ncol = 2, cex = 0.8)
            } #End if statement

            #Gap-filling
            dataset$filled <- dataset[, var]
            dataset$filled[pck] <- fit[pck, 1]
            # plot(my.series, type='l'); lines(fit[,1], col='red')
            # lines(dataset$filled, col='blue')
            return(dataset[, "filled"])
        } #End if statement
        return(dataset[, var])
    } #End if statement

    #If the data contains only NA values
    if(all(pck == TRUE) == TRUE){
        return(rep(NA, nrow(dataset)))
    } #End if statement
} #End fillMiss3 function
