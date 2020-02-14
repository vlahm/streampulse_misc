#' Filters metabolism data based on selection criteria
#' @description This funciton filters the combined metabolism data based on
#' the calculated diagnostics for each site-year.
#'
#' @param diag The site-year diagnostics calculated from diagnostics_fun
#' @param filters The filtering criteria. For example "num_days >= 180". Multiple
#' criteria can also be used, for example "num_days >= 180 & ER_K <= 0.5".
#' @return Returns a filtered set of metabolism data
#'
#' @export
#'
#===============================================================================
#Defining a function to filter the dataset
#===============================================================================
filter_metab <- function(diag, filters, metab_rds){
  #Catch if no filters are applied
    if(is.null(filters)){}

  #Filter the site-years (Update this so it is not deprecated)
    sy_filter <- dplyr::filter_(diag, filters)

  #Split the filtered sites by Site_ID
    site_split <- split(sy_filter, sy_filter[, "Site_ID"])

  #Helper function to retrieve the standardized data that meets filtering requirements
    retrieve_metab <- function(Site_ID){
      #Years that meet the filtering requirement
        years <- site_split[[Site_ID]][, "Year", drop=TRUE]

      #Get the standardized metabolism for the year of interest
        metab_SOI <- metab_rds[[Site_ID]]

      return(metab_SOI[metab_SOI[, "Year"] %in% years, ])

    } #End retrieve_metab function

  #Get the filtered metabolism data
    filtered <- lapply(names(site_split), retrieve_metab)
      names(filtered) <- names(site_split)

  return(filtered)

} #End filter_metab function
