#===============================================================================
#GPP-ER scatterplot by site
#Created 12/4/2018
#===============================================================================
gpp_er_scatter <- function(synthesis_metrics, FLUXNET_metrics, plot_type){
  #Cumulative annual GPP plot
    if(plot_type == "cumulative" ){
      plot(synthesis_metrics[, "ann_GPP_C"], synthesis_metrics[, "ann_ER_C"], 
        ylim = c(-5300, 0), xlim = c(0, 5300), pch = 20, col = "#bae5f1", 
        xlab = expression("mean annual cumulative GPP" ~ (gC ~ m^{-2} ~ y^{-1})),
        ylab = expression("mean annual cumulative ER" ~ (gC ~ m^{-2} ~ y^{-1}))
      )
    
      points(FLUXNET_metrics[, "ann_GPP"], FLUXNET_metrics[, "ann_ER"], pch = 20,
        col = "#c5982e")
    
      abline(0, -1)
      
      #Site numbers
        text(4500, -250, paste0(sum(!is.na(FLUXNET_metrics[, "ann_GPP"])), " sites"),
          col = "#c5982e")
      
        text(4500, -600, paste0(sum(!is.na(synthesis_metrics[, "ann_GPP_C"])), " sites"),
          col = "#bae5f1")
      
    } #End if statement
  
  #P95 GPP plot  
    if(plot_type == "upper"){
      plot(synthesis_metrics[, "upper_GPP_C"], synthesis_metrics[, "lower_ER_C"], 
        ylim = c(-30, 0), xlim = c(0, 30), pch = 20, col = "#bae5f1", 
        xlab = expression("P95 GPP" ~ (gC ~ m^{-2} ~ d^{-1})),
        ylab = expression("P05 ER" ~ (gC ~ m^{-2} ~ d^{-1}))
      )   
      
      points(FLUXNET_metrics[, "upper_GPP"], FLUXNET_metrics[, "lower_ER"], pch = 20, 
        col = "#c5982e")
    
      abline(0, -1)   
      
      #Site numbers
        text(4500, -250, paste0(sum(!is.na(FLUXNET_metrics[, "upper_GPP"])), " sites"),
          col = "#c5982e")
      
        text(4500, -600, paste0(sum(!is.na(synthesis_metrics[, "upper_GPP_C"])), " sites"),
          col = "#bae5f1")      
    } #End if statement

} #End gpp_er_scatter function 
