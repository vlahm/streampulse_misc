#===============================================================================
#Ordered GPP plot function
#Created 2/6/2020
#===============================================================================
ordered_gpp_plot <- function(synthesis_metrics, FLUXNET_metrics, plot_type){
  #Add some labels to the site metrics
    synthesis_metrics$Source <- "aquatic"
    synthesis_metrics$color <- "#bae5f1"
  
    FLUXNET_metrics$Source <- "terrestrial"
    FLUXNET_metrics$color <- "#c5982e"
    
  #Compile the data into a single table
    synthesis_sub <- setNames(synthesis_metrics[, c("Site_ID", "ann_GPP_C", "upper_GPP_C", 
      "ann_ER_C", "lower_ER_C", "Source", "color")], c("Site_ID", "ann_GPP", "upper_GPP", 
      "ann_ER", "lower_ER", "Source", "color"))
      
    bound <- rbind(synthesis_sub, FLUXNET_metrics)  
    
  #Cumulative annual GPP plot
    if(plot_type == "cumulative"){
      #Order the GPP
        gpp_ord <- bound[order(bound$ann_GPP, decreasing = TRUE), 
          c("Source", "color", "ann_GPP")]
        
      #Add in an x position
        gpp_ord$xloc <- c(nrow(gpp_ord):1)
       
      #Make the barplot 
        barplot(gpp_ord[, "ann_GPP"], col = gpp_ord[, "color"], xlab = "Site", 
          ylab = "Cumulative GPP")
        
      #Add text on number of sites
        text(300, 4000, paste0(sum(gpp_ord[, "Source"] == "terrestrial"), " sites"), 
          col = "#c5982e")
 
        text(300, 3750, paste0(sum(gpp_ord[, "Source"] == "aquatic"), " sites"), 
          col = "#bae5f1")
        
    } #End if statement 

  #P95 GPP plot  
    if(plot_type == "upper"){
      #Order the GPP
        gpp_ord <- bound[order(bound$upper_GPP, decreasing = TRUE), c("Source", "color", "upper_GPP")]
        
      #Add in an x position
        gpp_ord$xloc <- c(nrow(gpp_ord):1)
       
      #Make the barplot 
        barplot(gpp_ord[, "upper_GPP"], col = gpp_ord[, "color"], xlab = "Site", ylab = "P95 GPP")
        
      #Add text on number of sites
        text(300, 20, paste0(sum(gpp_ord[, "Source"] == "terrestrial"), " sites"), 
          col = "#c5982e")
 
        text(300, 18.75, paste0(sum(gpp_ord[, "Source"] == "aquatic"), " sites"), 
          col = "#bae5f1")        
    } #End if statement   
    
} #End ordered_gpp_plot function
