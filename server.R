#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(RColorBrewer)
library(MoBspatial)
require(gridExtra) # for multiple plots using grid.arrange()
require(sp) # creating a grid polygon
require(broom) # converting spatial polygon to dataframe for plotting with ggplot


# Define server logic for slider examples
shinyServer(function(input, output) {
  
  # Reactive expression to compose a data frame containing all of
  # the values
  # sliderValues <- reactive({
  # 
  #   # Compose data frame
  #   data.frame(
  #     Name = c("S",
  #              "N",
  #              "cv.abund",
  #              "spatagg",
  #              "Cell",
  #              "Resolution"),
  #     Value = as.character(c(input$S,
  #                            input$N,
  #                            input$cv.abund,
  #                            input$spatagg,
  #                            input$cell,
  #                            input$Resolution)),
  #     stringsAsFactors=FALSE)
  # })
  # 
  # # Show the values using an HTML table
  # output$values <- renderTable({
  #   sliderValues()
  # })
  ## plot theme
  theme_biodiviz <- function (base_size = 12, base_family = "") {
    theme_bw() +
      theme(axis.title = element_text(size = rel(0.6)),
            axis.text = element_text(size = rel(0.6)),
            legend.text = element_text(size = rel(0.6)),
            legend.position = "bottom",
            plot.title = element_text(size=rel(0.8)))
  }
  
  ##
  spatial.plot <- function(sim.com, print=T){
    if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
    n.species <- length(table(sim.com$census$species))
    # Classic palette, with n.study colors
    coul = brewer.pal(8, "Accent") 
    # add more tones to this palette :
    coul = colorRampPalette(coul)(n.species)
    
    spat.plot <- ggplot(sim.com$census, aes(x,y, color=species)) +
      geom_point() +
      guides(color="none") +
      xlab("") +
      ylab("") +
      theme_biodiviz()
    
    if(print==T) print(spat.plot)
    return(spat.plot)
  }
  
  ##
  SAD.plot <- function(sim.com, title="Species-abundance distribution",col="black", print=T){
    if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
    require(ggplot2)
    SAD <- data.frame(abundance=as.numeric(community_to_sad(sim.com)))
    
    SAD.plot <- ggplot(data=SAD, aes(abundance)) +
      geom_histogram(bins=ceiling(max(SAD$abundance/10)), fill=paste(col)) +
      ggtitle(paste(title)) +
      theme_biodiviz()
    
    if(print==T) print(SAD.plot)
    return(SAD.plot)
  }
  
  ##
  SAC.plot <- function(sim.com, title="Species-accumulation curve", logxy=T, col="black", lty="solid", size=1, print=T){
    if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
    
    require(ggplot2)
    SAC <- spec_sample_curve(sim.com, method="rarefaction")
    
    if(logxy==T) {
      SAC.plot <- ggplot(data=SAC,aes(x=n,y=spec_rarefied)) +
        geom_line(color=paste(col), linetype=paste(lty), size=size) +
        xlab("# individuals sampled") +
        ylab("Species richness") +
        scale_x_log10() +
        scale_y_log10() +
        ggtitle(paste(title)) +
        theme_biodiviz()
    }
    
    if(logxy==F){
      SAC.plot <- ggplot(data=SAC,aes(x=n,y=spec_rarefied)) +
        geom_line(color=paste(col), linetype=paste(lty), size=size) +
        xlab("# individuals sampled") +
        ylab("Species richness") +
        ggtitle(paste(title)) +
        theme_biodiviz()
    }
    
    if(print==T) print(SAC.plot)
    return(SAC.plot)
  }
  
  ##
  SAR.plot <- function(sim.com, title="Species-area relationship", logxy=T, col="black", lty="solid", size=1, add.ribbon=F, print=T){
    if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
    require(ggplot2)
    prop.a <- c(0.01,seq(0.1, 1, by = 0.1)) # size of area samples in proportions of total area
    n <- 30 # number of samples used for each size of area
    SAR <- data.frame(divar(sim.com, prop_area=prop.a, n_samples = n))
    
    if(logxy==T){
      SAR.plot <- ggplot(data=SAR,aes(x=prop_area,y=m_species)) +
        geom_line(color=paste(col), linetype=paste(lty), size=size) +
        xlab("sampled area/total area") +
        ylab("Species richness") +
        scale_x_log10() +
        scale_y_log10() +
        ggtitle(paste(title)) +
        theme_biodiviz()
    }
    if(logxy==F){
      SAR.plot <- ggplot(data=SAR,aes(x=prop_area,y=m_species)) +
        geom_line(color=paste(col), linetype=paste(lty), size=size) +
        xlab("sampled area/total area") +
        ylab("Species richness") +
        ggtitle(paste(title)) +
        theme_biodiviz()
    }
    
    if(add.ribbon==T) {
      SAR.plot <- SAR.plot + geom_ribbon(aes(ymin=m_species-1.96*sd_species,ymax=m_species+1.96*sd_species),alpha=0.3)
    }
    
    if(print==T) print(SAR.plot)
    
    return(SAR.plot)
  }
  
  ##
  add.to.plot <- function(sim.com, plot.to.add, plot.type, col="red", lty="solid", size=1, print=T){
    if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
    if (!is.element(plot.type, c("SAD","SAR","SAC")))
      stop("Unknown plot type. plot.type must be either 'SAD', 'SAC', or 'SAR'")
    require(ggplot2)
    if(plot.type=="SAD"){
      SAD <- data.frame(abundance=as.numeric(community_to_sad(sim.com)))
      p <- plot.to.add + geom_histogram(data=SAD, aes(abundance), bins=ceiling(max(SAD$abundance/20)), fill=paste(col))
    }
    if(plot.type=="SAC"){
      SAC <- spec_sample_curve(sim.com, method="rarefaction")
      p <- plot.to.add + geom_line(data=SAC,aes(x=n,y=spec_rarefied), color=paste(col), linetype=paste(lty), size=size)
    }
    if(plot.type=="SAR"){
      prop.a <- c(0.01,seq(0.1, 1, by = 0.1)) # size of area samples in proportions of total area
      n <- 30 # number of samples used for each size of area
      SAR <- data.frame(divar(sim.com, prop_area=prop.a, n_samples = n))
      p <- plot.to.add + geom_line(data=SAR,aes(x=prop_area,y=m_species),color=paste(col), linetype=paste(lty), size=size)
    }
    if(print==T)   print(p)
    return(p)
  }
  
  # Show the values using an HTML table
  output$InteractivePlot <- renderPlot({
    set.seed(229376)

  #   abund1 <- sim_sad(s_pool = input$S, n_sim = input$N, cv_abund = input$cv.abund)
  #   plot(abund1, method = "rank")
  # })
  #   
    sim.com <- sim_thomas_community(s_pool = input$S, n_sim = input$N, sigma=input$spatagg, sad_coef=list(cv_abund = input$cv.abund), fix_s_sim = T)
    
    # create the grid
    grid.rast <- raster(extent(c(0,1,0,1)), nrows=input$Resolution, ncols=input$Resolution)
    grid.poly <- as(grid.rast, "SpatialPolygons")

    # extract the focal cell
    cell.poly <- as(grid.poly, "SpatialPolygons")[input$cell]

    # extract only the part of sim.com that falls into the focal grid
    is.in <- over(SpatialPoints(sim.com$census[,1:2]), cell.poly)
    is.in[is.na(is.in)] <- 0
    is.in <- is.in == 1
    sim.com.cell <- sim.com
    sim.com.cell$census <- sim.com.cell$census[is.in,]
    sim.com.cell$census$species <- factor(sim.com.cell$census$species) # updating species factor levels of the community subset

    # prepare the grid and the cell for ggplot
    grid.poly.gg <- fortify(grid.poly)
    cell.poly.gg <- fortify(cell.poly)

    ### plot SAD
    SAD.plot.global <- SAD.plot(sim.com, title= "Species-abundance distribution 'global'", print=F)
    SAD.plot.cell <- SAD.plot(sim.com.cell, title= "Species-abundance distribution 'local'", col="red", print=F)

    ### plot SAC
    SAC.plot.global <- SAC.plot(sim.com, title= "Species-accumulation curve 'global'", print=F)
    SAC.plot.cell <- SAC.plot(sim.com.cell, title= "Species-accumulation curve 'local'", print=F)

    ### plot SAR
    SAR.plot.global <- SAR.plot(sim.com, title= "Species-area relationship 'global'", print=F)
    SAR.plot.cell <- SAR.plot(sim.com.cell, title= "Species-area relationship 'local'", print=F)

    # plot community
    spat.plot <- spatial.plot(sim.com, print=F) +
      geom_polygon(data=grid.poly.gg, aes(y=lat, x=long, group=group), colour="grey", fill=NA, size=1.2) +
      geom_polygon(data=cell.poly.gg, aes(y=lat, x=long, group=group), colour="red", fill=NA, size=2)

    # plot all at once
    lay <- rbind(c(1,2),
                 c(3,4),
                 c(5,6),
                 c(7,7))
    grid.arrange(SAD.plot.global, SAD.plot.cell,
                 SAC.plot.global, SAC.plot.cell,
                 SAR.plot.global, SAR.plot.cell,
                 spat.plot,
                 ncol=ncol(lay), nrow=nrow(lay), heights=c(1, 1, 1, 3),
                 layout_matrix = lay)
  })


})