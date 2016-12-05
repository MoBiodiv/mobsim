# -----------------------------------------------------------
#' Plotting basic sampling relationships of simulated community
#' @param comm \code{\link{community}} object which includes three items:
#' \enumerate{
#'    \item census - data.frame with three columns: X, Y and species names for each individual
#'    \item x_min_max - extent of the community in x-direction
#'    \item y_min_max - extent of the community in y-direction
#' }, output from \code{\link{Sim.Poisson.Community} or \code{\link{Sim.Thomas.Community}.
#' @param print - logical. If TRUE plot is plotted directly, otherwise is returned with an object to be plotted later, default is TRUE
#' @param col - the color to be used for plotting
#' @param add.ribbon - logical. If TRUE plots a 95% confidence interval around the species-area curve. Default is FALSE.
#' @param plot.to.add - output of ggplot
#' @param plot.type - character. Type of plot, either "SAD", "SAC", or "SAR".
#'
#' @examples
#' sim.com <- Sim.Thomas.Community(S = 10, N = 100, sigma=0.02, cv = 1)
#' spatial.plot(sim.com)
#' SAD.plot(sim.com)
#' SAC.plot(sim.com)
#' SAR.plot(sim.com)
#'
#' sim.com1 <- Sim.Thomas.Community(S = 10, N = 1000, sigma=0.02, cv = 1)
#' sim.com2 <- Sim.Thomas.Community(S = 100, N = 1000, sigma=0.02, cv = 1)
#' SAC.plot1 <- SAC.plot(sim.com1, print=F)
#' SAC.plot2 <- add.to.plot(sim.com2, plot.to.add=SAC.plot1, plot.type="SAC", col="red",print=T)

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
   require(ggplot2)
   spat.plot <- ggplot(sim.com$census, aes(X,Y, color=Species)) +
      geom_point() +
      guides(color="none") +
      xlab("") +
      ylab("") +
      theme_biodiviz()

   if(print==T) print(spat.plot)
   if(print==F) return(spat.plot)
}

##
SAD.plot <- function(sim.com, title="Species-abundance distribution",col="black", print=T){
   if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
   require(ggplot2)
   SAD <- data.frame(specID = names(table(sim.com$census$Species)),
                     abundance = as.integer(table(sim.com$census$Species)))

   SAD.plot <- ggplot(data=SAD, aes(abundance)) +
      geom_histogram(bins=ceiling(max(SAD$abundance/20)), fill=paste(col)) +
      ggtitle(paste(title)) +
      theme_biodiviz()

   if(print==T) print(SAD.plot)
   if(print==F) return(SAD.plot)
}

##
SAC.plot <- function(sim.com, title="Species-accumulation curve", col="black", print=T){
   if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
   require(ggplot2)
   SAD <- data.frame(specID = names(table(sim.com$census$Species)),
                     abundance = as.integer(table(sim.com$census$Species)))
   SAC <- data.frame(ind=1:sum(SAD$abundance),SR=SAC(sim.com))

   SAC.plot <- ggplot(data=SAC,aes(x=ind,y=SR)) +
      geom_line(color=paste(col)) +
      xlab("# individuals sampled") +
      ylab("Species richness") +
      ggtitle(paste(title)) +
      theme_biodiviz()

   if(print==T) print(SAC.plot)
   if(print==F) return(SAC.plot)
}

##
SAR.plot <- function(sim.com, title="Species-area relationship", col="black", add.ribbon=F, print=T){
   if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
   require(ggplot2)
   prop.a <- c(0.01,seq(0.1, 1, by = 0.1)) # size of area samples in proportions of total area
   n <- 50 # number of samples used for each size of area
   SAR <- data.frame(DivAR(sim.com, prop.A=prop.a, nsamples = n))

   SAR.plot <- ggplot(data=SAR,aes(x=propArea,y=meanSpec)) +
      geom_line(color=paste(col)) +
      xlab("sampled area/total area") +
      ylab("Species richness") +
      ggtitle(paste(title)) +
      theme_biodiviz()

   if(add.ribbon==T) {
      SAR.plot <- SAR.plot + geom_ribbon(aes(ymin=meanSpec-1.96*sdSpec,ymax=meanSpec+1.96*sdSpec),alpha=0.3)
      }

   if(print==T) print(SAR.plot)
   if(print==F) return(SAR.plot)
}

##
add.to.plot <- function(sim.com, plot.to.add, plot.type, col="red",print=T){
   if (class(sim.com)!="community")
      stop("sim.com is not of class 'community'.")
   if (!is.element(plot.type, c("SAD","SAR","SAC")))
      stop("Unknown plot type. plot.type must be either 'SAD', 'SAC', or 'SAR'")
   require(ggplot2)
   if(plot.type=="SAD"){
      SAD <- data.frame(specID = names(table(sim.com$census$Species)),
                        abundance = as.integer(table(sim.com$census$Species)))
      p <- plot.to.add + geom_histogram(data=SAD, aes(abundance), bins=ceiling(max(SAD$abundance/20)), color=paste(col))
   }
   if(plot.type=="SAC"){
      SAD <- data.frame(specID = names(table(sim.com$census$Species)),
                        abundance = as.integer(table(sim.com$census$Species)))
      SAC <- data.frame(ind=1:sum(SAD$abundance),SR=SAC(sim.com))
      p <- plot.to.add + geom_line(data=SAC,aes(x=ind,y=SR), color=paste(col))
   }
   if(plot.type=="SAR"){
      prop.a <- c(0.01,seq(0.1, 1, by = 0.1)) # size of area samples in proportions of total area
      n <- 50 # number of samples used for each size of area
      SAR <- data.frame(DivAR(sim.com, prop.A=prop.a, nsamples = n))
      p <- plot.to.add + geom_line(data=SAR,aes(x=propArea,y=meanSpec),color=paste(col))
   }
   if(print==T)   print(p)
   if(print==F)   return(p)
}

##
# -----------------------------------------------------------
#' Interactive plot of biodiversity patterns at multiple spatial scales
#'
#' Visualize multiple biodiversity patterns at multiple spatial scales. The tool
#' \enumerate{
#'    \item simulates locations of individuals of different species in a location (plot, area)
#'    \item plots biodiversity patterns such as species-abundance distributions (SAD), species accummulation curves (SAC), species-area relationships (SAR) and alike for a selected grid cell and for the entire region.
#' The function sim.Thomas.Community is used to simulate point pattern distributions of individual species.
#' The interactive interface relies on package `manipulate`.
#' }
#' Key parameters to vary using the sliders are:
#' \enumerate{
#'    \item S - total number of species,
#'    \item N - total number of individuals,
#'    \item eveness - evenness of the abundances,
#'    \item spat.agg - spatial aggregation of species,
#'    \item resolution - the number of rows or columns used to dissect the 'global' community (total number of cells = resolution²),
#'    \item cell - cell number for which 'local' relationship should be visualized. Counts start in the upper left corner rowwise.
#' }
#'
#' @examples
#' interactive.plot()
#'
interactive.plot <- function(){
   require(ggplot2) # for plotting
   require(gridExtra) # for multiple plots using grid.arrange()
   require(raster) # creating a grid polygon
   require(broom) # converting spatial polygon to dataframe for plotting with ggplot
   require(manipulate) # interactive plotting

   plot.sim.com.grid <- function(S.pool, S.max=s.max, N.pool, spat.agg, evenness, resolution, cell.id){
      set.seed(229376)
      sim.com <- Sim.Thomas.Community(S = S.pool, N = N.pool, sigma=spat.agg, cv = evenness)

      # create the grid
      grid.rast <- raster(extent(c(0,1,0,1)), nrows=resolution, ncols=resolution)
      grid.poly <- as(grid.rast, "SpatialPolygons")

      # extract the focal cell
      cell.poly <- as(grid.poly, "SpatialPolygons")[cell.id]

      # extract only the part of sim.com that falls into the focal grid
      is.in <- over(SpatialPoints(sim.com$census[,1:2]), cell.poly)
      is.in[is.na(is.in)] <- 0
      is.in <- is.in == 1
      sim.com.cell <- sim.com
      sim.com.cell$census <- sim.com.cell$census[is.in,]
      sim.com.cell$census$Species <- factor(sim.com.cell$census$Species) # updating species factor levels of the community subset

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
   }

   manipulate(
      plot.sim.com.grid(S.pool=S, S.max=500, N.pool=N, spat.agg=spat.agg, evenness=evenness,
                        resolution = res, cell.id=cell),
      S = slider(10,500,step=10),
      N = slider(5000,100000,step=1000),
      evenness = slider(0.5,3,step=0.5),
      spat.agg = slider(0.01,1),
      cell = slider(1, 400), # give error if cell > res^2 but manipulate doesn't work with dependent sliders, e.g cell=slider(1,res^2)
      res = slider(2, 20))
}
