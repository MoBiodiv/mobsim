#####################################################################
# Functions for species and endemics-area realtionships (SAR + EAR)
# and for the distance decay of similarity
# All functions work on point-mapped data, not on plot based surveys!
#
# Felix May, 19. 02. 2016
#####################################################################


# # -----------------------------------------------------------
# # calculates number of endemic species in subsquares of a plot
# endemics.tess <- function(prop.loss=0.125,community)
# {
#    require(spatstat)
#
#    x <- community[,1]
#    y <- community[,2]
#
#    xmin <- floor(min(x)); xmax <- ceiling(max(x))
#    ymin <- floor(min(y)); ymax <- ceiling(max(y))
#
#    dx <- xmax - xmin
#    dy <- ymax - ymin
#
#    pp1 <- ppp(x,y,marks=community[,3],window=owin(c(xmin,xmax),c(ymin,ymax)))
#
#    area.loss <- dx*dy*prop.loss
#    square.size <- sqrt(area.loss)
#
#    if (square.size <= min(c(dx,dy))){
#       grid1 <- tess(xgrid=seq(xmin,xmax,by=square.size),
#                     ygrid=seq(ymin,ymax,by=square.size))
#    } else
#    {
#       if (dx>=dy){
#          grid1 <- tess(xgrid=seq(xmin,xmax,by=prop.loss*dx),
#                        ygrid=seq(ymin,ymax,by=dy))
#       } else {
#          grid1 <- tess(xgrid=seq(xmin,xmax,by=dx),
#                        ygrid=seq(ymin,ymax,by=prop.loss*dy))
#       }
#    }
#
#    pp.grid <- split(pp1,grid1)
#
#    # generate species x sites presence absence matrix
#    pa.mat <- sapply(pp.grid,function(pp){table(marks(pp))>0})
#    if (ncol(pa.mat) > 1){
#       nEndemics <- colSums(pa.mat[rowSums(pa.mat)==1,])
#    } else
#    {
#       nEndemics <- sum(pa.mat)
#    }
#
#    return(c(meanEnd=mean(nEndemics),sdEnd=sd(nEndemics)))
# }


# -----------------------------------------------------------
# get the number of species and endemics in a rectangle with upper left corner in x0,y0
# size of the rectangle: xsize*ysize
nSpec.nEnd.rect <- function(x0,y0,xsize,ysize,community)
{
   x <- community[,1]
   y <- community[,2]
   id.spec <- community[,3]

   # logical vector which trees are in the sampling rectangle
   in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   spec.in <- unique(id.spec[in.rect])
   spec.out <- unique(id.spec[!in.rect])

   nSpecies <- length(spec.in)
   nEndemics <- length(spec.in[!spec.in %in% spec.out])

   return(c(nSpec=nSpecies,nEnd=nEndemics))
}

# -----------------------------------------------------------
# get the abundance dsitribution in a rectangle with upper left corner in x0,y0
# size of the rectangle: xsize*ysize
abund.rect <- function(x0,y0,xsize,ysize,community)
{
   x <- community[,1]
   y <- community[,2]
   id.spec <- community[,3]

   # logical vector which trees are in the sampling rectangle
   in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   abund <- table(id.spec[in.rect])

   return(abund)
}

# -----------------------------------------------------------
# calculates number of species and endemic species in randomly located squares or rectangles
# the area is provided as proportion of the total
# maximum coordinates of the landscape are defined in xext and yext (min and max)
nSpec.nEnd.rand.rect <- function(prop.A=0.25,community,nrect=100,xext=c(0,1),yext=c(0,1))
{
   x <- community[,1]
   y <- community[,2]

   dx.plot <- xext[2] - xext[1]
   dy.plot <- yext[2] - yext[1]

#    xmin <- floor(min(x)); xmax <- ceiling(max(x))
#    ymin <- floor(min(y)); ymax <- ceiling(max(y))
#
#    dx.plot <- xmax - xmin
#    dy.plot <- ymax - ymin

   area <- dx.plot*dy.plot*prop.A
   square.size <- sqrt(area)

   if (square.size <= min(c(dx.plot,dy.plot))){
      dx.rect <- square.size
      dy.rect <- square.size
   } else
   {
      if (dx.plot>=dy.plot){
         dx.rect <- dx.plot*prop.A
         dy.rect <- dy.plot
      } else {
         dx.rect <- dx.plot
         dy.rect <- dy.plot*prop.A
      }
   }

   xpos <- runif(nrect,min=0,max=xext[2]-dx.rect)
   ypos <- runif(nrect,min=0,max=yext[2]-dy.rect)

   system.time(nSpec.End <- mapply(nSpec.nEnd.rect,xpos,ypos,
                                   MoreArgs=list(xsize=dx.rect,ysize=dy.rect,community=community)))
   #
   #    system.time({
   #    Endemics <- numeric(nrect)
   #    id.spec <- community[,3]
   #    for (i in seq_along(xpos)){
   #       # logical vector which trees are in the sampling rectangle
   #       in.rect <- (x >= xpos[i] & x < (xpos[i]+dx.rect) & y >= ypos[i] & y < (ypos[i]+dy.rect))
   #
   #       spec.in <- unique(id.spec[in.rect])
   #       spec.out <- unique(id.spec[!in.rect])
   #
   #       nEndemics[i] <- length(spec.in[!spec.in %in% spec.out])
   #    }  }
   #    )

   return(c(meanSpec=mean(nSpec.End["nSpec",]),
            sdSpec=sd(nSpec.End["nSpec",]),
            meanEnd=mean(nSpec.End["nEnd",]),
            sdEnd=sd(nSpec.End["nEnd",])))
}


# -----------------------------------------------------------
# estimates SAR amd EAR based on random squares or rectangles
SAR.EAR.rand <- function(community,prop.A=seq(0.1,1,by=0.1),nsamples=100,xext=c(0,1),yext=c(0,1))
{
   nscales <- length(prop.A)
   SAR.EAR <- data.frame(Area = prop.A,
                         nspec.mean = numeric(nscales),
                         nspec.sd = numeric(nscales),
                         nend.mean = numeric(nscales),
                         nend.sd = numeric(nscales))

   sar.ear1 <- sapply(prop.A[-length(prop.A)],
                      nSpec.nEnd.rand.rect,
                      community=community,
                      nrect=nsamples,
                      xext=xext,
                      yext=yext)

   SAR.EAR$nspec.mean[1:(nscales-1)] <- sar.ear1["meanSpec",]
   SAR.EAR$nspec.sd[1:(nscales-1)] <- sar.ear1["sdSpec",]

   SAR.EAR$nend.mean[1:(nscales-1)] <- sar.ear1["meanEnd",]
   SAR.EAR$nend.sd[1:(nscales-1)] <- sar.ear1["sdEnd",]

   SAR.EAR$nspec.mean[nscales] <- length(unique(community[,3]))
   SAR.EAR$nspec.sd[nscales] <- 0

   SAR.EAR$nend.mean[nscales] <- length(unique(community[,3]))
   SAR.EAR$nend.sd[nscales] <- 0

   return(SAR.EAR)
}

# -----------------------------------------------------------
# Distance decay
distance.decay <- function(community, prop.A=0.05, nsamples=30, xext=c(0,1), yext=c(0,1),method="jaccard")
{
   dx.plot <- xext[2] - xext[1]
   dy.plot <- yext[2] - yext[1]

   area <- dx.plot*dy.plot*prop.A
   square.size <- sqrt(area)

   xpos <- runif(nsamples,min=0,max=xext[2]-square.size)
   ypos <- runif(nsamples,min=0,max=yext[2]-square.size)

   d <- dist(cbind(xpos,ypos))

   com.tab <- mapply(abund.rect,xpos,ypos,
                     MoreArgs=list(xsize=square.size,ysize=square.size,community=community))

   jaccard <- 1 - vegdist(t(com.tab),method=method,binary=T)
   dat.out <- data.frame(distance = as.numeric(d), similarity = as.numeric(jaccard))
   return(dat.out)
}



#--------------------------------------------------------------------------------------------
#example data set all trees >= 10 cm dbh in BCI 2010
#
# bci1 <- read.table("ExampleData/BCI2010.txt",header=T)
# head(bci1)
#
# pArea <- c(0.001,0.005,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1)
# SAR.EAR.bci <- SAR.EAR.rand(bci1,prop.A=pArea,xext=c(0,1000),yext=c(0,500))
#
# library(plotrix)
# plot(pArea*5e5,SAR.EAR.bci$nspec.mean,log="",ylim=c(1,250))
# plotCI(pArea*5e5,SAR.EAR.bci$nspec.mean,uiw=SAR.EAR.bci$nspec.sd,ylim=c(0,240),add=T)
# plotCI(pArea*5e5,SAR.EAR.bci$nend.mean,uiw=SAR.EAR.bci$nend.sd,add=T,col=2)
