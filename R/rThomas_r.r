#------------------------------------------------------------------------------
# simulate spatial ThomaS process
# the function has to be applied species by species.
# the function is an R reimplementation of Felix May's cpp efficient re-implementation of the rThomas function from the spatstat package

rThomas_r <- function(n_points,
                      n_mother_points=NA,
							 xmother=NA,
							 ymother=NA,
                      sigma,
                      mu,     # mean number of points per cluster, not used anymore
                      xmin = 0,
                      xmax = 1,
                      ymin = 0,
                      ymax = 1) {
							 
	stopifnot(!is.na(n_mother_points) | (!is.na(xmother) & !is.na(ymother)))
	
   if (n_mother_points > 0 | (length(na.omit(xmother))>0 & all(xmother!="no clustering"))) {	# or any(xmother=="no clustering"
		if(!is.na(n_mother_points) & all(is.na(xmother)) & all(is.na(ymother))) {	# if n_mother_points, xmother and ymother are given n_mother_points is overridden. This actually can't happen since this case is already controled by sim_thomas_coords()
			xmother = runif(n_mother_points, xmin, xmax)
			ymother = runif(n_mother_points, ymin, ymax)
		} else {
			n_mother_points <- length(xmother)
		}
		
		xydat <- matrix(NA, n_points, 2, dimnames=list(c(), c("x","y")))
      for (ipoint in 1:n_points){
            imother = sample(x=1:n_mother_points, size=1, replace=T)	# random selection of one of the n_mother_points mother points
				while(sum(is.na(xydat[ipoint, ])) != 0) {	# while + ifelse, there must be a more efficient way to control that simulated coordinates stay in the range.
					dxy = rnorm(2, 0.0, sigma)
					xydat[ipoint, 1] = ifelse((xmother[imother] + dxy[1]) >= xmin & (xmother[imother] + dxy[1]) <= xmax, xmother[imother] + dxy[1], NA)
					xydat[ipoint, 2] = ifelse((ymother[imother] + dxy[2]) >= ymin & (ymother[imother] + dxy[2]) <= ymax, ymother[imother] + dxy[2], NA)
				}
      }
   } # end if n_mother_points > 0

   else {
      xpoints = runif(n_points, xmin, xmax)
      ypoints = runif(n_points, ymin, ymax)
		xydat <- data.frame(x=xpoints, y=ypoints)
   }


   return(xydat)
}
