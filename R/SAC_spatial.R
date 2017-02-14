
# -----------------------------------------------------------
#' Sample species richness
#'
#' Expected species richness in a random sample of fixed size
#'
#' @param n single integer - sample size in number of individuals
#' @param abund_vec integer vector - species abundance distribution of the community
#'
#' @return expected number of species in n individuals
#'
#' @details The expected number of species is calculated after Coleman 1982.
#'
#' @references
#' Coleman, B. D. et al. (1982). Randomness, area, and species richness. Ecology 63, 1121-1133
#'
spec_sample <- function(n, abund_vec)
{
  n_spec <- length(abund_vec)
  n_ind <- sum(abund_vec)
  spec_n <- n_spec - sum((1 - n/ n_ind)^abund_vec)
  return(spec_n)
}


#' species rarefaction curves
#'
#' Estimate expected species richness as a function of sample size
#'
#' @param abund_vec integer vector with species abundance distribution
#' @param method available methods are \code{"coleman"} or \code{"hurlbert"}.
#' The species richness estimations of both methods quickly converge for larger
#' numbers of individuals
#'
#' @return Numeric Vector with expected species richness in samples of 1, 2, 3 ... n individuals
#'
#' @references
#' Coleman, B. D. et al. 1982. Randomness, area, and species richness. Ecology 63, 1121-1133
#'
#' Hurlbert, S.H. 1971. The nonconcept of species diversity: a critique and alternative parameters.
#' Ecology 52, 577-586.
#'
#' @examples
#' sim_com1 <- SAD.lognorm(100, 10000)
#' rc1 <- rare_curve(sim_com1$abund)
#' rc2 <- rare_curve(sim_com1$abund, method = "hurlbert")
#'
#' plot(rc1, type = "l", log = "xy", xlab = "Sample size",
#'      ylab = "Expected species richness")
#' lines(1:length(rc2), rc2, lty = 2, col = 2)
#'
rare_curve <- function(abund_vec, method = "coleman")
{
  n <- sum(abund_vec)
  n_vec <- 1:n

  if (method == "hurlbert"){
     require(vegan)
     rc <- as.numeric(rarefy(abund_vec, sample = n_vec))
  } else {
     rc <- sapply(n_vec, spec_sample, abund_vec = abund_vec)
  }

  return(rc)
}

# -----------------------------------------------------------
#' Spatially-explicit species accumulation curve (SAC)
#'
#' The SAC is similar to the rarefaction curve, but samples of defined size
#' always contain nearest-neighbour individuals. Accordingly the SAC is
#' influenced by intraspecific aggregation, while the rarefaction curve is only
#' driven by the species abundance distribution.
#'
#' @param comm \code{\link{community}} object
#'
#' @return integer vector with expected species number when n neighbouring individuals
#'          are sampled
#'
#' @examples
#' sim_com1 <- Sim.Thomas.Community(100, 1000)
#' sac1 <- SAC(sim_com1)
#' rare_curve1 <- rare_curve(table(sim_com1$census$Species))
#'
#' plot(sac1, type = "l", col = 2, xlab = "Sample size",
#'       ylab = "Expected species richness")
#' lines(1:length(rare_curve1), rare_curve1, col = 1)
#' legend("bottomrigh",c("Rarefaction curve","Species accumulation curve"),
#'        col = 1:2, lwd = 2)
accum_curve <- function(comm)
{
   if (class(comm) != "community")
      stop("accum_curve requires a community object as input. See ?community.")

   SAC <- sSAC1_C(comm$census$x, comm$census$y, as.integer(comm$census$species))
   return(SAC)
}



# # -----------------------------------------------------------
# # Function for the spatial SAC written by Xiao Xiao, Dan McGlinn and Nick Gotelli
# # This function computes the sSAC from one random individual in the community
#
# near_neigh_ind = function(data){
#    # The input data has three columns: x, y, and species ID for each individual.
#    data = data[sample(1:dim(data)[1]), ]
#    focal_row = sample(dim(data)[1], 1)
#    # Compute Euclidean distances
#    x_diff = data[, 1] - as.numeric(data[focal_row, 1])
#    y_diff = data[, 2] - as.numeric(data[focal_row, 2])
#    dist_row = sqrt(x_diff^2 + y_diff^2)
#    data_order = data[order(dist_row), ]
#    S = c()
#    #vec_list = lapply(1:dim(data_order)[1], seq)
#    #lapply(vec_list, length(unique(data_order[vec_list, 3])))
#    for (i in 1:dim(data_order)[1]){
#       sp_id_list = data_order[1:i, 3]
#       i_rich = length(unique(sp_id_list))
#       S = c(S, i_rich)
#    }
#    n = 1:dim(data_order)[1]
#    return(list(S = S, n = n))
# }
#
# # -----------------------------------------------------------
# # Average sSAC starting from n (=nsamples) random individuals in the local community
# sSAC.avg <- function(data, nsamples=20)
# {
#    n <- nrow(data)
#
#    sac.mat <- matrix(NA,nrow=nsamples,ncol=n)
#    for (i in 1:nsamples)
#       sac.mat[i,] <- near_neigh_ind(data)$S
#    return(colMeans(sac.mat))
# }
#



