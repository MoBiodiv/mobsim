#' -----------------------------------------------------------------------------
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
#' Hurlbert, S.H. 1971. The nonconcept of species diversity: a critique and alternative parameters.
#' Ecology 52, 577-586.
#'
#' @export
#'
spec_sample <- function(abund_vec, n)
{
   abund_vec <- abund_vec[abund_vec > 0]
   n_total <- sum(abund_vec)
   ldiv <- lchoose(n_total, n)
   p1 <- exp(lchoose(n_total - abund_vec, n) - ldiv)
   out <- sum(1 - p1)
   names(out) <- n

   return(out)
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
#' @export
#'
rare_curve <- function(abund_vec)
{
   abund_vec <- abund_vec[abund_vec > 0]
   n_total <- sum(abund_vec)
   n_vec <- 1:n_total

   rc <- sapply(n_vec, function(n_sample){spec_sample(abund_vec, n = n_sample)})
   names(rc) <- n_vec
   return(rc)
}

#' -----------------------------------------------------------------------------
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
#'        col = 1:2, lwd = 2)#
#'
#' @export
#'
spec_sample_curve <- function(comm, method = c("accumulation" ,"rarefaction"))
{
   if (class(comm) != "community")
      stop("spec_sample_curve requires a community object as input. See ?community.")

   out_dat <- data.frame(n = 1:nrow(comm$census))

   method <- match.arg(method, several.ok = TRUE)

   if ("accumulation" %in% method){
      out_dat$spec_accum <- sSAC1_C(comm$census$x, comm$census$y, as.integer(comm$census$species))
   }

   if ("rarefaction" %in% method){
      abund <- community_to_sad(comm)
      out_dat$spec_rarefied <- rare_curve(abund)

   }

   class(out_dat) <- c("spec_sample_curve", "data.frame")

   return(out_dat)
}


#' -----------------------------------------------------------------------------
#' @export
plot.spec_sample_curve <- function(curve)
{
   plot(curve[[1]], curve[[2]], type = "n", xlab = "No. of individuals sampled",
        ylab = "Expected no.of species", main = "Species sampling curves")

   if (ncol(curve) == 2){
      if (names(curve)[2] == "spec_accum"){
         legend_text <- c("Accumulation")
         line_col <- "red"
      }
      if (names(curve)[2] == "spec_rarefied"){
         legend_text <- c("Rarefaction")
         line_col <- "blue"
      }
      lines(curve[[1]], curve[[2]], col = line_col)
   }

   if (ncol(curve) == 3){
      legend_text <- c("Accumulation","Rarefaction")
      line_col <- c("red","blue")
      lines(curve[[1]], curve[[2]], col = line_col[1])
      lines(curve[[1]], curve[[3]], col = line_col[2])
   }

   legend("bottomright", legend = legend_text, lty = 1, col = line_col)
}

