# R code to reproduce Fig. 2 and Fig. S1 in the article
# mobsim: An R package for the simulation and measurement of biodiversity
# across spatial scales
# R code to reproduce Figs. 2 and S1 in the article
# mobsim: An R package for the simulation and measurement of
# biodiversity across spatial scales

# Felix May
# October 17 2017

# This script reproduces Fig. 2 in the main text.
# To reproduce Fig. S1 in lines 45 and 63 the argument
# fix_s_sim has to be set to fix_s_sim = TRUE


library(mobsim)

# set reference parameters

S0 <- 100  # no. of species
N0 <- 2000 # no. of individuals

nsim <- 100 # no. of replicate simulations (1000 in the article)

# define matrices to store simulation output for species rarefaction (src) and
# species accumulation (sac) curves

# reference
src_ref_mat <- matrix(NA, nrow = nsim, ncol = N0)
sac_ref_mat <- matrix(NA, nrow = nsim, ncol = N0)

# lower N
src_low_n_mat <- matrix(NA, nrow = nsim, ncol = N0/2)
sac_low_n_mat <- matrix(NA, nrow = nsim, ncol = N0/2)

# lower evenness
src_uneven_mat <- matrix(NA, nrow = nsim, ncol = N0)
sac_uneven_mat <- matrix(NA, nrow = nsim, ncol = N0)

# higher aggregation
src_agg_mat <- matrix(NA, nrow = nsim, ncol = N0)
sac_agg_mat <- matrix(NA, nrow = nsim, ncol = N0)

# run replicate simulations
for (i in 1:nsim){

   # simulate reference community
   sad_ref <- sim_sad(s_pool = S0, n_sim = N0, sad_coef = list(meanlog = 3, sdlog = 1),
                      fix_s_sim = FALSE)
   com_ref <- sim_poisson_coords(sad_ref)

   spec_curves_ref <- spec_sample_curve(com_ref)
   src_ref_mat[i, 1:nrow(spec_curves_ref)] <- spec_curves_ref$spec_rarefied
   sac_ref_mat[i, 1:nrow(spec_curves_ref)] <- spec_curves_ref$spec_accum

   #simulate community with lower number of individuals
   com_low_n <- com_ref
   com_low_n$census <- com_low_n$census[sample(1:N0, N0/2), ]

   spec_curves_low_n <- spec_sample_curve(com_low_n)
   src_low_n_mat[i, 1:nrow(spec_curves_low_n)] <- spec_curves_low_n$spec_rarefied
   sac_low_n_mat[i, 1:nrow(spec_curves_low_n)] <- spec_curves_low_n$spec_accum

   #simulate more uneven community
   com_uneven <- sim_poisson_community(s_pool = S0, n_sim = N0,
                                       sad_coef = list(meanlog = 3, sdlog = 1.5),
                                       fix_s_sim = FALSE)

   spec_curves_uneven <- spec_sample_curve(com_uneven)
   src_uneven_mat[i, 1:nrow(spec_curves_uneven)] <- spec_curves_uneven$spec_rarefied
   sac_uneven_mat[i, 1:nrow(spec_curves_uneven)] <- spec_curves_uneven$spec_accum

   #simulate more aggregated community
   com_agg <- sim_thomas_coords(sad_ref, sigma = 0.02)
   spec_curves_agg <- spec_sample_curve(com_agg)

   src_agg_mat[i, 1:nrow(spec_curves_agg)] <- spec_curves_agg$spec_rarefied
   sac_agg_mat[i, 1:nrow(spec_curves_agg)] <- spec_curves_agg$spec_accum
}

# calculate means and 95% confidence intervals from rarefaction and aggregation
# curves

# reference
src_ref_mean <- colMeans(src_ref_mat, na.rm = T)
src_ref_ci   <- apply(src_ref_mat, MARGIN = 2, FUN = "quantile",
                      probs = c(0.025, 0.975), na.rm = T)

sac_ref_mean <- colMeans(sac_ref_mat, na.rm = T)
sac_ref_ci   <- apply(sac_ref_mat, MARGIN = 2, FUN = "quantile",
                      probs = c(0.025, 0.975), na.rm = T)

# lower no. of individuals
src_low_n_mean <- colMeans(src_low_n_mat, na.rm = T)
src_low_n_ci   <- apply(src_low_n_mat, MARGIN = 2, FUN = "quantile",
                      probs = c(0.025, 0.975), na.rm = T)

sac_low_n_mean <- colMeans(sac_low_n_mat, na.rm = T)
sac_low_n_ci   <- apply(sac_low_n_mat, MARGIN = 2, FUN = "quantile",
                      probs = c(0.025, 0.975), na.rm = T)

# lower evenness
src_uneven_mean <- colMeans(src_uneven_mat, na.rm = T)
src_uneven_ci   <- apply(src_uneven_mat, MARGIN = 2, FUN = "quantile",
                        probs = c(0.025, 0.975), na.rm = T)

sac_uneven_mean <- colMeans(sac_uneven_mat, na.rm = T)
sac_uneven_ci   <- apply(sac_uneven_mat, MARGIN = 2, FUN = "quantile",
                        probs = c(0.025, 0.975), na.rm = T)
# higher aggregation
src_agg_mean <- colMeans(src_agg_mat, na.rm = T)
src_agg_ci   <- apply(src_agg_mat, MARGIN = 2, FUN = "quantile",
                         probs = c(0.025, 0.975), na.rm = T)

sac_agg_mean <- colMeans(sac_agg_mat, na.rm = T)
sac_agg_ci   <- apply(sac_agg_mat, MARGIN = 2, FUN = "quantile",
                         probs = c(0.025, 0.975), na.rm = T)

# Create plots
# png("Fig2_biodiv_components.png",width = 7.5, height = 5, units = "in", res = 200)

pdf("Figure_2.pdf", width = 7.5, height = 5)
par(mfcol = c(2,3), mar = c(1,4,3,1), las = 1, font.main = 1, oma = c(4,4,0,0),
    cex.lab = 1.2, cex.main = 1.4)

# ------------------------------------------------------------------------------
# Change in no. of individuals

# Rarefaction curve
label_x <- 100
label_y <- 96
label_cex <- 1.5

plot(1:N0, src_ref_mean, type = "l", ylim = c(1, S0),
     main = "Change in no. of individuals", xlab = "", ylab = "")
text(label_x, label_y, "(a)", cex = label_cex)
polygon(c(1:N0, N0:1), c(src_ref_ci[1,], rev(src_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:(N0/2), src_low_n_mean, type = "l", col = "red")
polygon(c(1:(N0/2), (N0/2):1), c(src_low_n_ci[1,], rev(src_low_n_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)

# Accumulation curves
plot(1:N0, sac_ref_mean, type = "l", ylim = c(1, S0),
     main = "",  xlab = "", ylab = "")
text(label_x, label_y, "(b)", cex = label_cex)
polygon(c(1:N0, N0:1), c(sac_ref_ci[1,], rev(sac_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:(N0/2), sac_low_n_mean, type = "l", col = "red")
polygon(c(1:(N0/2), (N0/2):1), c(sac_low_n_ci[1,], rev(sac_low_n_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)

# ------------------------------------------------------------------------------
# Change in evenness

# Rarefaction curve

plot(1:N0, src_ref_mean, type = "l", ylim = c(1, S0),
     main = "Change in evenness", xlab = "", ylab = "")
text(label_x, label_y, "(c)", cex = label_cex)
polygon(c(1:N0, N0:1), c(src_ref_ci[1,], rev(src_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:N0, src_uneven_mean, type = "l", col = "red")
polygon(c(1:N0, N0:1), c(src_uneven_ci[1,], rev(src_uneven_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)

# Accumulation curves
plot(1:N0, sac_ref_mean, type = "l", ylim = c(1, S0),
     main = "",  xlab = "No. of individuals sampled", ylab = "", xpd = NA)
text(label_x, label_y, "(d)", cex = label_cex)
polygon(c(1:N0, N0:1), c(sac_ref_ci[1,], rev(sac_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:N0, sac_uneven_mean, type = "l", col = "red")
polygon(c(1:N0, N0:1), c(sac_uneven_ci[1,], rev(sac_uneven_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)


# ------------------------------------------------------------------------------
# Change in aggregation

# Rarefaction curve

plot(1:N0, src_ref_mean, type = "l", ylim = c(1, S0),
     main = "Change in aggregation", xlab = "", ylab = "")
text(label_x, label_y, "(e)", cex = label_cex)
polygon(c(1:N0, N0:1), c(src_ref_ci[1,], rev(src_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:N0, src_agg_mean, type = "l", col = "red")
polygon(c(1:N0, N0:1), c(src_agg_ci[1,], rev(src_agg_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)

# Accumulation curves
plot(1:N0, sac_ref_mean, type = "l", ylim = c(1, S0),
     main = "", xlab = "", ylab = "")
text(label_x, label_y, "(f)", cex = label_cex)
polygon(c(1:N0, N0:1), c(sac_ref_ci[1,], rev(sac_ref_ci[2,])),
        col = adjustcolor("black",alpha.f = 0.1), border = NA)
lines(1:N0, sac_agg_mean, type = "l", col = "red")
polygon(c(1:N0, N0:1), c(sac_agg_ci[1,], rev(sac_agg_ci[2,])),
        col = adjustcolor("red",alpha.f = 0.1), border = NA)

mtext("No. of expected species", side = 2, line = -1, outer = T, las = 0,
      cex = 0.8, at = 0.45)
mtext("Rarefaction", side = 2, line = 1, outer = T, las = 0, at = 0.75, cex = 1.1)
mtext("Accumulation", side = 2, line = 1, outer = T, las = 0, at = 0.23, cex = 1.1)

dev.off()

