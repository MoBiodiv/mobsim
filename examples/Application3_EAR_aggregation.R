# R code to reproduce Fig. S2 in the article
# mobsim: An R package for the simulation and measurement of
# biodiversity across spatial scales

# Felix May
# October 17 2017

library(mobsim)

# set reference parameters
S0 <- 200
N0 <- 10000
cv0 <- 1

# parameters for scenarios with higher aggregation
sigma_high <- 0.05
sigma_low <- 0.01

# sampling areas for the evaluation of the endemics-area relationship (EAR)
p_area <- c(0.01,0.05,0.1,0.2,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0)

# number of replicate simulations (1000 in the article)
nsim <- 100

# define matrices to store simulation output for the EAR
ear_random_mat <- matrix(NA, nrow = nsim, ncol = length(p_area))
ear_sigma_high_mat <- matrix(NA, nrow = nsim, ncol = length(p_area))
ear_sigma_low_mat <- matrix(NA, nrow = nsim, ncol = length(p_area))
ear_clust1_mat <- matrix(NA, nrow = nsim, ncol = length(p_area))

# run replicate simulations for all scenarios
for (i in 1:nsim){

   # random
   com_random <- sim_poisson_community(s_pool = S0, n_sim = N0, sad_coef = list(meanlog = 3, sdlog = 1))
   divar1 <- divar(com_random, prop_area = p_area)
   ear_random_mat[i, ] <- divar1$m_endemics

   # sigma high
   com_sigma_high <- sim_thomas_community(s_pool = S0, n_sim = N0, sigma = sigma_high,
                                          sad_coef = list(meanlog = 3, sdlog = 1))
   divar1 <- divar(com_sigma_high, prop_area = p_area)
   ear_sigma_high_mat[i, ] <- divar1$m_endemics

   # sigma low
   com_sigma_low <- sim_thomas_community(s_pool = S0, n_sim = N0, sigma = sigma_low,
                                         sad_coef = list(meanlog = 3, sdlog = 1))
   divar1 <- divar(com_sigma_low, prop_area = p_area)
   ear_sigma_low_mat[i, ] <- divar1$m_endemics

   # one cluster per species
   com_clust1 <- sim_thomas_community(s_pool = S0, n_sim = N0, sigma = sigma_high,
                                      mother_points = 1,
                                      sad_coef = list(meanlog = 3, sdlog = 1))
   divar1 <- divar(com_clust1, prop_area = p_area)
   ear_clust1_mat[i, ] <- divar1$m_endemics
}

# calculate means and 95% confidence intervals
ear_random_mean <- colMeans(ear_random_mat, na.rm = T)
ear_random_ci   <- apply(ear_random_mat, MARGIN = 2, FUN = "quantile",
                         probs = c(0.025, 0.975), na.rm = T)
ear_sigma_high_mean <- colMeans(ear_sigma_high_mat, na.rm = T)
ear_sigma_high_ci   <- apply(ear_sigma_high_mat, MARGIN = 2, FUN = "quantile",
                         probs = c(0.025, 0.975), na.rm = T)

ear_sigma_low_mean <- colMeans(ear_sigma_low_mat, na.rm = T)
ear_sigma_low_ci   <- apply(ear_sigma_low_mat, MARGIN = 2, FUN = "quantile",
                            probs = c(0.025, 0.975), na.rm = T)

ear_clust1_mean <- colMeans(ear_clust1_mat, na.rm = T)
ear_clust1_ci   <- apply(ear_clust1_mat, MARGIN = 2, FUN = "quantile",
                                  probs = c(0.025, 0.975), na.rm = T)


# Create plot

png("FigS2_EAR.png", width = 5, height = 5, units = "in", res = 200)
par(las = 1, font.main = 1, cex.lab = 1.2, cex.main = 1.4)

plot(p_area, ear_random_mean, type = "n", las = 1,
     xlab = "Proportion of area lost", ylab = "No. of species lost",
     main = "Endemics area relationships")
polygon(c(p_area, rev(p_area)), c(ear_random_ci[1,], rev(ear_random_ci[2,])),
        col = adjustcolor(1, alpha.f = 0.1), border = NA)
lines(p_area, ear_random_mean)

polygon(c(p_area, rev(p_area)), c(ear_sigma_high_ci[1,], rev(ear_sigma_high_ci[2,])),
        col = adjustcolor(2, alpha.f = 0.1), border = NA)
lines(p_area, ear_sigma_high_mean, col = 2)

polygon(c(p_area, rev(p_area)), c(ear_sigma_low_ci[1,], rev(ear_sigma_low_ci[2,])),
        col = adjustcolor(3, alpha.f = 0.1), border = NA)
lines(p_area, ear_sigma_low_mean, col = 3)


polygon(c(p_area, rev(p_area)),
        c(ear_clust1_ci[1,], rev(ear_clust1_ci[2,])),
        col = adjustcolor(4,alpha.f = 0.1), border = NA)
lines(p_area, ear_clust1_mean, col = 4)

legend("topleft",legend = c("Random","Large clusters",
                             "Small clusters","One large cluster"),
       title = "Species distributions",
       lwd = 2, col = 1:4)

dev.off()
