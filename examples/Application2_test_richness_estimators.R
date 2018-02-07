# R code to reproduce a subset of Fig. 3 in the article
# mobsim: An R package for the simulation and measurement of
# biodiversity across spatial scales

# Felix May
# November 7 2017

# This script reproduces just a subset of the Fig. 3 in the article
# The full analysis in the article used many scenarios and large
# communities with 1,000,000 individuals.
# The simulations were done on a High Performance Computing Cluster

# This example illustrates the code and the general idea, but can be run
# in much shorter time than the full analysis

#------------------------------------------------------------------------------
# Function definitions

# Function to calculate richness estimator from community object
# for different sampling schemes
s_chao_sim <- function(comm, n_samples, area_total){

   require(mobsim)
   require(vegan)
   sample1 <- sample_quadrats(comm,
                              n_quadrats = n_samples,
                              quadrat_area = area_total/n_samples,
                              avoid_overlap = T, plot = F)
   sample_sad1 <- colSums(sample1$spec_dat)
   s_estimate <- estimateR(sample_sad1)

   s_chao1 = s_estimate["S.chao1"]

   return(s_chao1)
}

# Function for replicated richness estimations
rep_chao_sim <- function(pars, comm, nrep = 100){

   s_dat <- replicate(nrep, s_chao_sim(comm,
                                       n_samples = pars["n_quadrats"],
                                       area_total = pars["area_sample"])
   )
   return(s_dat)
}

# Function to calculate mean and 95% confidence interval from data vector
mean_ci95 <- function(x){
   ci <- quantile(x, prob = c(0.025, 0.975))
   out <- c(ci[1], mean(x), ci[2])
   names(out) <- c("ci_low","mean","ci_up")
   return(out)
}

#------------------------------------------------------------------------------
# Start of the analysis
library(mobsim)
library(vegan)
library(parallel)

sad1 <- sim_sad(s_pool = 1000, n_sim = 1e5, fix_s_sim = T)
# In the article we used s_pool = 1,000 and n_sim = 1e6

s_true <- length(sad1)

# Simulate scenario with one large cluster (Fig. 3d)
sim1 <- sim_thomas_coords(sad1, sigma = 0.05, mother_points = 1)

# Define sampling schemes
n_quadrats <- c(1, 10, 100) # n_quadrats <- c(1, 2, 5, 10, 20, 50, 100)
                            # in the full analysis

area_sample <- c(1e-4, 1e-3, 1e-2) # area_sample <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
                                   # in the full analysis

sample_pars <- expand.grid(n_quadrats = n_quadrats, area_sample = area_sample)

# Set up local computing cluster for parallel evaluation of the scenarios
n_cores <- detectCores()
print(n_cores)
clust1 <- makeCluster(n_cores)
clusterExport(cl = clust1, "s_chao_sim")

# Parallel execution - in the full analysis we used 1,000 replicates per parameter set
# On my machine (Intel® Core™ i5-4300U CPU @ 1.90GHz × 4, ubuntu 16.04 LTS) the calculations
# shown here took 200 seconds
system.time({
chao_out <- parRapply(cl = clust1, sample_pars,  FUN = rep_chao_sim, comm = sim1, nrep = 100)
})

stopCluster(cl = clust1)

# Calculate mean and confidence intervals and store in data frame
s_chao1_mat <-  matrix(chao_out, nrow = nrow(sample_pars), byrow = T)
s_chao_stats <- t(apply(s_chao1_mat, 1, mean_ci95))
chao_dat <- cbind(sample_pars, s_chao_stats)

#------------------------------------------------------------------------------
# Create plot
library(ggplot2)

chao_dat$n_quadrats <- factor(chao_dat$n_quadrats)

pd <- position_dodge(0.2) # move them .05 to the left and right
ggplot(data = chao_dat, aes(x = area_sample, y = mean, color = n_quadrats)) +
   geom_point(position = pd) +
   scale_x_log10() +
   geom_errorbar(aes(ymin = ci_low, ymax = ci_up),
                 width = 0.2, position = pd) +
   geom_hline(yintercept = s_true, linetype = 2) +
   xlab("Proportion of area sampled") +
   ylab("Estimated species richness") +
   scale_color_discrete(name="No. of quadrats") +
   theme_bw()



