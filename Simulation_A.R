################################################################################
#### Simulation A: Effect of target-distribution misspecification ####
################################################################################

# The simulation evaluates inverse probability of participation weighting (IPPW).
# IPPW aims to recover a population effect estimate from a selective sample. The
# sample is selected based on effect modifier Z. This simulation shows the 
# sensitivity of IPPW to noise in the assumed target population, under various 
# levels of effect modification.

################################################################################
#### Notation ####

# A: exposure
# Z: effect modifier
# Y: outcome

# beta_A: effect of A on Y
# beta_Z: effect of Z on Y
# beta_AZ: effect modification between A and Z on Y

# shift: level of noise in measurement of population distribution of Z

################################################################################
#### Load Libraries ####

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

################################################################################
#### Data generating process ####

# Population

simulate_population <- function(N, pZ,
                                beta0, betaA, betaZ, betaAZ) {
  
  Z <- rbinom(N, 1, pZ)
  A <- rbinom(N, 1, 0.5)
  
  lp <- beta0 +
    betaA * A +
    betaZ * Z +
    betaAZ * A * Z
  
  Y <- rbinom(N, 1, plogis(lp))
  
  data.frame(A = A, Z = Z, Y = Y)
}

# Biased sampling on Z - function depending on pZ_sample

sample_data <- function(dat, pZ_sample) {
  
  w <- ifelse(dat$Z == 1,
              pZ_sample,
              1 - pZ_sample)
  
  dplyr::slice_sample(dat,
                      n = floor(nrow(dat) * 0.2), # sample 20% of the population
                      weight_by = w)
}

################################################################################
#### Generating simulation ####

sim_run <- function(betaAZ, shift, seed) {
  
  set.seed(seed)
  
  pop <- simulate_population(
    N = 10000,
    pZ = 0.5, # proportion of Z=1 in population
    beta0 = 0.5,
    betaA = -0.2,
    betaZ = -0.1,
    betaAZ = betaAZ
  )
  
  # proportion of Z in sample misaligned with population
  samp <- sample_data(pop, pZ_sample = 0.3) 
  
  # Assumed target population proportion of Z = 1
  # true population proportion + measurement noise
  pZ_true <- rnorm(1, mean = 0.5, sd = shift)
  
  # truncation so not to have little/no values of either level of Z
  pZ_true <- min(max(pZ_true, 0.1), 0.9) 
  
  # computing participation weights
  w <- ifelse(
    samp$Z == 1,
    pZ_true / 0.3,
    (1 - pZ_true) / (1 - 0.3)
  )
  
  # truncating weights to avoid extreme weights
  w <- pmin(pmax(w, 0.1), 5)
  
  m_true <- glm(
    Y ~ A,
    data = pop,
    family = binomial()
  )
  
  m_naive <- glm(
    Y ~ A,
    data = samp,
    family = binomial()
  )
  
  m_w <- glm(
    Y ~ A,
    data = samp,
    weights = w,
    family = binomial()
  )
  
  extract_ci <- function(model){
    
    cf <- summary(model)$coefficients
    
    if(!"A" %in% rownames(cf)){
      return(c(
        OR = NA_real_,
        lo = NA_real_,
        hi = NA_real_
      ))
    }
    
    est <- cf["A","Estimate"]
    se  <- cf["A","Std. Error"]
    
    c(
      OR = exp(est),
      lo = exp(est - 1.96 * se),
      hi = exp(est + 1.96 * se)
    )
  }
  
  true_ci <- extract_ci(m_true)
  naive_ci <- extract_ci(m_naive)
  weighted_ci <- extract_ci(m_w)
  
  data.frame(
    betaAZ = betaAZ,
    shift = shift,
    
    true_OR = true_ci["OR"],
    true_lo = true_ci["lo"],
    true_hi = true_ci["hi"],
    
    naive_OR = naive_ci["OR"],
    naive_lo = naive_ci["lo"],
    naive_hi = naive_ci["hi"],
    
    weighted_OR = weighted_ci["OR"],
    weighted_lo = weighted_ci["lo"],
    weighted_hi = weighted_ci["hi"]
  )
}

################################################################################
#### Grid ####

grid <- expand.grid(
  betaAZ = c(0, -0.5, -1),
  shift = c(0, 0.2, 0.4, 0.6),
  rep = 1:50
)

results <- purrr::pmap_dfr(
  grid,
  function(betaAZ, shift, rep) {
    sim_run(betaAZ, shift, seed = rep)
  }
)

################################################################################
#### Results ####

summary_results <- results %>%
  group_by(betaAZ, shift) %>%
  summarise(
    true_mean = mean(true_OR),
    naive_mean = mean(naive_OR),
    weighted_mean = mean(weighted_OR),
    
    true_q025 = quantile(true_OR, 0.025),
    true_q975 = quantile(true_OR, 0.975),
    
    naive_q025 = quantile(naive_OR, 0.025),
    naive_q975 = quantile(naive_OR, 0.975),
    
    weighted_q025 = quantile(weighted_OR, 0.025),
    weighted_q975 = quantile(weighted_OR, 0.975),
    
    mean_true_CI_width =
      mean(true_hi - true_lo),
    
    mean_naive_CI_width =
      mean(naive_hi - naive_lo),
    
    mean_weighted_CI_width =
      mean(weighted_hi - weighted_lo),
    
    .groups = "drop"
  )