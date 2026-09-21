################################################################################
#### Simulation A: Effect of target-distribution misspecification ####
################################################################################

# Monte Carlo simulation study comparing:
#   - GEE on a simple random sample (SRS)
#   - GEE on an exposure-history-enriched sample (SEQ)
#   - Inverse probability of participation weighting GEE (IPPW-GEE) on the
#     enriched sample
# against the population-level GEE estimand.
#
#
# Inputs (see data/README.md):
#   data/Simulation_B_data.csv  - simulated population data
#   data/Sample_data.csv        - SEQUENCE-like enriched sample (real data;
#                                 not included in this repository)

################################################################################
#### Load Libraries ####

library(tidyr)
library(ggplot2)
library(readxl)
library(purrr)
library(synthpop)
library(combinat)
library(lubridate)
library(patchwork)
library(fitdistrplus)
library(dplyr)
library(openxlsx)
library(forcats)
library(tidyverse)
library(gee) 
library(ggpattern)
library(geepack)
library(furrr)
library(future)
library(ggh4x)

################################################################################
#### Load Data ####

pop_data <- read.csv("Simulation_B_data.csv") # generated population

pop_data_patient_list <- pop_data %>%
  distinct(patient_id, total_lines)

samp_data <- read.csv("Sample_data.csv") # real sample

n_workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = n_workers)

################################################################################
#### Target Estimand ####

# Fits the "gold standard" model on the full simulated population, giving
# the true odds ratios that sample-based estimators are compared against.

pop_gee <- gee(
  psarc ~ line_cat,
  data = pop_data,
  id = patient_id,
  family = binomial,
  corstr = "exchangeable"
)

# Population OR 

true_or_df <- data.frame(est = exp(summary(pop_gee)$coefficients[2:3,1]),
                         L = exp(summary(pop_gee)$coefficients[2:3,1]-
                                   1.96*summary(pop_gee)$coefficients[2:3,4]),
                         U = exp(summary(pop_gee)$coefficients[2:3,1]+
                                   1.96*summary(pop_gee)$coefficients[2:3,4]))

true_or <- data.frame(
  dataset = "Population",
  line = c("First", "4+"),
  OR  = true_or_df$est,
  LCI = true_or_df$L,
  UCI = true_or_df$U
)

################################################################################
#### Distribution of line of therapy in sample v population ####

seq_props <- samp_data %>%
  distinct(patient_id, total_lines) %>%  
  count(total_lines) %>%
  mutate(
    prop = n / sum(n),
    total_lines = as.numeric(as.character(total_lines)),
    source = "SEQUENCE"
  ) %>%
  dplyr::select(-n)


pop_props <- pop_data %>%
  distinct(patient_id, total_lines) %>%
  count(total_lines) %>%
  mutate(
    prop = n / sum(n),
    total_lines = as.numeric(as.character(total_lines)),
    source = "Population"
  ) %>%
  dplyr::select(-n)

distribution_of_lines <- bind_rows(
  pop_props,
  seq_props
) %>%
  mutate(total_lines = factor(total_lines, levels = 1:6),
         source = factor(source, levels = c("Population", "SEQUENCE")))

################################################################################
#### Monte Carlo Experiments ####

# For a given sample size, random seed, and enrichment design (SRS vs.
# SEQUENCE-like), draws a sample from the population, fits an unweighted
# GEE, and (for enriched samples) an IPPW-weighted GEE using weights
# derived from the ratio of population to sample line of therapy
# proportions.

run_all_gees <- function(n_sample, seed, enrich_level) {
  
  set.seed(seed)
  
  # Define sampling distributions 
  
  seq_props <- samp_data %>%
    distinct(patient_id, total_lines) %>%
    count(total_lines) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(total_lines)
  
  pop_props <- pop_data %>%
    distinct(patient_id, total_lines) %>%
    count(total_lines) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(total_lines)
  
  seq_props <- seq_props$prop / sum(seq_props$prop)
  pop_props <- pop_props$prop / sum(pop_props$prop)
  
  get_sel_props <- function(level) {
    if (level == "SRS") return(pop_props) # simple random sample
    if (level == "SEQ") return(seq_props) # SEQUENCE-like enrichment
  }
  
  sel_props_vec <- get_sel_props(enrich_level)
  
  n_per_line <- data.frame(
    total_lines = 1:6,
    prop_sel = sel_props_vec
  ) %>%
    mutate(n_to_sample = round(prop_sel * n_sample))
  
  # Draw sample 
  
  sampled_patients <- pop_data_patient_list %>%
    left_join(n_per_line, by = "total_lines") %>%
    group_by(total_lines) %>%
    group_modify(~ slice_sample(.x, n = unique(.x$n_to_sample))) %>%
    ungroup()
  
  sampled_data <- pop_data %>%
    semi_join(sampled_patients, by = "patient_id") %>%
    mutate(
      line = as.numeric(as.character(line)),
      line_cat = factor(cut(line, c(0, 1, 3, 10),
                            labels = c("First", "2-3", "4+")))
    )
  
  sampled_data$line_cat <- relevel(sampled_data$line_cat, ref = "2-3")
  
  # GEE (unweighted)
  
  gee_fit <- gee(
    psarc ~ line_cat,
    data = sampled_data,
    family = binomial,
    id = patient_id,
    corstr = "exchangeable"
  )
  
  OR_gee <- data.frame(
    dataset = paste0("GEE_", enrich_level),
    line = c("First", "4+"),
    OR  = exp(coef(gee_fit)[2:3]),
    SE  = summary(gee_fit)$coefficients[2:3, 4]
  ) %>%
    mutate(
      LCI = exp(log(OR) - 1.96 * SE),
      UCI = exp(log(OR) + 1.96 * SE)
    ) %>%
    dplyr::select(-SE)
  
  ## WGEE for enriched samples
  OR_wgee <- NULL
  
  if (enrich_level != "SRS") {
    
    target_props <- pop_data %>%
      distinct(patient_id, total_lines) %>%
      count(total_lines) %>%
      mutate(prop_pop = n / sum(n))
    
    sel_props <- sampled_data %>%
      distinct(patient_id, total_lines) %>%
      count(total_lines) %>%
      mutate(prop_sel = n/sum(n))
    
    ippw <- target_props %>%
      left_join(sel_props, by = "total_lines") %>%
      mutate(weight = prop_pop / prop_sel)
    
    sampled_w <- sampled_data %>%
      left_join(ippw, by = "total_lines")
    
    gee_w <- geeglm(
      psarc ~ line_cat,
      data = sampled_w,
      family = binomial,
      id = patient_id,
      weights = weight,
      corstr = "exchangeable"
    )
    
    OR_wgee <- data.frame(
      dataset = paste0("WGEE_", enrich_level),
      line = c("First", "4+"),
      OR  = exp(coef(gee_w)[2:3]),
      SE  = summary(gee_w)$coefficients[2:3, 2]
    ) %>%
      mutate(
        LCI = exp(log(OR) - 1.96 * SE),
        UCI = exp(log(OR) + 1.96 * SE)
      ) %>%
      dplyr::select(-SE)
    
    weights_out <- sampled_w %>%
      distinct(patient_id, total_lines, weight)
    
  } else {
    
    OR_wgee <- NULL
    
    weights_out <- sampled_data %>%
      distinct(patient_id, total_lines) %>%
      mutate(weight = 1)
  }
  
  
  # Output 
  
  return(tibble(
    OR_results = bind_rows(OR_gee, OR_wgee),
    weights = list(weights_out),
    sampled_data = list(sampled_data),
    n_sample = n_sample,
    enrich_level = enrich_level,
    seed = seed
  ))
}

sample_sizes <- c(500, 1000, 2000, 3000, 4000, 5000)
seeds <- 1:500
enrich_levels <- c("SRS", "SEQ")

mc_grid <- expand.grid(
  n_sample = sample_sizes,
  seed = seeds,
  enrich_level = enrich_levels
)

mc_results <- furrr::future_pmap_dfr( mc_grid, run_all_gees)
closeAllConnections()

future::plan(future::sequential)
closeAllConnections()

################################################################################
#### Results ####

# Weights

weights_all <- mc_results %>%
  dplyr::select(weights, n_sample, enrich_level, seed) %>%
  unnest(weights)

mean_weights_by_line_design <- weights_all %>%
  group_by(enrich_level, total_lines) %>%
  summarise(
    mean_weight = mean(weight),
    lower_CI = quantile(weight, 0.025),
    upper_CI = quantile(weight, 0.975),
    .groups = "drop"
  )

# Odds ratios

or_all <- mc_results %>%
  tidyr::unnest(OR_results)

or_summary <- or_all %>%
  group_by(dataset, n_sample, line) %>%
  summarise(
    OR = mean(OR),
    LCI = quantile(LCI, 0.025),
    UCI = quantile(UCI, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    true_or %>%
      dplyr::select(line, OR_true = OR),
    by = "line"
  ) %>%
  mutate(
    line = factor(line, levels = c("First", "2-3", "4+")),
    bias = OR - OR_true
  )

or_summary2 <- or_summary %>%
  mutate(
    line = factor(line, levels = c("First", "2-3", "4+")),
    
    method = ifelse(grepl("WGEE", dataset), "WGEE", "GEE"),
    design = case_when(
      grepl("SRS", dataset) ~ "SRS",
      grepl("SEQ", dataset) ~ "SEQ"
    ),
    
    group = factor(
      paste(design, method),
      levels = c(
        "SRS GEE",
        "SEQ GEE", 
        "SEQ WGEE"
      )
    ),
    
    col = case_when(
      line == "First" ~ "1 vs 2-3",
      line == "4+"    ~ "4+ vs 2-3"
    ),
    col = factor(col, levels = c("1 vs 2-3", "4+ vs 2-3"))
  )%>%
  mutate(group = fct_recode(group,
                            "SRS" = "SRS GEE",
                            "SEQ IPPW-GEE" = "SEQ WGEE"))


################################################################################
#### Performance ####

compute_mc_metrics <- function(or_data, true_or_vec) {
  
  or_data %>%
    left_join(
      true_or_vec %>% dplyr::select(line, true_OR = OR),
      by = "line"
    ) %>%
    
    mutate(
      error = OR - true_OR
    ) %>%
    
    group_by(dataset, n_sample, enrich_level, line) %>%
    
    summarise(
      Bias = mean(error, na.rm = TRUE), #bias	systematic error
      RMSE = sqrt(mean(error^2, na.rm = TRUE)), #rmse	total estimation error
      Coverage = mean(LCI <= true_OR & UCI >= true_OR, na.rm = TRUE), #coverage	CI validity
      
      lower = quantile(OR, 0.025, na.rm = TRUE),
      upper = quantile(OR, 0.975, na.rm = TRUE),
      Spread = upper - lower, #spread	empirical uncertainty (2.5–97.5%)
      
      .groups = "drop"
    )
}

mc_metrics <- compute_mc_metrics(or_all, true_or) %>%
  mutate(
    
    estimator = sub("_.*", "", dataset),
    design    = sub(".*_", "", dataset)
    
  )