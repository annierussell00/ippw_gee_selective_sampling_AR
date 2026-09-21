# Simulation B: Data generation
#
# Generates a synthetic PsA target population with drug-survival and
# response trajectories calibrated to the SEQUENCE study.
#
# IMPORTANT — data access:
#   This script requires two restricted inputs that are NOT included in
#   this repository and cannot be shared publicly:
#     - data/RNHRD Bio PsA data 2025 (no ID).xlsx  (target-population line
#       distribution)
#     - data/Sample_data.csv (SEQUENCE cohort; used here to draw real baseline dates 
#                             into the synthetic population)
#   It is included for provenance/transparency (to document exactly how
#   Simulation_B_data.csv was generated), not as something external readers
#   can execute themselves. 
#
#
# Output:
#   data/Simulation_B_data.csv
#
# Reproducibility note: this script has two separate seeding points
# (set.seed(1) before baseline-date assignment, set.seed(24) before PsARC/
# survival generation). Run top-to-bottom in this order to reproduce
# Simulation_B_data.csv exactly.


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
library(ggh4x)

################################################################################
#### Load & Prepare Data ####


target_pop_data <- read_excel("Data/RNHRD Bio PsA data 2025 (no ID).xlsx", sheet = 2) # Target population data

# Distribution of total exposures in target population

total_lines_target_pop_data <- target_pop_data %>%
  dplyr::select(CG_ID,"Max line")%>%
  rename(Patient_ID = "CG_ID",
         Total_lines = "Max line") %>%
  filter(Total_lines <= 6) %>%
  count(Total_lines) %>%
  mutate(prop = n/sum(n))

# Sample data

samp_data <- read.csv("Sample_data.csv")

################################################################################
#### Simulation Parameters 
max_line <- c("1", "2", "3", "4", "5", "6")
max_line_num <- 6
N_pop <- 10000


################################################################################
#### Generate Synthetic Population 

n_obs_line_population <- data.frame(total_lines = 1:6, n_obs_pop = round(N_pop*(total_lines_target_pop_data[1:6,]$prop),0))

population_data <- n_obs_line_population %>%
  uncount(n_obs_pop) %>%  
  mutate(patient_id = row_number(),
         total_lines = as.factor(total_lines))%>%
  dplyr::select(patient_id, total_lines)

# Assign each synthetic patient a baseline date (line-1 initiation date),
# drawn from real SEQUENCE baseline dates for patients with the same total
# number of treatment lines, so synthetic trajectories align with real
# calendar-time patterns in the SEQUENCE cohort.

baseline_samp_data <- samp_data[
  samp_data$total_lines %in% max_line,
  c(
    "patient_id",
    "total_lines",
    "date_at_line",
    "line"
  )
]

baseline_samp_data <- samp_data %>%
  filter(line == 1) %>%
  rename(date_baseline = date_at_line) %>%
  dplyr::select(-line)

baseline_dates_lookup <- baseline_samp_data %>%
  filter(total_lines %in% max_line) %>%
  group_by(total_lines) %>%
  summarise(baseline_dates = list(date_baseline),
            .groups = "drop")

set.seed(1)
population_data <- population_data %>%
  left_join(baseline_dates_lookup, by = c("total_lines" = "total_lines")) %>%
  rowwise() %>%
  mutate(
    date_baseline = sample(baseline_dates, size = 1)[[1]],
    D2T = case_when(
      total_lines %in% c("1") ~ "easy",
      total_lines %in% c("2", "3", "4") ~ "medium",
      total_lines %in% c("5", "6") ~ "difficult"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-baseline_dates)

final_date <- as.Date("2023-12-01")

# Total time in the study 

population_data <- population_data %>%
  mutate(total_time = interval(date_baseline, final_date) %/% months(1))

samp_data <- samp_data %>%
  mutate(total_time = interval(date_baseline, final_date) %/% months(1))

################################################################################
#### Parameters for drug survival generation 

alpha_beta <- 2.88
beta_beta  <- 1.03

beta_exp <- 0.2
M_D2T_medium <- 1
M_D2T_difficult <- 1.2

min_ds <- 3
max_ds <- 120

################################################################################
#### Drug Survival Generating Function

generate_drug_survival_lines <- function(date_start,
                                         total_lines,
                                         D2T,
                                         final_date,
                                         alpha_beta,
                                         beta_beta,
                                         beta_exp,
                                         M_D2T_medium,
                                         M_D2T_difficult,
                                         min_ds,
                                         max_ds) {
  
  # Compute total study time in months
  total_time <- interval(date_start, final_date) %/% months(1)
  
  # Choose M_D based on D2T
  M_D <- ifelse(D2T == "medium", M_D2T_medium, M_D2T_difficult)
  
  times <- numeric(total_lines)
  
  if (total_lines == 1) {
    
    times <- rbeta(1, shape1 = alpha_beta, shape2 = beta_beta) * total_time
    
  } else {
    
    # Decay per line
    h <- exp(beta_exp * (0:(total_lines-1)))
    
    # Baseline lambda
    lambda0 <- sum(1/h) / (0.8 * total_time)
    
    # Lambda per line
    lambda_lines <- lambda0 * M_D * h
    
    for (t in 1:total_lines) {
      repeat {
        x <- rexp(1, rate = lambda_lines[t])
        if (x >= min_ds && x <= max_ds) break
      }
      times[t] <- x
    }
  }
  
  if (sum(times) > total_time) {
    times <- times * (total_time / sum(times))
  }
  
  times <- pmin(pmax(times, min_ds), max_ds)
  
  return(times)
}

################################################################################
#### PsARC Generating Function 

generate_psarc_lines <- function(total_lines, drug_survival, total_time) {
  psarc <- numeric(total_lines)
  for (k in 1:total_lines) {
    prob <- plogis(-0.48 + 0.071 * drug_survival[k])
    psarc[k] <- rbinom(1, 1, prob)
  }
  return(psarc)
}

################################################################################
#### Generation of survival and PsARC 

generate_synthetic_ds <- function() {
  pop_data <- population_data %>%
    mutate(total_lines_num = as.integer(as.character(total_lines))) %>%
    rowwise() %>%
    mutate(
      survival_list = list(
        generate_drug_survival_lines(
          date_start = as.Date(date_baseline),
          final_date = final_date,
          total_lines = total_lines_num,
          D2T = D2T,
          alpha_beta = alpha_beta,
          beta_beta = beta_beta,
          beta_exp = beta_exp,
          M_D2T_medium = M_D2T_medium,
          M_D2T_difficult = M_D2T_difficult,
          min_ds = min_ds,
          max_ds = max_ds
        )
      ),
      psarc_list = list(
        generate_psarc_lines(
          total_lines = total_lines_num,
          drug_survival = survival_list
        )
      )
    ) %>%
    ungroup()
  
  drug_psarc_expand <- function(lst, prefix) {
    map_dfc(seq_len(max_line_num), function(j) {
      tibble(!!paste0(prefix, "_", j) :=
               map_dbl(lst, function(x) {
                 if (length(x) >= j) {
                   x[j]
                 } else {
                   NA_real_
                 }
               }))
    })
  }
  
  pop_data <- bind_cols(
    pop_data,
    drug_psarc_expand(pop_data$survival_list, "drug_survival"),
    drug_psarc_expand(pop_data$psarc_list, "psarc")
  ) %>%
    dplyr::select(-survival_list, -psarc_list)
  
  pop_data_long <- pop_data %>%
    rename_with(~ trimws(.x)) %>%
    pivot_longer(
      cols = matches("^(psarc|drug_survival).*_\\d+$", ignore.case = TRUE),
      names_to = c(".value", "line"),
      names_pattern = "^(psarc|drug_survival).*_(\\d+)$",
      values_drop_na = FALSE
    ) %>%
    mutate(
      patient_id = as.character(patient_id),
      scenario = "Synthetic"
    ) %>%
    filter(!is.na(psarc) | !is.na(drug_survival)) %>%
    dplyr::select(-total_lines_num) %>%
    arrange(patient_id, line)
  
  return(pop_data_long)
}

set.seed(24)
target_population <-  generate_synthetic_ds()
target_population <- target_population%>%
  mutate(
    line = as.numeric(as.character(line)),
    line_cat = factor(cut(line, c(0, 1, 3, 10),
                          labels = c("First", "2-3", "4+")))
  )
target_population$line_cat <- relevel(target_population$line_cat, ref = "2-3")

target_population_patients <- target_population %>%
  distinct(patient_id, total_lines)

write.csv(target_population, "Simulation_B_data.csv", row.names = FALSE)
