library(here)
library(magrittr)
library(dplyr)

source(here("multivariable_MR_simulation.R"))
source(here("multivariable_MR_analysis.R"))

# Variable A represents education and variable B cognitive ability. The parameters are set up so the B0 causes A1 and A1 causes B1.

# Set number of replicates 
iter <- 10

results_no_baseline_B0 <- replicate(n = iter, {
  sim_data <- multivariable_MR_simulation(n = 30000,B0_A0 = 0,B1_A1 = 1,A0_B1 = 0, B0_A1 = 1,A0_Y = 0,B0_Y = 0.5)$medi %>%
    multivariable_MR_analysis()
})  %>% t %>% apply(., 2, mean) %>% round(2)

results2_w_baseline_B0 <- replicate(n = iter, {
  sim_data <- multivariable_MR_simulation(n = 30000,B0_A0 = 0,B1_A1 = 1,A0_B1 = 0, B0_A1 = 1,A0_Y = 0,B0_Y = 0.5)$medi %>%
    multivariable_MR_analysis()
})  %>% t %>% apply(., 2, mean) %>% round(2)


