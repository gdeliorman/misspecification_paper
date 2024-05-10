library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(Matrix)

# Load data set containing the computed ICAs for all scenarios.
results_ICA_tbl = readRDS("results_ICA_tbl.rds")

# We compute the average ICA for each scenario; this value is used as reference
# value. Ideally, the ICAs lie closely to this reference value. 
reference_ICA_tbl = results_ICA_tbl %>%
  group_by(id, data_set, assumptions) %>%
  summarize(mean_ICA = mean(ICA))
# Add reference ICA to the original results tibble.
results_ICA_tbl = results_ICA_tbl %>%
  left_join(reference_ICA_tbl, by = c("id", "data_set", "assumptions"))

# Data Visualization ------------------------------------------------------
results_ICA_tbl %>%
  ggplot(aes(x = mean_ICA, y = ICA)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  facet_grid(data_set~assumptions)


# Reliability Analysis ----------------------------------------------------

# Fit mixed model for each data set-assumptions combination.
results_ICA_tbl %>%
  group_by(data_set, assumptions) %>%
  summarize(lme_fit = list(lmer(
    ICA ~ 1 + (1 | id) + (1 | copula_id), data = pick(everything())
  )))

# Extract Reliability coefficients.

