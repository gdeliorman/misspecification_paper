library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(Matrix)

# Load data set containing the computed ICAs for all scenarios.
results_ICA_tbl = readRDS("results_ICA_tbl.rds")
copula_id_tbl = readRDS("copula_id_tbl.rds")

# We compute the average ICA for each scenario; this value is used as reference
# value. Ideally, the ICAs lie closely to this reference value. In addition, we
# extract the ICA under the multivariate normal scenario as a second reference.
reference_ICA_tbl = results_ICA_tbl %>%
  group_by(id, data_set, assumptions) %>%
  summarize(mean_ICA = mean(ICA))
# Extract id of the set of copula families that corresponds to the multivariate
# normal distribution.
mvn_id = copula_id_tbl %>%
  filter(c23 == "gaussian",
         c13_2 == "gaussian",
         c24_3 == "gaussian",
         c14_23 == "gaussian") %>%
  `$`("copula_id")
mvn_ICA_tbl = results_ICA_tbl %>%
  filter(copula_id == mvn_id) %>%
  rename(mvn_ICA = ICA) %>%
  select(-copula_id)
# Add reference ICA to the original results tibble.
results_ICA_tbl = results_ICA_tbl %>%
  left_join(reference_ICA_tbl, by = c("id", "data_set", "assumptions")) %>%
  left_join(mvn_ICA_tbl, by = c("id", "assumptions", "data_set"))

# Data Visualization ------------------------------------------------------
results_ICA_tbl %>%
  ggplot(aes(x = mean_ICA, y = ICA)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  facet_grid(data_set~assumptions)

results_ICA_tbl %>%
  ggplot(aes(x = mvn_ICA, y = ICA)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  facet_grid(data_set~assumptions)


# Reliability Analysis ----------------------------------------------------

# Fit mixed model for each data set-assumptions combination.
results_ICA_tbl = results_ICA_tbl %>%
  group_by(data_set, assumptions) %>%
  summarize(lme_fit = list(lmer(
    ICA ~ 1 + (1 | id) + (1 | copula_id), data = pick(everything())
  )))

# Extract Reliability coefficients.
reliability_coef = function(lme_fit) {
  coefficients = as.data.frame(VarCorr(lme_fit))
  copula_id_var = coefficients$vcov[coefficients$grp == "copula_id"]
  id_var = coefficients$vcov[coefficients$grp == "id"]
  residual_var = coefficients$vcov[coefficients$grp == "Residual"]
  return(id_var / (residual_var + id_var + copula_id_var))
}
results_ICA_tbl = results_ICA_tbl %>%
  mutate(reliability = purrr::map_dbl(.x = lme_fit, .f = reliability_coef))

