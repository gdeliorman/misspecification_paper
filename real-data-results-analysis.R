library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(Matrix)

# Size parameter for saving plots to disk.
single_width = 9
double_width = 14
single_height = 8.2
double_height = 12.8
res = 600

# Set the theme for all plots.
theme_set(
  theme_bw()
)

# Load data set containing the computed ICAs for all scenarios.
results_ICA_tbl = readRDS("results_ICA_tbl.rds")
# Data set containing the parametric copulas corresponding to copula_id.
copula_id_tbl = readRDS("copula_id_tbl.rds")
# Data set containing the assumptions corresponding to each id.
id_tbl = readRDS("id_tbl.rds") %>%
  select(assumptions, id)

# Add assumption variable to results_ICA_tbl.
results_ICA_tbl = results_ICA_tbl %>%
  left_join(id_tbl, by = "id")

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
  mutate(
    assumptions = forcats::fct_recode(
      assumptions,
      "-" = "no",
      "PA" = "positive associations",
      "PA + CI" = "positive associations and conditional independence"
    )
  ) %>%
  ggplot(aes(x = mean_ICA, y = ICA)) +
  geom_point(alpha = 0.01, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(c(-0.1, 1)) +
  xlab(latex2exp::TeX("$\\bar{ICA}$")) +
  ylim(c(-0.1, 1)) +
  facet_grid(data_set ~ assumptions)
ggsave(filename = "figures/web-appendices/vine-copula-mean-reference.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_ICA_tbl %>%
  mutate(
    assumptions = forcats::fct_recode(
      assumptions,
      "-" = "no",
      "PA" = "positive associations",
      "PA + CI" = "positive associations and conditional independence"
    )
  ) %>%
  ggplot(aes(x = mvn_ICA, y = ICA)) +
  geom_point(alpha = 0.01, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(c(-0.1, 1)) +
  xlab(latex2exp::TeX("$ICA_{MVN}$")) +
  ylim(c(-0.1, 1)) +
  facet_grid(data_set~assumptions)
ggsave(filename = "figures/main-text/vine-copula-mvn-reference.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

# Reliability Analysis ----------------------------------------------------

# Fit mixed model for each data set-assumptions combination.
results_ICA_tbl = results_ICA_tbl %>%
  group_by(data_set, assumptions) %>%
  summarize(lme_fit = list(lmer(
    ICA ~ 1 + (1 | id) + (1 | copula_id), data = pick(everything())
  )))
# `data_set = "ARMD"` and `assumptions = "positive associations"` fails to
# converge.

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

# Print results to .txt file
sink(file = "tables/web-appendices/reliability-estimates.txt")
cat("Estimated reliability coefficients for each scenario.\n\n")
print(
  results_ICA_tbl %>%
    select(-lme_fit) %>%
    mutate(assumptions = forcats::fct_recode(assumptions, "no assumptions" = "no")) %>%
    pivot_wider(names_from = "assumptions", id_cols = "data_set", values_from = "reliability")
)
sink()

