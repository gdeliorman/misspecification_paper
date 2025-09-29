# Setup ------------------------------------------------------------

# Load required packages
library(tidyverse)
library(Surrogate)

# Directories to save results in.
dir_tables = "tables/main-text"
dir_figures = "figures/main-text"

# Set seed for reproducibility.
set.seed(123)

# Load results from the sensitivity analysis.
results_ICA_tbl = readRDS("results_ICA_tbl.rds")
copula_id_tbl = readRDS("copula_id_tbl.rds")
id_tbl = readRDS("id_tbl.rds")

# Join the copula id and the id table to the results table. The former two
# contain information about the parametric assumptions made in the sensitivity
# analysis and the assumptions about the unidentifiable parameters.
results_ICA_tbl = results_ICA_tbl %>%
  left_join(copula_id_tbl, by = "copula_id") %>%
  left_join(id_tbl, by = "id")

# Compute the intervals of ignorance for each combination of (i) parametric
# assumptions and (ii) assumptions about the unidentifiable parameters.
intervals_of_ignorance_tbl = results_ICA_tbl %>%
  group_by(copula_id, assumptions, data_set) %>%
  summarise(lower_bound = min(ICA), upper_bound = max(ICA)) %>%
  ungroup()

# Save Results ------------------------------------------------------

## Intervals of Ignorance ------------------------------------------------

# Save the intervals of ignorance in a text file. We first joint with the
# copula_id_tbl to ensure that information about the parametric
# assumptions is included.
sink(file.path(dir_tables, "intervals_of_ignorance-Schizo.txt"))
intervals_of_ignorance_tbl %>%
  filter(data_set == "Schizo") %>%
  left_join(copula_id_tbl, by = "copula_id")
sink()

sink(file.path(dir_tables, "intervals_of_ignorance-ARMD.txt"))
intervals_of_ignorance_tbl %>%
  filter(data_set == "ARMD") %>%
  left_join(copula_id_tbl, by = "copula_id")
sink()

sink(file.path(dir_tables, "intervals_of_ignorance-Ovarian.txt"))
intervals_of_ignorance_tbl %>%
  filter(data_set == "Ovarian") %>%
  left_join(copula_id_tbl, by = "copula_id")
sink()

# Random copula family ids to display in the plots below.
random_copula_ids = sample(unique(intervals_of_ignorance_tbl$copula_id), 30)

plot_tbl = intervals_of_ignorance_tbl %>%
  left_join(copula_id_tbl, by = "copula_id") %>%
  # We randomly select a subset of the copula ids to ensure that the plot is
  # readable.
  filter(copula_id %in% random_copula_ids) %>%
  mutate(
    assumptions = forcats::fct_recode(
      assumptions,
      "-" = "no",
      "PA" = "positive associations",
      "PA + CI" = "positive associations and conditional independence"
    )
  )

# Visualize the intervals of ignorance for Schizo data set.
plot_tbl %>%
  filter(data_set == "Schizo") %>%
  ggplot(aes(
    y = as.factor(copula_id),
    xmin = lower_bound,
    xmax = upper_bound
  )) +
  geom_errorbarh(height = 0.5) +
  scale_y_discrete(name = "Parametric assumptions", breaks = NULL) +
  facet_grid(. ~ assumptions) +
  # Add vertical dashed lines at 0 and 1.
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "red") +
  scale_x_continuous(name = "Estimated interval of ignorance for ICA") +
  theme_bw()
ggsave(
  file.path(dir_figures, "intervals_of_ignorance-Schizo.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

# For ARMD data set.
plot_tbl %>%
  filter(data_set == "ARMD") %>%
  ggplot(aes(
    y = as.factor(copula_id),
    xmin = lower_bound,
    xmax = upper_bound
  )) +
  geom_errorbarh(height = 0.5) +
  scale_y_discrete(name = "Parametric assumptions", breaks = NULL) +
  facet_grid(. ~ assumptions) +
  # Add vertical dashed lines at 0 and 1.
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "red") +
  scale_x_continuous(name = "Estimated interval of ignorance for ICA") +
  theme_bw()
ggsave(
  file.path(dir_figures, "intervals_of_ignorance-ARMD.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

# For Ovarian data set.
plot_tbl %>%
  filter(data_set == "Ovarian") %>%
  ggplot(aes(
    y = as.factor(copula_id),
    xmin = lower_bound,
    xmax = upper_bound
  )) +
  geom_errorbarh(height = 0.5) +
  scale_y_discrete(name = "Parametric assumptions", breaks = NULL) +
  facet_grid(. ~ assumptions) +
  # Add vertical dashed lines at 0 and 1.
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "red") +
  scale_x_continuous(name = "Estimated interval of ignorance for ICA") +
  theme_bw()
ggsave(
  file.path(dir_figures, "intervals_of_ignorance-Ovarian.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

## Distribution of ICA values --------------------------------------------

# Table for abbreviations for the copula family names.
copula_family_abbreviations = tibble::tibble(
  family = c(
    "gaussian",
    "clayton",
    "gumbel",
    "frank"
  ),
  abbreviation = c("Ga", "C", "Gu", "F")
)

# For each data set, we visualize the distribution of the ICA values (obtained
# in the MC sensitivity analysis) under the three different assumptions about
# the unidentifiable parameters (as columns). We have the multivariate normal
# distribution as first row (for reference) and then three randomly selected
# D-vine copulas on the other rows.
set.seed(123)
# Copula id for the multivariate normal distribution.
copula_id_mvn = copula_id_tbl %>%
  filter(c23 == "gaussian",
         c13_2 == "gaussian",
         c24_3 == "gaussian",
         c14_23 == "gaussian") %>%
  `$`("copula_id")
random_copula_ids = sample(
  results_ICA_tbl %>%
    filter(copula_id != copula_id_mvn) %>%
    pull(copula_id) %>%
    unique(),
  3
)

plot_tbl = results_ICA_tbl %>%
  filter(copula_id %in% c(copula_id_mvn, random_copula_ids)) %>%
  mutate(
    copula_id = as.factor(copula_id),
    # Variable that describes the parametric assumptions as the first letter of
    # the copula families in the following order: C23, C13|2, C24|3, C14|23. We
    # use the abbreviations in the copula_family_abbreviations table.
    copula_description = paste0(
      copula_family_abbreviations$abbreviation[match(c23, copula_family_abbreviations$family)],
      copula_family_abbreviations$abbreviation[match(c13_2, copula_family_abbreviations$family)],
      copula_family_abbreviations$abbreviation[match(c24_3, copula_family_abbreviations$family)],
      copula_family_abbreviations$abbreviation[match(c14_23, copula_family_abbreviations$family)]
    ),
    assumptions = forcats::fct_recode(
      assumptions,
      "-" = "no",
      "PA" = "positive associations",
      "PA + CI" = "positive associations and conditional independence"
    )
  ) %>%
  mutate(
    copula_description = factor(
      copula_description,
      levels = c(
        "GaGaGaGa",
        sort(unique(copula_description[copula_id != copula_id_mvn]))
      )
    )
  )

# Plot for the Schizo data set.
plot_tbl %>%
  filter(data_set == "Schizo") %>%
  ggplot(aes(x = ICA)) +
  geom_density() +
  xlab("ICA") +
  ylab("Frequency") +
  facet_grid(copula_description ~ assumptions) +
  theme_bw()
ggsave(
  file.path(dir_figures, "distribution_of_ICA_values-Schizo.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

# Plot for the ARMD data set.
plot_tbl %>%
  filter(data_set == "ARMD") %>%
  ggplot(aes(x = ICA)) +
  geom_density() +
  xlab("ICA") +
  ylab("Frequency") +
  facet_grid(copula_description ~ assumptions) +
  theme_bw()
ggsave(
  file.path(dir_figures, "distribution_of_ICA_values-ARMD.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

# Plot for the Ovarian data set.
plot_tbl %>%
  filter(data_set == "Ovarian") %>%
  ggplot(aes(x = ICA)) +
  geom_density() +
  xlab("ICA") +
  ylab("Frequency") +
  facet_grid(copula_description ~ assumptions) +
  theme_bw()
ggsave(
  file.path(dir_figures, "distribution_of_ICA_values-Ovarian.pdf"),
  width = 8,
  height = 6,
  device = "pdf"
)

