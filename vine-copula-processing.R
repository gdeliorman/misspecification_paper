# Load required packages
library(tidyverse)
library(Surrogate)
library(Rsurrogate)

# Directories to save results in.
dir_tables = "tables/main-text"
dir_figures = "figures/main-text"

# Set seed for reproducibility.
set.seed(123)

# Load results from the sensitivity analysis.
results_ICA_tbl = readRDS("results_ICA_tbl.rds")
copula_id_tbl = readRDS("copula_id_tbl.rds")
id_tbl = readRDS("id_tbl.rds")

# The above loaded tibbles contain the results for multiple data sets, we only
# need the results for the Schizo data here.
results_ICA_tbl = results_ICA_tbl %>%
  filter(data_set == "Schizo")
# Join the copula id and the id table to the results table. The former two
# contain information about the parametric assumptions made in the sensitivity
# analysis and the assumptions about the unidentifiable parameters.
results_ICA_tbl = results_ICA_tbl %>%
  left_join(copula_id_tbl, by = "copula_id") %>%
  left_join(id_tbl, by = "id")

# Compute the intervals of ignorance for each combination of (i) parametric
# assumptions and (ii) assumptions about the unidentifiable parameters.
intervals_of_ignorance_tbl = results_ICA_tbl %>%
  group_by(copula_id, assumptions) %>%
  summarise(lower_bound = min(ICA), upper_bound = max(ICA))

# Save the intervals of ignorance in a text file. We first joint with the
# copula_id_tbl to ensure that information about the parametric
# assumptions is included.
sink(file.path(dir_tables, "intervals_of_ignorance.txt"))
intervals_of_ignorance_tbl %>%
  left_join(copula_id_tbl, by = "copula_id")
sink()

# Visualize the intervals of ignorance.
intervals_of_ignorance_tbl %>%
  left_join(copula_id_tbl, by = "copula_id") %>%
  # We randomly select a subset of the copula ids to ensure that the plot is
  # readable.
  filter(copula_id %in% sample(unique(intervals_of_ignorance_tbl$copula_id), 30)) %>%
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
  file.path(dir_figures, "intervals_of_ignorance.png"),
  width = 8,
  height = 6
)

