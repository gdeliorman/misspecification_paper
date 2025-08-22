# Load required packages
library(tidyverse)
library(Surrogate)
library(Rsurrogate)
# Load the dataset
data("Schizo")

# Directories to save results in.
dir_tables = "tables/main-text"
dir_figures = "figures/main-text"

# Set seed for reproducibility. Some of the methods will use a Wild bootstrap.
set.seed(123)
# Estimate the proportion of treatment effect explained (PTE).

# It is convenient to split the Schizo data into a data set for the control and a data set for
# the active treatment groups.
Schizo_control = Schizo %>%
  filter(Treat == -1) %>%
  # Remove missing observations
  drop_na(BPRS, PANSS)
Schizo_active = Schizo %>%
  filter(Treat == 1) %>%
  # Remove missing observations
  drop_na(BPRS, PANSS)

# The PTE will be estimated using three different statistical methods: (i)
# robust, (ii) model, and (iii) Freedman. The robust approach is a
# non-parametric approach where nuisance functions are estimated through kernel
# regression. The model and Freedman approaches are model-based (i.e., they rely
# on parametric models); they are based on linear regression. The model approach
# uses a linear regression model with an interaction between treatment and the
# surrogate, the Freedman approach does not.
PTE_estimate_robust = R.s.estimate(
  sone = Schizo_active$BPRS,
  szero = Schizo_control$BPRS,
  yone = Schizo_active$PANSS,
  yzero = Schizo_control$PANSS,
  var = TRUE,
  conf.int = TRUE,
  type = "robust"
)

PTE_estimate_model = R.s.estimate(
  sone = Schizo_active$BPRS,
  szero = Schizo_control$BPRS,
  yone = Schizo_active$PANSS,
  yzero = Schizo_control$PANSS,
  var = TRUE,
  conf.int = TRUE,
  type = "model"
)

PTE_estimate_freedman = R.s.estimate(
  sone = Schizo_active$BPRS,
  szero = Schizo_control$BPRS,
  yone = Schizo_active$PANSS,
  yzero = Schizo_control$PANSS,
  var = TRUE,
  conf.int = TRUE,
  type = "freedman"
)

# Summarize the results and save them in a text file.
results <- tibble(
  Method = c("Robust", "Model", "Freedman"),
  PTE = c(PTE_estimate_robust$R.s, PTE_estimate_model$R.s, PTE_estimate_freedman$R.s),
  Lower_CI = c(PTE_estimate_robust$conf.int.quantile.R.s[1], PTE_estimate_model$conf.int.quantile.R.s[1], PTE_estimate_freedman$conf.int.quantile.R.s[1]),
  Upper_CI = c(PTE_estimate_robust$conf.int.quantile.R.s[2], PTE_estimate_model$conf.int.quantile.R.s[2], PTE_estimate_freedman$conf.int.quantile.R.s[2])
)
# Print the results
sink(file.path(dir_tables, "PTE_results.txt"))
print(results)
sink()
