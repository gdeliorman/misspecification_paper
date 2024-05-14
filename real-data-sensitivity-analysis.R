#!/usr/bin/env Rscript
.libPaths()
args = commandArgs(trailingOnly=TRUE)
options(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
Sys.setenv(TZ='Europe/Brussels')
n_cores = as.integer(args[1])
# Print R-version for better reproducibility.
print(version)
# Ensure to the state of packages is up-to-date.
renv::restore()
# Load the required packages.
library(Surrogate)
library(tidyr)
library(dplyr)

# Load data sets from the Surrogate package.
data("ARMD")
data("Ovarian")
data("Schizo")
Schizo = na.omit(Schizo)

# Number of unidentifiable copula parameter sets for each data set and
# additional assumptions.
n_sim = 1e2
# Number of Monte Carlo samples to compute the ICA.
n_prec = 1e4
# Seed for reproducibility.
set.seed(1)

# Estimating Bivariate Normal Margins -------------------------------------

# Fit multivariate normal distribution to the data sets. This boils down to
# estimating the correlation between S and T, and the corresponding means in
# each treatment group.
estimate_mvn = function(surr, true) {
  mean_surr = mean(surr)
  mean_true = mean(true)
  cor_surr_true = cor(surr, true)
  return(
    c(mean_surr, mean_true, cor_surr_true)
  )
}

mvn_est_tbl = tibble(
  data_set = c("ARMD", "Ovarian", "Schizo"),
  mvn_estimates = list(
    list(
      control = estimate_mvn(ARMD$Diff24[ARMD$Treat == -1], ARMD$Diff52[ARMD$Treat == -1]),
      treat = estimate_mvn(ARMD$Diff24[ARMD$Treat == 1], ARMD$Diff52[ARMD$Treat == 1])
    ),
    list(
      control = estimate_mvn(Ovarian$Pfs[Ovarian$Treat == 0], Ovarian$Surv[Ovarian$Treat == 0]),
      treat = estimate_mvn(Ovarian$Pfs[Ovarian$Treat == 1], Ovarian$Surv[Ovarian$Treat == 1])
    ),
    list(control = estimate_mvn(Schizo$CGI[Schizo$Treat == -1], Schizo$PANSS[Schizo$Treat == -1]),
         treat = estimate_mvn(Schizo$CGI[Schizo$Treat == 1], Schizo$PANSS[Schizo$Treat == 1]))
    )
)

# Save information about fitted multivariate normal distributions to a .txt
# file.
sink("tables/web-appendices/mvn-estimates.txt")
print(mvn_est_tbl %>%
        rowwise(data_set) %>%
        summarize(control_mean_S = mvn_estimates$control[1],
                  control_mean_T = mvn_estimates$control[2],
                  control_corr = mvn_estimates$control[3],
                  treat_mean_S = mvn_estimates$control[1],
                  treat_mean_T = mvn_estimates$control[2],
                  treat_corr = mvn_estimates$control[3]))
sink()

# Unidentifiable Copula Families ------------------------------------------

# Copula families under considerations.
copula_fam_vec = c("gaussian", "clayton", "frank", "gumbel")

# Tibble of all possible combinations of unidentifiable copula families.
copula_fam_df = expand_grid(
  c23 = copula_fam_vec,
  c13_2 = copula_fam_vec,
  c24_3 = copula_fam_vec,
  c14_23 = copula_fam_vec
)

# Unidentifiable Copula Parameters ----------------------------------------

## Sample Parameters on Spearman's rho Scale ------------------------------

# No additional assumptions.
sp_rho_df1 = tibble(
  sp_rho23 = runif(n_sim, -1, 1),
  sp_rho13_2 = runif(n_sim, -1, 1),
  sp_rho24_3 = runif(n_sim, -1, 1),
  sp_rho14_23 = runif(n_sim, -1, 1),
)

# Positive restricted associations only.
sp_rho_df2 = tibble(
  sp_rho23 = runif(n_sim, 0.2, 0.95),
  sp_rho13_2 = runif(n_sim, 0, 0.5),
  sp_rho24_3 = sp_rho13_2,
  sp_rho14_23 = runif(n_sim, 0.15, 0.80),
)

# Positive restricted associations and conditional independence.
sp_rho_df3 = tibble(
  sp_rho23 = runif(n_sim, 0.2, 0.95),
  sp_rho13_2 = 0,
  sp_rho24_3 = sp_rho13_2,
  sp_rho14_23 = runif(n_sim, 0.15, 0.80),
)

# Combine sampled parameters on Spearman's rho scale into a single data set and
# add a new variable that indicates the additional assumptions.
copula_sp_rho_tbl = bind_rows(
  sp_rho_df1 %>% mutate(assumptions = "no"),
  sp_rho_df2 %>% mutate(assumptions = "positive associations"),
  sp_rho_df3 %>% mutate(assumptions = "positive associations and conditional independence")
) %>%
  # Each set of sampled rho-parameters is identified by a unique id.
  mutate(id = row_number()) 



## Convert Sampled Parameters to Copula Scale ------------------------

# Function that converts the Spearman's rho parameter to the copula scale, given
# the parametric copula.
sp_rho_to_copula = function(sp_rho, copula_fam) {
  switch(
    copula_fam,
    frank = {
      copula_fun = copula::frankCopula()
      upper_limit = 35
    },
    gaussian = {
      copula_fun = copula::ellipCopula(family = "normal")
      upper_limit = 1
    },
    clayton = {
      copula_fun = copula::claytonCopula()
      upper_limit = 28
      # Clayton copula can only model positive associations. The corresponding
      # copula will be rotated such that it corresponds to sp_rho.
      sp_rho = abs(sp_rho)
    },
    gumbel = {
      copula_fun = copula::gumbelCopula()
      upper_limit = 50
      # Gumbel copula can only model positive associations. The corresponding
      # copula will be rotated such that it corresponds to sp_rho.
      sp_rho = abs(sp_rho)
    }
  )
  c_pm = copula::iRho(copula_fun, rho = sp_rho)
  if (copula_fam == "clayton")
    c_pm = ifelse(sp_rho < 5 * 1e-04, 1e-05, c_pm)
  c_pm = ifelse(c_pm > upper_limit | is.na(c_pm), upper_limit, c_pm)
  return(c_pm)
}

# Function that returns 90 (degrees rotation) if the Spearman's rho value is
# negative and the parametric copula is a Clayton or Gumbel Copula. Otherwise 0
# (rotation) is returned.
rotation = function(sp_rho, copula_fam) {
  if (copula_fam %in% c("clayton", "gumbel") & sp_rho < 0) {
    return(90)
  }
  else 
    return(0)
}

# Combine parametric copulas and sample parameters into the same data set. 
copula_pm_tbl = copula_fam_df %>%
  cross_join(
    copula_sp_rho_tbl
  )

# Convert the Spearman's rho value to the corresponding value on the copula
# parameter scale for all four copulas. These values are computed before
# constructing the tibble with all possible vine copulas. This avoids computing
# the same copula parameter value multiple times.
copula_pm_tbl = copula_sp_rho_tbl %>%
  cross_join(tibble(copula_fam = c("gaussian", "clayton", "frank", "gumbel"))) %>%
  mutate(
    c23_pm = purrr::map2_dbl(
      .x = sp_rho23,
      .y = copula_fam,
      .f = sp_rho_to_copula
    ),
    r23 = purrr::map2_dbl(.x = sp_rho23, .y = copula_fam, .f = rotation),
    c13_2_pm = purrr::map2_dbl(
      .x = sp_rho13_2,
      .y = copula_fam,
      .f = sp_rho_to_copula
    ),
    r13_2 = purrr::map2_dbl(.x = sp_rho13_2, .y = copula_fam, .f = rotation),
    c24_3_pm = purrr::map2_dbl(
      .x = sp_rho24_3,
      .y = copula_fam,
      .f = sp_rho_to_copula
    ),
    r24_3 = purrr::map2_dbl(.x = sp_rho24_3, .y = copula_fam, .f = rotation),
    c14_23_pm = purrr::map2_dbl(
      .x = sp_rho14_23,
      .y = copula_fam,
      .f = sp_rho_to_copula
    ),
    r14_23 = purrr::map2_dbl(.x = sp_rho14_23, .y = copula_fam, .f = rotation)
  )

# We consider all possible combinations of unidentifiable copulas. 
copula_pm_tbl = copula_pm_tbl %>%
  select(sp_rho23, copula_fam, c23_pm, r23, id, assumptions) %>%
  rename(c23 = copula_fam) %>%
  left_join(
    copula_pm_tbl %>%
      select(sp_rho13_2, copula_fam, c13_2_pm, r13_2, id, assumptions) %>%
      rename(c13_2 = copula_fam),
    by = c("id", "assumptions"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    copula_pm_tbl %>%
      select(sp_rho24_3, copula_fam, c24_3_pm, r24_3, id, assumptions) %>%
      rename(c24_3 = copula_fam),
    by = c("id", "assumptions"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    copula_pm_tbl %>%
      select(sp_rho14_23, copula_fam, c14_23_pm, r14_23, id, assumptions) %>%
      rename(c14_23 = copula_fam),
    by = c("id", "assumptions"),
    relationship = "many-to-many"
  )
# Add copula ID. This is useful later one to keep the data set smaller.
copula_id_tbl = expand_grid(
  c23 = copula_fam_vec,
  c13_2 = copula_fam_vec,
  c24_3 = copula_fam_vec,
  c14_23 = copula_fam_vec
) %>%
  mutate(copula_id = row_number())
copula_pm_tbl = copula_pm_tbl %>%
  left_join(copula_id_tbl)


# Compute ICA  ------------------------------------------------------------

# Helper function to compute the ICA given the bivariate normal margins.
compute_ICA = function(copula_fam, mvn_est, c_unid, r_unid) {
  mean_surr_control = mvn_est$control[1]
  mean_true_control = mvn_est$control[2]
  rho_12 = mvn_est$control[3]
  
  mean_surr_treat = mvn_est$treat[1]
  mean_true_treat = mvn_est$treat[2]
  rho_34 = mvn_est$treat[3]
  
  Surrogate:::compute_ICA_SurvSurv(
    copula_par = c(rho_12, c_unid[1], rho_34, c_unid[2:4]),
    rotation_par = c(0, r_unid[1], 0, r_unid[2:4]),
    copula_family1 = "gaussian",
    copula_family2 = copula_fam,
    n_prec = n_prec,
    q_S0 = function(p) {
      qnorm(p, mean = mean_surr_control)
    },
    q_S1 = function(p) {
      qnorm(p, mean = mean_surr_treat)
    },
    q_T0 = function(p) {
      qnorm(p, mean = mean_true_control)
    },
    q_T1 = function(p) {
      qnorm(p, mean = mean_true_treat)
    },
    composite = FALSE,
    marginal_sp_rho = FALSE,
    seed = 1
  )
}

# Compute ICA for every D-vine copula distribution in copula_pm_tbl and every
# estimated bivariate margin. Because this step is computer intensive, we use
# parallel computing.
copula_pm_tbl = copula_pm_tbl %>%
  cross_join(mvn_est_tbl) %>%
  rowwise(everything()) %>%
  summarise(arguments_list = list(list(
    copula_fam = c(c23, c13_2, c24_3, c14_23),
    c_unid = c(c23_pm, c13_2_pm, c24_3_pm, c14_23_pm),
    r_unid = c(r23, r13_2, r24_3, r14_23),
    mvn_estimates = mvn_estimates
  ))) %>%
  ungroup()

a = Sys.time()
cl = parallel::makeCluster(n_cores, type = "FORK")
parallel::clusterEvalQ(cl, library(Surrogate))
parallel::clusterEvalQ(cl, library(dplyr))
parallel::clusterExport(cl, c("compute_ICA", "n_prec"))
arguments_list = copula_pm_tbl %>%
  `$`("arguments_list")
ICA = parallel::clusterApply(
  x = arguments_list,
  fun = function(arguments_list) {
    compute_ICA(
      arguments_list$copula_fam,
      arguments_list$mvn_estimates,
      arguments_list$c_unid,
      arguments_list$r_unid
    )["ICA"]
  },
  cl = cl
) %>%
  as.numeric() # parallel::clusterApply() returns a list of numerics.
parallel::stopCluster(cl)
print(Sys.time() - a)

# To save space, we split the results_tbl into small blocks without losing
# information.
results_ICA_tbl = copula_pm_tbl %>%
  select(id, copula_id, assumptions, data_set) %>%
  mutate(ICA = ICA)
saveRDS(results_ICA_tbl, file = "results_ICA_tbl.rds")
id_tbl = copula_sp_rho_tbl
saveRDS(id_tbl, file = "id_tbl.rds")
copula_id_tbl
saveRDS(copula_id_tbl, file = "copula_id_tbl.rds")