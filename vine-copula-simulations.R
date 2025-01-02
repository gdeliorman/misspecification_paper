# Preparation -------------------------------------------------------------

library(tidyverse)
library(compositions)
library(rvinecopulib)
library(Surrogate)
ncores = 10
n_lognormal = 2e3
n_prec = 1e4

# Read in sampled parameters.
Sigma_df = readxl::read_excel("Log_normal_results_23december.xlsx", sheet = 3)
mu_df = readxl::read_excel("Log_normal_results_23december.xlsx", sheet = 2)
# Convert imported data to a tibble with matrices as elements.
parameters_tbl = Sigma_df %>%
  group_by(i) %>%
  summarize(Sigma = list(matrix(c(sigma1, sigma2, sigma3, sigma4), ncol = 4))) %>%
  left_join(mu_df)


# Helper Functions --------------------------------------------------------

# Function to generate data from the lognormal multivariate distribution.
sample_lnormal = function(n, mu, Sigma) {
  rlnorm.rplus(n, mu, Sigma)
}

# Function to fit the D-vine copula to the simulated data.
fit_dvine_copula = function(data, copula_family){
  # Convert data to pseudo-observations with uniform margins. 
  u = pseudo_obs(data)
  # Estimate copula parameters.
  fitted_vine = vinecop(
    data = u, 
    family_set = copula_family,
    structure = dvine_structure(order = 1:4)
  )
  # Return fitted vine copula object.
  return(fitted_vine)
}

# Extract required parameters from fitted D-vine copula model.
extract_params_dvine = function(fitted_object) {
  copula_par = summary(fitted_object)$parameters %>% as.numeric()
  rotation_par = summary(fitted_object)$rotation
  copula_family1 = summary(fitted_object)$family[c(1, 3)]
  copula_family2 = summary(fitted_object)$family[-c(1, 3)]
  
  return(
    list(
      copula_par = copula_par,
      rotation_par = rotation_par,
      copula_family1 = copula_family1,
      copula_family2 = copula_family2
    )
  )
}

# Function to fit (log)normal marginal distribution
fit_marginal = function(x, distr_family = "lognormal"){
  if (distr_family == "lognormal") {
    # Parameters of the lognormal distribution are estimated through maximum likelihood.
    mu = mean(log(x))
    sigma_sq = mean((log(x) - mean(log(x)))**2)
    est_params = c(
      mu = mu,
      sigma_sq = sigma_sq
    )
    est_q = function(p) {
      qlnorm(p, meanlog = mu, sdlog = sqrt(sigma_sq))
    }
  }
  else if (distr_family == "normal") {
    mu = mean(x)
    sigma_sq = var(x)
    est_params = c(
      mu = mu,
      sigma_sq = sigma_sq
    )
    est_q = function(p) {
      qnorm(p, mean = mu, sd = sqrt(sigma_sq))
    }
  }
  return(list(params = est_params,
              q = est_q))
}

# Function to estimate the ICA given the fitted D-vine copula model.
compute_ICA = function(copula_par,
                       rotation_par,
                       copula_family1,
                       copula_family2,
                       q_S0,
                       q_T0,
                       q_S1,
                       q_T1,
                       seed = 1,
                       n_prec = 1e4) {
  Surrogate:::compute_ICA_ContCont(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family1,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_T0 = q_T0,
    q_S1 = q_S1,
    q_T1 = q_T1,
    seed = seed,
    marginal_sp_rho = FALSE
  )
}

estimate_ICA = function(data, distr_family, copula_family, seed = 1, n_prec = 1e4) {
  # Estimate dvine copula model.
  fitted_dvine = fit_dvine_copula(data, copula_family)
  # Extract dvine parameters.
  dvine_params = extract_params_dvine(fitted_dvine)
  
  # Estimate marginal distributions.
  est_marg1 = fit_marginal(data[, 1], distr_family)
  est_marg2 = fit_marginal(data[, 2], distr_family)
  est_marg3 = fit_marginal(data[, 3], distr_family)
  est_marg4 = fit_marginal(data[, 4], distr_family)
  
  computed_ICA_row = compute_ICA(
    copula_par = dvine_params$copula_par,
    rotation_par = dvine_params$rotation_par,
    copula_family1 = dvine_params$copula_family1,
    copula_family2 = dvine_params$copula_family2,
    n_prec = n_prec,
    q_S0 = est_marg1$q,
    q_T0 = est_marg2$q,
    q_S1 = est_marg3$q,
    q_T1 = est_marg4$q, 
    seed = seed
  )
  
  return(computed_ICA_row["ICA"])
}

compute_ICA_lognormal = function(mu,
                                 Sigma,
                                 distr_family,
                                 copula_family,
                                 seed = 1,
                                 n_lognormal,
                                 n_prec) {
  # Simulate multivariate lognormal data.
  data = sample_lnormal(n_lognormal, mu, Sigma)
  # Estimate ICA.
  estimated_ICA = NA
  try({
    estimated_ICA = estimate_ICA(
      data = data,
      distr_family = distr_family,
      copula_family = copula_family,
      seed = seed,
      n_prec = n_prec
    )
  })
  
  return(estimated_ICA)
}


# Computations ------------------------------------------------------------

# For each set of simulation parameters, consider all combinations of parametric
# copulas and marginal distributions.
simulation_parameters_tbl = parameters_tbl %>%
  cross_join(expand_grid(
    copula_family = c("clayton", "gaussian", "gumbel", "frank"),
    distr_family = c("normal", "lognormal")
  ))
# Extract lists with information for computing the ICA.
list_mu = simulation_parameters_tbl %>%
  rowwise() %>%
  summarise(list(c(mu1, mu2, mu3, mu4))) %>%
  pull()
list_Sigma = simulation_parameters_tbl$Sigma
list_distr_family = simulation_parameters_tbl$distr_family
list_copula_family = simulation_parameters_tbl$copula_family
# Compute the ICA for all settings.
a = Sys.time()
cl = parallel::makeCluster(ncores)
parallel::clusterExport(
  cl = cl,
  varlist = c(
    "sample_lnormal",
    "fit_dvine_copula",
    "extract_params_dvine",
    "fit_marginal",
    "compute_ICA",
    "estimate_ICA",
    "compute_ICA_lognormal"
  )
)
parallel::clusterEvalQ(cl = cl, {
  library(tidyverse)
  library(compositions)
  library(rvinecopulib)
  library(Surrogate)
})
ICA_vec = parallel::clusterMap(
  cl = cl, 
  fun = compute_ICA_lognormal,
  mu = list_mu,
  Sigma = list_Sigma,
  distr_family = list_distr_family,
  copula_family = list_copula_family, 
  MoreArgs = list(
    seed = 1,
    n_lognormal = n_lognormal,
    n_prec = n_prec
  ) 
) %>% as.numeric()
parallel::stopCluster(cl)
print(Sys.time() - a)
# Add computed ICAs to the tibble.
simulation_dvine_results_tbl = simulation_parameters_tbl %>%
  mutate(ICA_vec)

# Save Results ------------------------------------------------------------
saveRDS(simulation_dvine_results_tbl,
        "simulation_dvine_results_tbl.rds")

