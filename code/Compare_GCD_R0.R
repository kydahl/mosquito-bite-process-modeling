# Code to generate the GCD and R0 data from ``Once bitten, twice shy: A modeling framework for incorporating heterogeneous mosquito biting into transmission models''

# Description: ----

# There are four model types:
# 1. Standard exponential
# 2. Exponential
# 2. Empirical
# 3. Phenomenological
# 4. Mechanistic

# For each model type, we
#   i) Calculate parameters from specified GCD values
#      - This differs for the mechanistic model where instead we vary specific parameters
#  ii) Calculate R0

# This is placed into a dataframe ("Full_df") containing values of GCD and R0 (and some auxiliary quantities) for each model type

#### Load libraries ####----
library(tidyverse)
library(matlib)

# Fixed Parameter values ----
# Oviposition and resting
gammaV = 1 / (5 * 1440) # exit rate from oviposition to resting, including bloodmeal digestion and site search (5 days)
gammaR = 1/ (2 * 1440) # exit rate from resting to return to blood-feeding (2 days)
# Transmission parameters
betaH = betaB = 0.5
eta = 1/(7 * 1440) # 6 days for infection to develop in vector
mu = 1/(20 * 1440) # 20 day lifespan for vector
gammaH = 1/(5 * 1440) # rate of recovery in hosts (5 days)
muH = 1/(365.25 * 65 * 1440) # host mortality rate (65 years)
KH = 1E3 # host population density

KJ = 0.75 * KH # Aquatic-stage mosquito population
rhoJ = 1 / (12 * 1440) # 12 day larval development period
muJ = 1 / (20 * 1440) # 62.5% probability of larval survival (20 / (20 + 12))
varPhi = 3 / 1440 # on average 3 eggs per female per day over the course of the oviposition period

### Generic functions: ----
# Functions used for several model types 

# Function: replace small numerical values with zeroes
zero_out = function(vector_in) {
  in_rows = dim(vector_in)[1]
  in_cols = dim(vector_in)[2]
  vector_out = vector_in
  vector_out[vector_out < .Machine$double.eps] <- 0
  vector_out = matrix(vector_out, nrow = in_rows, ncol = in_cols)
  return(vector_out)
}

# Function: Get stable mosquito population distributions
get_stable_pop <- function(A_matrix, v_alpha) {
  
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  
  if (A_dim == 1) {
    green_matrix = 1/(mu * diag(A_dim) - t(A_matrix))
  } else {
    green_matrix = inv(mu * diag(A_dim) - t(A_matrix))
  }
  
  tau = t(out_rates) %*% green_matrix %*% v_alpha %>% as.double()
  varrho = 1 - tau * (gammaR / (mu + gammaR)) * (gammaV / (mu + gammaV))
  nG = 1 / varrho
  N_offspring = (varPhi / (mu + gammaV)) * (rhoJ / (muJ + rhoJ)) * tau * nG
  
  if (N_offspring > 1) {
    J_star = KJ * (1 - 1 / N_offspring)
    v_B_star = rhoJ * J_star * green_matrix %*% v_alpha * nG
    V_star = rhoJ * J_star * (1 / (mu + gammaV)) * tau * nG
    R_star = rhoJ * J_star * (1 / (mu + gammaR)) * (gammaV / (mu + gammaV)) * tau * nG
  } else {
    J_star = 0
    v_B_star = matrix(rep(0, A_dim), ncol = 1)
    V_star = 0
    R_star = 0
  }
  
  return(
    c(
      tau = tau,
      nG = nG,
      N_offspring = N_offspring,
      J_star = J_star,
      v_B_star = list(v_B_star),
      V_star = V_star,
      R_star = R_star)
  )
}

# Function: Calculate theta, duration of the biting process
theta_calc <- function(A_matrix, v_alpha) {
  A_dim = dim(A_matrix)[1]
  green_A = inv(-(A_matrix)*1440) * 1440
  v_one = matrix(rep(1, A_dim), ncol = 1)
  theta = t(v_alpha) %*% green_A %*% v_one
}

# Function: Calculate GCD, gonotrophic cycle duration
GCD_calc <- function(A_matrix, v_alpha) {
  # Calculate gonotrophic cycle duration as a function of phase-type distribution parameters
  GCD = theta_calc(A_matrix, v_alpha) + (1/gammaV) + (1/gammaR)
}

# Function: Calculate R0, basic reproduction number
R0_calc <- function(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB) {
  # Calculates R0 as a function of the phase-type distribution parameters
  # and specific transmission formulation
  if (is.null(dim(A_matrix))) {
    A_matrix = as.matrix(A_matrix)
    v_alpha = as.matrix(v_alpha)
  }
  
  stable_pop = get_stable_pop(as.matrix(A_matrix), as.matrix(v_alpha))
  if (stable_pop$N_offspring < 1) {
    R0 = 0
  } else {
    
    KB = sum(stable_pop$v_B_star)
    A_dim = dim(A_matrix)[1]
    v_ones = matrix(rep(1, A_dim), ncol = 1)
    out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
    # Independent terms
    host_infectious_period = 1/(gammaH + muH)
    
    # Intermediate terms
    spec_mat = v_alpha %*% t(out_rates)
    
    if (dim(A_matrix)[1] == 1) {
      GammaI = 1/(mu * diag(A_dim) - t(A_matrix) - (gammaV / (mu + gammaV)) * (gammaR / (mu + gammaR)) * spec_mat)
      GammaE = 1/((eta + mu) * diag(A_dim) - t(A_matrix) - (gammaV / (mu + eta + gammaV)) * (gammaR / (mu + eta + gammaR)) * spec_mat)
    } else {
      GammaI = inv(mu * diag(A_dim) - t(A_matrix) - (gammaV / (mu + gammaV)) * (gammaR / (mu + gammaR)) * spec_mat)
      GammaE = inv((eta + mu) * diag(A_dim) - t(A_matrix) - (gammaV / (mu + eta + gammaV)) * (gammaR / (mu + eta + gammaR)) * spec_mat)
    }
    
    complicated_probability = (gammaV/(mu+gammaV)) * ((eta / (mu+gammaV+eta)) * (gammaR/(mu+gammaR+eta)) + (eta/(mu+gammaR+eta) * (gammaR/(mu+gammaR))))
    
    tauE = (eta * diag(A_dim) + complicated_probability * spec_mat) %*% GammaE
    
    R02 = host_infectious_period * t(v_ones) %*% BetaH %*% LambdaH %*% GammaI %*% tauE %*% BetaB %*% LambdaB %*% stable_pop$v_B_star  / KB
    R0 = sqrt(as.double(R02))
    
  }
  
  return(R0)
}

### Standard case functions: ----
# Function: Get stable population distribution
get_standard_stable_pop <- function(b) {
  
  check = (mu / (1440 *varPhi * b)) * ((rhoJ + muJ)/ rhoJ)
  
  if (check > 1) {
    J_star = 0
    B_star = 0
  } else {
    J_star = KJ * (1 - check)
    B_star = rhoJ * J_star / mu
  }
  
  return(
    c(
      J_star = J_star,
      B_star = B_star)
  )
}
# Get R0 value
standard_R0_calc <- function(b) {
  B_star = as.list(get_standard_stable_pop(b))$B_star
  host_infectious_period = (1 / (muH + gammaH))
  standard_beta = 0.175 * betaB 
  # this is approximately = betaB * 1440 / sqrt(host_infectious_period * (eta / (mu + eta)) * (1 / mu) * KJ / KH)
  
  R0 = b * sqrt(host_infectious_period * standard_beta * (eta / (mu + eta)) * (1 / mu) * standard_beta * B_star / KH)
}

### Exponential case functions: ----
# Function: Get epidemiological terms (LambdaH and LambdaB and BetaH and BetaB)
exp_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = -as.matrix(A_matrix) %*% v_ones
  
  LambdaB = diag(as.vector(out_rates), nrow = length(out_rates)) %>%
    list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% list()
  
  BetaH = (betaH * diag(A_dim)) %>%
    list()
  BetaB = (betaB * BetaH[[1]] / betaH) %>% list()
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

### Empirical case functions: ----
# Function: transient rate matrix A from given theta value
get_emp_A <- function(theta) {
  ## Generate 100 samples from a lognormal distribution with mean theta and variance 1...
  sample_num = 100
  
  # Correct values of the lognormal parameters to match R's implementation
  lognorm_sigma <- sqrt(log(1 + (1 / theta^2)))  # Corrected sigma
  lognorm_mu <- log(theta) - 0.5 * lognorm_sigma^2  # Corrected mu
  
  samples <- tibble(theta = theta) %>%
    expand_grid(tibble(n = 1:sample_num)) %>%
    rowwise() %>%
    mutate(sample = rlnorm(1, meanlog = lognorm_mu, sdlog = lognorm_sigma))
  
  # Set dimension of PH distribution to be fit
  emp_A_dim = 8
  
  ## Run EM algorithm to obtain best-fit A_matrix and v_alpha...
  library(mapfit)
  sample_in = samples$sample
  # set.seed(90210) # if you'd like to ensure same initial state
  out = phfit.point(ph = ph(emp_A_dim),
                    x = sample_in)
  v_alpha = as.matrix(out$alpha)
  
  return(c(
    v_alpha = list(v_alpha),
    A_matrix = list(as.matrix(out$Q)))
  )
  
}

# Function: Get epidemiological terms (LambdaH and LambdaB and BetaH and BetaB)
emp_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = zero_out(-as.matrix(A_matrix) %*% v_ones)
  
  LambdaB = diag(as.vector(out_rates)) %>%
    zero_out() %>%  list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% zero_out() %>%  list()
  
  BetaH = diag(A_dim) %>%
    zero_out() %>%  list()
  BetaB = BetaH
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

### Phenomenological case functions: ----
# Function: get transient rate matrix from theta
get_phenom_A <- function(theta) {
  b = 1 / theta
  
  v_alpha = matrix(c(0.53, 0.42, 0, 0.05, 0, 0), ncol = 1)
  
  temp_matrix = matrix(rep(0, 36), nrow = 6, ncol = 6 )
  diag(temp_matrix) = c(-b, -2*b, -2*b, -3*b, -3*b, -3*b)
  temp_matrix[3,2] = 2*b
  temp_matrix[5,4] = 3*b
  temp_matrix[6,5] = 3*b
  
  A_matrix = t(temp_matrix)
  
  return(c(
    v_alpha = list(v_alpha),
    A_matrix = list(A_matrix)
  ))
  
}

# Function: Get epidemiological terms (LambdaH and LambdaB and BetaH and BetaB)
phenom_EpiTerms <- function(A_matrix, v_alpha) {
  KB = sum(get_stable_pop(A_matrix, v_alpha)$v_B_star)
  A_dim = dim(A_matrix)[1]
  v_ones = matrix(rep(1, A_dim), ncol = 1)
  out_rates = -as.matrix(A_matrix) %*% v_ones
  
  LambdaB = -diag(diag(A_matrix)) %>%
    list()
  LambdaH = (LambdaB[[1]] * KB / KH) %>% list()
  
  BetaH = (betaH * diag(A_dim)) %>%
    list()
  BetaB = (BetaH[[1]] * betaB / betaH) %>% list()
  
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

### Mechanistic case functions: ----
# Function: build a dataframe of the mechanistic parameters with 'varied_parameter' varied throughout it's largest biologically resonable range
vary_parameter_function = function(df_in, varied_parameter) {
  # Collect and keep track of the baseline value of the chosen parameter
  baseline_values = df_in %>%
    select(mosquito_type, any_of(varied_parameter)) %>%
    filter(!is.na(!!as.symbol(varied_parameter))) %>%
    unique() %>%
    mutate(parameter_type = "baseline") %>%
    right_join(df_in %>%
                 select(-c(any_of(varied_parameter))),
               by = "mosquito_type") %>%
    distinct() %>%
    mutate(varied_parameter = varied_parameter)
  
  # Set up variation in probability parameters
  if (varied_parameter %in% c("pQ", "pL", "pP", "pG", "sigma")) {
    # The probability parameters (pL, pP, pG, and sigma) will vary from 0 to 1
    prob_vec = seq(0, 1, length.out = variation_resolution+1)
    prob_vec = prob_vec[-1] # get rid of the zero entry
    
    # Dataframe with varied values
    vary_df = baseline_values %>%
      cross_join(tibble(prob = prob_vec)) %>%
      mutate(across(any_of(as.name(varied_parameter)), ~  prob))%>%
      select(-prob) %>%
      mutate(parameter_type = "varied") %>%
      rbind(baseline_values)
    # Set up variation in rate parameters
  } else {
    rate_vec = 10^seq(-6, 2, length.out = variation_resolution)
    # Dataframe with varied values
    vary_df = baseline_values %>%
      cross_join(tibble(rate = rate_vec)) %>%
      mutate(across(any_of(as.name(varied_parameter)), ~  rate * .x)) %>%
      select(-rate) %>%
      mutate(parameter_type = "varied") %>%
      rbind(baseline_values)
  }
  
  return(vary_df)
}

# Function: Get epidemiological terms (LambdaH and LambdaB and BetaH and BetaB)
mech_EpiTerms <- function(v_B_star, lP, lG, varied_parameter) {
  
  if (varied_parameter == "pG") {
    correction_for_plot = 0.175
  } else {
    correction_for_plot = 1
  }
  
  KB = sum(unlist(v_B_star))
  v_ones = matrix(rep(1, 4), ncol = 1)
  
  LambdaH = matrix(c(0, 0, 0,            0,
                     0, 0, 0,            0,
                     0, 0, lP * KB / KH, 0,
                     0, 0, 0,            0
  ),
  nrow = 4) %>% list()
  LambdaB = matrix(c(0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, lG
  ),
  nrow = 4) %>% list()
  
  BetaH = matrix(c(0, 0, 0,                           0,
                   0, 0, 0,                           0,
                   0, 0, correction_for_plot * betaH, 0,
                   0, 0, 0,                           0
  ),
  nrow = 4) %>% list()
  BetaB = matrix(c(0, 0, 0, 0,     
                   0, 0, 0, 0,     
                   0, 0, 0, 0,     
                   0, 0, 0, correction_for_plot * betaB
                   ),
  nrow = 4) %>% list()
  return(
    c(
      LambdaB = LambdaB,
      LambdaH = LambdaH,
      BetaB = BetaB,
      BetaH = BetaH)
  )
}

# GCD range ----
# Set range of GCD values to consider
# Note: GCD = biting duration + oviposition duration + resting duration
#       Only theta = biting duration will be varied
resolution = 501 # resolution(-ish) of the GCD VECTOR

theta_min = 0 # minimum biting state duration of zero (will be removed later)
theta_max = 1 / mu # maximum biting state duration equivalent to mosquito lifespan
theta_vec = seq(theta_min, theta_max, length.out = resolution)
theta_vec = theta_vec[-1]
inv_theta_vec = seq(1e-9, 3, length.out = resolution) # we want evenly spaced values for 1/theta as well
# Combine linear theta and 1/theta vectors
theta_vec = sort(unique(c(theta_vec, 1/inv_theta_vec, (1/8) * 1440, (1/4) * 1440, (1/3) * 1440, (1/2) * 1440, 1 * 1440, 2 * 1440, 3 * 1440)))
# Make vector of GCD values
GCD_vec = theta_vec + (1/gammaR) + (1/gammaV)

# 1. Standard exponential model ----
# Use an alternative vector of GCD values for the standard model
Standard_vec = seq(1/(3 * 1440), min(GCD_vec), length.out = resolution)
Standard_vec = c(Standard_vec, GCD_vec, (1/2) * 1440, 1440, 2 * 1440)

Standard_df <- tibble(theta = Standard_vec) %>%
  mutate(b = 1/theta) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_standard_stable_pop(b))) %>% unnest_wider(stable_pops) %>%
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = standard_R0_calc(b))

# 2. Exponential model with resting and ovipositing ----
Exponential_df <- tibble(theta = theta_vec) %>%
  # get A and alpha
  mutate(
    A_matrix = matrix(-1/theta),
    v_alpha = matrix(1)
  ) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(exp_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  mutate(
    LambdaB = as.list(LambdaB),
    LambdaH = as.list(LambdaH),
    BetaB = as.list(BetaB),
    BetaH = as.list(BetaH),
  ) %>% 
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# 3. Empirical model ----
Empirical_df <- tibble(theta = theta_vec) %>%
  # get A and alpha
  rowwise() %>%
  mutate(params = list(get_emp_A(theta))) %>% unnest_wider(params) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(emp_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  # basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# 4. Phenomenological model ----
Phenom_df <- tibble(theta = theta_vec) %>%
  rowwise() %>%
  mutate(mats = list(get_phenom_A(theta))) %>% unnest_wider(mats) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # Epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(phenom_EpiTerms(A_matrix, v_alpha))) %>% unnest_wider(EpiTerms) %>%
  # Basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# 5. Mechanistic model ----
# Set up parameter ranges

# Define "flighty" and "persistent" parameter sets
flighty_parameters = tibble(
  mosquito_type = "flighty",
  # Questing
  pQ = 1,
  lQ = 1 / 480, # 8 hours = 480 minutes
  # Landing
  pL =  0.5,
  lL = 0.1, # 10 minutes
  # Probing
  pP = 0.5,
  lP = 0.2, # 5 minutes
  # Ingesting
  pG = 0.5,
  lG = 1, # 1 minutes
  # Fleeing
  sigma = 1 - 0.9
)
persistent_parameters = tibble(
  mosquito_type = "persistent",
  # Questing
  pQ = 1,
  lQ = 1 / 480,  # 8 hours = 480 minutes
  # Landing
  pL =  0.7,
  lL = 0.1, # 10 minutes
  # Probing
  pP = 0.8,
  lP = 0.2, # 5 minutes
  # Ingesting
  pG = 0.9,
  lG = 1, # 1 minutes
  # Fleeing
  sigma = 1 - 0.66
)
base_parameters = rbind(flighty_parameters, persistent_parameters)

# Set up parameter variation
variation_resolution = resolution

# Initialize data frame
full_variation_df = as_tibble(matrix(
  nrow = 0,
  ncol = length(colnames(base_parameters))+2,
  dimnames = list(NULL, c("parameter_type", "varied_parameter", colnames(base_parameters)))
))

# Name the parameters we'll vary
parameter_characters = c("lQ","pQ", "pL", "lL", "pP", "lP", "pG", "lG", "sigma")

# Append all the varied parameter sets
for (parameter_name in parameter_characters) {
  new_df = vary_parameter_function(base_parameters, parameter_name)
  full_variation_df = rbind(full_variation_df, new_df)
}

# Build data frame
Mech_df <- full_variation_df %>%
  relocate(mosquito_type, parameter_type, varied_parameter) %>%
  # Transform parameters into A matrix
  rowwise() %>%
  # Calculate transient rate matrix directly
  mutate(A_matrix = list(matrix(c(
    -lQ+(1-pQ)*lQ,                                 pQ * lQ,       0,       0,
    (1 - sigma) * (1- pL) * lL, -lL + sigma * (1- pL) * lL, pL * lL,       0,
    (1 - sigma) * (1 - pP) * lP,     sigma * (1 - pP) * lP,     -lP, pP * lP,
    (1 - sigma) * (1 - pG) * lG,     sigma * (1 - pG) * lG,       0,      -lG
  ),
  ncol = 4, byrow = TRUE) )) %>%
  mutate(v_alpha = list(matrix(c(1, 0, 0, 0), ncol = 1))) %>%
  # Get the theta value for each parameter set
  rowwise() %>%
  mutate(theta = theta_calc(A_matrix, v_alpha)) %>%
  # stable population terms
  rowwise() %>%
  mutate(stable_pops = list(get_stable_pop(A_matrix, v_alpha))) %>% unnest_wider(stable_pops) %>%
  # Epidemiological terms
  rowwise() %>%
  mutate(EpiTerms = list(mech_EpiTerms(v_B_star, lP, lG, varied_parameter))) %>% unnest_wider(EpiTerms) %>%
  # Basic reproduction number
  rowwise() %>%
  mutate(R0 = R0_calc(A_matrix, v_alpha, LambdaH, LambdaB, BetaH, BetaB))

# Save Mechanistic table separately (to use for an independent figure)
saveRDS(Mech_df, "data/Mechanistic_results.rds")

# Join all data ----
Full_df <- bind_rows(
  # 1. Standard model
  Standard_df %>% 
    mutate(`Model type` = "Standard",
           A_matrix = as.list(-b),
           v_alpha = as.list(c(1)),
           B_star = as.list(B_star),
           N_offspring = 1/((mu / (1440 *varPhi * b)) * ((rhoJ + muJ)/ rhoJ))) %>% 
    rename(v_B_star = B_star),
  # 2. Exponential model
  Exponential_df %>% 
    mutate(`Model type` = "Exponential",
           v_B_star = as.list(v_B_star),
           A_matrix = as.list(A_matrix),
           v_alpha = as.list(v_alpha)),
  # 3. Empirical model
  Empirical_df %>% 
    mutate(`Model type` = "Empirical"),
  # 4. Phenomenological model
  Phenom_df %>% 
    mutate(`Model type` = "Phenomenological"),
  # 5. Mechanistic model (pG, lQ, and pP)
  Mech_df %>%
    filter(
      mosquito_type == "persistent",
      varied_parameter %in% c("pG", "lQ", "pP")) %>%
    mutate(`Model type` = "Mechanistic")
) %>%
  rowwise() %>%
  mutate(GCD = if_else(
    `Model type` == "Standard",
    theta - (1/gammaR) - (1/gammaV),
    theta))

# Save data
saveRDS(Full_df, "data/GCD_R0_data.rds")
