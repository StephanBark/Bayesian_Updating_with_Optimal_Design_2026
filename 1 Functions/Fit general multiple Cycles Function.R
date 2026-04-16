#### "fit_general_multiple_cycles" function ############################################
#
# Idea: This function fits a flexible Bayesian hierarchical mixed model to multiple
# historical cycles using the Stan framework. It automatically splits data into
# cycles, allows customizable random and fixed effect structures with flexible
# prior specifications, estimates posterior distributions of model parameters, and 
# updates priors based on posterior estimates from previous cycles for sequential 
# Bayesian updating.
#
# Input: MET_data: A data frame containing multi-environment trial data with variables
#                  such as yield, year, trial, rep, genotype, environment, etc.
#        fixed: A list specifying fixed effects. Default is list("beta_0") for 
#               intercept only. Categorical fixed effects use reference coding.
#        random: A list specifying random effects. Empty list defaults to standard 
#                MET random effects structure: year, trial, trial:year, trial:rep:year,
#                genotype:year, genotype:trial. Interaction effects use ":" separator.
#        priors: A list specifying prior distributions for variance components. Can use
#                R/G structure (GLMM_MCMC style) or direct effect naming. Supports
#                inverse_gamma and inverse_wishart priors. Default values are 1 for
#                inverse gamma parameters.
#        cycles: A list specifying years per cycle. Empty list defaults to 3 balanced
#                cycles from the most recent years in the data.
#        iter: A numeric vector of MCMC iterations per cycle, or single value applied
#              to all cycles. Default is 2000 iterations per cycle.
#        warmup: A numeric vector of warmup iterations per cycle, or single value
#                applied to all cycles. Default is 60% of iter per cycle.
#        chains: Number of MCMC chains. Default is 4.
#        thin: Thinning factor for MCMC sampling. Default is 2.
#        delta: A numeric vector of adaptation delta hyperparameters for HMC sampling
#               per cycle, or single value applied to all cycles. Default is 0.8.
#        stepsize: A numeric vector of stepsize parameters for HMC per cycle, or
#                  single value applied to all cycles. Default is 0.1.
#        max_treedepth: A numeric vector of maximum tree depth for HMC per cycle, or
#                       single value applied to all cycles. Default is 10.
#        seed: Random seed for reproducibility. Default is 123.
#        init: Initialization method for MCMC ('random' or numeric). Default is 'random'.
#
# Output: A list containing multiple elements:
#         $cycles: A list of fitted Stan model objects and posterior samples
#                  for each cycle, with comprehensive diagnostics
#         $parameters: Initial priors, final posteriors, parameter evolution,
#                      and posterior means from the final cycle
#         $diagnostics: Stan convergence information, divergent
#                       transitions, Rhat values, effective sample sizes
#         $model_spec: Complete model specification and MCMC settings
#         $rerun_info: Information for pipeline continuation and reruns
#         $session_info: Call details, timestamps, and package versions

fit_general_multiple_cycles <- function(MET_data, 
                                        fixed = list("beta_0"),
                                        random = list(),
                                        priors = list(),
                                        cycles = list(),
                                        iter = NULL, warmup = NULL, delta = NULL, 
                                        stepsize = NULL, max_treedepth = NULL, 
                                        chains = 4, thin = 2, seed = 123, 
                                        init = 'random') {
  
  ## Load packages -----------------------------------------------------------
  
  library(showtext)
  library(data.table)
  library(lubridate)
  library(OptimalDesign)
  library(fitdistrplus)
  library(matrixcalc)
  library(reshape2)
  library(viridis)
  library(writexl)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(rstan)
  library(MCMCpack)
  library(MASS)
  library(coda)
  library(tidyr)
  library(ggmcmc)
  library(GGally)
  library(SpATS)
  library(fitdistrplus)
  
  ## Section 1: HELPER FUNCTIONS ---------------------------------------------
  
  # Helper function to parse model effects
  parse_effects <- function(effect_list) {
    if (length(effect_list) == 0) {
      # Default effects structure
      return(list(
        "beta_0" = list(type = "simple", vars = "beta_0")
      ))
    }
    
    parsed <- list()
    for (effect in effect_list) {
      if (grepl(":", effect)) {
        vars <- strsplit(effect, ":")[[1]]
        parsed[[effect]] <- list(type = "interaction", vars = vars)
      } else {
        parsed[[effect]] <- list(type = "simple", vars = effect)
      }
    }
    return(parsed)
  }
  
  # Helper function to parse prior specifications
  parse_priors <- function(priors_list, random_effects) {
    # Check for R/G structure (GLMM_MCMC style)
    if ("R" %in% names(priors_list) || "G" %in% names(priors_list)) {
      # R/G structure
      named_priors <- list()
      
      # Handle residual (R component)
      if ("R" %in% names(priors_list)) {
        named_priors[["residual"]] <- priors_list[["R"]]
      }
      
      # Handle random effects (G component)
      if ("G" %in% names(priors_list)) {
        g_priors <- priors_list[["G"]]
        effect_names <- names(random_effects)
        
        if (is.null(names(g_priors)) || all(names(g_priors) == "")) {
          # Unnamed G list - map by position
          for (i in 1:min(length(g_priors), length(effect_names))) {
            named_priors[[effect_names[i]]] <- g_priors[[i]]
          }
        } else {
          # Named G list - map by name
          for (effect_name in effect_names) {
            if (effect_name %in% names(g_priors)) {
              named_priors[[effect_name]] <- g_priors[[effect_name]]
            }
          }
        }
      }
    } else {
      # Original structure - check if named or unnamed
      if (is.null(names(priors_list)) || all(names(priors_list) == "")) {
        # Unnamed list - create mapping (old way)
        if (length(priors_list) == 0) {
          named_priors <- list()
        } else {
          effect_names <- c("residual", names(random_effects))
          if (length(priors_list) > length(effect_names)) {
            stop("Too many prior specifications provided for available effects")
          }
          
          named_priors <- list()
          for (i in 1:length(priors_list)) {
            if (i <= length(effect_names)) {
              named_priors[[effect_names[i]]] <- priors_list[[i]]
            }
          }
        }
      } else {
        # Named list - use as is
        named_priors <- priors_list
      }
    }
    
    parsed_priors <- list()
    
    # Handle residual variance prior
    if ("residual" %in% names(named_priors)) {
      residual_spec <- named_priors[["residual"]]
      
      # Validate residual prior specification
      if (!"type" %in% names(residual_spec)) {
        stop("Prior type not specified for residual variance")
      }
      
      if (residual_spec$type == "inverse_gamma") {
        # Set default parameters if not provided
        if (!"alpha" %in% names(residual_spec)) {
          residual_spec$alpha <- 1
        }
        if (!"beta" %in% names(residual_spec)) {
          residual_spec$beta <- 1
        }
      } else {
        stop(paste("Unsupported prior type for residual variance:", residual_spec$type))
      }
      
      parsed_priors[["residual"]] <- residual_spec
    } else {
      # Default residual prior
      parsed_priors[["residual"]] <- list(type = "inverse_gamma", alpha = 1, beta = 1)
    }
    
    for (effect_name in names(random_effects)) {
      if (effect_name %in% names(named_priors)) {
        prior_spec <- named_priors[[effect_name]]
        
        # Validate prior specification
        if (!"type" %in% names(prior_spec)) {
          stop(paste("Prior type not specified for effect:", effect_name))
        }
        
        if (prior_spec$type == "inverse_wishart") {
          # Validate inverse Wishart specification
          if (!"blocks" %in% names(prior_spec)) {
            stop(paste("Block structure not specified for inverse Wishart prior on effect:", effect_name))
          }
          
          # Set default parameters if not provided
          if (!"nu" %in% names(prior_spec)) {
            prior_spec$nu <- "default"
          }
          if (!"psi" %in% names(prior_spec)) {
            prior_spec$psi <- "identity"
          }
        } else if (prior_spec$type == "inverse_gamma") {
          # Set default parameters if not provided
          if (!"alpha" %in% names(prior_spec)) {
            prior_spec$alpha <- 1
          }
          if (!"beta" %in% names(prior_spec)) {
            prior_spec$beta <- 1
          }
        } else {
          stop(paste("Unsupported prior type:", prior_spec$type, "for effect:", effect_name))
        }
        
        parsed_priors[[effect_name]] <- prior_spec
      } else {
        # Default to inverse gamma
        parsed_priors[[effect_name]] <- list(type = "inverse_gamma", alpha = 1, beta = 1)
      }
    }
    
    return(parsed_priors)
  }
  
  # Helper function to get all variables used in the model
  get_model_variables <- function(random_effects, fixed_effects) {
    all_vars <- c("yield") # Always need the response variable
    
    # Add variables from random effects
    for (effect in random_effects) {
      all_vars <- c(all_vars, effect$vars)
    }
    
    # Add variables from fixed effects (excluding intercept)
    for (effect in fixed_effects) {
      if (is.character(effect) && length(effect) == 1 && effect != "beta_0") {
        all_vars <- c(all_vars, effect)
      }
    }
    
    return(unique(all_vars))
  }
  
  # Helper function to process cycles
  process_cycles <- function(data, cycles_spec) {
    data$year_nr <- as.numeric(data$year)
    unique_years <- sort(unique(data$year_nr))
    
    if (length(cycles_spec) == 0) {
      # Default: 3 balanced cycles
      n_years <- length(unique_years)
      years_per_cycle <- floor(n_years / 3)
      if (years_per_cycle == 0) stop("Not enough years for 3 cycles")
      
      # Use only the most recent years that can be evenly divided
      years_to_use <- years_per_cycle * 3
      start_year <- unique_years[n_years - years_to_use + 1]
      unique_years <- unique_years[unique_years >= start_year]
      
      cycles_spec <- rep(years_per_cycle, 3)
    }
    
    if (sum(cycles_spec) > length(unique_years)) {
      stop("Specified cycles require more years than available in data - If correct check if year variable is written out small!")
    }
    
    cycles_data <- list()
    year_idx <- 1
    
    for (i in 1:length(cycles_spec)) {
      years_in_cycle <- unique_years[year_idx:(year_idx + cycles_spec[i] - 1)]
      cycle_data <- subset(data, year_nr %in% years_in_cycle)
      
      # Adjust years to start from 1 within each cycle
      cycle_data$year_nr <- as.integer(factor(cycle_data$year_nr))
      
      cycles_data[[i]] <- cycle_data
      year_idx <- year_idx + cycles_spec[i]
    }
    
    return(cycles_data)
  }
  
  # Helper function to determine the structure of inverse Wishart effects
  get_wishart_structure <- function(effect_name, effect_info, prior_spec, cycle_data) {
    blocks_var <- prior_spec$blocks
    covariance_vars <- setdiff(effect_info$vars, blocks_var)
    
    if (length(covariance_vars) != 1) {
      stop(paste("Inverse Wishart prior for", effect_name, "must have exactly one covariance dimension"))
    }
    
    n_blocks <- length(unique(cycle_data[[blocks_var]]))
    cov_dimension <- length(unique(cycle_data[[covariance_vars[1]]]))
    
    return(list(
      blocks_var = blocks_var,
      covariance_var = covariance_vars[1],
      n_blocks = n_blocks,
      cov_dimension = cov_dimension
    ))
  }
  
  # Helper function to compute total size for interaction
  compute_interaction_size <- function(effect_vars) {
    size_terms <- paste0("n_", effect_vars, "s")
    return(paste(size_terms, collapse = " * "))
  }
  
  # Helper function to compute linear index for interaction
  compute_interaction_index <- function(effect_vars) {
    if (length(effect_vars) == 1) {
      return(paste0(effect_vars[1], "[i]"))
    } else if (length(effect_vars) == 2) {
      return(paste0("(", effect_vars[1], "[i] - 1) * n_", effect_vars[2], "s + ", effect_vars[2], "[i]"))
    } else {
      # For 3+ way interactions, build nested index calculation
      index_parts <- c()
      multiplier <- 1
      
      for (j in length(effect_vars):2) {
        if (j == length(effect_vars)) {
          index_parts <- c(paste0("(", effect_vars[j], "[i] - 1)"), index_parts)
        } else {
          index_parts <- c(paste0("(", effect_vars[j], "[i] - 1) * n_", effect_vars[j+1], "s"), index_parts)
        }
      }
      
      # Build the full index expression
      index_expr <- paste(index_parts, collapse = " + ")
      index_expr <- paste0(index_expr, " + ", effect_vars[1], "[i]")
      
      return(index_expr)
    }
  }
  
  # Helper function to sanitize effect name notation
  sanitize_effect_name <- function(effect_name) {
    gsub(":", "_", effect_name)
  }
  
  # Helper function to generate Stan model code
  generate_stan_code <- function(random_effects, fixed_effects, model_vars, parsed_priors, cycle_data) {
    
    #---------------- DATA BLOCK ----------------
    
    data_block <- "data {\n\n  int<lower=0> N;\n  vector[N] y;\n\n"
    
    for (var in model_vars) {
      if (var != "yield") {
        data_block <- paste0(data_block,
                             "  int<lower=1> ", var, "[N];\n",
                             "  int<lower=1> n_", var, "s;\n\n")
      }
    }
    
    data_block <- paste0(data_block,
                         "  int<lower=1> environment[N];\n",
                         "  int<lower=1> n_environments;\n\n")
    
    for (effect_name in names(random_effects)) {
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_gamma") {
        data_block <- paste0(data_block,
                             "  real<lower=0> alpha_var_", stan_effect_name, ";\n",
                             "  real<lower=0> beta_var_", stan_effect_name, ";\n\n")
      } else if (prior_spec$type == "inverse_wishart") {
        wishart_struct <- get_wishart_structure(effect_name, random_effects[[effect_name]], prior_spec, cycle_data)
        data_block <- paste0(data_block,
                             "  cov_matrix[", wishart_struct$cov_dimension, "] psi_", stan_effect_name, ";\n",
                             "  real<lower=", wishart_struct$cov_dimension, "> nu_", stan_effect_name, ";\n\n")
      }
    }
    
    data_block <- paste0(data_block,
                         "  real<lower=0> alpha_var_resid_env;\n",
                         "  real<lower=0> beta_var_resid_env;\n",
                         "}\n\n")
    
    #---------------- PARAMETERS BLOCK ----------------
    
    params_block <- "parameters {\n\n"
    
    if (length(fixed) == 0) {
      
      params_block <- paste0(params_block, "  real beta_0;\n\n")
      
    } else {
      
      params_block <- paste0(params_block, "  real beta_0;\n\n")
      
      for (effect in fixed_effects) {
        
        vars <- effect$vars
        
        if (length(vars) == 1) {
          params_block <- paste0(params_block,
                                 "  vector[n_", vars[1], "s - 1] beta_", vars[1], ";\n\n")
        } else {
          dim_sizes <- paste0("n_", vars, "s - 1", collapse = ", ")
          name <- paste(vars, collapse = "_")
          
          params_block <- paste0(params_block,
                                 "  array[", dim_sizes, "] real beta_", name, ";\n\n")
        }
        
      }
      
    }
    
    params_block <- paste0(params_block, "\n")
    
    for (effect_name in names(random_effects)) {
      
      effect <- random_effects[[effect_name]]
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_wishart") {
        
        wishart_struct <- get_wishart_structure(effect_name, effect, prior_spec, cycle_data)
        
        params_block <- paste0(params_block,
                               "  matrix[n_", wishart_struct$blocks_var, "s, n_", wishart_struct$covariance_var, "s] u_", stan_effect_name, ";\n",
                               "  cov_matrix[n_", wishart_struct$covariance_var, "s] Sigma_", stan_effect_name, ";\n\n")
        
      } else {
        
        vars <- effect$vars
        
        if (effect$type == "simple") {
          
          params_block <- paste0(params_block,
                                 "  vector[n_", vars[1], "s] u_", stan_effect_name, ";\n\n")
          
        } else if (length(vars) == 2) {
          
          params_block <- paste0(params_block,
                                 "  matrix[n_", vars[1], "s, n_", vars[2], "s] u_", stan_effect_name, ";\n\n")
          
        } else {
          
          last_var <- vars[length(vars)]
          other_vars <- vars[-length(vars)]
          
          row_size <- paste0("(", paste0("n_", other_vars, "s", collapse = " * "), ")")
          col_size <- paste0("n_", last_var, "s")
          
          params_block <- paste0(params_block,
                                 "  matrix[", row_size, ", ", col_size, "] u_", stan_effect_name, ";\n\n")
        }
      }
    }
    
    for (effect_name in names(random_effects)) {
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_gamma") {
        params_block <- paste0(params_block,
                               "  real<lower=0> var_", stan_effect_name, ";\n")
      }
    }
    
    params_block <- paste0(params_block,
                           "\n  vector<lower=0>[n_environments] var_resid_env;\n",
                           "}\n\n")
    
    #---------------- TRANSFORMED PARAMETERS BLOCK ----------------
    
    tp_block <- "transformed parameters {\n\n  vector[N] mu;\n  vector[n_environments] sigma_resid_env = sqrt(var_resid_env);\n\n"
    
    tp_block <- paste0(tp_block, "  for (i in 1:N) {\n    mu[i] =\n\n")
    
    # Fixed effects
    fixed_terms <- c()
    
    if (length(fixed) == 0) {
      
      fixed_terms <- c(fixed_terms, "beta_0")
      
    } else {
      
      fixed_terms <- c(fixed_terms, "beta_0")
      
      for (effect in fixed_effects) {
        
        vars <- effect$vars
        
        if (length(vars) == 1) {
          term <- paste0("(", vars[1], "[i] == 1 ? 0 : beta_", vars[1], "[", vars[1], "[i] - 1])")
        } else {
          name <- paste(vars, collapse = "_")
          zero_condition <- paste0(vars, "[i] == 1", collapse = " || ")
          index_terms <- paste0(vars, "[i] - 1", collapse = ", ")
          
          term <- paste0("(", zero_condition, " ? 0 : beta_", name, "[", index_terms, "])")
        }
        
        fixed_terms <- c(fixed_terms, term)
        
      }
      
    }
    
    if (length(fixed_terms) > 0) {
      tp_block <- paste0(tp_block, "      // Fixed effects\n")
      for (term in fixed_terms) {
        tp_block <- paste0(tp_block, "      ", term, " +\n")
      }
      tp_block <- paste0(tp_block, "\n\n")
    }
    
    # Random effects
    random_terms <- c()
    
    for (effect_name in names(random_effects)) {
      
      effect <- random_effects[[effect_name]]
      vars <- effect$vars
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (effect$type == "simple") {
        term <- paste0("u_", stan_effect_name, "[", vars[1], "[i]]")
      } else if (length(vars) == 2) {
        term <- paste0("u_", stan_effect_name, "[", vars[1], "[i], ", vars[2], "[i]]")
      } else {
        last_var <- vars[length(vars)]
        other_vars <- vars[-length(vars)]
        
        row_index <- paste0("(", other_vars[1], "[i] - 1)")
        if (length(other_vars) > 1) {
          for (j in 2:length(other_vars)) {
            row_index <- paste0("(", row_index, ") * n_", other_vars[j], "s + (", other_vars[j], "[i] - 1)")
          }
        }
        row_index <- paste0("(", row_index, ") + 1")
        
        term <- paste0("u_", stan_effect_name, "[", row_index, ", ", last_var, "[i]]")
      }
      
      random_terms <- c(random_terms, term)
    }
    
    if (length(random_terms) > 0) {
      tp_block <- paste0(tp_block, "      // Random effects\n")
      for (j in seq_along(random_terms)) {
        tp_block <- paste0(tp_block,
                           "      ", random_terms[j],
                           ifelse(j < length(random_terms), " +\n", ";\n"))
      }
    }
    
    tp_block <- paste0(tp_block, "  }\n}\n\n")
    
    #---------------- MODEL BLOCK ----------------
    
    model_block <- "model {\n\n"
    
    # Fixed effects priors - not used but potentially more stable!
    #
    #if (length(fixed) == 0) {
    #  
    #  model_block <- paste0(model_block, "  beta_0 ~ normal(0, 2);\n")
    #  
    #} else {
    #  
    #  model_block <- paste0(model_block, "  beta_0 ~ normal(0, 2);\n")
    #  
    #  for (effect in fixed_effects) {
    #    
    #    vars <- effect$vars
    #    name <- paste(vars, collapse = "_")
    #    model_block <- paste0(model_block,
    #                          "  to_array_1d(beta_", name, ") ~ normal(0, 2);\n")
    #    
    #  }
    #  
    #}
    
    # Random effects priors
    for (effect_name in names(random_effects)) {
      
      effect <- random_effects[[effect_name]]
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_wishart") {
        model_block <- paste0(model_block,
                              "  for (i in 1:rows(u_", stan_effect_name, "))\n",
                              "    u_", stan_effect_name, "[i] ~ multi_normal(rep_vector(0, cols(u_", stan_effect_name, ")), Sigma_", stan_effect_name, ");\n")
      } else {
        model_block <- paste0(model_block,
                              "  to_vector(u_", stan_effect_name, ") ~ normal(0, sqrt(var_", stan_effect_name, "));\n")
      }
    }
    
    model_block <- paste0(model_block, "\n")
    
    # Variance priors
    for (effect_name in names(random_effects)) {
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_gamma") {
        model_block <- paste0(model_block,
                              "  var_", stan_effect_name, " ~ inv_gamma(alpha_var_", stan_effect_name, ", beta_var_", stan_effect_name, ");\n")
      } else if (prior_spec$type == "inverse_wishart") {
        model_block <- paste0(model_block,
                              "  Sigma_", stan_effect_name, " ~ inv_wishart(nu_", stan_effect_name, ", psi_", stan_effect_name, ");\n")
      }
    }
    
    model_block <- paste0(model_block,
                          "\n  var_resid_env ~ inv_gamma(alpha_var_resid_env, beta_var_resid_env);\n\n")
    
    model_block <- paste0(model_block,
                          "  y ~ normal(mu, sigma_resid_env[environment]);\n")
    
    model_block <- paste0(model_block, "}\n")
    
    return(paste0(data_block, params_block, tp_block, model_block))
  }
  
  
  # Helper Function to extract Stan diagnostics
  extract_stan_diagnostics <- function(stan_fit, max_treedepth_val = 10) {
    if (is.null(stan_fit)) {
      return(list(
        converged = FALSE,
        rhat_max = NA,
        ess_min = NA,
        divergent_transitions = NA,
        max_treedepth_hits = NA,
        ebfmi_min = NA,
        warnings = "Stan fit is NULL"
      ))
    }
    
    tryCatch({
      # Extract basic diagnostics
      summary_fit <- summary(stan_fit)$summary
      rhat_vals <- summary_fit[, "Rhat"]
      ess_vals <- summary_fit[, "n_eff"]
      
      # Get sampler diagnostics
      sampler_params <- get_sampler_params(stan_fit, inc_warmup = FALSE)
      divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
      max_td_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= max_treedepth_val)))
      
      # Energy diagnostics
      ebfmi_vals <- sapply(sampler_params, function(x) {
        energy <- x[, "energy__"]
        numer <- sum(diff(energy)^2) / length(energy)
        denom <- var(energy)
        numer / denom
      })
      
      list(
        converged = all(rhat_vals < 1.1, na.rm = TRUE) && all(ess_vals > 400, na.rm = TRUE) && divergent == 0,
        rhat_max = max(rhat_vals, na.rm = TRUE),
        ess_min = min(ess_vals, na.rm = TRUE),
        divergent_transitions = divergent,
        max_treedepth_hits = max_td_hits,
        ebfmi_min = min(ebfmi_vals),
        warnings = if(divergent > 0) "Divergent transitions detected" else NULL
      )
    }, error = function(e) {
      list(
        converged = FALSE,
        rhat_max = NA,
        ess_min = NA,
        divergent_transitions = NA,
        max_treedepth_hits = NA,
        ebfmi_min = NA,
        warnings = paste("Error extracting diagnostics:", e$message)
      )
    })
  }
  
  # Helper Function to extract posterior means for rerun starting values
  extract_posterior_means <- function(post_samples) {
    if (is.null(post_samples)) return(NULL)
    
    means <- list()
    for (param_name in names(post_samples)) {
      param_data <- post_samples[[param_name]]
      if (is.matrix(param_data) || is.array(param_data)) {
        # Take mean over first dimension (iterations), keep remaining dimensions
        n_dims <- length(dim(param_data))
        if (n_dims > 1) {
          margin_dims <- 2:n_dims  # All dimensions except the first
          means[[param_name]] <- apply(param_data, margin_dims, mean, na.rm = TRUE)
        } else {
          means[[param_name]] <- mean(param_data, na.rm = TRUE)
        }
      } else {
        means[[param_name]] <- mean(param_data, na.rm = TRUE)
      }
    }
    return(means)
  }
  
  # Helper function to fit model on a cycle
  fit_cycle_model <- function(cycle_data, cycle_num, updated_priors = NULL) {
    # Prepare Stan data
    stan_data <- list(
      N = nrow(cycle_data),
      y = cycle_data$yield,
      environment = as.integer(cycle_data$environment),
      n_environments = length(unique(cycle_data$environment))
    )
    
    # Add model variables to stan_data
    for (var in model_vars) {
      if (var != "yield") {
        stan_data[[var]] <- as.integer(cycle_data[[var]])
        stan_data[[paste0("n_", var, "s")]] <- length(unique(cycle_data[[var]]))
      }
    }
    
    # Set priors based on prior type
    for (effect_name in names(random_effects)) {
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_gamma") {
        if (!is.null(updated_priors) && paste0("alpha_var_", stan_effect_name) %in% names(updated_priors)) {
          stan_data[[paste0("alpha_var_", stan_effect_name)]] <- updated_priors[[paste0("alpha_var_", stan_effect_name)]]
          stan_data[[paste0("beta_var_", stan_effect_name)]] <- updated_priors[[paste0("beta_var_", stan_effect_name)]]
        } else {
          # Use values from parsed priors
          stan_data[[paste0("alpha_var_", stan_effect_name)]] <- prior_spec$alpha
          stan_data[[paste0("beta_var_", stan_effect_name)]] <- prior_spec$beta
        }
      } else if (prior_spec$type == "inverse_wishart") {
        if (!is.null(updated_priors) && paste0("nu_", stan_effect_name) %in% names(updated_priors)) {
          stan_data[[paste0("nu_", stan_effect_name)]] <- updated_priors[[paste0("nu_", stan_effect_name)]]
          stan_data[[paste0("psi_", stan_effect_name)]] <- updated_priors[[paste0("psi_", stan_effect_name)]]
        } else {
          # Set default values
          wishart_struct <- get_wishart_structure(effect_name, random_effects[[effect_name]], prior_spec, cycle_data)
          
          if (prior_spec$nu == "default") {
            stan_data[[paste0("nu_", stan_effect_name)]] <- wishart_struct$cov_dimension + 1
          } else {
            stan_data[[paste0("nu_", stan_effect_name)]] <- prior_spec$nu
          }
          
          if (is.character(prior_spec$psi) && identical(prior_spec$psi, "identity")) {
            psi_matrix <- diag(1, wishart_struct$cov_dimension)
          } else if (is.character(prior_spec$psi) && identical(prior_spec$psi, "corr")) {
            psi_matrix <- matrix(0.4, wishart_struct$cov_dimension, wishart_struct$cov_dimension)
            diag(psi_matrix) <- 0.5
          } else {
            psi_matrix <- prior_spec$psi
          }
          stan_data[[paste0("psi_", stan_effect_name)]] <- psi_matrix
        }
      }
    }
    
    # Residual environment priors
    residual_spec <- parsed_priors[["residual"]]
    if (!is.null(updated_priors)) {
      stan_data$alpha_var_resid_env <- updated_priors$alpha_var_resid_env_mean
      stan_data$beta_var_resid_env <- updated_priors$beta_var_resid_env_mean
    } else {
      stan_data$alpha_var_resid_env <- residual_spec$alpha
      stan_data$beta_var_resid_env <- residual_spec$beta
    }
    
    # Fit the model
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    set.seed(seed)
    
    # Get cycle-specific MCMC parameters
    cycle_iter <- iter[cycle_num]
    cycle_warmup <- warmup[cycle_num]
    cycle_delta <- delta[cycle_num]
    cycle_stepsize <- stepsize[cycle_num]
    cycle_max_treedepth <- max_treedepth[cycle_num]
    
    # Trace computing time
    start_time <- Sys.time()
    
    cycle_warnings <- list()
    
    stan_fit <- withCallingHandlers(
      sampling(stan_model, data = stan_data,
               iter = cycle_iter,
               chains = chains,
               warmup = cycle_warmup,
               thin = thin,
               refresh = 1000,
               control = list(adapt_delta = cycle_delta,
                              stepsize = cycle_stepsize,
                              max_treedepth = cycle_max_treedepth),
               cores = getOption("mc.cores", 1L),
               seed = seed,
               init = init),
      #open_progress = TRUE),
      warning = function(w) {
        cycle_warnings <<- c(cycle_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")  # prevents printing to console
      }
    )
    
    # Trace computing time
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    
    # Extract posterior samples
    post_samples <- rstan::extract(stan_fit)
    
    return(list(fit = stan_fit, 
                post_samples = post_samples,
                time = elapsed_time,
                warnings = cycle_warnings
    ))
  }
  
  # Helper function to compute posterior parameters
  compute_params <- function(post_sample, parsed_priors, cycle_num = 1) {
    
    # Inverse Gamma fitting function
    neg_log_likelihood_invgamma_transformed <- function(theta, data) {
      alpha <- 2 + exp(theta[1])   # ensures finite variance
      beta <- exp(theta[2])
      
      n <- length(data)
      term1 <- n * (alpha * log(beta) - lgamma(alpha))
      term2 <- - (alpha + 1) * sum(log(data))
      term3 <- - beta * sum(1 / data)
      
      return(- (term1 + term2 + term3))
    }
    
    # Precomputation function
    precompute_iwishart_data <- function(data_array) {
      n <- dim(data_array)[1]
      p <- dim(data_array)[2]
      
      log_dets <- numeric(n)
      inv_mats <- vector("list", n)
      
      for (i in 1:n) {
        X <- data_array[i,,]
        
        chol_X <- tryCatch(chol(X), error = function(e) return(NULL))
        if (is.null(chol_X)) return(NULL)
        
        log_dets[i] <- 2 * sum(log(diag(chol_X)))
        inv_mats[[i]] <- chol2inv(chol_X)
      }
      
      list(log_dets = log_dets, inv_mats = inv_mats, n = n, p = p)
    }
    
    
    # Inverse Wishart fitting function
    log_likelihood_iwishart_cholesky <- function(params, data) {
      
      # Extract precomputed values
      log_dets <- data$log_dets
      inv_mats <- data$inv_mats
      n <- data$n
      p <- data$p
      
      # Degrees of freedom (nu)
      nu <- (p + 3) + exp(params[1])
      
      # Reconstruct the Cholesky factor L
      L <- matrix(0, nrow = p, ncol = p)
      L[lower.tri(L, diag = TRUE)] <- params[-1]
      diag(L) <- exp(diag(L))
      
      # Scale matrix psi = LL^T
      psi <- L %*% t(L)
      
      # Regularize psi
      epsilon <- 1e-6
      psi <- psi + epsilon * diag(p)
      
      # Multivariate gamma
      log_gamma_p <- function(a, p) {
        if (a <= 0) return(1e10)
        sum(lgamma(a - (0:(p-1)) / 2)) + (p * (p - 1) / 4) * log(pi)
      }
      
      # Log-determinant of psi
      log_det_psi <- 2 * sum(log(diag(L)))
      if (!is.finite(log_det_psi)) return(1e10)
      
      # Base likelihood
      log_likelihood <- n * (nu / 2 * log_det_psi -
                               (nu * p / 2) * log(2) -
                               log_gamma_p(nu / 2, p))
      
      # Vectorized likelihood contributions
      log_likelihood <- log_likelihood -
        ((nu + p + 1) / 2) * sum(log_dets)
      
      # Trace term (vectorized)
      trace_term <- sum(sapply(inv_mats, function(Xi) sum(psi * Xi)))
      
      log_likelihood <- log_likelihood - 0.5 * trace_term
      
      return(-log_likelihood)
    }
    
    result <- list()
    
    # Trace computing time
    start_time <- Sys.time()
    
    # Process each random effect based on its prior type
    for (effect_name in names(parsed_priors)) {
      if (effect_name == "residual") next  # Skip residual, handle separately
      
      prior_spec <- parsed_priors[[effect_name]]
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      if (prior_spec$type == "inverse_gamma") {
        var_name <- paste0("var_", stan_effect_name)
        sample_var <- post_sample[[var_name]]
        
        ## Using two stage approach:
        
        n_bins <- max(10, as.integer(length(sample_var)/100))
        
        ## Idea 1: ML on equal-probability histogram bins
        
        # Compute quantile-based breaks
        breaks <- quantile(sample_var, probs = seq(0, 1, length.out = n_bins + 1))
        # Plot histogram with these breaks
        simple_var <- (hist(sample_var, breaks = breaks, plot=FALSE))$mids
        
        # Initial guess
        start_points <- list(
          c(log(0.2), log(0.1)),
          c(log(0.3), log(0.2)),
          c(log(0.5), log(0.3)),
          c(log(1.0), log(0.5)),
          c(log(2.0), log(1.0)),
          c(log(3.0), log(1.5)),
          c(log(4.0), log(2.0))
        )
        
        best_fit <- NULL
        best_loglik <- Inf
        
        # Bounds
        lower_bounds <- c(log(1e-5), log(1e-5))
        upper_bounds <- c(log(1.5e2), log(1.5e1))
        
        for (start in start_points) {
          fit_var <- tryCatch({
            optim(par = start,
                  fn = neg_log_likelihood_invgamma_transformed,
                  data = simple_var,
                  method = "L-BFGS-B",
                  lower = lower_bounds,
                  upper = upper_bounds)
          }, error = function(e) NULL)
          
          if (!is.null(fit_var) && fit_var$value < best_loglik) {
            best_fit <- fit_var
            best_loglik <- fit_var$value
          }
        }
        
        alpha_est <- 2 + exp(best_fit$par[1])
        beta_est <- exp(best_fit$par[2])
        
        # Slow down learning across cycles (Prior tempering) - not used but potentially more stable!
        #w_cycle = 1 / sqrt(cycle_num)
        #alpha_est =  (w_cycle*alpha_est) + (1-w_cycle)*parsed_priors[[effect_name]]$alpha
        #beta_est = (w_cycle*beta_est) + (1-w_cycle)*parsed_priors[[effect_name]]$beta
        
        result[[paste0("alpha_var_", stan_effect_name)]] <- alpha_est
        result[[paste0("beta_var_", stan_effect_name)]] <- beta_est
        
        cat("Estimated shape (alpha_var_", stan_effect_name, "):", alpha_est, "\n")
        cat("Estimated scale (beta_var_", stan_effect_name, "):", beta_est, "\n")
        
      } else if (prior_spec$type == "inverse_wishart") {
        sigma_name <- paste0("Sigma_", stan_effect_name)
        sample_sigma <- post_sample[[sigma_name]]
        
        p <- nrow(sample_sigma[1,,]) 
        
        # Precompute once
        precomp <- precompute_iwishart_data(sample_sigma)
        if (is.null(precomp)) {
          cat("Precomputation failed: non-PD matrix detected\n")
          return(NULL)
        }
        
        # Moment-based initials and bounds
        n_scale_params <- p * (p + 1) / 2
        nu_init <- p + 5
        params_init <- c(log(nu_init - (p + 3)), rep(0.1, n_scale_params))
        
        lower_bounds <- c(log(1e-6), rep(-10, n_scale_params))
        upper_bounds <- c(log(200), rep(10, n_scale_params))
        
        fit_sigma <- tryCatch({
          optim(par = params_init, 
                fn = log_likelihood_iwishart_cholesky, 
                data = precomp,
                method = "L-BFGS-B", 
                lower = lower_bounds,
                upper = upper_bounds,
                control = list(fnscale = 1, maxit = 2000))
        }, error = function(e) {
          cat("Inv Wishart Optimization failed:", e$message, "\n")
          return(list(par = params_init, convergence = 1, value = Inf))
        })
        
        nu_est <- (p + 3) + exp(fit_sigma$par[1])
        
        L_est <- matrix(0, nrow = p, ncol = p)
        L_est[lower.tri(L_est, diag = TRUE)] <- fit_sigma$par[-1]
        diag(L_est) <- exp(diag(L_est))
        
        psi_est <- L_est %*% t(L_est)
        
        # Slow down learning across cycles (Prior tempering) - not used but potentially more stable!
        #w_cycle = 1 / sqrt(cycle_num)
        #nu_est =  (w_cycle*nu_est) + (1-w_cycle)*parsed_priors[[effect_name]]$nu
        #psi_est = (w_cycle*psi_est) + (1-w_cycle)*parsed_priors[[effect_name]]$psi
        
        result[[paste0("nu_", stan_effect_name)]] <- nu_est
        result[[paste0("psi_", stan_effect_name)]] <- psi_est
        
        cat("Estimated degrees of freedom (nu_", stan_effect_name, "):", nu_est, "\n")
        cat("Estimated scale matrix (psi_", stan_effect_name, "):\n")
        print(psi_est)
      }
    }
    
    # Process residual environment variance separately as inverse gammas
    sample_var_resid_env <- post_sample$var_resid_env
    
    # Peparing lists
    alpha_est_resid_env <- vector("list", ncol(sample_var_resid_env))
    beta_est_resid_env <- vector("list", ncol(sample_var_resid_env))
    breaks <- list()
    simple_var_resid_env <- list()
    
    for (i in 1:ncol(sample_var_resid_env)) {
      
      start_points <- list(
        c(log(0.2), log(0.1)),
        c(log(0.3), log(0.2)),
        c(log(0.5), log(0.3)),
        c(log(1.0), log(0.5)),
        c(log(2.0), log(1.0)),
        c(log(3.0), log(1.5)),
        c(log(4.0), log(2.0))
      )
      
      best_fit <- NULL
      best_loglik <- Inf
      
      lower_bounds <- c(log(1e-5), log(1e-5))
      upper_bounds <- c(log(1.5e2), log(1.5e1))
      
      ## Using two stage approach:
      
      n_bins <- max(10L, as.integer(length(sample_var_resid_env[, i]) / 100L))
      
      ## Idea 1: ML on equal-probability histogram bins
      
      # Compute quantile-based breaks
      breaks[[i]] <- quantile(sample_var_resid_env[,i], probs = seq(0, 1, length.out = n_bins + 1))
      
      # Plot histogram with these breaks
      simple_var_resid_env[[i]] <- (hist(sample_var_resid_env[,i], breaks = breaks[[i]], plot=FALSE))$mids
      
      for (start in start_points) {
        fit_var <- tryCatch({
          optim(par = start,
                fn = neg_log_likelihood_invgamma_transformed,
                data = simple_var_resid_env[[i]],
                method = "L-BFGS-B",
                lower = lower_bounds,
                upper = upper_bounds)
        }, error = function(e) NULL)
        
        if (!is.null(fit_var) && fit_var$value < best_loglik) {
          best_fit <- fit_var
          best_loglik <- fit_var$value
        }
      }
      
      alpha_est_resid_env[i] <- 2 + exp(best_fit$par[1])
      beta_est_resid_env[i] <- exp(best_fit$par[2])
      
    }
    
    # Slow down learning across cycles (Prior tempering) - not used but potentially more stable!
    #w_cycle = 1 / sqrt(cycle_num)
    #alpha_var_resid_env_mean =  (w_cycle*mean(unlist(alpha_est_resid_env))) + (1-w_cycle)*parsed_priors[["residual"]]$alpha
    #beta_var_resid_env_mean = (w_cycle*mean(unlist(beta_est_resid_env))) + (1-w_cycle)*parsed_priors[["residual"]]$beta
    
    result$alpha_var_resid_env_mean <- mean(unlist(alpha_est_resid_env)) #alpha_var_resid_env_mean
    result$beta_var_resid_env_mean <- mean(unlist(beta_est_resid_env)) #beta_var_resid_env_mean
    
    # Trace computing time
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    
    cat("Estimated shape (alpha_var_resid_env_mean):", result$alpha_var_resid_env_mean, "\n")
    cat("Estimated scale (beta_var_resid_env_mean):", result$beta_var_resid_env_mean, "\n")
    
    return(list(result = result,
                time = elapsed_time))
  }
  
  
  ## Section 2: BAYESIAN UPDATING --------------------------------------------
  
  # Parse main function arguments
  random_effects <- parse_effects(random)
  parsed_priors <- parse_priors(priors, random_effects)
  fixed_effects <- parse_effects(fixed)
  model_vars <- get_model_variables(random_effects, fixed_effects)
  
  # Check if all required variables exist in data
  missing_vars <- setdiff(model_vars, names(MET_data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables in MET_data:", paste(missing_vars, collapse = ", ")))
  }
  
  #Yield target variable
  MET_data$yield <- as.numeric(MET_data$yield)
  
  # Process historical data cycles
  cycles_data <- process_cycles(MET_data, cycles)
  n_cycles <- length(cycles_data)
  
  # Set default or validate MCMC parameters per cycle
  if (is.null(iter)) {
    iter <- rep(2000, n_cycles)
  } else if (length(iter) == 1) {
    iter <- rep(iter, n_cycles)
  } else if (length(iter) != n_cycles) {
    stop("Length of 'iter' must be 1 or match the number of cycles")
  }
  
  if (is.null(warmup)) {
    warmup <- round(iter * 0.6)  # Default to 60% of iterations
  } else if (length(warmup) == 1) {
    warmup <- rep(warmup, n_cycles)
  } else if (length(warmup) != n_cycles) {
    stop("Length of 'warmup' must be 1 or match the number of cycles")
  }
  
  if (is.null(delta)) {
    delta <- rep(0.8, n_cycles)
  } else if (length(delta) == 1) {
    delta <- rep(delta, n_cycles)
  } else if (length(delta) != n_cycles) {
    stop("Length of 'delta' must be 1 or match the number of cycles")
  }
  
  if (is.null(stepsize)) {
    stepsize <- rep(0.1, n_cycles)
  } else if (length(stepsize) == 1) {
    stepsize <- rep(stepsize, n_cycles)
  } else if (length(stepsize) != n_cycles) {
    stop("Length of 'stepsize' must be 1 or match the number of cycles")
  }
  
  if (is.null(max_treedepth)) {
    max_treedepth <- rep(10, n_cycles)
  } else if (length(max_treedepth) == 1) {
    max_treedepth <- rep(max_treedepth, n_cycles)
  } else if (length(max_treedepth) != n_cycles) {
    stop("Length of 'max_treedepth' must be 1 or match the number of cycles")
  }
  
  # Process each cycle data
  for (i in 1:length(cycles_data)) {
    cat("Processing cycle", i, "data...\n")
    
    # Only keep needed variables plus environment
    needed_vars <- c(model_vars, "environment")
    missing_vars_in_cycle <- setdiff(needed_vars, names(cycles_data[[i]]))
    if (length(missing_vars_in_cycle) > 0) {
      stop(paste("Cycle", i, "is missing variables:", paste(missing_vars_in_cycle, collapse = ", ")))
    }
    
    cycles_data[[i]] <- cycles_data[[i]][, needed_vars, drop = FALSE]
    
    # Check for empty data
    if (nrow(cycles_data[[i]]) == 0) {
      stop(paste("Cycle", i, "has no data after subsetting"))
    }
    
    # Encode categorical variables (exclude year_nr and yield)
    for (var in model_vars) {
      if (var != "yield" && var != "year_nr") {
        original_length <- length(cycles_data[[i]][[var]])
        cycles_data[[i]][[var]] <- as.integer(factor(droplevels(cycles_data[[i]][[var]])))
        new_length <- length(cycles_data[[i]][[var]])
        
        # Check for dimension consistency
        if (original_length != new_length) {
          stop(paste("Dimension mismatch in variable", var, "for cycle", i))
        }
        
        # Check for valid encoding
        if (any(is.na(cycles_data[[i]][[var]])) || any(cycles_data[[i]][[var]] < 1)) {
          stop(paste("Invalid encoding for variable", var, "in cycle", i))
        }
        
        cat("Variable", var, "encoded with", max(cycles_data[[i]][[var]]), "levels\n")
      }
    }
    
    # Encode environment variable
    original_env_length <- length(cycles_data[[i]]$environment)
    cycles_data[[i]]$environment <- as.integer(factor(cycles_data[[i]]$environment))
    
    if (length(cycles_data[[i]]$environment) != original_env_length) {
      stop(paste("Environment encoding failed for cycle", i))
    }
    
    # Reorder data - fix the ordering issue
    if (ncol(cycles_data[[i]]) > 0) {
      first_var <- names(cycles_data[[i]])[1]
      cycles_data[[i]] <- cycles_data[[i]][order(cycles_data[[i]][[first_var]]), , drop = FALSE]
    }
    
    cat("Cycle", i, "processed:", nrow(cycles_data[[i]]), "observations\n")
  }
  
  # Generate Stan model code (pass first cycle data for structure determination)
  stan_code <- generate_stan_code(random_effects, fixed_effects, model_vars, parsed_priors, cycles_data[[1]])
  
  # Compile Stan model
  stan_model <- stan_model(model_code = stan_code)
  
  # Initialize results storage with better structure
  cycle_results <- list()
  convergence_summary <- list()
  parameter_evolution <- list()
  
  # Iterate over cycles (modified to collect more diagnostics)
  for (i in 1:length(cycles_data)) {
    cat("Fitting cycle", i, "of", length(cycles_data), "...\n")
    flush.console()
    
    if (i == 1) {
      current_fit <- fit_cycle_model(cycles_data[[i]], cycle_num = i)
    } else {
      previous_post <- cycle_results[[i - 1]]$posterior_samples
      updated_priors <- compute_params(previous_post, parsed_priors, cycle_num = i)
      current_fit <- fit_cycle_model(cycles_data[[i]], cycle_num = i, updated_priors$result)
    }
    
    # Extract diagnostics
    stan_diagnostics <- extract_stan_diagnostics(current_fit$fit, max_treedepth_val = max_treedepth[i])
    posterior_means <- extract_posterior_means(current_fit$post_samples)
    
    # Store cycle results
    if (i == 1) {
      cycle_results[[i]] <- list(
        stan_fit = current_fit$fit,
        posterior_samples = current_fit$post_samples,
        posterior_means = posterior_means,
        cycle_time = current_fit$time,
        priors_used = if(i == 1) parsed_priors else updated_priors$result,
        diagnostics = stan_diagnostics,
        mcmc_settings = list(
          iter = iter[i],
          warmup = warmup[i],
          delta = delta[i],
          stepsize = stepsize[i],
          max_treedepth = max_treedepth[i]
        ),
        data_summary = list(
          n_obs = nrow(cycles_data[[i]]),
          n_environments = length(unique(cycles_data[[i]]$environment)),
          years_included = if("year" %in% names(cycles_data[[i]])) 
            unique(cycles_data[[i]]$year) else NULL
        )
      )
    } else {
      cycle_results[[i]] <- list(
        stan_fit = current_fit$fit,
        posterior_samples = current_fit$post_samples,
        posterior_means = posterior_means,
        cycle_time = current_fit$time,
        post_parm_time = updated_priors$time,
        priors_used = if(i == 1) parsed_priors else updated_priors$result,
        diagnostics = stan_diagnostics,
        mcmc_settings = list(
          iter = iter[i],
          warmup = warmup[i],
          delta = delta[i],
          stepsize = stepsize[i],
          max_treedepth = max_treedepth[i]
        ),
        data_summary = list(
          n_obs = nrow(cycles_data[[i]]),
          n_environments = length(unique(cycles_data[[i]]$environment)),
          years_included = if("year" %in% names(cycles_data[[i]])) 
            unique(cycles_data[[i]]$year) else NULL
        )
      )
    }
    
    # Track convergence
    convergence_summary[[i]] <- list(
      cycle = i,
      stan_converged = stan_diagnostics$converged,
      rhat_max = stan_diagnostics$rhat_max,
      ess_min = stan_diagnostics$ess_min,
      divergent_transitions = stan_diagnostics$divergent_transitions,
      our_warnings = stan_diagnostics$warnings,
      stan_warnings = current_fit$warnings
    )
    
    # Track parameter evolution
    if (i > 1 && !is.null(updated_priors$result)) {
      parameter_evolution[[i-1]] <- updated_priors$result
    }
    
    cat("Cycle", i, "complete! Converged:", stan_diagnostics$converged, "\n")
    if (!is.null(stan_diagnostics$warnings)) {
      cat("Warnings:", stan_diagnostics$warnings, "\n")
    }
  }
  
  # Compute final posterior parameters
  final_post_parms <- compute_params(cycle_results[[length(cycles_data)]]$posterior_samples, parsed_priors, cycle_num = length(cycles_data))
  
  
  ## Section 3: PROCESS FINAL RESULTS ----------------------------------------
  
  # Add %||% operator if not defined
  `%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x
  
  # Structure the final output
  result <- list(
    # Main Results
    cycles = cycle_results,
    
    # Parameter Information
    parameters = list(
      initial_priors = parsed_priors,
      final_posteriors = final_post_parms,
      evolution = parameter_evolution,
      posterior_means_final = extract_posterior_means(cycle_results[[length(cycles_data)]]$posterior_samples)
    ),
    
    # Convergence and Diagnostics
    diagnostics = list(
      stan_convergence = convergence_summary,
      overall_converged = all(sapply(convergence_summary, function(x) x$stan_converged)),
      total_divergent = sum(sapply(convergence_summary, function(x) x$divergent_transitions %||% 0)),
      max_rhat = max(sapply(convergence_summary, function(x) x$rhat_max %||% 1), na.rm = TRUE),
      min_ess = min(sapply(convergence_summary, function(x) x$ess_min %||% Inf), na.rm = TRUE)
    ),
    
    # Model Specification
    model_spec = list(
      fixed_effects = fixed_effects,
      random_effects = random_effects,
      stan_model = stan_code,
      cycles_specification = cycles,
      n_cycles_fitted = length(cycles_data),
      total_observations = sum(sapply(cycles_data, nrow)),
      mcmc_settings = data.frame(
        cycle = 1:length(cycles_data),
        iter = iter,
        warmup = warmup,
        delta = delta,
        stepsize = stepsize,
        max_treedepth = max_treedepth
      )
    ),
    
    # Rerun Information (for pipeline continuation)
    rerun_info = list(
      successful_cycles = which(sapply(convergence_summary, function(x) x$stan_converged)),
      failed_cycles = which(!sapply(convergence_summary, function(x) x$stan_converged)),
      suggested_starting_values = final_post_parms,
      suggested_mcmc_adjustments = {
        failed_idx <- which(!sapply(convergence_summary, function(x) x$stan_converged))
        if (length(failed_idx) > 0) {
          list(
            increase_adapt_delta = any(sapply(convergence_summary[failed_idx], function(x) (x$divergent_transitions %||% 0) > 0)),
            increase_max_treedepth = any(sapply(cycle_results[failed_idx], function(x) (x$diagnostics$max_treedepth_hits %||% 0) > 0)),
            increase_iterations = any(sapply(convergence_summary[failed_idx], function(x) (x$ess_min %||% Inf) < 400))
          )
        } else {
          NULL
        }
      }
    ),
    
    # Session Information
    session_info = list(
      call = match.call(),
      timestamp = Sys.time(),
      seed = seed,
      r_version = R.version.string,
      package_versions = list(
        rstan = as.character(packageVersion("rstan")),
        asreml = if(requireNamespace("asreml", quietly = TRUE)) as.character(packageVersion("asreml")) else "not available"
      )
    )
  )
  
  class(result) <- c("fit_general_multiple_cycles", "list")
  return(result)
  
}

## APPENDIX WITH TO HELPER FUNCTION TO PROCESS FINAL RESULT OBJECT -------------

# Add a simple print method for better output display
print.fit_general_multiple_cycles <- function(x, ...) {
  cat("Multi-Cycle Bayesian Analysis Results\n")
  cat("=====================================\n")
  cat("Cycles fitted:", x$model_spec$n_cycles_fitted, "\n")
  cat("Total observations:", x$model_spec$total_observations, "\n")
  cat("Overall convergence:", x$diagnostics$overall_converged, "\n")
  if (isTRUE(x$diagnostics$total_divergent > 0)) {
    cat("Warning:", x$diagnostics$total_divergent, "divergent transitions detected\n")
  }
  cat("\nUse summary() for detailed results\n")
}

# Add a summary method
summary.fit_general_multiple_cycles <- function(object, ...) {
  cat("Multi-Cycle Bayesian Analysis Summary\n")
  cat("====================================\n\n")
  
  # Model specification
  cat("Model Specification:\n")
  cat("- Fixed effects:\n", paste(" ", names(object$model_spec$fixed_effects), collapse = "\n"), "\n")
  cat("- Random effects:\n", paste(" ", names(object$model_spec$random_effects), collapse = "\n"), "\n")
  cat("- Cycles:", object$model_spec$n_cycles_fitted, "\n\n")
  
  # Convergence summary
  cat("Convergence Summary:\n")
  for (i in 1:length(object$diagnostics$stan_convergence)) {
    conv <- object$diagnostics$stan_convergence[[i]]
    cat(sprintf("Cycle %d: %s (Rhat max: %.3f, ESS min: %.0f, Divergent: %d)\n",
                i, ifelse(conv$stan_converged, "CONVERGED", "FAILED"),
                conv$rhat_max %||% NA, conv$ess_min %||% NA, conv$divergent_transitions %||% 0))
  }
  cat("\n")
  
  # Rerun suggestions
  if (!object$diagnostics$overall_converged && !is.null(object$rerun_info$suggested_mcmc_adjustments)) {
    cat("Rerun Suggestions:\n")
    adj <- object$rerun_info$suggested_mcmc_adjustments
    if (adj$increase_adapt_delta) cat("- Increase adapt_delta (current max divergent cycles)\n")
    if (adj$increase_max_treedepth) cat("- Increase max_treedepth\n")
    if (adj$increase_iterations) cat("- Increase iterations\n")
  }
}
