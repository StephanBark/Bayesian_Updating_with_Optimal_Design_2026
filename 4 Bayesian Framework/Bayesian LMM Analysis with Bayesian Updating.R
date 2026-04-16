####### Bayesian LMM and Optimal Design Analysis with Bayesian Updating ########


### Bayesian Updating via MCMC Methods -----------------------------------------

## NOTE: The user is invited to adjust the model specification and hyper-parameters in the scope of the function call below.
##       Individual HTML reports can be generated after running the sampler.
##       The current setup represents the exact as presented in our paper.

cycle_fits_US_zones_model <- fit_general_multiple_cycles(
  
  ## Data:
  yield_winter_medium,
  
  ## Bayesian linear mixed model specification:
  # Fixed effects
  fixed = c("Zone"),
  # Random effects
  random = c("year",
            "Zone:year",
            "Zone:Location:year", 
             "Zone:Location:Rep:year",
             "Genotype:year", 
             "Genotype:Zone", 
             "Genotype:Zone:year", 
             "Genotype:Zone:Location:year"),
  # Priors for variance components:
  # NOTE: For the residual variance (R) and each random effect variance (G),
  #       either an inverse-gamma prior or an inverse-Wishart prior for covariance matrices can be specified.
  #       G must be specified in the same order as random effects vector above.
  priors = list("R" = list(type = "inverse_gamma", alpha = 5, beta = 1),
                "G" = list(list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_wishart",
                                blocks = "Genotype",
                                nu = 10,
                                psi = `diag<-`(matrix(0.9, 4, 4), 1)
                                ),
                           list(type = "inverse_gamma", alpha = 5, beta = 1),
                           list(type = "inverse_gamma", alpha = 5, beta = 1))),
  ## Other hyper-parameters:
  # Sequence of multi-year cycles (windows) to be fitted, e.g., c(1, 2, 3) for 3 cycles with 1 year added in each cycle
  cycles = c(8, 5, 3, 3, 3),
  # MCMC total number of iterations per window
  iter = c(37500, 27500, 27500, 27500, 27500),
  # MCMC warmup iterations per window (must be less than total iter)
  warmup = c(30000, 20000, 20000, 20000, 20000),
  # Target acceptance probability for the NUTS sampler
  delta = c(0.95, 0.9, 0.9, 0.9, 0.9),
  # Initial step size for the NUTS sampler
  stepsize = c(0.1, 0.1, 0.1, 0.1, 0.1),
  # Maximum tree depth for the NUTS sampler
  max_treedepth = c(12, 10, 10, 10, 10),
  # Number of MCMC chains to run in parallel
  chains = 4,
  # Number of MCMC samples to save per chain (after thinning)
  thin = 3,
  # Random seed for reproducibility
  seed = 123,
  # Initialization method for MCMC chains
  init = 'random'
  
  )

#writeLines(cycle_fits_US_zones_model$model_spec$stan_model, con = "cycle_fits_US_zones_model.stan")
#save(cycle_fits_US_zones_model, file = "5_cycle_fits_37500_27500_27500_27500_27500_2_090_corr09_random_init_coverged_13_04_26.rda")

### HTML Diagnosis and Result Report ------------------------------------------

# REML estimates from equivalent Frequentist LMM Analysis
reference_values <- c(
  var_year = 0.0287,
  var_Zone_year = 0.0000,
  var_Zone_Location_year = 0.5610,
  var_Zone_Location_Rep_year = 0.0046,
  var_Genotype_year = 0.0246,
  var_Genotype_Zone_year = 0.0000,
  var_Genotype_Zone_Location_year = 0.2504,
  var_resid_env_mean = 0.3249,
  `Sigma_Genotype_Zone[1,1]` = 0.3121,
  `Sigma_Genotype_Zone[1,2]` = 0.2478,
  `Sigma_Genotype_Zone[1,3]` = 0.2943,
  `Sigma_Genotype_Zone[1,4]` = 0.1713,
  `Sigma_Genotype_Zone[2,2]` = 0.2157,
  `Sigma_Genotype_Zone[2,3]` = 0.2388,
  `Sigma_Genotype_Zone[2,4]` = 0.1285,
  `Sigma_Genotype_Zone[3,3]` = 0.2955,
  `Sigma_Genotype_Zone[3,4]` = 0.1559,
  `Sigma_Genotype_Zone[4,4]` = 0.1031
)

render_fit_general_multiple_cycles_report(
  fit_object = cycle_fits_US_zones_model,
  output_file = "cycle_fits_US_zones_model_report.html",
  reference_values = reference_values,
  report_title = "US zones Bayesian mixed-model report"
)


### Result Grid Plots - Posterior Histogram with MLEs --------------------------

## NOTE: This section contains the exact plots from our paper.

# Use paper font
font_add("CMU Serif", "cmunrm.ttf")
showtext_auto()

# Helper: extract posterior inverse gamma parameters for a given cycle
get_cycle_posterior_params <- function(fit_obj, cycle_index) {
  n_cycles <- length(fit_obj$cycles)
  if (cycle_index < n_cycles) {
    fit_obj$parameters$evolution[[cycle_index]]
  } else {
    fit_obj$parameters$final_posteriors$result
  }
}

# Inverse gamma density functions
dinvgamma <- function(x, alpha, beta) {
  (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x)
}

dinvgamma_stable <- function(x, alpha, beta) {
  ifelse(x <= 0, 0,
         exp(alpha*log(beta) - lgamma(alpha) -
               (alpha+1)*log(x) - beta/x))
}

# Mapping from sample index to posterior parameter name stems
var_param_names <- c(
  "year",
  "Zone_year",
  "Zone_Location_year",
  "Zone_Location_Rep_year",
  "Genotype_year",
  "Genotype_Zone_year",
  "Genotype_Zone_Location_year"
)

## All Variance Components else than Genotype x Zone ---------------------------

# Access individual fits
cycle_1_fit <- cycle_fits_US_zones_model$cycles[[1]]
# for the 2th cycle:
cycle_2_fit <- cycle_fits_US_zones_model$cycles[[2]]
# for the 3th cycle:
cycle_3_fit <- cycle_fits_US_zones_model$cycles[[3]]
# for the 4th cycle:
cycle_4_fit <- cycle_fits_US_zones_model$cycles[[4]]
# for the 5th cycle:
cycle_5_fit <- cycle_fits_US_zones_model$cycles[[5]]

# Initialize an empty list to store all samples
all_samples <- list()

# Loop through all cycles and construct each sample list dynamically
for (i in seq_along(cycle_fits_US_zones_model$cycles)) {
  cycle_fit <- get(paste0("cycle_", i, "_fit"))  # Dynamically get the cycle fit object
  
  all_samples[[i]] <- list(
    cycle_fit$posterior_samples$var_year,
    cycle_fit$posterior_samples$var_Zone_year,
    cycle_fit$posterior_samples$var_Zone_Location_year,
    cycle_fit$posterior_samples$var_Zone_Location_Rep_year,
    cycle_fit$posterior_samples$var_Genotype_year,
    cycle_fit$posterior_samples$var_Genotype_Zone_year,
    cycle_fit$posterior_samples$var_Genotype_Zone_Location_year
    
  )
}

# Define corresponding variance component names of interest
titles <- c(
  "Year Variance",
  "Year × Zone Variance",
  "Year × Zone × Location Variance",
  "Year × Zone × Location × Rep Variance",
  "Genotype × Year Variance",
  "Genotype × Zone × Year Variance",
  "Genotype × Zone × Location × Year Variance"
)

## MLEs namend
MLEs_named <- c(
  var_year                           = 0.0287,
  var_Zone_year                      = 0.0000,
  var_Zone_Location_year             = 0.5610,
  var_Zone_Location_Rep_year         = 0.0046,
  var_Genotype_year                  = 0.0246,
  var_Genotype_Zone_year             = 0.0000,
  var_Genotype_Zone_Location_year    = 0.2504,
  var_resid_env_mean                 = 0.3249
)

# Open PDF device
cairo_pdf("5_cycle_fit_37500_27500_27500_27500_27500_ALL.pdf", width = 15, height = 15)

# Initialize vector to store x-axis and y-axis limits for each variance component
x_max_values <- rep(NA, length(titles))
y_max_values <- rep(NA, length(titles))

# Loop through all cycles
for (cycle_index in seq_along(all_samples)) {
  
  samples <- all_samples[[cycle_index]]  # Extract sample list for the current cycle
  plot_list <- vector("list", length(samples))
  
  for (i in seq_along(samples)) {
    sample_data <- samples[[i]]
    
    # Extract posterior inverse gamma parameters from the fit object
    cycle_params <- get_cycle_posterior_params(cycle_fits_US_zones_model, cycle_index)
    alpha_est <- cycle_params[[paste0("alpha_var_", var_param_names[i])]]
    beta_est <- cycle_params[[paste0("beta_var_", var_param_names[i])]]
    
    # Determine x-axis limits (only set during the first cycle)
    if (cycle_index == 1) {
      x_max_values[i] <- min(10, mean(sample_data) * 2)  # Limit x_max to at most 10
    }
    
    # Compute y-axis limit based on histogram and fitted curve
    hist_data <- hist(sample_data, breaks = 10, plot = FALSE)  # Get histogram data
    hist_densities <- hist_data$counts / sum(hist_data$counts) / diff(hist_data$breaks)  # Compute density manually
    
    # Compute fitted curve density values manually (no unwanted plot!)
    x_vals <- seq(0, x_max_values[i], length.out = 1000)  # Generate x-values
    y_vals <- dinvgamma(x_vals, alpha_est, beta_est)  # Compute corresponding densities
    
    if (cycle_index == 1) {
      y_max_values[i] <- max(hist_densities, y_vals, na.rm = TRUE) * 3 # Add margin
    } 
    
    dinvgamma <- function(x, alpha, beta) {
      (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x)
    }
    
    # Build x-grid for the fitted PDF curve
    x_grid <- seq(0, x_max_values[i], length.out = 2000)
    pdf_df <- data.frame(
      x   = x_grid,
      pdf = dinvgamma_stable(x_grid, alpha_est, beta_est)
    )
    
    # Build plot
    plot_list[[i]] <- ggplot(data.frame(sample_data = sample_data), aes(x = sample_data)) +
      # Histogram scaled to density (freq = FALSE equivalent)
      geom_histogram(aes(y = after_stat(density)),
                     breaks = seq(min(sample_data), max(sample_data), length.out = 15 + 1),
                     fill = "grey80",
                     color = "white") +
      # Fitted inverse-gamma PDF
      geom_line(data = pdf_df, aes(x = x, y = pdf),
                color = "blue", linewidth = 0.7) +
      geom_vline(xintercept = MLEs_named[i],
                 color = "red", linewidth = 1.3) +
      labs(
        title = paste0(titles[i]),
        x = NULL,
        y = NULL
      ) +
      coord_cartesian(xlim = c(0, x_max_values[i]),
                      ylim = c(0, y_max_values[i])) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 20),
        axis.text.y  = element_text(size = 20),
        text = element_text(family = "CMU Serif")
      )
  }
  
  combined <- wrap_plots(plot_list, nrow = 4)
  
  # Convert patchwork to a grob so we can add shared axis titles
  g <- patchwork::patchworkGrob(combined)
  
  gridExtra::grid.arrange(
    g,
    bottom = grid::textGrob(
      "Variance Value of Yield [t/ha]",
      gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
    ),
    left = grid::textGrob(
      "Inverse Gamma Density",
      rot = 90,
      gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
    ),
    newpage = cycle_index > 1
  )
  
}

# Close PDF device
dev.off()


## Genotype x Year Variance Component ------------------------------------------

# Access individual fits, e.g., for the first cycle:
cycle_1_fit <- cycle_fits_US_zones_model$cycles[[1]]
# for the 2th cycle:
cycle_2_fit <- cycle_fits_US_zones_model$cycles[[2]]
# for the 3th cycle:
cycle_3_fit <- cycle_fits_US_zones_model$cycles[[3]]
# for the 4th cycle:
cycle_4_fit <- cycle_fits_US_zones_model$cycles[[4]]
# for the 5th cycle:
cycle_5_fit <- cycle_fits_US_zones_model$cycles[[5]]

# Initialize an empty list to store all samples
all_samples <- list()

# Loop through all cycles and construct each sample list dynamically
for (i in seq_along(cycle_fits_US_zones_model$cycles)) {
  cycle_fit <- get(paste0("cycle_", i, "_fit"))  # Dynamically get the cycle fit object
  
  all_samples[[i]] <- list(
    cycle_fit$posterior_samples$var_Genotype_year
    
  )
}

# Define corresponding variance component names of interest
titles <- c(
  "Genotype × Year Variance"
)


# Vector of MLE estimates of original NON Bayesian unstructured Mixed model
MLE <- c(0.0246)

# Open PDF device
cairo_pdf("5_cycle_fit_37500_27500_27500_27500_27500_GENOTYPE_YEAR.pdf", width = 15, height = 5)

# Initialize vector to store x-axis and y-axis limits for each variance component
x_max_values <- NA
y_max_values <- NA

plot_list <- vector("list", length(all_samples))

# Loop through all cycles
for (cycle_index in seq_along(all_samples)) {
  
  samples <- all_samples[[cycle_index]]  # Extract sample list for the current cycle
  
  sample_data <- samples[[1]]
  
  # Extract posterior inverse gamma parameters from the fit object
  cycle_params <- get_cycle_posterior_params(cycle_fits_US_zones_model, cycle_index)
  alpha_est <- cycle_params[["alpha_var_Genotype_year"]]
  beta_est <- cycle_params[["beta_var_Genotype_year"]]
  
  # Determine x-axis limits (only set during the first cycle)
  if (cycle_index == 1) {
    x_max_values <- min(10, mean(sample_data) * 2)  # Limit x_max to at most 10
  }
  
  # Compute y-axis limit based on histogram and fitted curve
  hist_data <- hist(sample_data, breaks = 10, plot = FALSE)  # Get histogram data
  hist_densities <- hist_data$counts / sum(hist_data$counts) / diff(hist_data$breaks)  # Compute density manually
  
  # Compute fitted curve density values manually
  x_vals <- seq(0, x_max_values, length.out = 1000)  # Generate x-values
  y_vals <- dinvgamma(x_vals, alpha_est, beta_est)  # Compute corresponding densities
  
  if (cycle_index == 1) {
    y_max_values <- max(hist_densities, y_vals, na.rm = TRUE) * 3 # Add margin
  } 
  
  dinvgamma <- function(x, alpha, beta) {
    (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x)
  }
  
  # Build x-grid for the fitted PDF curve
  x_grid <- seq(0, x_max_values, length.out = 2000)
  pdf_df <- data.frame(
    x   = x_grid,
    pdf = dinvgamma_stable(x_grid, alpha_est, beta_est)
  )
  
  # Build plot
  plot_list[[cycle_index]] <- ggplot(data.frame(sample_data = sample_data), aes(x = sample_data)) +
    # Histogram scaled to density (freq = FALSE equivalent)
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(min(sample_data), max(sample_data), length.out = 15 + 1),
                   fill = "grey80",
                   color = "white") +
    # Fitted inverse-gamma PDF
    geom_line(data = pdf_df, aes(x = x, y = pdf),
              color = "blue", linewidth = 0.7) +
    geom_vline(xintercept = MLE[1],
               color = "red", linewidth = 1.0) +
    labs(
      title = paste0(" Window ", cycle_index),
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(xlim = c(0, x_max_values),
                    ylim = c(0, y_max_values)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 12),
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20),
      text = element_text(family = "CMU Serif")
    )
}

combined <- wrap_plots(plot_list, ncol = 5)

# Convert patchwork to a grob so we can add shared axis titles
g <- patchwork::patchworkGrob(combined)

gridExtra::grid.arrange(
  g,
  top = grid::textGrob(
    "Year × Zone Variance",
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  bottom = grid::textGrob(
    "Variance Value of Yield [t/ha]",
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  left = grid::textGrob(
    "Inverse Gamma Density",
    rot = 90,
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  newpage = FALSE
)


# Close PDF device
dev.off()


## Year x Zone Variance Component ----------------------------------------------

# Access individual fits, e.g., for the first cycle:
cycle_1_fit <- cycle_fits_US_zones_model$cycles[[1]]
# for the 2th cycle:
cycle_2_fit <- cycle_fits_US_zones_model$cycles[[2]]
# for the 3th cycle:
cycle_3_fit <- cycle_fits_US_zones_model$cycles[[3]]
# for the 4th cycle:
cycle_4_fit <- cycle_fits_US_zones_model$cycles[[4]]
# for the 5th cycle:
cycle_5_fit <- cycle_fits_US_zones_model$cycles[[5]]

# Initialize an empty list to store all samples
all_samples <- list()

# Loop through all cycles and construct each sample list dynamically
for (i in seq_along(cycle_fits_US_zones_model$cycles)) {
  cycle_fit <- get(paste0("cycle_", i, "_fit"))  # Dynamically get the cycle fit object
  
  all_samples[[i]] <- list(
    cycle_fit$posterior_samples$var_Zone_year
    
  )
}

# Define corresponding variance component names of interest
titles <- c(
  "Year × Zone Variance"
)


# Vector of MLE estimates of original NON Bayesian unstructured Mixed model
MLE <- c(0.0000)

# Open PDF device
cairo_pdf("5_cycle_fit_37500_27500_27500_27500_27500_YEAR_ZONE.pdf", width = 15, height = 5)

# Initialize vector to store x-axis and y-axis limits for each variance component
x_max_values <- NA
y_max_values <- NA

plot_list <- vector("list", length(all_samples))

# Loop through all cycles
for (cycle_index in seq_along(all_samples)) {
  
  samples <- all_samples[[cycle_index]]  # Extract sample list for the current cycle
  
  sample_data <- samples[[1]]
  
  # Extract posterior inverse gamma parameters from the fit object
  cycle_params <- get_cycle_posterior_params(cycle_fits_US_zones_model, cycle_index)
  alpha_est <- cycle_params[["alpha_var_Zone_year"]]
  beta_est <- cycle_params[["beta_var_Zone_year"]]
  
  # Determine x-axis limits (only set during the first cycle)
  if (cycle_index == 1) {
    x_max_values <- min(10, mean(sample_data) * 2)  # Limit x_max to at most 10
  }
  
  # Compute y-axis limit based on histogram and fitted curve
  hist_data <- hist(sample_data, breaks = 10, plot = FALSE)  # Get histogram data
  hist_densities <- hist_data$counts / sum(hist_data$counts) / diff(hist_data$breaks)  # Compute density manually
  
  # Compute fitted curve density values manually
  x_vals <- seq(0, x_max_values, length.out = 1000)  # Generate x-values
  y_vals <- dinvgamma(x_vals, alpha_est, beta_est)  # Compute corresponding densities
  
  if (cycle_index == 1) {
    y_max_values <- max(hist_densities, y_vals, na.rm = TRUE) * 1.5 # Add margin
  } 
  
  dinvgamma <- function(x, alpha, beta) {
    (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x)
  }
  
  # Build x-grid for the fitted PDF curve
  x_grid <- seq(0, x_max_values, length.out = 2000)
  pdf_df <- data.frame(
    x   = x_grid,
    pdf = dinvgamma_stable(x_grid, alpha_est, beta_est)
  )
  
  # Build plot
  plot_list[[cycle_index]] <- ggplot(data.frame(sample_data = sample_data), aes(x = sample_data)) +
    # Histogram scaled to density (freq = FALSE equivalent)
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(min(sample_data), max(sample_data), length.out = 15 + 1),
                   fill = "grey80",
                   color = "white") +
    # Fitted inverse-gamma PDF
    geom_line(data = pdf_df, aes(x = x, y = pdf),
              color = "blue", linewidth = 0.7) +
    geom_vline(xintercept = MLE[1],
               color = "red", linewidth = 1.0) +
    labs(
      title = paste0(" Window ", cycle_index),
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(xlim = c(0, x_max_values),
                    ylim = c(0, y_max_values)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 12),
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20),
      text = element_text(family = "CMU Serif")
    )
}

combined <- wrap_plots(plot_list, ncol = 5)

# Convert patchwork to a grob so we can add shared axis titles
g <- patchwork::patchworkGrob(combined)

gridExtra::grid.arrange(
  g,
  top = grid::textGrob(
    "Year × Zone Variance",
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  bottom = grid::textGrob(
    "Variance Value of Yield [t/ha]",
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  left = grid::textGrob(
    "Inverse Gamma Density",
    rot = 90,
    gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
  ),
  newpage = FALSE
)


# Close PDF device
dev.off()


## Genotype x Zone Variance Component ------------------------------------------

# Access individual fits, e.g., for the first cycle:
cycle_1_fit <- cycle_fits_US_zones_model$cycles[[1]]
# for the 2th cycle:
cycle_2_fit <- cycle_fits_US_zones_model$cycles[[2]]
# for the 3th cycle:
cycle_3_fit <- cycle_fits_US_zones_model$cycles[[3]]
# for the 4th cycle:
cycle_4_fit <- cycle_fits_US_zones_model$cycles[[4]]
# for the 5th cycle:
cycle_5_fit <- cycle_fits_US_zones_model$cycles[[5]]

# Initialize an empty list to store all samples
all_samples <- list()

# Loop through all cycles and construct each sample list dynamically
for (i in seq_along(cycle_fits_US_zones_model$cycles)) {
  cycle_fit <- get(paste0("cycle_", i, "_fit"))  # Dynamically get the cycle fit object
  
  all_samples[[i]] <- list(
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,1,1],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,1,2],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,1,3],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,1,4],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,2,2],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,2,3],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,2,4],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,3,3],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,3,4],
    cycle_fit$posterior_samples$Sigma_Genotype_Zone[,4,4]
  )
}

# Define corresponding variance component names of interest
titles <- c("Genotype × Zone Covariance [1,1]", 
            "Genotype × Zone Covariance [1,2]", 
            "Genotype × Zone Covariance [1,3]",
            "Genotype × Zone Covariance [1,4]",
            "Genotype × Zone Covariance [2,2]",
            "Genotype × Zone Covariance [2,3]", 
            "Genotype × Zone Covariance [2,4]", 
            "Genotype × Zone Covariance [3,3]",
            "Genotype × Zone Covariance [3,4]",
            "Genotype × Zone Covariance [4,4]")

# Vector of REML estimates of original NON Bayesian unstructured Mixed model  
MLEs_geno_zone <- c(var_Genotype_Zone_1_1 = 0.3121,
                    cov_Genotype_Zone_1_2 = 0.2478,
                    cov_Genotype_Zone_1_3 = 0.2943,
                    cov_Genotype_Zone_1_4 = 0.1713,
                    var_Genotype_Zone_2_2 = 0.2157,
                    cov_Genotype_Zone_2_3 = 0.2388,
                    cov_Genotype_Zone_2_4 = 0.1285,
                    var_Genotype_Zone_3_3 = 0.2955,
                    cov_Genotype_Zone_3_4 = 0.1559,
                    var_Genotype_Zone_4_4 = 0.1031
                    )

# Open PDF device
cairo_pdf("5_cycle_fit_37500_27500_27500_27500_27500_GEN_ZONE.pdf", width = 15, height = 15)


# Loop through all cycles
for (cycle_index in seq_along(all_samples)) {
  
  samples <- all_samples[[cycle_index]]  # Extract sample list for the current cycle
  plot_list <- vector("list", length(samples))
  
  for (i in seq_along(samples)) {
    
    sample_data <- samples[[i]]
    
    # Build plot
    plot_list[[i]] <- ggplot(data.frame(sample_data = sample_data), aes(x = sample_data)) +
      # Histogram scaled to density (freq = FALSE equivalent)
      geom_histogram(breaks = seq(min(sample_data), 
                                  max(sample_data), 
                                  length.out = 15 + 1),
                     fill = "grey80",
                     color = "white") +
      # Vertical line at MLE
      geom_vline(xintercept = MLEs_geno_zone[i],
                 color = "red", linewidth = 1.3) +
      # Title + axis labels
      labs(
        title = paste0(titles[i]),
        x = NULL,
        y = NULL
      ) +
      coord_cartesian(xlim = c(0, 0.5),
                      ylim = c(0, 5000)) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 20),
        axis.text.y  = element_text(size = 20),
        text = element_text(family = "CMU Serif")
      )
  }
  
  combined <- wrap_plots(plot_list, nrow = 5)
  
  # Convert patchwork to a grob so we can add shared axis titles
  g <- patchwork::patchworkGrob(combined)
  
  gridExtra::grid.arrange(
    g,
    bottom = grid::textGrob(
      "Variance Value of Yield [t/ha]",
      gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
    ),
    left = grid::textGrob(
      "Histogram",
      rot = 90,
      gp = grid::gpar(fontfamily = "CMU Serif", fontsize = 25)
    ),
    newpage = cycle_index > 1
  )
  
  
}

# Close PDF device
dev.off()


### Grid of optimal design output ----------------------------------------------

## NOTE: Gurobi optimization solver container friendly license required!

grid_bayes_design_US <- grid_bayes_design(cycle_fits_US_zones_model, 
                                          zone_nr = 4, sampsize = 100, seed = 123)

# Average optimal allocations of trials and standard deviations
sampsize = 100
grid_bayes_design_US[[sampsize+1]]

#save(grid_bayes_design_US, file = "grid_bayes_design_US.rda")
#write_xlsx(grid_bayes_design_US, "grid_bayes_design_US.xlsx")
#write.csv(grid_bayes_design_US[[101]], "grid_bayes_design_US.csv")
