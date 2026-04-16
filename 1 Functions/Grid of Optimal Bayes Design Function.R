#### "grid_bayes_design" function ##############################################
#
#Idea: This function calculates a grid of standard optimal design means and standard deviations
#      based on the paper of Prus and Piepho 2021 or Prus and Piepho 2024 but using posterior samples of designs. 
#      This function makes use of the "optimal_bayes_design" and "create_bayes_variance_info" function 
#      to calculate multiple design means and standard deviations
#      for the specified grid of error variances, location and/or years.
#
#Input: cycle_fits: A list containing the fitted Bayesian mixed model results.
#       zone_nr: Number of Zones
#       sampsize: Number of samples to draw in the "create_bayes_variance_info" function
#       type: Decide of the approach of Prus and Piepho 2021 ("Without Year") or the approach
#             of Prus and Piepho 2024 ("With Year") should be used
#       seed: Random seed for reproducibility.
#
#Output: design_df_list_all: Approximated -/ and Exact Design posterior mean and standard deviation 
#                            with related efficiencies and MSE trace evaluation criteria as a list 
#                            for a specified grid of error variances / Total Number of Location /
#                            Number of Years dependent on the chosen option type


grid_bayes_design <- function(cycle_fits, zone_nr = 4, sampsize = 10, seed = 123) {

  require_gurobi_for_bayes_design <- function() {
    if (!requireNamespace("gurobi", quietly = TRUE)) {
      stop(
        paste(
          "grid_bayes_design() requires the R package 'gurobi' because this",
          "workflow uses OptimalDesign::od_MISOCP() for both approximate and",
          "exact designs.",
          "\n\n",
          "The current container was built without Gurobi support, so the",
          "package is not available.",
          "\n\n",
          "To enable Bayesian optimal design in this project, rebuild the",
          "image after adding the Linux Gurobi installer archive to the",
          "project root so the Dockerfile can install the R interface, and",
          "provide a valid licence file for runtime use."
        ),
        call. = FALSE
      )
    }
  }
  
  ## Section 1: HELPER FUNCTIONS -----------------------------------------------
  
  ### "create_bayes_variance_info" function
  #
  #Idea: This function processes the variance components posterior sample information of a fitted 
  #      Bayesian mixed model via "fit_multiple_cycles" function to be valid as 
  #      input for the "optimal_bayes_design" function.
  #
  # Input: cycle_fits: A list containing the fitted Bayesian mixed model results.
  #        zone_nr1: Number of Zones.
  #        sampsize: Number of samples to draw
  #        seed: Random seed for reproducibility.
  #        
  # Output: cycle_fits: A list containing the drawn variance components based on the 
  #                     approximated inverse gamma and inverse Wishart assumptions
  #                     which are estimated based on the variance component 
  #                     posterior samples extracted from the "fit_multiple_cycles" function output.
  
  create_bayes_variance_info <- function(cycle_fits, zone_nr1 = 4, sampsize = 10, seed = 123){
    
    # Set seed for reproducibility
    set.seed(seed)
    
    # Function to draw samples from the inverse gamma distribution
    draw_invgamma <- function(n, shape, scale) {
      inv_gamma_sample <- rinvgamma(n, shape, scale = scale)
      return(inv_gamma_sample)
    }
    
    # Function to draw samples from the inverse Wishart distribution
    draw_invwishart <- function(n, df, scale_matrix) {
      inv_wishart_sample <- list()
      for (i in 1:n) {
        sample <- riwish(df, scale_matrix)
        inv_wishart_sample[[i]] <- sample
      }
      return(inv_wishart_sample)
    }
    
    # List to store the samples
    sample_list <- list()
    
    # Loop through the parameter pairs
    for (i in seq_along(cycle_fits$parameters$final_posteriors$result)) {
      if (startsWith(names(cycle_fits$parameters$final_posteriors$result)[i], "alpha_var")) {
        # Extract the pair of shape and scale
        pair_name <- sub("alpha_var", "", names(cycle_fits$parameters$final_posteriors$result)[i])
        shape <- cycle_fits$parameters$final_posteriors$result[[paste0("alpha_var", pair_name)]]
        scale <- cycle_fits$parameters$final_posteriors$result[[paste0("beta_var", pair_name)]]
        
        # Draw samples
        samples <- draw_invgamma(n = sampsize, shape = shape, scale = scale)
        
        # Store the samples in the list
        sample_list[[paste0("samples_", pair_name)]] <- samples
      }
    }
    
    
    # Check if scale matrix is positive definite
    if (!is.positive.definite(cycle_fits$parameters$final_posteriors$result$psi_Genotype_Zone[,])) {
      stop("Scale matrix must be positive definite.")
    }
    
    # Draw samples from the inverse Wishart distribution
    inv_wishart_samples <- draw_invwishart(n = sampsize,  cycle_fits$parameters$final_posteriors$result$nu_Genotype_Zone, cycle_fits$parameters$final_posteriors$result$psi_Genotype_Zone[,])
    inv_wishart_samples1 <- array(unlist(inv_wishart_samples), dim = c(zone_nr1, zone_nr1, sampsize))
    # Permute dimensions to dim = c(num_matrices, nrow, ncol)
    inv_wishart_samples1 <- aperm(inv_wishart_samples1, c(3, 1, 2))
    
    
    # Store results in lists
    sample_list[[paste0("samples_gen_zone")]] <- inv_wishart_samples1
    
    
    return(sample_list)
    
  }
  
  ### "optimal_bayes_design" function
  #
  #Idea: This function calculates the optimal design as given in the paper of
  #      Prus and Piepho 2021 or Prus and Piepho 2024.The only main difference 
  #      is that this function was written to be usable based on fitted ONE STAGE mixed models!
  #
  #Input: bayes_cov_info: The calculated variance component samples as generated by the
  #                       "create_bayes_variance_info" function
  #
  #       p: Number of Zones
  #       J: Total Number of Location
  #       H: Number of Years considered for the Optimal Design Computation
  #       sig: error variance corresponding to the model specification of 
  #            Prus and Piepho 21, default is 0.1
  #       rep: Number of Replicates in the underlying data set
  #       type: Decide if "Standard" or "Weighted" Design will be computed
  #       type2: Decide of the approach of Prus and Piepho 2021 ("Without Year") or the approach
  #              of Prus and Piepho 2024 ("With Year") should be used
  #       
  #Output: Design_Info: List of:
  #                     Approximated -/ and Exact Design Sample wrt. sample of variance parameters 
  #                     by "create_bayes_variance_info" function
  #                     Estimated standard deviations corresponding to components of 
  #                     calculated Approximated -/ and Exact Design samples 
  #                     Efficiency and MSE trace evaluation criteria as a list
  
  
  optimal_bayes_design <- function(bayes_cov_info, p = 4, J = 10, H = 10, rep = 3, 
                                   type = c("Standard", "Weighted")){

    require_gurobi_for_bayes_design()
    
    Design_Info <- list()
    
    for (i in 1:length(bayes_cov_info$samples__resid_env)) {
      
      ## Known constants
      
      s1<-bayes_cov_info$samples__Genotype_year[i]/H
      
      s2<-bayes_cov_info$samples__Genotype_Zone_year[i]/H
      
      s3<-(bayes_cov_info$samples__Genotype_Zone_Location_year[i] +
             bayes_cov_info$samples__resid_env_mean[i])/H
      
      #Design matrix (not the same as in paper)
      F<-diag(p)
      
      if (length(bayes_cov_info$samples__resid_env_mean) == 1) {
        #Model Genotype x Zone covariance information 
        V <- bayes_cov_info$samples_gen_zone[,]
      }
      else {
        #Model Genotype x Zone covariance information 
        V <- bayes_cov_info$samples_gen_zone[i,,]
      }
      
      #Standard (Bayesian) A-criterion
      if (type == "Standard") {
        L <- diag(p)
      }
      
      #Matrix of weights for weighted A-criterion
      else if (type == "Weighted") {
        L <- rbind(c(813685,0,0,0,0),c(0,432716,0,0,0),c(0,0,477365,0,0),
                   c(0,0,0,995298,0),c(0,0,0,0,1174818)) #subregion related area
      }
      
      L<-L/sum(L)
      
      ## Constrains
      
      #Constraints needed for Bayesian criterion only
      b<-c(t(c(rep(1,p))),J)
      A<-cbind(matrix(0,nrow=p,ncol=p),diag(p))
      A<-rbind(A,cbind(t(c(rep(1,p))),t(c(rep(0,p)))))
      
      #Constraints needed for J_i larger than 0 for all i (at least one location in each sub-region)
      bJ<-c(rep(1,p))
      AJ<-cbind(diag(p),matrix(0,nrow=p,ncol=p))
      
      ## Optimal Design calculation
      
      #Matrices B and Q (Adjusted Variance-covariance matrix)
      B <- V+s2*diag(p)+s1*matrix(1,nrow=p,ncol=p)
      Q <- 1/s3*B
      
      #Define matrices of approximate and exact designs 
      DesA <- matrix(0, ncol=2*p, nrow=1)
      DesE <- matrix(0, ncol=2*p, nrow=1)
      
      #Prepare main input for the od_MISOCP function (depended on cov-comps and design weighting) 
      W<-solve(B)%*%V%*%L%*%V%*%solve(B)
      C<-chol(W)
      C<-t(C)
      U<-solve(C)%*%solve(Q)%*%solve(t(C))
      res<-eigen(U)
      Z<-diag(sqrt(res$values))%*%t(res$vector)
      O<-Z
      F1t <- F%*%solve(t(C))
      FN<-rbind(F1t,O)
      
      #Computing approximate and exact designs 
      res1 <- od_MISOCP(FN, b2=bJ, A2=AJ, b3=b, A3=A, crit="A",type="approximate", gap=0, t.max=Inf)
      #res1 <- od_MISOCP(FN, b3=b, A3=A, crit="A",type="approximate", gap=0, t.max=Inf) #Without constrain
      DesA<-res1$w.best
      res <- od_MISOCP(FN, b2=bJ, A2=AJ, b3=b, A3=A, crit="A",type="exact", gap=0, t.max=Inf)
      #res <- od_MISOCP(FN, b3=b, A3=A, crit="A",type="exact", gap=0, t.max=Inf) #Without constrain
      DesE <- res$w.best
      
      DesAin9<-DesA[1:p]
      DesEin9<-DesE[1:p]
      DesAin9w<-DesAin9/J
      
      ## Efficiency / Trace of MSE Calculation based on A-Criteria ratios
      
      #Standard A-criterion
      if (type == "Standard") {
        
        #Specify information matrices
        MA<-diag(DesAin9)
        ME<-diag(DesEin9)
        Mb<-diag(p)/p*J
        
        #Calculate A-Criteria
        PhiA<-sum(diag(W%*%solve(MA+solve(Q))))
        PhiE<-sum(diag(W%*%solve(ME+solve(Q))))
        Phib<-sum(diag(W%*%solve(Mb+solve(Q))))
        
        #Calculate Efficiency
        EffAb<-PhiA/Phib
        EffEb<-PhiE/Phib
        
        Efficiency <- list(EffAb, EffEb)
        names(Efficiency) <- c("Approximate vs Balanced", 
                               "Exact vs Balanced")
        
        #Trace of MSE (adjusted, as in paper)
        W1<-solve(B)%*%V%*%V%*%solve(B)
        MSEtr<-sum(diag(W1%*%solve(MA+solve(Q))))*s3
        
        MSE_trace <- MSEtr
      }
      
      #Weighted A-criterion
      else if (type == "Weighted") {
        
        #Specify information matrices
        MA<-diag(DesAin9)
        ME<-diag(DesEin9)
        Mb<-diag(p)/p*J
        Ml<-1/sum(diag(L))*L*J
        
        #Calculate A-Criteria
        PhiA<-sum(diag(W%*%solve(MA+solve(Q))))
        Phib<-sum(diag(W%*%solve(Mb+solve(Q))))
        Phil<-sum(diag(W%*%solve(Ml+solve(Q))))
        PhiE<-sum(diag(W%*%solve(ME+solve(Q))))
        
        #Calculate Efficiency
        EffAb<-PhiA/Phib
        EffEb<-PhiE/Phib
        EffAl<-PhiA/Phil
        
        Efficiency <- list(EffAb, EffEb, EffAl)
        names(Efficiency) <- c("Approximate vs Balanced", 
                               "Exact vs Balanced",
                               "Approximate vs Weighted")
        
        #Trace of MSE (adjusted, as in paper)
        W1<-solve(B)%*%V%*%V%*%solve(B)
        MSEtr<-sum(diag(W1%*%solve(MA+solve(Q))))*s3
        
        MSE_trace <- MSEtr
      }
      
      ## Output
      
      Design_Info[[i]] <- list(DesAin9w, DesEin9, Efficiency, MSE_trace)
      
      names(Design_Info[[i]]) <- c("Approximated Design", "Exact Design", 
                                   "Efficiency", "Trace of MSE Matrix")
      
    }
    
    return(Design_Info)
    
  }
  
  ## Section 2: OPTIMAL DESIGN ANALYSES ----------------------------------------
  
  ## Create List with Design relevant Variance info
  bayes_cov_info <- create_bayes_variance_info(cycle_fits, zone_nr1 = zone_nr, 
                                               sampsize = sampsize, seed = seed)
    
  # Initialize an empty list to store data frames for each bayes_cov_info
  bayes_design_df_list_all <- list()
  
  # Loop over each bayes_cov_info type
  for (i in 1:length(bayes_cov_info$samples__resid_env)) {
    
    # Initialize an empty list to store each row
    design_df_list <- list()
    
    # Define the possible values for J and sig
    J_values <- c(10, 20, 40, 100, 200)
    H_values <- c(3)
    
    # Counter to keep track of the row index
    row_index <- 1
    
    # Loop through each J value
    for (H in H_values) {
      
      # Loop through each sig value within each J value
      for (J in J_values) {
        
        # Extract current sample component
        bayes_cov_info1 <- lapply(bayes_cov_info, function(x) {
          if (is.array(x) && length(dim(x)) == 3) {
            x[i,,] # Extract [i,,] if it is a 3D array
          } else if (is.vector(x) || is.array(x)) {
            x[i] # Extract the first element for 1D arrays or vectors
          } else {
            NULL # Handle unexpected cases
          }
        })
        
        # Run the optimal_design function for current J and sig
        design <- optimal_bayes_design(bayes_cov_info1, J = J, H = H, type = "Standard")
        
        # Create a data frame row for the current J and sig
        design_row <- data.frame(
          H = H,
          J = J,
          round(t(design[[1]]$`Approximated Design`), 2),
          Eff_a = round(design[[1]]$Efficiency$`Approximate vs Balanced`, 2),
          t(design[[1]]$`Exact Design`),
          Eff_e = round(design[[1]]$Efficiency$`Exact vs Balanced`, 2),
          MSE_Tr = round(design[[1]]$`Trace of MSE Matrix`, 4)
        )
        
        # Store the row in the results list
        design_df_list[[row_index]] <- design_row
        row_index <- row_index + 1
      }
    }
    
    # Combine all rows into a single data frame
    bayes_design_df_list_all[[i]] <- do.call(rbind, design_df_list)
    names(bayes_design_df_list_all[[i]]) <- c("H", "J", "w_1", "w_2", "w_3", "w_4", "Eff_a", 
                                               "J_1", "J_2", "J_3", "J_4", "Eff_e", "MSE_Tr")
    
  }
    
  ## Section 3: PROCESS FINAL RESULTS ------------------------------------------
    
  # Create a summary data frame for each J and sigma2 combination
  all_data_combined <- do.call(rbind, bayes_design_df_list_all)
  summary_df <- all_data_combined %>%
    group_by(H, J) %>%
    summarise(
      w_1 = paste0(round(mean(w_1, na.rm = TRUE), 2), " (", round(sd(w_1, na.rm = TRUE), 2), ")"),
      w_2 = paste0(round(mean(w_2, na.rm = TRUE), 2), " (", round(sd(w_2, na.rm = TRUE), 2), ")"),
      w_3 = paste0(round(mean(w_3, na.rm = TRUE), 2), " (", round(sd(w_3, na.rm = TRUE), 2), ")"),
      w_4 = paste0(round(mean(w_4, na.rm = TRUE), 2), " (", round(sd(w_4, na.rm = TRUE), 2), ")"),
      Eff_a = paste0(round(mean(Eff_a, na.rm = TRUE), 2), " (", round(sd(Eff_a, na.rm = TRUE), 3), ")"),
      J_1 = paste0(round(mean(J_1, na.rm = TRUE), 2), " (", round(sd(J_1, na.rm = TRUE), 2), ")"),
      J_2 = paste0(round(mean(J_2, na.rm = TRUE), 2), " (", round(sd(J_2, na.rm = TRUE), 2), ")"),
      J_3 = paste0(round(mean(J_3, na.rm = TRUE), 2), " (", round(sd(J_3, na.rm = TRUE), 2), ")"),
      J_4 = paste0(round(mean(J_4, na.rm = TRUE), 2), " (", round(sd(J_4, na.rm = TRUE), 2), ")"),
      Eff_e = paste0(round(mean(Eff_e, na.rm = TRUE), 2), " (", round(sd(Eff_e, na.rm = TRUE), 3), ")"),
      MSE_Tr = paste0(round(mean(MSE_Tr, na.rm = TRUE), 3), " (", round(sd(MSE_Tr, na.rm = TRUE), 4), ")"),
      .groups = "drop"
    )
  
  # Append the summary data frame to the list
  bayes_design_df_list_all[[length(bayes_design_df_list_all) + 1]] <- as.data.frame(summary_df)
  
  return(bayes_design_df_list_all)
    
}
  
