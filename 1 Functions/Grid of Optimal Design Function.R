#### "grid_design" function ####################################################
#
#Idea:This function calculates a grid of standard optimal designs as given in the paper of
#     Prus and Piepho 2021 or Prus and Piepho 2024. This function makes use of the 
#     "optimal_design" function to calculate multiple designs for the specified grid of 
#     error variances, location and/or years.
#
#Input: varcomp: The variance component results of an asreml object (summary(US_zones_model)$varcomp)
#                 corresponding to the Frequentist LMM which will than be transformed by the 
#                 "create_variance_info" function here
#       zone_nr: Number of Zones
#       type: Decide of the approach of Prus and Piepho 2021 ("Without year") or the approach
#             of Prus and Piepho 2024 ("With year") should be used
#
#Output: design_df_list_all: Approximated -/ and Exact Design with related
#                            efficiencies and MSE trace evaluation criteria as a list 
#                            for a specified grid of error variances / Total Number of Location /
#                            Number of years dependent on the chosen option type

grid_design <- function(varcomp, zone_nr = 4) {

  require_gurobi_for_grid_design <- function() {
    if (!requireNamespace("gurobi", quietly = TRUE)) {
      stop(
        paste(
          "grid_design() requires the R package 'gurobi' because this workflow",
          "uses OptimalDesign::od_MISOCP() for both approximate and exact",
          "designs.",
          "\n\n",
          "The current container was built without Gurobi support, so the",
          "package is not available.",
          "\n\n",
          "To enable optimal design in this project, rebuild the image after",
          "adding the Linux Gurobi installer archive to the project root so",
          "the Dockerfile can install the R interface, and provide a valid",
          "licence file for runtime use.",
          "\n\n",
          "If you only need the mixed-model fit, stop before calling",
          "grid_design()."
        ),
        call. = FALSE
      )
    }
  }
  
  ## Section 1: HELPER FUNCTIONS -----------------------------------------------
  
  ### "create_variance_info" function
  #
  # Idea: This function processes the variance components information of a fitted 
  #       single stage mixed model to be valid as input for the "optimal_design" function 
  #
  # Input: model: The summary(model)$varcomp asreml-package mixed model object of interest
  #
  # Output: cov_info: The calculated variance covariance matrices of US structure
  #                   from the fitted model corresponding to the random 
  #                   Genotype x Zone interaction effect, together with all other 
  #                   corresponding variance component information of the models
  #                   as a list
  
  create_variance_info <- function(varcomp, zone_nr1 = 4){
    
    #Build Genotype x Zone Covariance Structure
    cov_geno_US <- varcomp[grep("Genotype:Zone!", rownames(varcomp)),"component", drop = FALSE]
    
    # Extract row and column indices from the row names column
    cov_geno_US$row <- as.numeric(gsub(".*Zone_([0-9]+):.*", "\\1", row.names(cov_geno_US)))
    cov_geno_US$col <- as.numeric(gsub(".*Zone_.*:([0-9]+)", "\\1", row.names(cov_geno_US)))
    
    # Initialize an empty matrix with NA values
    cov_geno_US_mat <- matrix(NA, nrow = zone_nr1, ncol = zone_nr1)
    
    # Fill in the matrix using values from df
    for (i in 1:nrow(cov_geno_US)) {
      row <- cov_geno_US$row[i]
      col <- cov_geno_US$col[i]
      component <- cov_geno_US$component[i]
      cov_geno_US_mat[row, col] <- component
      cov_geno_US_mat[col, row] <- component  # Mirror the value for symmetry
    }
    
    cov_geno_US <- cov_geno_US_mat
    
    #Process final function output
    varcomp <- varcomp[!grepl("Genotype:Zone!", rownames(varcomp)), ]
    varcomp <- varcomp[, "component", drop = FALSE]
    
    #Filter rows containing "environment" in the row name and calculate mean
    environment_rows <- grepl("environment", rownames(varcomp))
    environment_resid_mean <- mean(varcomp[environment_rows, "component"])
    
    #Remove the "environment" rows from the data frame
    varcomp <- varcomp[!environment_rows, , drop = FALSE]
    
    #Add a new row with the mean to the bottom of the data frame
    varcomp <- rbind(varcomp, environment_resid_mean = environment_resid_mean)
    
    # Convert the data frame to a list with each row as a separate component
    varcomp_list <- lapply(1:nrow(varcomp), function(i) varcomp[i, 1])
    names(varcomp_list) <- rownames(varcomp)
    
    varcomp_list[["Genotype:Zone"]] <- cov_geno_US
    
    return(varcomp_list)
    varcomp_list
    
  }
  
  ### "optimal_design" function
  #
  #Idea:This function calculates the optimal design as given in the paper of 
  #     Prus and Piepho 2024.The only main difference 
  #     is that this function was written to be usable based on fitted ONE STAGE mixed models
  #
  #Input: varcomp_processed: The calculated variance covariance matrix from 
  #                 the fitted model corresponding to the random 
  #                 Genotype x Zone interaction effect, together with all other 
  #                 variance component information as a list created by the 
  #                 "create_variance_info" function
  #       p: Number of Zones
  #       J: Total Number of Location
  #       H: Number of years considered for the Optimal Design Computation
  #       rep: Number of Replicates in the underlying data set
  #       type: Decide if "Standard" or "Weighted" Design will be computed
  #       
  #Output: Design_Info: Approximated -/ and Exact Design 
  #                     Efficiency and MSE trace evaluation criteria as a list
  
  
  optimal_design <- function(varcomp_processed, p = 4, J = 10, H = 10, rep = 3, 
                             type = c("Standard", "Weighted")){

    require_gurobi_for_grid_design()
    
    ## Known constants
    
    s1<-varcomp_processed$`Genotype:year`/H
    
    s2<-varcomp_processed$`Genotype:Zone:year`/H
    
    s3<-(varcomp_processed$`Genotype:Zone:Location:year` +
           varcomp_processed$environment_resid_mean)/H
    
    #Design matrix (not the same as in paper)
    F<-diag(p)
    
    #Model Genotype x Zone covariance information 
    V <- varcomp_processed$`Genotype:Zone`
    
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
    
    Design_Info_with_year <- list(DesAin9w, DesEin9, Efficiency, MSE_trace)
    
    names(Design_Info_with_year) <- c("Approximated Design", "Exact Design", 
                                      "Efficiency", "Trace of MSE Matrix")
    
    return(Design_Info_with_year)
    Design_Info_with_year
    
  }
  
  ## Section 2: OPTIMAL DESIGN ANALYSES ----------------------------------------
  
  # Initialize an empty list to store data frames for each varcomp
  design_df_list_all <- list()

  # Generate the appropriate varcomp using create_variance_info for each cov_type
  varcomp_processed <- create_variance_info(varcomp_US_zones_model, zone_nr1 = zone_nr)
  
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
      
      # Run the optimal_design function for current J and sig
      design <- optimal_design(varcomp_processed, J = J, H = H, type = "Standard")
      
      # Create a data frame row for the current J and sig
      design_row <- data.frame(
        H = H,
        J = J,
        round(t(design$`Approximated Design`), 2),
        Eff_a = round(design$Efficiency$`Approximate vs Balanced`, 2),
        t(design$`Exact Design`),
        Eff_e = round(design$Efficiency$`Exact vs Balanced`, 2),
        MSE_Tr = round(design$`Trace of MSE Matrix`, 4)
      )
      
      # Store the row in the results list
      design_df_list[[row_index]] <- design_row
      row_index <- row_index + 1
    }
  }
    
  ## Section 3: PROCESS FINAL RESULTS ------------------------------------------
    
  # Combine all rows into a single data frame
  design_df_list_all <- do.call(rbind, design_df_list)
  names(design_df_list_all) <- c("H", "J", "w_1", "w_2", "w_3", "w_4", "Eff_a", 
                                             "J_1", "J_2", "J_3", "J_4", "Eff_e", "MSE_Tr")

  # View the final data frame
  print(design_df_list_all)
    
  return(design_df_list_all)
  
}
