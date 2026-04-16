### "my_residuals" function ####################################################
#
# Idea: This function stores the marginal and conditional residuals of an asreml model object
#
# Input: model: asreml model object
#        n: Number of observations, default is 4272 (sammple size of genotype group Medium,
#                                                    Winter season)
# 
# Output: my_residuals: Marginal -/ and conditional residual vector in a list

my_residuals <- function(model, n = 4272){
  
    ### Process Data -----------------------------------------------------------
    
    #Define vectors of fixed effects: Intercept and Zone
    intercept <- rep(as.vector(model$coefficients$fixed)[1], times = n)
    
    #Create a data frame with labeled raw conditional residuals
    daten <- data.frame(Raw_cond_residuals = c(model$residuals),
                        environment = c(model$mf$environment))
    
    #Create a vector of estimated standard errors of residuals
    Cond_Residual_Var <- tapply(daten$Raw_cond_residuals, daten$environment, var, na.rm = TRUE)
    daten2 <- data.frame(Cond_Residual_Var = c(Cond_Residual_Var),
                         environment = c(unique(model$mf$environment)))
    daten2 <- merge(daten, daten2, by.x = "environment", by.y = "environment", all.x = TRUE)
    
    #Create a data frame with labeled environment specific conditional response 
    #standard errors (residual standard errors) 
    vparms <- as.data.frame(model$vparameters)
    resid_vparms <- vparms[grep("environment", rownames(vparms)),, drop = FALSE]
    resid_vparms <- resid_vparms$`model$vparameters`
    daten3 <- data.frame(Cond_Response_Var = c(resid_vparms),
                         environment = c(unique(model$mf$environment)))
    
    ### Store Residuals --------------------------------------------------------
    
    #Store Raw Marginal Residuals
    Raw_marg_residuals <- model$mf$yield - intercept
    Raw_marg_residuals <- as.vector(Raw_marg_residuals)
    
    #Store Raw Conditional Residuals
    Raw_cond_residuals <- model$residuals
    Raw_cond_residuals <- as.vector(Raw_cond_residuals)
    
    #Store Student Conditional Residuals
    Student_cond_residuals <- daten2$Raw_cond_residuals / sqrt(daten2$Cond_Residual_Var)
    
    #Store Pearson Conditional Residuals
    daten4 <- merge(daten, daten3, by.x = "environment", by.y = "environment", all.x = TRUE)
    daten4$Cond_Residual_Var <- daten4$Cond_Residual_Var
    daten4$Pearson_cond_residuals <- daten4$Raw_cond_residuals / sqrt(daten4$Cond_Response_Var)
    Pearson_cond_residuals <- daten4$Pearson_cond_residuals
    
    ### Store Output -----------------------------------------------------------
    
    my_residuals <- list(Raw_marg_residuals, Raw_cond_residuals, 
                         Student_cond_residuals, Pearson_cond_residuals)
    names(my_residuals) <- c("Raw Marginal Residuals", "Raw Conditional Residuals", 
                             "Student Conditional Residuals", "Pearson Conditional Residuals")
    
    return(my_residuals)
    
}
