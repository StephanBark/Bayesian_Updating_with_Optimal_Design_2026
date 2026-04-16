################# Frequentist LMM and Optimal Design Analysis ##################


### One stage winter medium group - unstructured zones in genotype and zone random effect

asreml.options(ai.sing = TRUE, workspace = "8gb", maxit = 50, trace = TRUE)

yield_winter_medium$yield <- as.numeric(yield_winter_medium$yield)
yield_winter_medium_copy <- copy(yield_winter_medium)
setorder(yield_winter_medium_copy, Group, environment, ID)

US_zones_model <- asreml(#Block model:
                         fixed = yield ~ Zone, 
                         random = ~  year +
                                     Zone:year +
                                     Zone:Location:year +
                                     Zone:Location:Rep:year +
                         #Treatment model:  
                                     Genotype:year + 
                                     Genotype:us(Zone) +
                                     Genotype:Zone:year +
                                     Genotype:Zone:Location:year, 
                         #Residual Model:
                         residual = ~ dsum(~ ID|(environment)),
                         data = yield_winter_medium_copy)

summary(US_zones_model)

# Extract variance information
varcomp_US_zones_model <- summary(US_zones_model)$varcomp

# Extract REML estimates
geno_zone_US <- varcomp_US_zones_model[grep("Genotype:Zone!", 
                                            rownames(varcomp_US_zones_model)), , drop = FALSE]
MLEs_geno_zone <- geno_zone_US$component

MLEs_geno_zone <- c(var_Genotype_Zone_1_1 = MLEs_geno_zone[1],  # 1,1
                    cov_Genotype_Zone_1_2 = MLEs_geno_zone[2],  # 2,1 -> 1,2
                    cov_Genotype_Zone_1_3 = MLEs_geno_zone[4],  # 3,1 -> 1,3
                    cov_Genotype_Zone_1_4 = MLEs_geno_zone[7],  # 4,1 -> 1,4
                    var_Genotype_Zone_2_2 = MLEs_geno_zone[3],  # 2,2
                    cov_Genotype_Zone_2_3 = MLEs_geno_zone[5],  # 3,2 -> 2,3
                    cov_Genotype_Zone_2_4 = MLEs_geno_zone[8],  # 4,2 -> 2,4
                    var_Genotype_Zone_3_3 = MLEs_geno_zone[6],  # 3,3
                    cov_Genotype_Zone_3_4 = MLEs_geno_zone[9],  # 4,3 -> 3,4
                    var_Genotype_Zone_4_4 = MLEs_geno_zone[10]  # 4,4
                    )
MLEs_geno_zone

resid_env <- varcomp_US_zones_model[grep("environment", 
                                         rownames(varcomp_US_zones_model)), , drop = FALSE]
var_resid_env_mean <- mean(resid_env$component)
var_resid_env_mean

# Other MLEs named
MLEs_named <- c(
  var_year                           = 0.0287,
  var_Zone_year                      = 0.0000,
  var_Zone_Location_year             = 0.5813,
  var_Zone_Location_Rep_year         = 0.0046,
  var_Genotype_year                  = 0.0246,
  var_Genotype_Zone_year             = 0.0000,
  var_Genotype_Zone_Location_year    = 0.2502,
  var_resid_env_mean                 = 0.3243
)

### LMM Assumptions via Residual Plots ----------------------------------------

# Use paper font
font_add("CMU Serif", "cmunrm.ttf")
showtext_auto() 

# Set font globally
par(family = "CMU Serif")

## Preprocessing --------------------------------------------------------------

# Get fits and residuals
daten <- as.data.table(US_zones_model$mf)
daten$fitted <- as.vector(US_zones_model$linear.predictors)

daten$raw_marg_resid <- my_residuals(US_zones_model, 
                                     n = nrow(daten))$`Raw Marginal Residuals`
daten$raw_cond_resid <- my_residuals(US_zones_model, 
                                     n = nrow(daten))$`Raw Conditional Residuals`
daten$student_cond_resid <- my_residuals(US_zones_model, 
                                         n = nrow(daten))$`Student Conditional Residuals`
daten$pearson_cond_resid <- my_residuals(US_zones_model, 
                                         n = nrow(daten))$`Pearson Conditional Residuals`


## Basic Yield based cond residual plots --------------------------------------

a <- ggplot(daten, aes(x = fitted, y = student_cond_resid)) + 
            geom_point(size = 0.6) +
            labs(x = "Fitted Yield", 
                 y = "Conditional student residuals") +
            ggtitle("Residuals vs. Fitted Yield") +
            theme_bw() +
            theme(text = element_text(family = "CMU Serif"), 
                  axis.title = element_text(size = 19), 
                  axis.text = element_text(size = 15),
                  plot.title = element_text(size = 20))

b <- ggplot(daten, aes(x = student_cond_resid)) + 
            geom_histogram(color="black", fill="white") + 
            labs(x = "Contidional student residuals", 
                 y = "Absolute frequency") + 
            ggtitle("Histogram of residuals") +
            theme_bw() +
            theme(text = element_text(family = "CMU Serif"), 
                  axis.title = element_text(size = 19), 
                  axis.text = element_text(size = 15),
                  plot.title = element_text(size = 20))

c <- ggplot(daten, aes(sample = student_cond_resid)) + 
            stat_qq(size = 0.6) + 
            stat_qq_line() + 
            labs(x = "Theoretical quantiles of normal distribution", 
                 y = "Observed quantiles of residuals") +
            ggtitle("QQ-Plot of residuals") +
            theme_bw() +
            theme(text = element_text(family = "CMU Serif"), 
                  axis.title = element_text(size = 19), 
                  axis.text = element_text(size = 15),
                  plot.title = element_text(size = 20))

grid.arrange(a, b, c, ncol = 3)


## Ordered marginal yield residual plots --------------------------------------

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$year, daten$Location)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_marg_resid, na.rm = TRUE)

ggplot(daten) + aes(y = raw_marg_resid, x = index, colour = year) +
  geom_point(size = 1.5) +
  labs(x = "Observation count with coloured location ordered by year / location",
       y = "Marginal residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(10, "mm"), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$Location, daten$Genotype)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_marg_resid, na.rm = TRUE)

ggplot(daten) + 
  aes(y = raw_marg_resid, x = index, colour = Location) +
  geom_point(size = 1.5) + 
  labs(x = "Observation count with coloured location ordered by location",
       y = "Marginal residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$year, daten$Location)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_marg_resid, na.rm = TRUE)

ggplot(daten) + 
  aes(y = raw_marg_resid, x = index, colour = Location) +
  geom_point(size = 1.5) + 
  labs(x = "Observation count with coloured location ordered by year / location",
       y = "Marginal residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))


## Ordered conditional Yield residual plots -----------------------------------

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$year, daten$Genotype)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_cond_resid, na.rm = TRUE)

ggplot(daten) + 
  aes(y = raw_cond_resid, x = index, colour = year) +
  geom_point(size = 1.5) + 
  labs(x = "Observation count with coloured year ordered by year",
       y = "Conditional residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(10, "mm"), 
        axis.text.x = element_blank(),
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$Location, daten$Genotype)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_cond_resid, na.rm = TRUE)

ggplot(daten) + aes(y = raw_cond_resid, x = index, colour = Location) +
  geom_point(size = 1.5) + 
  labs(x = "Observation count with coloured location ordered by location",
       y = "Conditional residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))

##------------------------------------------------------------------------------
# Ordering
daten <- daten[order(daten$year, daten$Location)]

# Add an index column for plotting
daten$index <- seq_len(nrow(daten))

# Check range of yield
range(daten$raw_cond_resid, na.rm = TRUE)

ggplot(daten) + aes(y = raw_cond_resid, x = index, colour = Location) +
  geom_point(size = 1.5) + 
  labs(x = "Observation count with coloured location ordered by year / location",
       y = "Conditional residuals of yield in tons per hectare [t/ha]") +
  coord_cartesian(ylim = c(-3.6, 3.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text = element_text(size = 20)) + 
  scale_colour_viridis_d(option = "H") +
  guides(colour = guide_legend(override.aes = list(size = 6)))


## Check Heteroscedasticity of Variance ---------------------------------------

ggplot(daten, aes(x = fitted, y = (student_cond_resid)^2)) + 
  geom_point(size = 0.8, shape = 1) +
  geom_smooth(se = TRUE) +
  labs(x = "Fitted values of yield", 
       y = "Squared conditional student residuals") +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))

ggplot(daten, aes(x = environment, y = student_cond_resid)) + 
  geom_point(size = 0.8, shape = 1) +
  geom_smooth(se = TRUE) +
  labs(x = "Unique year / location combination", 
       y = "Conditional studentized residuals") +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))


## Check Residuals close to zero ----------------------------------------------

ggplot(daten, aes(x = fitted, y = student_cond_resid)) + 
  geom_point(size = 0.8, shape = 1) +
  labs(x = "Fitted values of yield", 
       y = "Squared conditional residuals") +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme_bw() +
  theme(text = element_text(family = "CMU Serif"), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))


### Grid of optimal design output ----------------------------------------------

## NOTE: Gurobi optimization solver container friendly license required!

## Load varcomp_US_zones_model.rda file if previous section is not fordable because of missing asreml license
#load("./3 Frequentist Framework/varcomp_US_zones_model.rda")

grid_design_US <- grid_design(varcomp_US_zones_model, zone_nr = 4)
#write.csv(grid_design_US, "grid_design_US.csv")
