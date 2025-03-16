# Load required libraries
library(ggplot2)
library(reshape2)
library(gridExtra)


set.seed(100)

# Function to simulate data from the endemic-epidemic model with Xu's dispersion specification
simulate_ee_ar_xu <- function(T = 552, lambda_start = 100, theta = 0.5, 
                              alpha_nu = 5, beta_nu = 0.1, gamma_nu = 0.1,
                              alpha_phi = 0.2, beta_phi = 0.1, gamma_phi = 0.1,
                              gamma = 0.5, S = 52) {
  # Initialize variables
  Xt <- numeric(T)  # Infection counts per week
  lambda_t <- numeric(T)  # Conditional mean
  nu_t <- numeric(T)  # Seasonal components
  phi_t_values <- numeric(T)  # Epidemic components
  gamma_values <- numeric(T)  # Autoregressive components
  
  # Initial values
  Xt[1] <- rnbinom(1, size = lambda_start/theta, mu = lambda_start)
  lambda_t[1] <- lambda_start
  nu_t[1] <- exp(alpha_nu + beta_nu * sin(2 * pi * 1 / S) + gamma_nu * cos(2 * pi * 1 / S))
  phi_t_values[1] <- 0  # No epidemic force at t = 1
  gamma_values[1] <- 0  # No autoregression at t = 1
  
  # Generate values for t = 2 to T
  for (t in 2:T) {
    # Compute components
    nu_t[t] <- exp(alpha_nu + beta_nu * sin(2 * pi * t / S) + gamma_nu * cos(2 * pi * t / S))
    phi_t_values[t] <- exp((alpha_phi + beta_phi * sin(2 * pi * t / S) + gamma_phi * cos(2 * pi * t / S))) * Xt[t-1]
    gamma_values[t] <- gamma * lambda_t[t-1]
    
    # Compute conditional mean
    lambda_t[t] <- nu_t[t] + phi_t_values[t] + gamma_values[t]
    
    # Sample new infections using Xu's dispersion model
    Xt[t] <- rnbinom(1, size = lambda_t[t]/theta, mu = lambda_t[t])
  }
  
  return(list(Xt = Xt[-(1:S)], lambda_t = lambda_t[-(1:S)], nu_t = nu_t[-(1:S)], 
              phi_t_values = phi_t_values[-(1:S)], gamma_values = gamma_values[-(1:S)])) 
}

scenarios <- list(
  # Original Scenarios
  list(gamma = 0.3, alpha_phi = -1.0, beta_phi = 0.3, gamma_phi = 0.3,
       alpha_nu = 1.3, beta_nu = 0.3, gamma_nu = 0.3, theta = 0.5, 
       title = "Low Autoregression, Low Epidemic, High Seasonality"),
  
  list(gamma = 0.5, alpha_phi = -1.2, beta_phi = 0.4, gamma_phi = 0.4,
       alpha_nu = 0.7, beta_nu = 0.4, gamma_nu = 0.4,  theta = 6,
       title = "Medium Autoregression, Medium Epidemic, Medium Seasonality"),
  
  list(gamma = 0.7, alpha_phi = -1.6, beta_phi = 0.3, gamma_phi = 0.2,
       alpha_nu = 2, beta_nu = 1, gamma_nu = 0.5, theta =  100,  
       title = "High Autoregression, High Epidemic, Low Seasonality")
)

# Simulate all scenarios
results_simulation <- lapply(scenarios, function(scenario) {
  simulate_ee_ar_xu(gamma = scenario$gamma, 
                    alpha_phi = scenario$alpha_phi, 
                    beta_phi = scenario$beta_phi, 
                    gamma_phi = scenario$gamma_phi, 
                    alpha_nu = scenario$alpha_nu, 
                    beta_nu = scenario$beta_nu, 
                    gamma_nu = scenario$gamma_nu,
                    theta = scenario$theta)
})

# Plot function for scenarios
plot_scenarios <- function(results_simulation, scenarios) {
  time_series_plots <- list()
  for (i in 1:length(results_simulation)) {
    df <- data.frame(Time = 1:length(results_simulation[[i]]$Xt), Cases = results_simulation[[i]]$Xt)
    time_series_plots[[i]] <- ggplot(df, aes(x = Time, y = Cases)) +
      geom_line(color = "blue") +
      theme_minimal() +
      labs(title = paste("Scenario:", i), x = "Time", y = "Cases")
  }
  return(time_series_plots)
}

# Generate plots
plots_scenarios <- plot_scenarios(results_simulation, scenarios)

# Arrange plots in a grid
grid.arrange(grobs = plots_scenarios, ncol = 3, nrow = 1)

# 6) Plot Decomposition of seasonal component nu_t, epidemic component phi_t * X_t 
#    and autoregressive component gamma * lambda_t

# Function to plot decomposition for all scenarios
plot_decompositions <- function(results_simulation, scenarios) {
  
  decomposition_plots <- list()
  
  for (i in seq_along(results_simulation)) {
    result <- results_simulation[[i]]  # Select scenario results
    
    df_components <- data.frame(
      Time = 1:length(result$Xt),
      "Endemic" = result$nu_t,
      "Epidemic" = result$phi_t_values,
      "Autoregressive" = result$gamma_values
    )
    
    # Reshape for ggplot
    df_long <- melt(df_components, id.vars = "Time")
    
    # Create decomposition plot without legend
    p <- ggplot(df_long, aes(x = Time, y = value, color = variable)) +
      geom_line() +
      theme_minimal() +
      labs(title = paste("Scenario", i), 
           x = "Time", y = "Component Contribution") +
      scale_color_manual(values = c("red", "blue", "green")) +
      theme(legend.position = "none")  # Remove legend from individual plots
    
    decomposition_plots[[i]] <- p  # Store each plot in a list
  }
  
  # Extract legend from the first plot
  p_legend <- ggplot(df_long, aes(x = Time, y = value, color = variable)) +
    geom_line() +
    theme_minimal() +
    labs(color = "Component") +  # Ensure the legend appears
    scale_color_manual(values = c("red", "blue", "green")) 
  
  legend <- g_legend(p_legend)  # Extract the legend
  
  return(list(plots = decomposition_plots, legend = legend))
}

# Function to extract the legend from a ggplot object
g_legend <- function(a_gplot){
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Generate decomposition plots and extract the legend
decomposition_result <- plot_decompositions(results_simulation, scenarios)

# Arrange plots with a single legend on the right
grid.arrange(
  grobs = c(decomposition_result$plots, list(decomposition_result$legend)), 
  ncol = 4,  # 3 plots + 1 legend column
  widths = c(1, 1, 1, 0.3)  # Adjust legend width
)

# See what was computed for each scenario
results_1 <- as.data.frame(results_simulation[[1]])
results_2 <- as.data.frame(results_simulation[[2]])
results_3 <- as.data.frame(results_simulation[[3]])

# 7) Run simulation over multiple scenarios
set.seed(123)

start_time <- Sys.time()

num_simulations <- 1000

results_simulation <- lapply(scenarios, function(scenario) {
  replicate(num_simulations,   simulate_ee_ar_xu(gamma = scenario$gamma, 
                                                  alpha_phi = scenario$alpha_phi, 
                                                  beta_phi = scenario$beta_phi, 
                                                  gamma_phi = scenario$gamma_phi, 
                                                  alpha_nu = scenario$alpha_nu, 
                                                  beta_nu = scenario$beta_nu, 
                                                  gamma_nu = scenario$gamma_nu,
                                                  theta = scenario$theta),
            simplify = FALSE)
})

# Negative Log-Likelihood Function for Xu's Model
neg_log_likelihood_xu <- function(params, Xt, S = 52) {
  alpha_nu <- params[1]
  beta_nu <- params[2]
  gamma_nu <- params[3]
  alpha_phi <- params[4]
  beta_phi <- params[5]
  gamma_phi <- params[6]
  gamma <- params[7]
  theta <- params[8]
  
  T <- length(Xt)
  lambda_t <- numeric(T)
  lambda_t[1] <- Xt[1]
  
  for (t in 2:T) {
    nu_t <- exp(alpha_nu + beta_nu * sin(2 * pi * t / S) + gamma_nu * cos(2 * pi * t / S))
    phi_t <- exp(alpha_phi + beta_phi * sin(2 * pi * t / S) + gamma_phi * cos(2 * pi * t / S))
    lambda_t[t] <- nu_t + phi_t * Xt[t-1] + gamma * lambda_t[t-1]
    
    if (lambda_t[t] <= 0) return(1e10)
  }
  
  ll <- sum(dnbinom(Xt[2:T], size = lambda_t[2:T]/theta, mu = lambda_t[2:T], log = TRUE), na.rm = TRUE)
  return(-ll)
}

estimate_mle <- function(Xt) {
  # Define initial parameters and bounds
  init_params <- c(5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5)
  lower_bounds <- c(0, 0, 0, -3, 0, 0, 0, 0.001)
  upper_bounds <- c(5, 3, 1, 3, 1, 1, 1, 1000)
  
  # Optimization using L-BFGS-B
  res_lbfgsb <- optim(par = init_params,
                      fn = neg_log_likelihood_xu,
                      Xt = Xt,
                      method = "L-BFGS-B",
                      lower = lower_bounds,
                      upper = upper_bounds)
  
  # Optimization using nlminb
  res_nlminb <- nlminb(start = init_params,
                       objective = function(params) neg_log_likelihood_xu(params, Xt = Xt),
                       lower = lower_bounds,
                       upper = upper_bounds)
  
  # Choose the result with the lower negative log-likelihood
  if (res_lbfgsb$value < res_nlminb$objective) {
    chosen <- list(method = "L-BFGS-B", par = res_lbfgsb$par, value = res_lbfgsb$value)
  } else {
    chosen <- list(method = "nlminb", par = res_nlminb$par, value = res_nlminb$objective)
  }
  
  return(chosen$par)
}

# Extract true values from scenarios
true_values <- lapply(scenarios, function(scenario) {
  c(scenario$alpha_nu, scenario$beta_nu, scenario$gamma_nu,
    scenario$alpha_phi, scenario$beta_phi, scenario$gamma_phi,
    scenario$gamma, scenario$theta)  
})

scenario_names <- c("Low Autoregression, Low Epidemic, High Seasonality", 
                    "Medium Autoregression, Medium Epidemic, Medium Seasonality", 
                    "High Autoregression, High Epidemic, Low Seasonality")

print(true_values)

# Apply MLE estimation to each simulated dataset (1000 per scenario)
mle_results <- lapply(results_simulation, function(scenario_simulations) {
  lapply(scenario_simulations, function(sim) estimate_mle(sim$Xt))
})

# 9) Compute comparison true values vs. estimated values

# Convert mle results and true values into a matrix n_sim by parameters
mle_matrix <- do.call(rbind, unlist(mle_results, recursive = FALSE)) 
true_values_expanded <- do.call(rbind, lapply(true_values, function(x) {
  matrix(rep(x, times = num_simulations), ncol = length(x), byrow = TRUE)
})) 

param_means_1 <- colMeans(mle_matrix[1:num_simulations,], na.rm = TRUE)
param_means_2 <- colMeans(mle_matrix[(num_simulations + 1):(num_simulations*2),], na.rm = TRUE) 
param_means_3 <- colMeans(mle_matrix[(num_simulations*2 + 1):(3 * num_simulations),], na.rm = TRUE) 

param_sd_1 <- apply(mle_matrix[1:num_simulations,], 2, sd, na.rm = TRUE)
param_sd_2 <- apply(mle_matrix[(num_simulations + 1):(num_simulations*2),], 2, sd, na.rm = TRUE) 
param_sd_3<- apply(mle_matrix[(num_simulations*2 + 1):(3 * num_simulations),], 2, sd, na.rm = TRUE) 

# Compute Root Mean Squared Error (RMSE)
rmse_values_1 <- sqrt(colMeans((mle_matrix[1:num_simulations,] - 
                                  true_values_expanded[1:num_simulations,])^2, na.rm = TRUE))
rmse_values_2 <- sqrt(colMeans((mle_matrix[(num_simulations + 1):(num_simulations*2),] - 
                                  true_values_expanded[(num_simulations + 1):(num_simulations*2),])^2, na.rm = TRUE))
rmse_values_3 <- sqrt(colMeans((mle_matrix[(num_simulations*2 + 1):(3 * num_simulations),] - 
                                  true_values_expanded[(num_simulations*2 + 1):(3 * num_simulations),] )^2, na.rm = TRUE))

# Compute Relative Standard Error (RSE)
rse_values_1 <- rmse_values_1 / colMeans(true_values_expanded[1:num_simulations,], na.rm = TRUE)
rse_values_2 <- rmse_values_2 / colMeans(true_values_expanded[(num_simulations + 1):(num_simulations*2),], na.rm = TRUE)
rse_values_3 <- rmse_values_3 / colMeans(true_values_expanded[(num_simulations*2 + 1):(3 * num_simulations),], na.rm = TRUE)

true = as.numeric(as.character(unlist(true_values)))
estimated = as.numeric(as.character(c(param_means_1, param_means_2, param_means_3)))

# Combine parameter estimates into one long format using rbind()
comparison_df <- data.frame(
  Parameter = rep(c("Alpha_Nu", "Beta_Nu", "Gamma_Nu", 
                    "Alpha_Phi", "Beta_Phi", "Gamma_Phi", "Gamma", "Psi"), 3),
  True = true,
  Estimated = estimated,
  Difference = abs(true - estimated),
  Relative_Diff = abs(true - estimated)/abs(true),
  Sd = c(param_sd_1,param_sd_2,param_sd_3),
  RMSE = c(rmse_values_1, rmse_values_2, rmse_values_3),
  RSE = c(rse_values_1, rse_values_2, rse_values_3),
  Scenario = rep(c(scenario_names[1], scenario_names[2], scenario_names[3]), each = 8) 
)

# Compute Execution Time
end_time <- Sys.time()  
execution_time <- end_time - start_time  
print(paste("Total computation time:", execution_time))


# Function to format numeric columns with comma as decimal separator
format_decimal <- function(df) {
  # Identify numeric columns
  num_cols <- sapply(df, is.numeric)
  
  # Apply formatting to numeric columns
  df[, num_cols] <- lapply(df[, num_cols], function(x) 
    gsub("\\.", ",", format(round(x, 3), nsmall = 3)))
  
  return(df)
}

# Initialize a data frame to store results
pooled_stats <- data.frame(
  Scenario = scenario_names,
  Pooled_Mean = NA,
  Pooled_Variance = NA,
  Average_Mean = NA,
  Average_Variance = NA,
  stringsAsFactors = FALSE
)

# Iterate over each scenario
for (i in seq_along(results_simulation)) {
  sim_list <- results_simulation[[i]]
  
  # Pool all Xt values from every simulation run into one vector
  pooled_Xt <- unlist(lapply(sim_list, function(sim) sim$Xt))
  
  # Calculate pooled mean and variance
  pooled_mean <- mean(pooled_Xt)
  pooled_var  <- var(pooled_Xt)
  
  # For comparison: compute mean and variance for each simulation run
  individual_means <- sapply(sim_list, function(sim) mean(sim$Xt))
  individual_vars  <- sapply(sim_list, function(sim) var(sim$Xt))
  
  # Average the individual estimates
  average_mean <- mean(individual_means)
  average_var  <- mean(individual_vars)
  
  # Store results in the data frame
  pooled_stats[i, "Pooled_Mean"] <- pooled_mean
  pooled_stats[i, "Pooled_Variance"] <- pooled_var
  pooled_stats[i, "Average_Mean"] <- average_mean
  pooled_stats[i, "Average_Variance"] <- average_var
}

# Apply the function to `comparison_df`
comparison_df <- format_decimal(comparison_df)
pooled_stats <- format_decimal(pooled_stats)
# Print to check
print(comparison_df)
print(pooled_stats)

