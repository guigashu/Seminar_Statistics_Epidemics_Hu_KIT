library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(grid)
library(Metrics)
library(forecast)
library(lmtest)
library(skedastic)


# 1) Data Preparation
data_set <- read_csv("Downloads/National-level data_SINGAPORE_20010101_20221224.csv", 
                     col_types = cols(adm_0_name = col_skip(), 
                                      adm_1_name = col_skip(), adm_2_name = col_skip(), 
                                      full_name = col_skip(), ISO_A0 = col_skip(), 
                                      FAO_GAUL_code = col_skip(), RNE_iso_code = col_skip(), 
                                      IBGE_code = col_skip(), calendar_start_date = col_character(), 
                                      calendar_end_date = col_character(), 
                                      case_definition_standardised = col_skip(), 
                                      S_res = col_skip(), 
                                      UUID = col_skip()))

# Filter data set for only weekly values
data_set <- data_set[data_set$T_res == "Week", ]

data_set <- data_set %>%
  mutate(
    calendar_start_date = ymd(calendar_start_date),
    calendar_end_date   = ymd(calendar_end_date)
  )

# Ensure your date is in Date format 
data_set$calendar_end_date <- as.Date(data_set$calendar_end_date, format = "%Y%m%d")

# Separate data for the analysis
data_set <- data_set[data_set$Year %in% 2006:2017, ]

# Split into training and testing 

n <- nrow(data_set)
split_index <- floor(0.7 * n)
train_set <- data_set[1:split_index, ]
test_set  <- data_set[(split_index + 1):n, ]
Xt_data_training <- train_set$dengue_total
Xt_data_test <- test_set$dengue_total

# 2) Model Estimation via MLE with CLS initial values

# Define the three negative log-likelihood functions:

neg_log_likelihood_zhu <- function(params, Xt, S = 52) {
  # Zhu model:
  #   X_t | past ~ NegBin(mean = lambda_t, size = 1/psi)
  #   => Var(X_t | past) = lambda_t + psi * lambda_t^2
  endemic_intercept   <- params[1]
  endemic_sine_coef   <- params[2]
  endemic_cosine_coef <- params[3]
  epidemic_intercept  <- params[4]
  epidemic_sine_coef  <- params[5]
  epidemic_cosine_coef<- params[6]
  gamma               <- params[7]
  psi                 <- params[8]
  
  T <- length(Xt)
  lambda_t <- numeric(T)
  lambda_t[1] <- Xt[1]
  
  for (t in 2:T) {
    nu_t  <- exp(endemic_intercept + endemic_sine_coef * sin(2 * pi * t / S) +
                   endemic_cosine_coef * cos(2 * pi * t / S))
    phi_t <- exp(epidemic_intercept + epidemic_sine_coef * sin(2 * pi * t / S) +
                   epidemic_cosine_coef * cos(2 * pi * t / S))
    lambda_t[t] <- nu_t + phi_t * Xt[t-1] +gamma* lambda_t[t-1]
    if (is.na(lambda_t[t]) || lambda_t[t] <= 0) return(1e10)
  }
  # Using constant size = 1/psi
  ll <- sum(dnbinom(Xt[2:T], size = 1/psi, mu = lambda_t[2:T], log = TRUE), na.rm = TRUE)
  return(-ll)
}

neg_log_likelihood_xu <- function(params, Xt, S = 52) {
  # Xu model:
  #   X_t | past ~ NegBin(mean = lambda_t, size = lambda_t/theta)
  #   => Var(X_t | past) = (1 + theta) * lambda_t
  endemic_intercept   <- params[1]
  endemic_sine_coef   <- params[2]
  endemic_cosine_coef <- params[3]
  epidemic_intercept  <- params[4]
  epidemic_sine_coef  <- params[5]
  epidemic_cosine_coef<- params[6]
  gamma               <- params[7]
  theta               <- params[8]
  
  T <- length(Xt)
  lambda_t <- numeric(T)
  lambda_t[1] <- Xt[1]
  
  for (t in 2:T) {
    nu_t  <- exp(endemic_intercept + endemic_sine_coef * sin(2 * pi * t / S) +
                   endemic_cosine_coef * cos(2 * pi * t / S))
    phi_t <- exp(epidemic_intercept + epidemic_sine_coef * sin(2 * pi * t / S) +
                   epidemic_cosine_coef * cos(2 * pi * t / S))
    lambda_t[t] <- nu_t + phi_t * Xt[t-1] +gamma* lambda_t[t-1]
    if (is.na(lambda_t[t]) || lambda_t[t] <= 0) return(1e10)
    
  }
  # For Xu: size = lambda_t/theta
  ll <- sum(dnbinom(Xt[2:T], size = lambda_t[2:T] / theta, mu = lambda_t[2:T], log = TRUE), na.rm = TRUE)
  return(-ll)
}

neg_log_likelihood_generalized <- function(params, Xt, S = 52) {
  endemic_intercept   <- params[1]
  endemic_sine_coef   <- params[2]
  endemic_cosine_coef <- params[3]
  epidemic_intercept  <- params[4]
  epidemic_sine_coef  <- params[5]
  epidemic_cosine_coef<- params[6]
  gamma               <- params[7]
  psi                 <- params[8]
  theta               <- params[9]
  
  T <- length(Xt)
  lambda_t <- numeric(T)
  lambda_t[1] <- Xt[1]
  
  for (t in 2:T) {
    nu_t <- exp(endemic_intercept + endemic_sine_coef * sin(2 * pi * t / S) + endemic_cosine_coef * cos(2 * pi * t / S))
    phi_t <- exp(epidemic_intercept + epidemic_sine_coef * sin(2 * pi * t / S) + epidemic_cosine_coef * cos(2 * pi * t / S))
    lambda_t[t] <- nu_t + phi_t * Xt[t-1] + gamma * lambda_t[t-1]
    
    if (is.na(lambda_t[t]) || lambda_t[t] <= 0) return(1e10) # Prevent invalid values
  }
  
  dispersion_t <- psi + (theta / lambda_t[2:T])
  ll <- sum(dnbinom(Xt[2:T], size = 1 / dispersion_t, mu = lambda_t[2:T], log = TRUE), na.rm = TRUE)
  return(-ll)  # Negative log-likelihood
}


# Define estimation function

estimate_NBINGARCH_model <- function(Xt, S = 52, model = c("zhu", "xu", "generalized"), 
                                     n_starts = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  
  # ---------------------------------------------
  # Step 1: Obtain CLS Estimates for Baseline Dynamics
  # ---------------------------------------------
  # Use a simple regression of Xt[t] on Xt[t-1]
  cls_model <- lm(Xt[-1] ~ Xt[-length(Xt)])
  a0_cls <- coef(cls_model)[1]   # baseline intercept
  a1_cls <- coef(cls_model)[2]   # lag coefficient
  
  # ---------------------------------------------
  # Step 2: Initialize Seasonal Parameters
  # ---------------------------------------------
  # We initialize the endemic and epidemic intercepts by log-transforming the CLS estimates
  endemic_intercept_init  <- log(max(a0_cls, 1e-6))
  epidemic_intercept_init <- log(max(a1_cls, 1e-6))
  
  # Initialize seasonal sine and cosine coefficients to zero
  endemic_sine_coef_init    <- 0
  endemic_cosine_coef_init  <- 0
  epidemic_sine_coef_init   <- 0
  epidemic_cosine_coef_init <- 0
  
  # Initialize gamma (autoregressive effect on lambda) to a default value, e.g., 0.5
  gamma_init <- 0.5
  
  # ---------------------------------------------
  # Step 3: Initialize Dispersion Parameters Based on Model
  # ---------------------------------------------
  if (model %in% c("zhu", "xu")) {
    # For Zhu and Xu models: 8 parameters (last parameter: psi for Zhu or theta for Xu)
    dispersion_init <- 0.1
    base_init <- c(endemic_intercept_init,
                   endemic_sine_coef_init,
                   endemic_cosine_coef_init,
                   epidemic_intercept_init,
                   epidemic_sine_coef_init,
                   epidemic_cosine_coef_init,
                   gamma_init,
                   dispersion_init)
    lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0, 1e-6)
    upper_bounds <- c( Inf,  Inf,  Inf,  Inf,  Inf,  Inf, 0.99, Inf)
  } else if (model == "generalized") {
    # For the generalized model: 9 parameters (dispersion parameters: psi and theta)
    psi_init <- 0.1
    theta_init <- 0.1
    base_init <- c(endemic_intercept_init,
                   endemic_sine_coef_init,
                   endemic_cosine_coef_init,
                   epidemic_intercept_init,
                   epidemic_sine_coef_init,
                   epidemic_cosine_coef_init,
                   gamma_init,
                   psi_init,
                   theta_init)
    lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0, 1e-6, 1e-6)
    upper_bounds <- c( Inf,  Inf,  Inf,  Inf,  Inf,  Inf, 0.99, Inf, Inf)
  }
  
  # ---------------------------------------------
  # Step 4: Select the Appropriate Negative Log-Likelihood Function
  # ---------------------------------------------
  if (model == "zhu") {
    neg_log_like <- neg_log_likelihood_zhu
  } else if (model == "xu") {
    neg_log_like <- neg_log_likelihood_xu
  } else if (model == "generalized") {
    neg_log_like <- neg_log_likelihood_generalized
  }
  
  # ---------------------------------------------
  # Step 5: Multi-Start Optimization
  # ---------------------------------------------
  best_val <- Inf
  best_result <- NULL
  best_params <- NULL
  
  for (i in 1:n_starts) {
    # Generate a random perturbation around the base initial values
    perturbation <- rnorm(length(base_init), mean = 0, sd = 0.1)
    start_params <- base_init + perturbation
    # Ensure starting parameters satisfy lower bounds
    start_params <- pmax(start_params, lower_bounds)
    
    result <- optim(
      par = start_params,
      fn = neg_log_like,
      Xt = Xt,
      S = S,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds
    )
    
    if (result$value < best_val) {
      best_val <- result$value
      best_result <- result
      best_params <- result$par
    }
  }
  
  # Assign names to parameters for output
  if (model %in% c("zhu", "xu")) {
    if (model == "zhu") {
      names(best_params) <- c("endemic_intercept", "endemic_sine_coef", "endemic_cosine_coef",
                              "epidemic_intercept", "epidemic_sine_coef", "epidemic_cosine_coef",
                              "gamma", "psi")
    } else {
      names(best_params) <- c("endemic_intercept", "endemic_sine_coef", "endemic_cosine_coef",
                              "epidemic_intercept", "epidemic_sine_coef", "epidemic_cosine_coef",
                              "gamma", "theta")
    }
  } else {
    names(best_params) <- c("endemic_intercept", "endemic_sine_coef", "endemic_cosine_coef",
                            "epidemic_intercept", "epidemic_sine_coef", "epidemic_cosine_coef",
                            "gamma", "psi", "theta")
  }
  
  return(list(estimated_parameters = best_params,
              best_value = best_val,
              best_result = best_result,
              model = model))
}


# Estimate parameters
result_zhu <- estimate_NBINGARCH_model(Xt_data_training, model = "zhu", n_starts = 20, seed = 123)
result_xu  <- estimate_NBINGARCH_model(Xt_data_training, model = "xu", n_starts = 20, seed = 123)
result_gen <- estimate_NBINGARCH_model(Xt_data_training, model = "generalized", n_starts = 20, seed = 123)

print(result_zhu$estimated_parameters)
print(result_xu$estimated_parameters)
print(result_gen$estimated_parameters)

params_zhu<- result_zhu$estimated_parameters
params_xu <- result_xu$estimated_parameters
params_gen <- result_gen$estimated_parameters

# 3) Validate Models

# Compute Recursive Conditional Means lambda_t
compute_lambda <- function(Xt, S = 52, params) {
  # Extract common parameters:
  endemic_intercept   <- params["endemic_intercept"]
  endemic_sine_coef   <- params["endemic_sine_coef"]
  endemic_cosine_coef <- params["endemic_cosine_coef"]
  epidemic_intercept  <- params["epidemic_intercept"]
  epidemic_sine_coef  <- params["epidemic_sine_coef"]
  epidemic_cosine_coef<- params["epidemic_cosine_coef"]
  gamma               <- params["gamma"]
  
  T <- length(Xt)
  lambda <- numeric(T)
  # Initialize: here we simply set lambda[1] = Xt[1]
  lambda[1] <- Xt[1]
  
  for (t in 2:T) {
    nu_t  <- exp(endemic_intercept +
                   endemic_sine_coef * sin(2 * pi * t / S) +
                   endemic_cosine_coef * cos(2 * pi * t / S))
    phi_t <- exp(epidemic_intercept +
                   epidemic_sine_coef * sin(2 * pi * t / S) +
                   epidemic_cosine_coef * cos(2 * pi * t / S))
    lambda[t] <- nu_t + phi_t * Xt[t-1] + gamma * lambda[t-1]
    # safeguard: if lambda becomes nonpositive, force a small positive value
    if(lambda[t] <= 0) lambda[t] <- 1e-6
  }
  return(lambda)
}

# Fit X_t
simulate_fitted_counts <- function(lambda_t, model, params) {
  T <- length(lambda_t)
  fitted_sim <- numeric(T)
  fitted_sim[1] <- lambda_t[1]  # For t=1, just use the initial value as an anchor
  
  # Identify the dispersion parameter(s) depending on the model
  if(model == "zhu") {
    # size = 1 / psi, mean = lambda_t[t]
    psi <- params["psi"]
    for(t in 2:T) {
      # Draw 1000 samples from NB
      draws <- rnbinom(1000, size = 1/psi, mu = lambda_t[t])
      fitted_sim[t] <- mean(draws)
    }
  } else if(model == "xu") {
    # size = lambda_t[t] / theta, mean = lambda_t[t]
    theta <- params["theta"]
    for(t in 2:T) {
      size_val <- lambda_t[t] / theta
      # Avoid negative or zero size
      if(size_val <= 0) size_val <- 1e-6
      draws <- rnbinom(1000, size = size_val, mu = lambda_t[t])
      fitted_sim[t] <- mean(draws)
    }
  } else if(model == "generalized") {
    # size = 1 / (psi + theta / lambda_t[t]), mean = lambda_t[t]
    psi <- params["psi"]
    theta <- params["theta"]
    for(t in 2:T) {
      denom <- psi + theta / lambda_t[t]
      if(denom <= 0) denom <- 1e-6
      size_val <- 1 / denom
      draws <- rnbinom(1000, size = size_val, mu = lambda_t[t])
      fitted_sim[t] <- mean(draws)
    }
  }
  
  return(fitted_sim)
}

# Validate the model
validate_model <- function(estimated_params, best_neglogL, Xt, S = 52, model = c("zhu", "xu", "generalized")) {
  model <- match.arg(model)
  T <- length(Xt)
  n_params <- length(estimated_params)
  
  # Compute conditional means once (λₜ)
  lambda_t <- compute_lambda(Xt, S, estimated_params)
  
  # Simulate fitted counts from the NB distribution (using median)
  fitted_sim <- simulate_fitted_counts(lambda_t, model, estimated_params)
  
  # Compute Information Criteria (using best_neglogL)
  n_eff <- T - 1  
  logL <- -best_neglogL
  AIC <- 2 * n_params - 2 * logL
  BIC <- n_params * log(n_eff) - 2 * logL
  
  # Set up the NB 'size' parameter vector for t = 2:T (used in PIT and residuals)
  if (model == "zhu") {
    size_vec <- rep(1 / estimated_params["psi"], T - 1)
  } else if (model == "xu") {
    size_vec <- lambda_t[2:T] / estimated_params["theta"]
  } else if (model == "generalized") {
    size_vec <- 1 / (estimated_params["psi"] + estimated_params["theta"] / lambda_t[2:T])
  }
  
  # Vectorized PIT calculations for t = 2:T
  cdf_lower <- pnbinom(Xt[2:T] - 1, size = size_vec, mu = lambda_t[2:T])
  cdf_upper <- pnbinom(Xt[2:T], size = size_vec, mu = lambda_t[2:T])
  u <- runif(T - 1)
  PIT <- cdf_lower + u * (cdf_upper - cdf_lower)
  
  # Vectorized Pearson residuals for t = 2:T
  if (model == "zhu") {
    var_vec <- lambda_t[2:T] + estimated_params["psi"] * lambda_t[2:T]^2
  } else if (model == "xu") {
    var_vec <- lambda_t[2:T] + estimated_params["theta"] * lambda_t[2:T]
  } else if (model == "generalized") {
    var_vec <- lambda_t[2:T] + (estimated_params["psi"] + estimated_params["theta"] / lambda_t[2:T]) * lambda_t[2:T]^2
  }
  pearson_res <- (Xt[2:T] - lambda_t[2:T]) / sqrt(var_vec)
  
  # Vectorized Tail Probabilities: P(X >= observed) for t = 2:T
  tail_probs <- 1 - pnbinom(Xt[2:T] - 1, size = size_vec, mu = lambda_t[2:T])
  
  # Perform Ljung-Box test on Pearson residuals
  lb_lag <- min(10, length(pearson_res) - 1)
  lb_test <- Box.test(pearson_res, lag = lb_lag, type = "Ljung-Box")
  
  # Plotting functions:
  
  # Observed vs Fitted (using λₜ)
  obs_vs_fitted_plot <- function() {
    plot(lambda_t, Xt,
         xlab = "Fitted values (lambda_t)",
         ylab = "Observed counts (X_t)",
         main = "Observed vs Fitted Values (λₜ)",
         pch = 19, col = "blue")
    abline(0, 1, col = "red", lwd = 2, lty = 2)
  }
  
  # Observed vs Simulated Fitted (fitted_sim)
  obs_vs_simulated_plot <- function() {
    plot(fitted_sim, Xt,
         xlab = "Simulated Fitted Counts",
         ylab = "Observed counts (X_t)",
         main = "Observed vs Simulated Fitted Values",
         pch = 19, col = "purple")
    abline(0, 1, col = "red", lwd = 2, lty = 2)
  }
  
  # Pearson Residuals vs Fitted (using λₜ)
  residuals_vs_fitted_plot <- function() {
    plot(lambda_t[-1], pearson_res,
         xlab = "Fitted values (lambda_t)",
         ylab = "Pearson Residuals",
         main = "Pearson Residuals vs Fitted Values",
         pch = 19, col = "blue")
    abline(h = 0, col = "red", lty = 2)
  }
  
  # QQ Plot and ACF/PACF for residuals
  plot_QQ_residuals <- function() {
    qqnorm(pearson_res, main = "QQ Plot of Pearson Residuals")
    qqline(pearson_res, col = "red", lwd = 2)
  }
  
  plot_residual_acf_pacf <- function() {
    par(mfrow = c(1,2))
    acf(pearson_res, main = "ACF of Pearson Residuals")
    pacf(pearson_res, main = "PACF of Pearson Residuals")
    par(mfrow = c(1,1))
  }
  
  # Generate the plots
  obs_vs_fitted_plot()
  obs_vs_simulated_plot()
  residuals_vs_fitted_plot()
  plot_QQ_residuals()
  plot_residual_acf_pacf()
  
  par(mfrow=c(2,2))
  acf(pearson_res, main = "ACF of Pearson Residuals")
  pacf(pearson_res, main = "PACF of Pearson Residuals")
  residuals_vs_fitted_plot()  # calls the function
  plot_QQ_residuals()         # calls the function
  par(mfrow=c(1,1))
  
  
  # Calculate error metrics comparing both fitted series to observed (for t = 2:T)
  mae_lambda <- mean(abs(Xt[2:T] - lambda_t[2:T]))
  mae_simulated <- mean(abs(Xt[2:T] - fitted_sim[2:T]))
  rmse_lambda <- sqrt(mean((Xt[2:T] - lambda_t[2:T])^2))
  rmse_simulated <- sqrt(mean((Xt[2:T] - fitted_sim[2:T])^2))
  error_table <- data.frame(
    Fit_Type = c("Conditional Mean", "Simulated Fitted"),
    MAE = c(mae_lambda, mae_simulated),
    RMSE = c(rmse_lambda, rmse_simulated)
  )
  print(error_table)
  
  return(list(AIC = AIC,
              BIC = BIC,
              PIT = PIT,
              PearsonResiduals = pearson_res,
              LjungBox = lb_test,
              TailProbabilities = tail_probs,
              lambda_t = lambda_t,
              fitted_sim = fitted_sim,
              obs_vs_fitted_plot = obs_vs_fitted_plot,
              obs_vs_simulated_plot = obs_vs_simulated_plot,
              residuals_vs_fitted_plot = residuals_vs_fitted_plot,
              plot_QQ_residuals = plot_QQ_residuals,
              plot_residual_acf_pacf = plot_residual_acf_pacf,
              ErrorMetrics = error_table))
}

validation_results_zhu <- validate_model(estimated_params = params_zhu, 
                            best_neglogL = result_zhu$best_value, Xt = Xt_data_training, 
                            model = "zhu")
validation_results_xu <- validate_model(estimated_params = params_xu, 
                           best_neglogL = result_xu$best_value, Xt = Xt_data_training, 
                           model = "xu")
validation_results_gen <- validate_model(estimated_params = params_gen, 
                            best_neglogL = result_gen$best_value, Xt = Xt_data_training, 
                            model = "generalized")

# Collect AIC and BIC values in a data frame
comparison_table_AIC_BIC <- data.frame(
  Model = c("Zhu", "Xu", "Generalized"),
  AIC = c(validation_results_zhu$AIC,
          validation_results_xu$AIC,
          validation_results_gen$AIC),
  BIC = c(validation_results_zhu$BIC,
          validation_results_xu$BIC,
          validation_results_gen$BIC)
)
print(comparison_table_AIC_BIC)


# Evaluation for the absence of autocorrelation in the Pearson residuals 
validation_results_zhu$LjungBox
validation_results_xu$LjungBox
validation_results_gen$LjungBox

par(mfrow=c(1,3))
# Evaluation of overall calibration of the predictive distribution via PIT
hist(validation_results_zhu$PIT, breaks = 20, main = "PIT Histogram Zhu", xlab = "PIT")
abline(h = length(validation_results_zhu$PIT) / 20, col = "red", lwd = 2, lty = 2)

hist(validation_results_xu$PIT, breaks = 20, main = "PIT Histogram Xu", xlab = "PIT")
abline(h = length(validation_results_xu$PIT) / 20, col = "red", lwd = 2, lty = 2)

hist(validation_results_gen$PIT, breaks = 20, main = "PIT Histogram Generalized", xlab = "PIT")
abline(h = length(validation_results_gen$PIT) / 20, col = "red", lwd = 2, lty = 2)

ks_test_zhu <- ks.test(validation_results_zhu$PIT, "punif")
print(ks_test_zhu)
ks_test_xu <- ks.test(validation_results_xu$PIT, "punif")
print(ks_test_xu)
ks_test_gen <- ks.test(validation_results_gen$PIT, "punif")
print(ks_test_gen)

# Evaluation tail probabilities
hist(validation_results_zhu$TailProbabilities, breaks = 20, 
     main = "Tail Probability Histogram", xlab = "Tail Probabilities", col = "lightgreen")
summary(validation_results_zhu$TailProbabilities)

hist(validation_results_xu$TailProbabilities, breaks = 20, 
     main = "Tail Probability Histogram", xlab = "Tail Probabilities", col = "lightgreen")
summary(validation_results_xu$TailProbabilities)

hist(validation_results_gen$TailProbabilities, breaks = 20, 
     main = "Tail Probability Histogram", xlab = "Tail Probabilities", col = "lightgreen")
summary(validation_results_gen$TailProbabilities)


# Ensure the train_set aligns with lambda_t (which starts at t=2)
train_set_adj <- train_set[-1, ]

# Compute quantile-based bins
quantiles <- quantile(validation_results_zhu$lambda_t, probs = c(0, 0.33, 0.67, 1))

# Define function to categorize based on data-driven thresholds
categorize_fitted_values <- function(fitted_values) {
  cut(fitted_values,
      breaks = c(-Inf, quantiles[2], quantiles[3], Inf),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE)
}

# Categorize fitted values
train_set_adj$fitted_category_zhu <- categorize_fitted_values(validation_results_zhu$lambda_t[-1])
train_set_adj$fitted_category_xu  <- categorize_fitted_values(validation_results_xu$lambda_t[-1])
train_set_adj$fitted_category_gen <- categorize_fitted_values(validation_results_gen$lambda_t[-1])

print(quantiles)

# Function to generate stratified QQ plots
plot_stratified_qq <- function(residuals, fitted_category, model_name) {
  library(ggplot2)
  
  qq_data <- data.frame(residuals = residuals, category = fitted_category)
  
  ggplot(qq_data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    facet_wrap(~ category, scales = "free") +
    labs(title = paste("Stratified QQ-Plot of Pearson Residuals -", model_name),
         x = "Standardized Normal Quantiles",
         y = "Sample Quantiles") +
    theme_minimal()
}

# Generate QQ plots for each model
qq_plot_zhu <- plot_stratified_qq(validation_results_zhu$PearsonResiduals, 
                                  train_set_adj$fitted_category_zhu, "Zhu Model")

qq_plot_xu <- plot_stratified_qq(validation_results_xu$PearsonResiduals, 
                                 train_set_adj$fitted_category_xu, "Xu Model")

qq_plot_gen <- plot_stratified_qq(validation_results_gen$PearsonResiduals, 
                                  train_set_adj$fitted_category_gen, "Generalized Model")


# Display plots
print(qq_plot_zhu)
print(qq_plot_xu)
print(qq_plot_gen)

library(patchwork)

# Combine the three QQ plots in a 3x3 vertical layout
combined_qq_plots <- qq_plot_zhu / qq_plot_xu / qq_plot_gen

# Print the combined plot
print(combined_qq_plots)


# Compute variance of Pearson residuals for each category
compute_residual_variance <- function(residuals, fitted_category) {
  data.frame(residuals = residuals, category = fitted_category) %>%
    group_by(category) %>%
    summarise(
      Mean_Residuals = mean(residuals),
      Variance_Residuals = var(residuals),
      Count = n()
    )
}

# Compute residual variance for each model
residual_variance_zhu <- compute_residual_variance(validation_results_zhu$PearsonResiduals, train_set_adj$fitted_category_zhu)
residual_variance_xu  <- compute_residual_variance(validation_results_xu$PearsonResiduals, train_set_adj$fitted_category_xu)
residual_variance_gen <- compute_residual_variance(validation_results_gen$PearsonResiduals, train_set_adj$fitted_category_gen)

# Print variance summaries
print(residual_variance_zhu)
print(residual_variance_xu)
print(residual_variance_gen)

# 4) Rolling window forecast

# Function to perform a rolling one‐step‐ahead forecast
rolling_forecast_t1 <- function(train_data, test_data, params, S = 52) {
  n_train <- length(train_data)
  n_test  <- length(test_data)
  forecast_values <- numeric(n_test)
  
  # Compute conditional means for the training data using your compute_lambda function
  lambda_train <- compute_lambda(train_data, S, params)
  # Initialize with the last training conditional mean and observed value
  lambda_last  <- lambda_train[n_train]
  last_observed <- train_data[n_train]
  
  # Rolling forecast loop for each time point in the test set
  for (i in 1:n_test) {
    t_index <- n_train + i  # overall time index
    
    # Compute seasonal components at forecast time t_index
    nu_t <- exp(params["endemic_intercept"] +
                  params["endemic_sine_coef"] * sin(2 * pi * t_index / S) +
                  params["endemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    phi_t <- exp(params["epidemic_intercept"] +
                   params["epidemic_sine_coef"] * sin(2 * pi * t_index / S) +
                   params["epidemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    # One-step ahead forecast using the model: λₜ = νₜ + φₜ * Xₜ₋₁ + γ * λₜ₋₁
    forecast_t <- nu_t + phi_t * last_observed + params["gamma"] * lambda_last
    forecast_values[i] <- forecast_t
    
    # Update the history: use the actual test observation and the forecasted λ for the next step
    last_observed <- test_data[i]
    lambda_last   <- forecast_t
  }
  return(forecast_values)
}

# Run the rolling forecast using the test set and Zhu model estimated parameters
rolling_forecasts_zhu <- rolling_forecast_t1(Xt_data_training, Xt_data_test, params_zhu, S = 52)
rolling_forecasts_xu <- rolling_forecast_t1(Xt_data_training, Xt_data_test, params_xu, S = 52)
rolling_forecasts_gen <- rolling_forecast_t1(Xt_data_training, Xt_data_test, params_gen, S = 52)

# Create a comparison table of forecasts vs. observed values
forecast_comparison_zhu <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_zhu
)

forecast_comparison_zhu <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_xu
)

forecast_comparison_zhu <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_gen
)

# Optionally, compute forecast error metrics (e.g., MAE and RMSE)
mae_forecast_zhu  <- mean(abs(rolling_forecasts_zhu - Xt_data_test))
rmse_forecast_zhu <- sqrt(mean((rolling_forecasts_zhu - Xt_data_test)^2))
mae_forecast_xu  <- mean(abs(rolling_forecasts_xu - Xt_data_test))
rmse_forecast_xu <- sqrt(mean((rolling_forecasts_xu - Xt_data_test)^2))
mae_forecast_gen  <- mean(abs(rolling_forecasts_gen - Xt_data_test))
rmse_forecast_gen <- sqrt(mean((rolling_forecasts_gen - Xt_data_test)^2))

cat("Forecast MAE Zhu:", mae_forecast_zhu, "\n")
cat("Forecast RMSE Zhu:", rmse_forecast_zhu, "\n")
cat("Forecast MAE Xu:", mae_forecast_xu, "\n")
cat("Forecast RMSE Xu:", rmse_forecast_xu, "\n")
cat("Forecast MAE Gen:", mae_forecast_gen, "\n")
cat("Forecast RMSE Gen:", rmse_forecast_gen, "\n")

# Rolling forecast function using simulated counts (Zhu model example)
rolling_forecast_t1_simulated <- function(train_data, test_data, params, S = 52) {
  n_train <- length(train_data)
  n_test  <- length(test_data)
  forecast_values <- numeric(n_test)
  
  # Compute conditional means for training data using your compute_lambda function
  lambda_train <- compute_lambda(train_data, S, params)
  
  # Initialize with the last training observation and its conditional mean
  lambda_last  <- lambda_train[n_train]
  last_observed <- train_data[n_train]
  
  for (i in 1:n_test) {
    t_index <- n_train + i  # overall time index
    
    # Compute seasonal components for time t_index
    nu_t <- exp(params["endemic_intercept"] +
                  params["endemic_sine_coef"] * sin(2 * pi * t_index / S) +
                  params["endemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    phi_t <- exp(params["epidemic_intercept"] +
                   params["epidemic_sine_coef"] * sin(2 * pi * t_index / S) +
                   params["epidemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    # Compute the updated conditional mean (λₜ)
    forecast_lambda <- nu_t + phi_t * last_observed + params["gamma"] * lambda_last
    
    # Instead of using forecast_lambda directly,
    # simulate 1,000 draws from NB (size = 1/psi, mean = forecast_lambda)
    draws <- rnbinom(1000, size = 1/params["psi"], mu = forecast_lambda)
    simulated_forecast <- mean(draws)
    
    forecast_values[i] <- simulated_forecast
    
    # For a one-step-ahead rolling forecast, update using the actual observed test value
    last_observed <- test_data[i]
    lambda_last   <- forecast_lambda
  }
  
  return(forecast_values)
}

# Run the rolling forecast using the test set and Zhu model estimated parameters
rolling_forecasts_zhu_simulated <- rolling_forecast_t1_simulated(Xt_data_training, Xt_data_test, params_zhu, S = 52)
rolling_forecasts_xu_simulated <- rolling_forecast_t1_simulated(Xt_data_training, Xt_data_test, params_xu, S = 52)
rolling_forecasts_gen_simulated <- rolling_forecast_t1_simulated(Xt_data_training, Xt_data_test, params_gen, S = 52)

# Create a comparison table of forecasts vs. observed values
forecast_comparison_zhu_simulated <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_zhu
)

forecast_comparison_xu_simulated <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_xu
)

forecast_comparison_gen_simulated <- data.frame(
  Time     = 1:length(Xt_data_test),
  Observed = Xt_data_test,
  Forecast = rolling_forecasts_gen
)

# Optionally, compute forecast error metrics (e.g., MAE and RMSE)
mae_forecast_zhu_simulated  <- mean(abs(rolling_forecasts_zhu_simulated - Xt_data_test))
rmse_forecast_zhu_simulated <- sqrt(mean((rolling_forecasts_zhu_simulated - Xt_data_test)^2))
mae_forecast_xu_simulated  <- mean(abs(rolling_forecasts_xu_simulated - Xt_data_test))
rmse_forecast_xu_simulated <- sqrt(mean((rolling_forecasts_xu_simulated - Xt_data_test)^2))
mae_forecast_gen_simulated  <- mean(abs(rolling_forecasts_gen_simulated - Xt_data_test))
rmse_forecast_gen_simulated <- sqrt(mean((rolling_forecasts_gen_simulated - Xt_data_test)^2))

cat("Forecast MAE Zhu:", mae_forecast_zhu_simulated, "\n")
cat("Forecast RMSE Zhu:", rmse_forecast_zhu_simulated, "\n")
cat("Forecast MAE Xu:", mae_forecast_xu_simulated, "\n")
cat("Forecast RMSE Xu:", rmse_forecast_xu_simulated, "\n")
cat("Forecast MAE Gen:", mae_forecast_gen_simulated, "\n")
cat("Forecast RMSE Gen:", rmse_forecast_gen_simulated, "\n")

rolling_forecast_t1_simulated_interval <- function(train_data, test_data, params, S = 52, n_sim = 1000, alpha = 0.05, model = c("zhu", "xu", "generalized")) {
  model <- match.arg(model)
  n_train <- length(train_data)
  n_test  <- length(test_data)
  forecast_mean <- numeric(n_test)
  lower_bound   <- numeric(n_test)
  upper_bound   <- numeric(n_test)
  
  # Compute conditional means for training data
  lambda_train <- compute_lambda(train_data, S, params)
  lambda_last  <- lambda_train[n_train]
  last_observed <- train_data[n_train]
  
  for (i in 1:n_test) {
    t_index <- n_train + i  # overall time index
    
    # Compute seasonal components at t_index
    nu_t <- exp(params["endemic_intercept"] +
                  params["endemic_sine_coef"] * sin(2 * pi * t_index / S) +
                  params["endemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    phi_t <- exp(params["epidemic_intercept"] +
                   params["epidemic_sine_coef"] * sin(2 * pi * t_index / S) +
                   params["epidemic_cosine_coef"] * cos(2 * pi * t_index / S))
    
    # Forecast the conditional mean λₜ
    forecast_lambda <- nu_t + phi_t * last_observed + params["gamma"] * lambda_last
    
    # Use the correct simulation for each model
    if(model == "zhu") {
      draws <- rnbinom(n_sim, size = 1 / params["psi"], mu = forecast_lambda)
    } else if(model == "xu") {
      draws <- rnbinom(n_sim, size = forecast_lambda / params["theta"], mu = forecast_lambda)
    } else if(model == "generalized") {
      draws <- rnbinom(n_sim, size = 1 / (params["psi"] + params["theta"] / forecast_lambda), mu = forecast_lambda)
    }
    
    # Compute the 95% prediction interval
    forecast_mean[i] <- mean(draws, na.rm = TRUE)
    lower_bound[i] <- quantile(draws, probs = alpha / 2, na.rm = TRUE)
    upper_bound[i] <- quantile(draws, probs = 1 - alpha / 2, na.rm = TRUE)
    
    # Update state for next forecast
    last_observed <- test_data[i]
    lambda_last   <- forecast_lambda
  }
  
  return(list(mean = forecast_mean, lower = lower_bound, upper = upper_bound))
}

compute_interval_score <- function(observed, lower, upper, alpha = 0.05) {
  n <- length(observed)
  score <- numeric(n)
  out_of_bounds <- numeric(n)
  
  for (i in 1:n) {
    width <- upper[i] - lower[i]
    penalty_lower <- ifelse(observed[i] < lower[i], (2 / alpha) * (lower[i] - observed[i]), 0)
    penalty_upper <- ifelse(observed[i] > upper[i], (2 / alpha) * (observed[i] - upper[i]), 0)
    score[i] <- width + penalty_lower + penalty_upper
    out_of_bounds[i] <- as.numeric(observed[i] < lower[i] || observed[i] > upper[i])
  }
  
  total_out_of_bounds <- sum(out_of_bounds)
  return(list(score = score, total_out_of_bounds = total_out_of_bounds))
}

interval_forecast_zhu <- rolling_forecast_t1_simulated_interval(Xt_data_training, Xt_data_test, params_zhu, S = 52, model = "zhu")
interval_forecast_xu  <- rolling_forecast_t1_simulated_interval(Xt_data_training, Xt_data_test, params_xu, S = 52, model = "xu")
interval_forecast_gen <- rolling_forecast_t1_simulated_interval(Xt_data_training, Xt_data_test, params_gen, S = 52, model = "generalized")

score_zhu <- compute_interval_score(Xt_data_test, interval_forecast_zhu$lower, interval_forecast_zhu$upper)
score_xu  <- compute_interval_score(Xt_data_test, interval_forecast_xu$lower, interval_forecast_xu$upper)
score_gen <- compute_interval_score(Xt_data_test, interval_forecast_gen$lower, interval_forecast_gen$upper)

avg_score_zhu <- mean(score_zhu$score)
avg_score_xu  <- mean(score_xu$score)
avg_score_gen <- mean(score_gen$score)

cat("Zhu Model: Average Interval Score =", avg_score_zhu, 
    ", Out-of-bound count =", score_zhu$total_out_of_bounds, "\n")
cat("Xu Model: Average Interval Score =", avg_score_xu, 
    ", Out-of-bound count =", score_xu$total_out_of_bounds, "\n")
cat("Generalized Model: Average Interval Score =", avg_score_gen, 
    ", Out-of-bound count =", score_gen$total_out_of_bounds, "\n")

