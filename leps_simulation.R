# -------- Set Working Directory --------
setwd("C:\\Users\\1012 G2\\Documents\\R files\\EPS")

# -------- Load Required Libraries --------
library(survival) # Still useful for Surv object, though not for fitting directly
library(ggplot2)
library(dplyr)
library(tidyr)

# -------- Define True Parameters --------
true_params <- list(
  c = 0.7,
  sigma = 1.0,
  beta = c(0.5, -0.4, 0.3, -0.2, 0.1, -0.1, 0.05)  # beta0 to beta6
)
# Make sure beta is named for easier access
names(true_params$beta) <- paste0("beta", 0:6)

# -------- LEPS Quantile Function (for data generation) --------
rLEPS <- function(n, c, sigma, eta) {
  u <- runif(n)
  log_arg <- 1 - u^(1 / c)
  # Handle potential issues with log_arg to avoid NaNs
  if (any(log_arg <= 0 | !is.finite(log_arg))) {
    warning("Invalid arguments in rLEPS, some quantiles might be NaN/Inf.")
    t_val <- rep(NaN, n)
    valid_idx <- which(log_arg > 0 & is.finite(log_arg))
    t_val[valid_idx] <- -log(log_arg[valid_idx]) / (sigma * exp(eta[valid_idx]))
    return(t_val)
  }
  t <- -log(log_arg) / (sigma * exp(eta))
  return(t)
}

# -------- LEPS PDF and Survival Function (for likelihood) --------
# t: time
# c_param, sigma_param: LEPS parameters (renamed to avoid conflict with function args)
# eta: linear predictor (beta0 + X %*% betas)
# log: return log-density/log-survival?
dLEPS <- function(t, c_param, sigma_param, eta, log = FALSE) {
  # Ensure all parameters are positive for a valid distribution
  if (c_param <= 0 || sigma_param <= 0) {
    return(rep(-Inf, length(t))) # Return -Inf log-likelihood for invalid parameters
  }
  
  exp_term <- exp(-t * sigma_param * exp(eta))
  
  # Calculate PDF
  pdf_val <- c_param * (sigma_param * exp(eta)) * exp_term * (1 - exp_term)^(c_param - 1)
  
  if (log) {
    # Handle log(0) for pdf_val which can happen for extreme t or parameters
    log_pdf_val <- log(pdf_val)
    log_pdf_val[pdf_val <= 0] <- -Inf # Assign -Inf for non-positive pdf
    return(log_pdf_val)
  } else {
    return(pdf_val)
  }
}

pLEPS <- function(q, c_param, sigma_param, eta) {
  if (c_param <= 0 || sigma_param <= 0) {
    return(rep(NaN, length(q))) 
  }
  return((1 - exp(-q * sigma_param * exp(eta)))^c_param)
}

sLEPS <- function(t, c_param, sigma_param, eta, log = FALSE) {
  # S(t) = 1 - F(t)
  if (c_param <= 0 || sigma_param <= 0) {
    return(rep(NaN, length(t))) 
  }
  survival_val <- 1 - (1 - exp(-t * sigma_param * exp(eta)))^c_param
  
  if (log) {
    log_survival_val <- log(survival_val)
    log_survival_val[survival_val <= 0] <- -Inf # Assign -Inf for non-positive survival
    return(log_survival_val)
  } else {
    return(survival_val)
  }
}


# -------- Negative Log-Likelihood Function for LEPS Regression --------
# This function will be minimized by optim
neg_log_likelihood_leps <- function(params, time, status, X) {
  # params vector: c, sigma, beta0, beta1, ..., beta6
  c_param <- params[1]
  sigma_param <- params[2]
  beta0 <- params[3]
  betas_cov <- params[4:length(params)] # beta1 to beta6
  
  # Calculate linear predictor eta
  eta <- beta0 + X %*% betas_cov
  
  # Initialize log-likelihood contributions
  loglik_contributions <- numeric(length(time))
  
  # Contributions from uncensored observations (status == 1)
  idx_event <- which(status == 1)
  if (length(idx_event) > 0) {
    loglik_contributions[idx_event] <- dLEPS(time[idx_event], c_param, sigma_param, eta[idx_event], log = TRUE)
  }
  
  # Contributions from censored observations (status == 0)
  idx_censored <- which(status == 0)
  if (length(idx_censored) > 0) {
    loglik_contributions[idx_censored] <- sLEPS(time[idx_censored], c_param, sigma_param, eta[idx_censored], log = TRUE)
  }
  
  # Sum the log-likelihood contributions
  total_loglik <- sum(loglik_contributions)
  
  # optim minimizes, so we return negative log-likelihood
  # Handle cases where log-likelihood might be -Inf (e.g., invalid parameters)
  if (!is.finite(total_loglik)) {
    return(.Machine$double.xmax) # Return a very large number for bad parameters
  }
  
  return(-total_loglik)
}

# -------- Simulation Function (updated for custom MLE) --------
simulate_leps_mle <- function(n, true_params, n_sim = 100) {
  # 9 parameters to estimate: c, sigma, beta0, beta1, ..., beta6
  all_estimates <- matrix(NA, nrow = n_sim, ncol = 9)
  colnames(all_estimates) <- c("c_est", "sigma_est", paste0("beta", 0:6, "_est"))
  
  failure_count <- 0 # Initialize failure counter
  
  # Initial values for optimization (crucial for convergence)
  # These are often close to the true values for simulation, or small positive values
  # for scale/shape and 0 for betas.
  initial_values <- c(
    c = 0.5,     # Initial guess for c
    sigma = 0.8, # Initial guess for sigma
    beta0 = 0,   # Initial guess for beta0
    beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0 # Initial guess for other betas
  )
  
  # Lower and upper bounds for parameters (for L-BFGS-B method)
  # c and sigma must be > 0. Betas can be any real number.
  lower_bounds <- c(1e-6, 1e-6, rep(-Inf, 7)) # c, sigma > 0
  upper_bounds <- rep(Inf, 9)
  
  for (i in 1:n_sim) {
    # Generate covariates
    X <- matrix(rnorm(n * 6), ncol = 6) # X1 to X6
    
    # Calculate linear predictor (eta = beta0 + beta1*X1 + ... + beta6*X6)
    eta_true <- true_params$beta["beta0"] + X %*% true_params$beta[paste0("beta", 1:6)]
    
    # Generate true event times from LEPS distribution
    T <- rLEPS(n, true_params$c, true_params$sigma, eta_true)
    
    # Check for invalid T values (NaN/Inf from rLEPS if input was problematic)
    if (any(!is.finite(T))) {
      failure_count <- failure_count + 1
      next # Skip this simulation if generated times are invalid
    }
    
    # Generate censoring times (exponential, ~40% censoring)
    rate_val <- 1 / (quantile(T, 0.6, na.rm = TRUE) + .Machine$double.eps)
    C <- rexp(n, rate = rate_val)
    
    # Determine observed time and status
    time <- pmin(T, C)
    status <- as.numeric(T <= C) # 1 if event, 0 if censored
    
    # Handle cases where all events are censored or other issues making fitting impossible
    if (sum(status) == 0 || length(unique(time)) < length(initial_values) / 2) {
      failure_count <- failure_count + 1
      next
    }
    
    # Attempt to fit the LEPS model using optim
    # Using "L-BFGS-B" for box constraints (lower/upper bounds)
    fit <- try(optim(
      par = initial_values,
      fn = neg_log_likelihood_leps,
      time = time,
      status = status,
      X = X,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(fnscale = 1) # fnscale = 1 means we are minimizing, not maximizing
    ), silent = TRUE)
    
    # Check for fitting errors or non-convergence
    if (inherits(fit, "try-error") || fit$convergence != 0 || any(is.na(fit$par)) || length(fit$par) != 9) {
      failure_count <- failure_count + 1
      next # Skip this iteration if fit failed or coefficients are invalid
    }
    
    # Store estimated parameters
    all_estimates[i, ] <- fit$par
  }
  
  # Filter out rows with NA estimates (from failed simulations)
  valid_rows_idx <- apply(all_estimates, 1, function(x) all(is.finite(x)))
  clean_estimates <- all_estimates[valid_rows_idx, , drop = FALSE]
  
  # True parameter values for comparison (all 9 parameters)
  true_all_params <- c(true_params$c, true_params$sigma, true_params$beta)
  
  # Calculate AE, Bias, MSE for all 9 parameters
  AE_all <- colMeans(clean_estimates, na.rm = TRUE)
  Bias_all <- AE_all - true_all_params
  MSE_all <- colMeans((clean_estimates - matrix(rep(true_all_params, each = nrow(clean_estimates)), ncol = 9))^2, na.rm = TRUE)
  
  return(list(
    AE = AE_all,
    Bias = Bias_all,
    MSE = MSE_all,
    failures = failure_count # Return the total failure count
  ))
}

# -------- Run Simulations for Sample Sizes --------
sample_sizes <- c(50, 100, 300, 600, 1000)
# Update param_names to match the order of 'params' in neg_log_likelihood_leps
param_names <- c("c", "sigma", paste0("beta", 0:6)) # Total 9 parameters
table_data <- data.frame()

for (n in sample_sizes) {
  message(paste("Running simulation for n =", n, "..."))
  # Changed to the new simulation function
  sim_result <- simulate_leps_mle(n, true_params, n_sim = 100)
  
  current_failures <- sim_result$failures
  
  for (j in seq_along(param_names)) {
    table_data <- rbind(table_data, data.frame(
      Sample_Size = n,
      Parameter = param_names[j],
      AE = round(sim_result$AE[j], 4),
      Bias = round(sim_result$Bias[j], 4),
      MSE = round(sim_result$MSE[j], 6),
      Failures = ifelse(j == 1, current_failures, NA) # Only assign to the first row for each n
    ))
  }
}

# Fill down the Failures column for display purposes
table_data <- table_data %>%
  group_by(Sample_Size) %>%
  fill(Failures, .direction = "downup") %>%
  ungroup()

# -------- Save Results as CSV --------
write.csv(table_data, "Simulation_Results_LEPS_MLE.csv", row.names = FALSE)
print("Simulation_Results_LEPS_MLE.csv created.")

# -------- Prepare Data for Plotting --------
output_for_plots <- table_data # Now all parameters will have AE/Bias/MSE

output_for_plots$Parameter <- factor(output_for_plots$Parameter, levels = param_names)

# -------- Create Expression Labels for Plotting --------
param_labels <- c(
  expression(c),
  expression(sigma),
  expression(beta[0]),
  expression(beta[1]),
  expression(beta[2]),
  expression(beta[3]),
  expression(beta[4]),
  expression(beta[5]),
  expression(beta[6])
)
names(param_labels) <- param_names

# Ensure labels only include the parameters present in the filtered data (which is all now)
plot_param_labels <- param_labels[names(param_labels) %in% unique(output_for_plots$Parameter)]

# -------- Plot and Save Bias --------
pdf("Bias_Plot_LEPS_MLE.pdf", width = 9, height = 6)
ggplot(output_for_plots, aes(x = as.factor(Sample_Size), y = Bias, group = Parameter, color = Parameter)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scales::hue_pal()(length(plot_param_labels)), labels = plot_param_labels) +
  labs(title = "Bias of LEPS MLEs", x = "Sample Size", y = "Bias", color = "Parameter") +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
dev.off()
print("Bias_Plot_LEPS_MLE.pdf created.")

# -------- Plot and Save MSE --------
pdf("MSE_Plot_LEPS_MLE.pdf", width = 9, height = 6)
ggplot(output_for_plots, aes(x = as.factor(Sample_Size), y = MSE, group = Parameter, color = Parameter)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scales::hue_pal()(length(plot_param_labels)), labels = plot_param_labels) +
  labs(title = "MSE of LEPS MLEs", x = "Sample Size", y = "MSE", color = "Parameter") +
  theme_minimal(base_size = 14)
dev.off()
print("MSE_Plot_LEPS_MLE.pdf created.")