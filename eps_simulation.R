# Set working directory (adjust the path as needed)
setwd("C:\\Users\\1012 G2\\Documents\\R files\\EPS")

# --- Parameters for EPS distribution (Global for convenience in this script) ---
# These are the true parameters used for data generation
TRUE_c     = 2.0
TRUE_theta = 1.5
TRUE_alpha = 0.5

# --- EPS distribution PDF (f_eps) ---
# This is the PDF of the EPS distribution based on your initial definition.
# It corresponds to equation (3.6) or the form used in your normalization.
f_eps = function(x, c_val = TRUE_c, theta_val = TRUE_theta, alpha_val = TRUE_alpha) {
  # Return 0 for non-positive x, as the distribution is defined for x > 0
  if (any(x <= 0)) {
    # If any x is non-positive, return 0 for those, and NA for others if they are valid
    # This ensures the function can handle vectors with mixed valid/invalid inputs
    return(ifelse(x > 0, NA, 0))
  }
  
  part1 = c_val * alpha_val * theta_val^2 / (theta_val^2 + 1)
  part2 = (theta_val + x^alpha_val) * x^(alpha_val - 1) * exp(-theta_val * x^alpha_val)
  
  # This part directly reflects the {1 - [...]}^(c-1) from your PDF definition
  inner_bracket = (1 + (theta_val * x^alpha_val) / (theta_val^2 + 1)) * exp(-theta_val * x^alpha_val)
  part3_base = 1 - inner_bracket
  
  # Handle cases where part3_base might be negative or zero.
  # If it's non-positive, the PDF is ill-defined or 0 in that region.
  # Using ifelse to handle vector inputs correctly.
  result = ifelse(part3_base > 0, part1 * part2 * part3_base^(c_val - 1), 0)
  
  # Ensure no NaNs or Infs propagate from other calculations within the result
  result[is.nan(result) | is.infinite(result) | result < 0] = 0
  
  return(result)
}

# --- g(x) PDF (for acceptance-rejection method) ---
# This is the PDF of the Exponential Power distribution, used as a proposal.
g_eps = function(x, theta_val = TRUE_theta, alpha_val = TRUE_alpha) {
  if (any(x <= 0)) return(0) # Return 0 for non-positive x
  return(theta_val * alpha_val * x^(alpha_val - 1) * exp(-theta_val * x^alpha_val))
}

# --- Quantile function for g(x) (inverse of g_eps) ---
# This generates random numbers from the g(x) distribution.
rw_eps = function(n, theta_val = TRUE_theta, alpha_val = TRUE_alpha) {
  u = runif(n)
  # Inverse CDF for G(x) = 1 - exp(-theta * x^alpha) is x = (-log(1 - u) / theta)^(1 / alpha)
  Q = (-log(1 - u) / theta_val)^(1 / alpha_val)
  return(Q)
}

# --- Acceptance-rejection sampling for EPS distribution ---
# Generates 'n' random variates from f_eps using g_eps as a proposal.
rchrisjerry_eps = function(f, g, M_val, n, theta_params = list(c = TRUE_c, theta = TRUE_theta, alpha = TRUE_alpha)) {
  y = numeric(n)
  n_accepts = 0
  n_proposals = 0 # To track efficiency
  
  # Ensure M_val is valid
  if (is.infinite(M_val) || is.na(M_val) || M_val <= 0) {
    stop("Invalid M value for acceptance-rejection sampling. Check f_eps/g_eps ratio.")
  }
  
  while (n_accepts < n) {
    x <- rw_eps(1, theta_params$theta, theta_params$alpha) # Propose from g(x)
    u <- runif(1) # Uniform random number for acceptance
    
    # Calculate w(x) = f(x) / (M * g(x))
    # Pass all required parameters to f() and g()
    fx = f(x, theta_params$c, theta_params$theta, theta_params$alpha)
    gx = g(x, theta_params$theta, theta_params$alpha)
    
    # Handle potential issues with PDF values
    if (is.na(fx) || is.infinite(fx) || fx < 0 || is.na(gx) || is.infinite(gx) || gx <= 0) {
      n_proposals <- n_proposals + 1 # Still count as a proposal attempt
      next # Skip this proposed value if PDF is invalid
    }
    
    w <- fx / (M_val * gx)
    
    # Check for NA/NaN/Inf in w
    if (!is.finite(w)) {
      n_proposals <- n_proposals + 1
      next
    }
    
    if (u <= w) {
      n_accepts <- n_accepts + 1
      y[n_accepts] <- x
    }
    n_proposals <- n_proposals + 1
  }
  # print(paste("Acceptance rate:", n / n_proposals)) # Uncomment to check efficiency
  return(y)
}


# --- Log-likelihood function for Maximum Likelihood Estimation ---
# It takes a vector 'theta_par' containing c, theta, alpha in that order,
# and 'n_obs' which is the number of observations (nobs from MCS_eps).
logLikk_eps = function(theta_par, n_obs) { # 'n_obs' added as an argument
  # Extract parameters from the input vector
  c_val     = theta_par[1]
  theta_val = theta_par[2]
  alpha_val = theta_par[3]
  
  # --- Parameter Constraints ---
  # Ensure parameters are positive. Returning Inf for invalid parameters.
  if (c_val <= 0 || theta_val <= 0 || alpha_val <= 0) {
    return(Inf)
  }
  
  # --- Log-likelihood Calculation (Element-wise for summation later) ---
  # y_eps is the data (global variable, assigned in MCS_eps)
  
  # Term 1: n * [log(c) + log(alpha) + 2*log(theta) - log(theta^2 + 1)]
  # 'n_obs' is used here for the count of observations.
  log_const_part = n_obs * (log(c_val) + log(alpha_val) + 2 * log(theta_val) - log(theta_val^2 + 1))
  
  # Term 2: -theta * sum(x^alpha) => -theta_val * y_eps^alpha_val for each observation
  log_exp_term = -theta_val * y_eps^alpha_val
  
  # Term 3: sum(log(theta + x^alpha)) => log(theta_val + y_eps^alpha_val) for each observation
  log_sum_theta_x_alpha = log(theta_val + y_eps^alpha_val)
  
  # Term 4: (alpha - 1) * sum(log(x_i)) => (alpha_val - 1) * log(y_eps) for each observation
  log_alpha_minus_1_x_log_x = (alpha_val - 1) * log(y_eps)
  
  # Term 5: (c - 1) * sum(log{1 - [1 + (theta * x^alpha) / (theta^2 + 1)] * exp(-theta * x^alpha)})
  inner_frac_term = (theta_val * y_eps^alpha_val) / (theta_val^2 + 1)
  exp_term = exp(-theta_val * y_eps^alpha_val)
  log_base_arg = 1 - (1 + inner_frac_term) * exp_term
  
  # --- Numerical Stability Checks for Logarithms ---
  # If any argument to log() is non-positive, the log-likelihood is undefined.
  # Return Inf to prevent 'optim' from exploring these invalid parameter regions.
  if (any(log_base_arg <= 0) || any(y_eps <= 0) || any(theta_val + y_eps^alpha_val <= 0)) {
    return(Inf)
  }
  
  log_c_minus_1_part = (c_val - 1) * log(log_base_arg)
  
  # Sum all the element-wise log-likelihood contributions
  total_log_likelihood_elements = log_exp_term + log_sum_theta_x_alpha + log_alpha_minus_1_x_log_x + log_c_minus_1_part
  
  # Check for non-finite values in the element-wise log-likelihoods
  if (any(!is.finite(total_log_likelihood_elements))) {
    return(Inf)
  }
  
  # Combine constant part with summed components
  total_log_likelihood = log_const_part + sum(total_log_likelihood_elements)
  
  # Final check for non-finite total log-likelihood
  if (!is.finite(total_log_likelihood)) {
    return(Inf)
  }
  
  # Return the negative log-likelihood for minimization (optim seeks minimum)
  return(-total_log_likelihood)
}


# --- Monte Carlo Simulation Function ---
# Simulates data, estimates parameters, and calculates performance metrics.
MCS_eps = function(nrep = 1000, nobs, seed = 1986,
                   c_true = TRUE_c, theta_true = TRUE_theta, alpha_true = TRUE_alpha) {
  
  set.seed(seed)
  
  # Initialize vectors to store MLEs from each replication
  mlec_eps   = numeric(nrep)
  mletheta   = numeric(nrep)
  mlealpha   = numeric(nrep)
  
  faultcounter = 0 # Counter for non-convergent optimization runs
  
  # Pre-calculate M for acceptance-rejection.
  # This M value should ideally be the supremum of f(x)/g(x) over the domain.
  # Using a sufficiently large sample from g_eps to estimate M.
  temp_x_for_M = rw_eps(10000, theta_true, alpha_true)
  ratios = f_eps(temp_x_for_M, c_true, theta_true, alpha_true) / g_eps(temp_x_for_M, theta_true, alpha_true)
  
  # Remove non-finite ratios before finding max
  ratios_finite = ratios[is.finite(ratios) & ratios > 0] 
  
  M_val = max(ratios_finite, na.rm = TRUE) 
  
  # Add a small buffer to M for safety, and handle problematic M values
  if (M_val <= 0 || !is.finite(M_val)) {
    warning("Calculated M for acceptance-rejection is problematic. Setting to a default value.")
    M_val = 10 # Fallback, adjust if needed for your specific distribution parameters
  } else {
    M_val = M_val * 1.05 # Add a small buffer
  }
  
  for (i in 1:nrep) {
    # Generate data using the acceptance-rejection method
    # y_eps <<- makes it available in the global environment for logLikk_eps
    y_eps <<- rchrisjerry_eps(f_eps, g_eps, M_val, nobs,
                              theta_params = list(c = c_true, theta = theta_true, alpha = alpha_true))
    
    # --- Initial guess for optimization ---
    # Using a slightly perturbed version of the true parameters for initial guess.
    # This often helps 'optim' find the optimum, especially with complex likelihoods.
    initial_guess = c(c_true + runif(1, -0.5, 0.5),
                      theta_true + runif(1, -0.5, 0.5),
                      alpha_true + runif(1, -0.5, 0.5))
    # Ensure initial guess is positive, as parameters must be > 0.
    initial_guess[initial_guess <= 0] = 0.1 
    
    # Perform optimization using Nelder-Mead (derivative-free method)
    # Pass 'nobs' to logLikk_eps using the '...' argument in optim
    ir = optim(initial_guess, fn = logLikk_eps,
               n_obs = nobs, # <--- Passing nobs here to the log-likelihood function
               method = "Nelder-Mead",
               control = list(maxit = 10000, reltol = 1e-8)) # Increased iterations and tolerance for better convergence
    
    # --- Store results or count failures ---
    # Assign NA to estimates if optimization failed (convergence code != 0).
    if (ir$convergence == 0) { # convergence == 0 means successful convergence
      mlec_eps[i] = ir$par[1]
      mletheta[i] = ir$par[2]
      mlealpha[i] = ir$par[3]
    } else {
      faultcounter = faultcounter + 1
      mlec_eps[i] = NA # Mark as NA if not converged
      mletheta[i] = NA
      mlealpha[i] = NA
    }
  }
  
  # --- Calculate simulation metrics, ignoring NA values (failed runs) ---
  # Average Estimates (mean of successful estimates)
  c_est_mean     = mean(mlec_eps, na.rm = TRUE)
  theta_est_mean = mean(mletheta, na.rm = TRUE)
  alpha_est_mean = mean(mlealpha, na.rm = TRUE)
  
  # Bias (Average Estimate - True Value)
  bias_c     = c_est_mean - c_true
  bias_theta = theta_est_mean - theta_true
  bias_alpha = alpha_est_mean - alpha_true
  
  # Mean Squared Error (MSE)
  mse_c     = mean((mlec_eps - c_true)^2, na.rm = TRUE)
  mse_theta = mean((mletheta - theta_true)^2, na.rm = TRUE)
  mse_alpha = mean((mlealpha - alpha_true)^2, na.rm = TRUE)
  
  # Average Absolute Error (AE / Absolute Bias) - optional, as Bias usually suffices
  abs_err_c     = mean(abs(mlec_eps - c_true), na.rm = TRUE)
  abs_err_theta = mean(abs(mletheta - theta_true), na.rm = TRUE)
  abs_err_alpha = mean(abs(mlealpha - alpha_true), na.rm = TRUE)
  
  # Return all results in a data frame
  return(data.frame(
    nobs = nobs,
    nrep = nrep,
    true_c = c_true,
    true_theta = theta_true,
    true_alpha = alpha_true,
    c_avg_est = c_est_mean,
    theta_avg_est = theta_est_mean,
    alpha_avg_est = alpha_est_mean,
    bias_c = bias_c,
    bias_theta = bias_theta,
    bias_alpha = bias_alpha,
    mse_c = mse_c,
    mse_theta = mse_theta,
    mse_alpha = mse_alpha,
    abs_err_c = abs_err_c,
    abs_err_theta = abs_err_theta,
    abs_err_alpha = abs_err_alpha,
    failures = faultcounter
  ))
}

# --- Run simulations for different sample sizes ---
# It's good practice to store the results of each run.
results_50 = MCS_eps(nobs = 50, c_true = TRUE_c, theta_true = TRUE_theta, alpha_true = TRUE_alpha)
print("Simulation for n=50 completed.")
print(results_50)

results_100 = MCS_eps(nobs = 100, c_true = TRUE_c, theta_true = TRUE_theta, alpha_true = TRUE_alpha)
print("Simulation for n=100 completed.")
print(results_100)

results_300 = MCS_eps(nobs = 300, c_true = TRUE_c, theta_true = TRUE_theta, alpha_true = TRUE_alpha)
print("Simulation for n=300 completed.")
print(results_300)

results_600 = MCS_eps(nobs = 600, c_true = TRUE_c, theta_true = TRUE_theta, alpha_true = TRUE_alpha)
print("Simulation for n=600 completed.")
print(results_600)

# You can combine these results into a single data frame if needed for a final table:
# all_results = rbind(results_50, results_100, results_300, results_600)
# print(all_results)