library(rBeta2009)
library(glasso)
library(mvtnorm)

# K-means initialization function
initialize_fuc <- function(data, K, n.start = 100) {
  n <- nrow(data)
  p <- ncol(data)
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, centers = K, nstart = n.start)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size / n
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  
  for (k in 1:K) {
    Mu[k, ] <- colMeans(data[memb == k, , drop = FALSE])
    S[, , k] <- cov(data[memb == k, , drop = FALSE]) + diag(1e-6, p)  # Adding small value to avoid singularity
    Theta[, , k] <- solve(S[, , k])
  }
  
  list(prob = prob, Mu = Mu, Theta = Theta, S = S, memb = memb)
}

# EM algorithm with penalized complete log-likelihood and K-means initialization
fit_cluster_model_with_reg <- function(data, 
                                       G, 
                                       lambda, 
                                       rho, 
                                       max_iter = 250, 
                                       tol = 1e-3,
                                       verbose = FALSE) {
  n <- nrow(data)
  p <- ncol(data)
  
  # Initialize parameters using K-means
  init <- initialize_fuc(data, G, n.start = 10)
  pi_k <- init$prob           # Initial proportions
  mu_k <- init$Mu             # Initial means
  Theta_k <- init$Theta       # Initial precision matrices
  Sigma_k <- init$S           # Initial covariance matrices
  
  # Precompute regularization parameter for glasso
  rho_k <- 2 * rho * n / G
  log_likelihood <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    log_resp <- matrix(0, n, G)
    for (k in 1:G) {
      log_resp[, k] <- log(pi_k[k]) + dmvnorm(data, mean = mu_k[k, ], 
                                              sigma = Sigma_k[, , k], log = TRUE)
    }
  
    max_log_resp <- apply(log_resp, 1, max)
    log_resp <- sweep(log_resp, 1, max_log_resp, FUN = "-")
    resp <- exp(log_resp)

    resp <- resp / rowSums(resp)
    t_ik <- resp
    
    ## M-step:
    # Update proportions
    n_k <- colSums(t_ik)
    pi_k <- n_k / n
    
    # Update means and covariance for each cluster
    for (k in 1:G) {
      # Update means
      mu_k[k, ] <- colSums(t_ik[, k] * data) / n_k[k]
      
      # Center the data around the updated mean
      centered_data <- sweep(data, 2, mu_k[k, ], FUN = "-")
      
      # Update means with penalization using the current precision matrix
      for (j in 1:p) {
        # Calculate the threshold condition for variable j in cluster k
        cross_terms <- numeric(n)
        
        for (v in 1:p) {
          if (v != j) {
            cross_terms <- cross_terms + centered_data[, v] * Theta_k[j, v, k]
          }
        }
        
        # Calculate the threshold value
        threshold_val <- abs(sum(t_ik[, k] * cross_terms) + sum(t_ik[, k] * centered_data[, j] * Theta_k[j, j, k]))
        
        if (threshold_val <= lambda) {
          # Set mean to 0 if threshold condition is met
          mu_k[k, j] <- 0
        } else {
          # Otherwise, update mean using penalization with precision matrix
          numerator <- sum(t_ik[, k] * data[, j] * Theta_k[j, j, k]) - 
            sum(sapply(1:p, function(v) {
              if (v != j) {
                sum(t_ik[, k] * centered_data[, v] * Theta_k[j, v, k])
              } else {
                0
              }
            })) + lambda * sign(mu_k[k, j])
          denominator <- sum(t_ik[, k] * Theta_k[j, j, k])
          
          mu_k[k, j] <- numerator / denominator
        }
      }
      
      # Compute the weighted covariance matrix with the updated means
      centered_data <- sweep(data, 2, mu_k[k, ], FUN = "-")
      S_k <- t(centered_data) %*% (centered_data * t_ik[, k]) / n_k[k]
      
      # Apply graphical lasso to get precision matrix (Theta_k) and update covariance matrix
      glasso_fit <- glasso::glasso(S_k, rho_k)
      Theta_k[, , k] <- glasso_fit$wi  # Precision matrix
      Sigma_k[, , k] <- solve(Theta_k[, , k])  # Covariance matrix for output
    }
    # Compute penalized log-likelihood for convergence check
    penalized_loglik <- 0
    for (i in 1:n) {
      for (k in 1:G) {
        penalized_loglik <- penalized_loglik + t_ik[i, k] * (log(pi_k[k]) + dmvnorm(data[i, ], mean = mu_k[k, ], sigma = Sigma_k[, , k], log = TRUE))
      }
    }
    penalized_loglik <- penalized_loglik - lambda * sum(abs(mu_k)) - rho * sum(sapply(1:G, function(k) sum(abs(Theta_k[, , k]))))
    log_likelihood[iter] <- penalized_loglik
    
    # Check for convergence
    if (iter > 1 && abs(log_likelihood[iter] - log_likelihood[iter - 1]) <= tol) {
      if (verbose) {
        cat("Converged after", iter, "iterations.\n")
      }
      break
    }
  }
  
  # Return the final parameters and log-likelihood progression
  list(pi_k = pi_k, mu_k = mu_k, Sigma_k = Sigma_k, log_likelihood = log_likelihood[1:iter])
}


# Function to generate true parameters and synthetic data
generate_test_data <- function(n, p, G, seed = 123) {
  set.seed(seed)
  
  # Generate true parameters
  true_params <- list()
  
  # True mixing proportions
  true_params$pi_k <- rep(1/G, G)
  
  # True means - well-separated clusters
  true_params$mu_k <- matrix(0, G, p)
  for(k in 1:G) {
    true_params$mu_k[k,] <- rnorm(p, mean = k*2, sd = 0.5)
  }
  
  # True precision and covariance matrices
  true_params$Theta_k <- array(0, dim = c(p, p, G))
  true_params$Sigma_k <- array(0, dim = c(p, p, G))
  
  for(k in 1:G) {
    # Generate random positive definite matrix
    temp <- matrix(rnorm(p*p), p, p)
    temp <- t(temp) %*% temp + diag(p)  # Ensure positive definiteness
    true_params$Sigma_k[,,k] <- temp / max(abs(temp)) # Scale entries
    true_params$Theta_k[,,k] <- solve(true_params$Sigma_k[,,k])
  }
  
  # Generate synthetic data
  data <- matrix(0, n, p)
  true_labels <- numeric(n)
  
  for(i in 1:n) {
    # Sample cluster
    cluster <- sample(1:G, 1, prob = true_params$pi_k)
    true_labels[i] <- cluster
    
    # Generate observation from corresponding multivariate normal
    data[i,] <- rmvnorm(1, 
                        mean = true_params$mu_k[cluster,], 
                        sigma = true_params$Sigma_k[,,cluster])
  }
  
  list(data = data, 
       true_labels = true_labels, 
       true_params = true_params)
}

# Modified initialization function that adds fixed small value to true parameters
initialize_with_truth <- function(true_params, epsilon = 1e-3) {
  G <- length(true_params$pi_k)
  p <- ncol(true_params$mu_k)
  
  init_params <- list()
  
  # Initialize proportions: add epsilon while ensuring sum = 1
  raw_props <- true_params$pi_k + epsilon
  init_params$prob <- raw_props / sum(raw_props)
  
  # Initialize means: add epsilon to all elements
  init_params$Mu <- true_params$mu_k + epsilon
  
  # Initialize precision and covariance matrices
  init_params$Theta <- array(0, dim = dim(true_params$Theta_k))
  init_params$S <- array(0, dim = dim(true_params$Sigma_k))
  
  for(k in 1:G) {
    # Add epsilon to diagonal elements of covariance matrix to ensure positive definiteness
    init_params$S[,,k] <- true_params$Sigma_k[,,k] + diag(epsilon, p)
    init_params$Theta[,,k] <- solve(init_params$S[,,k])
  }
  
  init_params
}

# Function to compute parameter estimation error
compute_estimation_error <- function(true_params, estimated_params) {
  errors <- list()
  
  # Error in mixing proportions
  errors$pi_error <- mean(abs(true_params$pi_k - estimated_params$pi_k))
  
  # Error in means
  errors$mu_error <- mean(abs(true_params$mu_k - estimated_params$mu_k))
  
  # Error in covariance matrices
  G <- length(true_params$pi_k)
  cov_errors <- numeric(G)
  for(k in 1:G) {
    cov_errors[k] <- norm(true_params$Sigma_k[,,k] - estimated_params$Sigma_k[,,k], "F")
  }
  errors$sigma_error <- mean(cov_errors)
  
  errors
}

# Test function
test_em_convergence <- function(n = 500, p = 3, G = 2, lambda = 0.1, rho = 0.1, 
                                epsilon = 1e-3, seed = 123) {
  # Generate test data
  test_data <- generate_test_data(n, p, G, seed)
  
  # Create initialization close to true parameters
  init_params <- initialize_with_truth(test_data$true_params, epsilon)
  
  # Run EM algorithm with the modified initialization
  result <- fit_cluster_model_with_reg(test_data$data, G, lambda, rho, 
                                       max_iter = 250, tol = 1e-6, verbose = TRUE)
  
  # Compute estimation errors
  errors <- compute_estimation_error(test_data$true_params, result)
  
  # Return results
  list(
    errors = errors,
    log_likelihood = result$log_likelihood,
    true_params = test_data$true_params,
    estimated_params = result,
    data = test_data$data,
    true_labels = test_data$true_labels,
    init_params = init_params  # Added to check initialization
  )
}

# Example usage
run_test_both_epsilon <- function() {
  # Test with epsilon = 1e-3
  cat("Testing with epsilon = 1e-3\n")
  result1 <- test_em_convergence(epsilon = 1e-3)
  print(result1$errors)
  
  # Test with epsilon = 1e-4
  cat("\nTesting with epsilon = 1e-4\n")
  result2 <- test_em_convergence(epsilon = 1e-4)
  print(result2$errors)
  
  list(result1 = result1, result2 = result2)
}

result <- run_test_both_epsilon()
