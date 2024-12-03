library(rBeta2009)
library(glasso)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(gridExtra)

# K-means initialization function (keeping your original implementation)
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
    S[, , k] <- cov(data[memb == k, , drop = FALSE]) + diag(1e-6, p)
    Theta[, , k] <- solve(S[, , k])
  }
  
  list(prob = prob, Mu = Mu, Theta = Theta, S = S, memb = memb)
}

# Modified EM algorithm with monitoring and diagnostics
fit_cluster_model_with_reg <- function(data, 
                                       G, 
                                       lambda, 
                                       rho, 
                                       max_iter = 100, 
                                       tol = 1e-6,
                                       verbose = TRUE) {
  n <- nrow(data)
  p <- ncol(data)
  
  # Initialize parameters using K-means
  init <- initialize_fuc(data, G, n.start = 10)
  pi_k <- init$prob
  mu_k <- init$Mu
  Theta_k <- init$Theta
  Sigma_k <- init$S
  
  # Storage for monitoring
  log_likelihood <- numeric(max_iter)
  means_history <- array(NA, dim = c(max_iter, G, p))
  pi_history <- matrix(NA, max_iter, G)
  
  # Precompute regularization parameter
  rho_k <- 2 * rho * n / G
  
  for (iter in 1:max_iter) {
    ## Store current state for monitoring
    means_history[iter,,] <- mu_k
    pi_history[iter,] <- pi_k
    
    ## E-step with numerical stability
    log_resp <- matrix(0, n, G)
    for (k in 1:G) {
      log_resp[, k] <- log(pi_k[k]) + dmvnorm(data, mean = mu_k[k, ], sigma = Sigma_k[, , k], log = TRUE)
    }
    
    # Numerical stable normalization
    max_log_resp <- apply(log_resp, 1, max)
    resp <- exp(log_resp - max_log_resp)
    resp <- resp / rowSums(resp)
    t_ik <- pmax(resp, 1e-10)  # Add small constant for numerical stability
    
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
      # Add monitoring of condition numbers
      if(verbose && iter %% 10 == 0) {
        cond_num <- kappa(Sigma_k[,,k])
        if(cond_num > 1e10) {
          warning(sprintf("High condition number (%.2e) in cluster %d", cond_num, k))
        }
      }
    }
    
    # Compute penalized log-likelihood with safeguards
    penalized_loglik <- try({
      ll <- 0
      for (i in 1:n) {
        for (k in 1:G) {
          ll <- ll + t_ik[i, k] * (log(pi_k[k]) + 
                                     dmvnorm(data[i, ], mean = mu_k[k, ], 
                                             sigma = Sigma_k[, , k], log = TRUE))
        }
      }
      ll - lambda * sum(abs(mu_k)) - rho * sum(sapply(1:G, function(k) sum(abs(Theta_k[, , k]))))
    }, silent = TRUE)
    
    if(inherits(penalized_loglik, "try-error")) {
      warning("Error in log-likelihood computation, using previous value")
      penalized_loglik <- if(iter > 1) log_likelihood[iter-1] else -Inf
    }
    
    log_likelihood[iter] <- penalized_loglik
    
    # Check convergence with relative change
    if (iter > 1) {
      rel_change <- abs(log_likelihood[iter] - log_likelihood[iter-1]) / 
        abs(log_likelihood[iter-1])
      
      if(verbose && iter %% 10 == 0) {
        cat(sprintf("Iteration %d: log-likelihood = %.4f, relative change = %.6f\n",
                    iter, log_likelihood[iter], rel_change))
      }
      
      if(rel_change < tol) {
        if(verbose) cat("Converged after", iter, "iterations\n")
        break
      }
    }
  }
  
  # Calculate clustering metrics
  final_clusters <- max.col(resp)
  silhouette_avg <- cluster::silhouette(final_clusters, dist(data))[, 3]
  
  # Create diagnostic plots if verbose
  if(verbose) {
    # Plot 1: Log-likelihood progression
    p1 <- ggplot(data.frame(iteration = 1:iter, ll = log_likelihood[1:iter]), 
                 aes(x = iteration, y = ll)) +
      geom_line() +
      theme_minimal() +
      labs(title = "Log-likelihood progression")
    
    # Plot 2: Means convergence
    means_df <- melt(means_history[1:iter,,])
    colnames(means_df) <- c("iteration", "cluster", "dimension", "value")
    p2 <- ggplot(means_df, aes(x = iteration, y = value, color = factor(cluster))) +
      geom_line() +
      facet_wrap(~dimension) +
      theme_minimal() +
      labs(title = "Means convergence")
    
    # Plot 3: Mixing proportions
    pi_df <- melt(pi_history[1:iter,])
    colnames(pi_df) <- c("iteration", "cluster", "value")
    p3 <- ggplot(pi_df, aes(x = iteration, y = value, color = factor(cluster))) +
      geom_line() +
      theme_minimal() +
      labs(title = "Mixing proportions convergence")
    
    grid.arrange(p1, p2, p3, ncol = 2)
  }
  
  # Return enriched results
  list(
    pi_k = pi_k,
    mu_k = mu_k,
    Sigma_k = Sigma_k,
    Theta_k = Theta_k,
    log_likelihood = log_likelihood[1:iter],
    responsibilities = resp,
    clustering = final_clusters,
    silhouette = mean(silhouette_avg),
    BIC = -2 * log_likelihood[iter] + log(n) * (G * p + G * p * (p-1)/2),
    AIC = -2 * log_likelihood[iter] + 2 * (G * p + G * p * (p-1)/2),
    iterations = iter,
    means_history = means_history[1:iter,,],
    pi_history = pi_history[1:iter,]
  )
}

# Function to run multiple trials with different parameters
tune_parameters <- function(data, G, lambda_grid, rho_grid, n_trials = 5) {
  results <- list()
  for(lambda in lambda_grid) {
    for(rho in rho_grid) {
      cat(sprintf("\nTrying lambda = %.3f, rho = %.3f\n", lambda, rho))
      
      # Run multiple trials
      trial_results <- replicate(n_trials, {
        fit <- fit_cluster_model_with_reg(data, G, lambda, rho, verbose = FALSE)
        c(BIC = fit$BIC, 
          silhouette = fit$silhouette, 
          loglik = tail(fit$log_likelihood, 1))
      })
      
      # Store average results
      results[[sprintf("lambda%.3f_rho%.3f", lambda, rho)]] <- 
        list(
          lambda = lambda,
          rho = rho,
          BIC_mean = mean(trial_results["BIC",]),
          BIC_sd = sd(trial_results["BIC",]),
          silhouette_mean = mean(trial_results["silhouette",]),
          loglik_mean = mean(trial_results["loglik",])
        )
    }
  }
  results
}


# Example usage with parameter tuning
if(TRUE) {  # Set to TRUE to run
  # Generate synthetic data (as in your original code)
  set.seed(123)
  n <- 300
  p <- 5
  G <- 3
  
  # Generate more separated clusters for testing
  true_mu <- list(
    c(3, 3, 0, 0, 0),
    c(-3, -3, 0, 0, 0),
    c(0, 0, 3, -3, 0)
  )
  
  true_sigma <- list(
    diag(0.5, p),
    diag(0.5, p),
    diag(0.5, p)
  )
  
  # Generate data
  true_pi <- c(0.3, 0.4, 0.3)
  data <- matrix(0, n, p)
  cluster_labels <- numeric(n)
  start <- 1
  for(k in 1:G) {
    n_k <- round(n * true_pi[k])
    end <- start + n_k - 1
    data[start:end,] <- rmvnorm(n_k, true_mu[[k]], true_sigma[[k]])
    cluster_labels[start:end] <- k
    start <- end + 1
  }
  
  # Scale data
  scaled_data <- scale(data)
  
  # Tune parameters
  lambda_grid <- c(0.01, 0.05, 0.1)
  rho_grid <- c(0.01, 0.05, 0.1)
  tuning_results <- tune_parameters(scaled_data, G, lambda_grid, rho_grid)
  
  # Find best parameters
  BICs <- sapply(tuning_results, function(x) x$BIC_mean)
  best_params <- tuning_results[[which.min(BICs)]]
  
  # Fit final model with best parameters
  final_model <- fit_cluster_model_with_reg(scaled_data, G, 
                                            best_params$lambda, 
                                            best_params$rho,
                                            verbose = TRUE)
  
  # Compare with true clusters
  table(final_model$clustering, cluster_labels)
}