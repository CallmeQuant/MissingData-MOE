library(mvtnorm)  
library(glasso)
# K-means initialization function
initialize_func <- function(data, K, n.start = 100) {
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

# Main EM algorithm function with regularization
fit_cluster_model_with_reg <- function(data, 
                                       G, 
                                       lambda, 
                                       rho, 
                                       max_iter = 250, 
                                       tol = 1e-3,
                                       verbose = FALSE) {
  # Initialize dimensions
  n <- nrow(data)
  p <- ncol(data)
  
  # Initialize parameters using K-means
  init <- initialize_func(data, G, n.start = 10)
  pi_k <- init$prob           # Initial proportions
  mu_k <- init$Mu             # Initial means
  Theta_k <- init$Theta       # Initial precision matrices
  Sigma_k <- init$S          # Initial covariance matrices
  
  # Storage for log-likelihood
  log_likelihood <- numeric(max_iter)
  
  # Main EM loop
  for (iter in 1:max_iter) {
    # E-step: Calculate conditional probabilities
    log_resp <- matrix(0, n, G)
    for (k in 1:G) {
      # Ensure symmetry in covariance matrix
      Sigma_k[, , k] <- (Sigma_k[, , k] + t(Sigma_k[, , k])) / 2
      
      # Calculate log probabilities
      det_result <- determinant(Sigma_k[, , k], logarithm = TRUE)
      log_resp[, k] <- log(pi_k[k]) + dmvnorm(data, 
                                              mean = mu_k[k, ],
                                              sigma = Sigma_k[, , k],
                                              log = TRUE)
    }
    
    # Numerical stability adjustment
    max_log_resp <- apply(log_resp, 1, max)
    log_resp <- sweep(log_resp, 1, max_log_resp, FUN = "-")
    resp <- exp(log_resp)
    resp <- resp / rowSums(resp)
    t_ik <- resp
    
    # M-step:
    # Update mixing proportions
    n_k <- colSums(t_ik)
    pi_k <- n_k / n
    
    # Update means and covariance for each cluster
    for (k in 1:G) {
      W <- Theta_k[, , k]          # Precision matrix for cluster k
      MuPrec <- mu_k[k, ]          # Previous mean vector
      MuNew <- numeric(p)          # New mean vector
      
      # Update means with penalization
      for (j in 1:p) {
        # Calculate total absolute value for thresholding
        Tabs <- sum(t_ik[, k] * (rowSums(sweep(data, 2, MuPrec, FUN = "-") * W[, j]) + 
                                   MuPrec[j] * W[j, j]))
        
        if (abs(Tabs) <= lambda) {
          MuNew[j] <- 0  # Threshold condition met
        } else {
          # Calculate penalized update
          T <- sum(t_ik[, k] * data[, j] * W[j, j]) - 
            sum(t_ik[, k] * rowSums(sweep(data, 2, MuPrec, FUN = "-") * W[j, -j])) + 
            lambda * sign(MuPrec[j])
          MuNew[j] <- T / (n_k[k] * W[j, j])
        }
      }
      mu_k[k, ] <- MuNew
      
      # Calculate empirical covariance matrix
      centered_data <- sweep(data, 2, MuNew, FUN = "-")
      S_k <- crossprod(centered_data * sqrt(t_ik[, k])) / n_k[k]
      
      # Update precision and covariance matrices using graphical lasso
      rhotilde <- (2 * rho) / n_k[k]
      glasso_fit <- glasso::glasso(S_k, 
                                   rho = rhotilde,
                                   penalize.diagonal = FALSE,
                                   thr = 0.001,
                                   maxit = 1000)
      
      Theta_k[, , k] <- glasso_fit$wi      # Precision matrix
      Sigma_k[, , k] <- glasso_fit$w       # Covariance matrix
      
      # Enforce symmetry
      Sigma_k[, , k] <- (Sigma_k[, , k] + t(Sigma_k[, , k])) / 2
    }
    
    # Compute penalized log-likelihood
    penalized_loglik <- sum(log(rowSums(resp))) - 
      (lambda * sum(abs(mu_k))) - 
      (rho * sum(sapply(1:G, function(k) sum(abs(Theta_k[, , k])))))
    
    log_likelihood[iter] <- penalized_loglik
    
    # Check convergence
    if (iter > 1) {
      rel_change <- abs(log_likelihood[iter] - log_likelihood[iter - 1]) / 
        abs(log_likelihood[iter - 1])
      abs_change <- abs(log_likelihood[iter] - log_likelihood[iter - 1])
      if (verbose) {
        cat(sprintf("Iteration %d: log-likelihood = %.4f, relative change = %.6f\n",
                    iter, log_likelihood[iter], rel_change))
      }
      if (abs_change < tol) {
        log_likelihood <- log_likelihood[1:iter]
        break
      }
    }
  }
  
  # Return results
  list(
    proportions = pi_k,
    means = mu_k,
    precision_matrices = Theta_k,
    covariance_matrices = Sigma_k,
    responsibilities = t_ik,
    log_likelihood = log_likelihood,
    n_iterations = iter,
    cluster_assignments = apply(t_ik, 1, which.max)
  )
}

match_labels <- function(estimated_means, true_means) {
  G <- nrow(estimated_means)
  best_perm <- 1:G
  min_dist <- Inf
  
  # For 2 clusters, we only need to check 2 permutations
  perms <- list(
    1:G,
    c(2:G, 1)
  )
  
  for (perm in perms) {
    dist <- sum((estimated_means[perm,] - true_means)^2)
    if (dist < min_dist) {
      min_dist <- dist
      best_perm <- perm
    }
  }
  
  return(best_perm)
}

# Test function with true parameters and convergence metrics
test_clustering_convergence <- function() {
  # Set random seed for reproducibility
  set.seed(123)
  
  # Define true parameters
  n <- 200  # Total sample size
  p <- 5    # Number of dimensions
  G <- 2    # Number of clusters
  
  # True mixing proportions
  true_pi <- c(0.4, 0.6)
  
  # True means
  true_mu <- matrix(c(
    0, 0, 0, 0, 0,    # Cluster 1 means
    3, 3, 3, 3, 3     # Cluster 2 means
  ), nrow = G, byrow = TRUE)
  
  # True covariance matrices
  true_Sigma <- array(0, dim = c(p, p, G))
  # Cluster 1: AR(1) structure with rho = 0.7
  rho1 <- 0.7
  true_Sigma[,,1] <- rho1^abs(outer(1:p, 1:p, "-"))
  # Cluster 2: Compound symmetry with variance 1 and correlation 0.3
  rho2 <- 0.3
  true_Sigma[,,2] <- matrix(rho2, p, p) + diag(1-rho2, p)
  
  # True precision matrices
  true_Theta <- array(0, dim = c(p, p, G))
  true_Theta[,,1] <- solve(true_Sigma[,,1])
  true_Theta[,,2] <- solve(true_Sigma[,,2])
  
  # Generate data
  data <- matrix(0, n, p)
  true_labels <- numeric(n)
  
  # Sample cluster assignments
  z <- sample(1:G, n, replace = TRUE, prob = true_pi)
  
  # Generate observations
  for (i in 1:n) {
    cluster <- z[i]
    data[i,] <- MASS::mvrnorm(1, 
                              mu = true_mu[cluster,], 
                              Sigma = true_Sigma[,,cluster])
    true_labels[i] <- cluster
  }
  
  # Fit model with different penalty parameters
  penalties <- list(
    weak = list(lambda = 0.01, rho = 0.01),
    medium = list(lambda = 0.1, rho = 0.1),
    strong = list(lambda = 0.5, rho = 0.5)
  )
  
  results <- list()
  for (penalty_name in names(penalties)) {
    cat(sprintf("\nFitting model with %s penalties...\n", penalty_name))
    
    results[[penalty_name]] <- fit_cluster_model_with_reg(
      data = data,
      G = G,
      lambda = penalties[[penalty_name]]$lambda,
      rho = penalties[[penalty_name]]$rho,
      max_iter = 100,
      tol = 1e-4,
      verbose = TRUE
    )
  }
  
  # Compute evaluation metrics
  evaluate_results <- function(result, penalty_name) {
    # Match labels based on means
    estimated_means <- result$means
    matching <- match_labels(estimated_means, true_mu)
    
    # Reorder estimated parameters
    result$means <- result$means[matching,]
    result$covariance_matrices <- result$covariance_matrices[,,matching]
    result$precision_matrices <- result$precision_matrices[,,matching]
    result$proportions <- result$proportions[matching]
    
    # Adjust cluster assignments
    old_to_new <- numeric(G)
    old_to_new[matching] <- 1:G
    result$cluster_assignments <- old_to_new[result$cluster_assignments]
    
    # Compute metrics
    metrics <- list(
      # Mean estimation error
      mean_error = mean(abs(result$means - true_mu)),
      
      # Proportion estimation error
      prop_error = mean(abs(result$proportions - true_pi)),
      
      # Covariance matrix error (Frobenius norm)
      cov_error = mean(sapply(1:G, function(k) {
        norm(result$covariance_matrices[,,k] - true_Sigma[,,k], "F")
      })),
      
      # Precision matrix error (Frobenius norm)
      prec_error = mean(sapply(1:G, function(k) {
        norm(result$precision_matrices[,,k] - true_Theta[,,k], "F")
      })),
      
      # Classification error rate
      class_error = mean(result$cluster_assignments != true_labels),
      
      # Number of iterations to convergence
      n_iter = length(result$log_likelihood),
      
      # Final log-likelihood
      final_loglik = tail(result$log_likelihood, 1)
    )
    
    # Print results
    cat(sprintf("\nResults for %s penalties:\n", penalty_name))
    cat("Mean estimation error:", metrics$mean_error, "\n")
    cat("Proportion estimation error:", metrics$prop_error, "\n")
    cat("Covariance matrix error:", metrics$cov_error, "\n")
    cat("Precision matrix error:", metrics$prec_error, "\n")
    cat("Classification error rate:", metrics$class_error, "\n")
    cat("Number of iterations:", metrics$n_iter, "\n")
    cat("Final log-likelihood:", metrics$final_loglik, "\n")
    
    return(metrics)
  }
  
  # Evaluate all results
  evaluation <- lapply(names(results), function(penalty_name) {
    evaluate_results(results[[penalty_name]], penalty_name)
  })
  names(evaluation) <- names(results)
  
  # Return full results
  return(list(
    true_parameters = list(
      proportions = true_pi,
      means = true_mu,
      covariance_matrices = true_Sigma,
      precision_matrices = true_Theta,
      labels = true_labels
    ),
    fitted_models = results,
    evaluation = evaluation,
    data = data
  ))
}

# Visualization function for convergence
plot_convergence <- function(test_results) {
  # Extract log-likelihoods
  log_liks <- lapply(test_results$fitted_models, function(x) x$log_likelihood)
  max_iter <- max(sapply(log_liks, length))
  
  # Create plot
  plot(1:max_iter, type = "n", 
       ylim = range(unlist(log_liks)),
       xlab = "Iteration", 
       ylab = "Penalized Log-likelihood",
       main = "Convergence of EM Algorithm")
  
  colors <- c("weak" = "blue", "medium" = "red", "strong" = "green")
  
  for (penalty_name in names(log_liks)) {
    lines(1:length(log_liks[[penalty_name]]), 
          log_liks[[penalty_name]], 
          col = colors[penalty_name],
          lwd = 2)
  }
  
  legend("bottomright", 
         legend = names(log_liks),
         col = colors,
         lwd = 2,
         title = "Penalty Strength")
}

# Run test
results <- test_clustering_convergence()

# Plot convergence
plot_convergence(results)