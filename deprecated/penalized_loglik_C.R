library(rBeta2009)
library(glasso)
library(mvtnorm)

# K-means initialization function remains unchanged
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

# Adjusted EM algorithm with penalized complete log-likelihood and K-means initialization
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
  
  log_likelihood <- numeric(max_iter)

  for (iter in 1:max_iter) {
    ## E-step: Calculate conditional probabilities (posterior probabilities)
    log_resp <- matrix(0, n, G)
    for (k in 1:G) {
      log_resp[, k] <- log(pi_k[k]) + dmvnorm(data, mean = mu_k[k, ], 
                                              sigma = (Sigma_k[, , k] + t(Sigma_k[, , k])) / 2, 
                                              log = TRUE)
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
      W <- Theta_k[, , k]          # Precision matrix for cluster k
      MuPrec <- mu_k[k, ]          # Previous mean vector
      MuNew <- numeric(p)          # New mean vector
      
      # Update means with penalization
      for (j in 1:p) {
        Tabs <- sum(t_ik[, k] * (rowSums(sweep(data, 2, MuPrec, FUN = "-") * W[, j]) + MuPrec[j] * W[j, j]))
        
        if (abs(Tabs) <= lambda) {
          MuNew[j] <- 0  # Threshold condition met, set mean to 0
        } else {
          # Update using penalization with precision matrix
          T <- sum(t_ik[, k] * data[, j] * W[j, j]) - 
            sum(t_ik[, k] * rowSums(sweep(data, 2, MuPrec, FUN = "-") * W[j, -j])) + lambda * sign(MuPrec[j])
          MuNew[j] <- T / (n_k[k] * W[j, j])
        }
      }
      mu_k[k, ] <- MuNew  # Update mean for cluster k
      
      # Calculate empirical covariance matrix S_k with updated means
      centered_data <- sweep(data, 2, MuNew, FUN = "-")
      S_k <- t(centered_data) %*% (centered_data * t_ik[, k]) / n_k[k]
      
      # Update precision and covariance matrices using graphical lasso
      rhotilde <- (2 * rho) / n_k[k]
      glasso_fit <- glasso::glasso(S_k, rho = rhotilde, 
                                   penalize.diagonal = FALSE, 
                                   thr = 0.001, maxit = 1000)
      Theta_k[, , k] <- glasso_fit$wi      # Precision matrix
      Sigma_k[, , k] <- solve(glasso_fit$wi)  # Covariance matrix
      
      # Enforce symmetry in Sigma_k to avoid numerical asymmetry issues
      Sigma_k[, , k] <- (Sigma_k[, , k] + t(Sigma_k[, , k])) / 2
    }
    
    ## Compute penalized log-likelihood for convergence check
    penalized_loglik <- sum(log(rowSums(resp))) - 
      (lambda * sum(abs(mu_k))) - 
      (rho * sum(sapply(1:G, function(k) sum(abs(Theta_k[, , k])))))
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
