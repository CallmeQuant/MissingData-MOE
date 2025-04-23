simulate_Maugis_data2012 <- function(
    n = 2000, 
    mixture_probs = c(0.25, 0.25, 0.2, 0.3),
    mixture_means = rbind(
      c(0, 0, 0),
      c(-6, 6, 0),
      c(0, 0, 6),
      c(-6, 6, 6)
    ),
    seed = 123,
    sigma_scale = 1,
    noise_vars = 2
) {
  set.seed(seed)
  
  A <- rotation_matrix_3d("z", pi/6) %*% rotation_matrix_3d("x", pi/3)
  
  Sigma <- A %*% diag(c(6*sqrt(2), 1, 2) * sigma_scale) %*% t(A)
  diag_vals <- diag(Sigma)
  Sigma <- diag(x = diag_vals)
  
  component <- sample(1:nrow(mixture_means), n, replace = TRUE, prob = mixture_probs)
  
  X <- matrix(0, nrow = n, ncol = 3)
  
  # Data for each mixture component
  for (k in 1:nrow(mixture_means)) {
    idx <- which(component == k)
    X[idx,] <- mvrnorm(length(idx), mu = mixture_means[k,], Sigma = Sigma)
  }
  
  # Additional correlated variables
  epsilon <- mvrnorm(n, mu = c(0, 0), 
                     Sigma = rotation_matrix_2d(pi/6) %*% 
                       diag(c(1, 3)) %*% 
                       t(rotation_matrix_2d(pi/6)))
  
  Y45 <- cbind(X[,1:2] %*% matrix(c(0.5, 1, 2, 0), nrow = 2) + 
                 epsilon + c(-1, 2))
  
  # Additional noise variables
  noise <- matrix(rnorm(n * noise_vars), ncol = noise_vars)
  
  data <- cbind(X, Y45, noise)
  
  colnames(data) <- paste0("V", 1:(3 + 2 + noise_vars))
  
  return(list(data=data, class=component))
}

simulate_Maugis_data2019 <- function(
    n = 2000, 
    seed = 123,
    scenario = 1,
    all_scenarios = FALSE
) {
  set.seed(seed)
  
  if (all_scenarios) {
    scenario_list <- lapply(1:8, function(s) {
      generate_multi_scenario_data(n = n, seed = seed, scenario = s, all_scenarios = FALSE)
    })
    names(scenario_list) <- paste0("Scenario_", 1:8)
    return(scenario_list)
  }

  # Cluster means for the first two variables
  mu_list <- list(
    c(0, 0),
    c(4, 0),
    c(0, 2),
    c(4, 2)
  )
  
  # Cluster assignments (equiprobable)
  cluster_assignments <- sample(1:4, size = n, replace = TRUE)
  
  # y1_2 based on cluster assignments
  Sigma_cluster <- diag(c(0.5, 0.5))
  
  y1_2 <- matrix(0, nrow = n, ncol = 2)
  for (k in 1:4) {
    idx <- which(cluster_assignments == k)
    nk <- length(idx)
    if (nk > 0) {
      y1_2[idx, ] <- mvrnorm(nk, mu = mu_list[[k]], Sigma = Sigma_cluster)
    }
  }
  
  # tilde_alpha 
  tilde_alpha <- c(0, 0, seq(0.4, 4, by = 0.4))
  
  # Scenario-specific parameter generation
  if (scenario == 1) {
    # Scenario 1
    tilde_beta <- matrix(0, nrow = 2, ncol = 12)
    tilde_Omega_diag <- rep(1, 12)
    
  } else if (scenario == 2) {
    # Scenario 2
    tilde_beta <- cbind(matrix(c(3, 0), nrow = 2), matrix(0, nrow = 2, ncol = 11))
    tilde_Omega_diag <- c(0.5, rep(1, 11))
    
  } else if (scenario == 3) {
    # Scenario 3
    tilde_beta <- cbind(matrix(c(0.5, 1), nrow = 2), matrix(0, nrow = 2, ncol = 11))
    tilde_Omega_diag <- rep(1, 12)
    
  } else if (scenario == 4) {
    # Scenario 4
    beta_1 <- matrix(c(0.5, 1, 2, 0), nrow = 2, ncol = 2)
    tilde_beta <- cbind(beta_1, matrix(0, nrow = 2, ncol = 10))
    tilde_Omega_diag <- rep(1, 12)
    
  } else if (scenario == 5) {
    # Scenario 5
    beta_1 <- matrix(c(0.5, 1, 2, 0), nrow = 2, ncol = 2)
    beta_2 <- matrix(c(0, 3, -1, 2, 2, -4), nrow = 2, ncol = 3)
    tilde_beta <- cbind(beta_1, beta_2, matrix(0, nrow = 2, ncol = 7))
    tilde_Omega_diag <- c(rep(1, 3), rep(0.5, 5), rep(1, 4))
    
  } else if (scenario == 6) {
    # Scenario 6
    beta_1 <- matrix(c(0.5, 1, 2, 0), nrow = 2)
    beta_2 <- matrix(c(0, 3, -1, 2, 2, -4), nrow = 2, ncol = 3)
    beta_3 <- matrix(c(0.5, 0, 4, 0.5, 3, 0, 2, 1), nrow = 2, ncol = 4)
    tilde_beta <- cbind(beta_1, beta_2, beta_3, matrix(0, nrow = 2, ncol = 3))
    tilde_Omega_diag <- c(rep(1, 3), rep(0.5, 2), rep(1, 7))
    
  } else if (scenario == 7) {
    # Scenario 7
    beta_1 <- matrix(c(0.5, 1, 2, 0), nrow = 2)
    beta_2 <- matrix(c(0, 3, -1, 2, 2, -4), nrow = 2, ncol = 3)
    beta_3 <- matrix(c(0.5, 0, 4, 0.5, 3, 0, 2, 1), nrow = 2, ncol = 4)
    additional_betas <- matrix(c(-1, -2, 0, 0.5, 1, 1), nrow = 2, ncol = 3)
    tilde_beta <- cbind(beta_1, beta_2, beta_3, additional_betas)
    tilde_Omega_diag <- c(rep(1, 3), rep(0.5, 2), rep(1, 7))
    
  } else if (scenario == 8) {
    # Scenario 8
    intercept_8 <- c(0, 0, seq(0.4, by = 0.4, length.out = 7))
    b_8 <- matrix(c(
      0.5,  1,
      2,    0,
      0,    3,
      -1,   2,
      2,   -4,
      0.5,  0,
      4,    0.5,
      3,    0,
      2,    1
    ), nrow = 2, byrow = FALSE)
    Omega_8_diag <- c(rep(1,3), rep(0.5,2), 2.5, 1.5, 3.0, 5.0)
    epsilon_8 <- matrix(rnorm(n * 9), nrow = n, ncol = 9) %*% diag(sqrt(Omega_8_diag))
    y3_11 <- matrix(rep(intercept_8, each = n), n, byrow = FALSE) + (y1_2 %*% b_8) + epsilon_8
    
    # Add 3 noise variables 12:14
    y12_14 <- matrix(rnorm(n * 3), n, 3)
    
    data <- cbind(y1_2, y3_11, y12_14)
    colnames(data) <- paste0("V", 1:14)
    return(list(data=data, class=cluster_assignments))
  }
  
  # For scenarios 1 to 7:
  # Generate epsilon with the chosen diagonal covariance
  epsilon <- matrix(rnorm(n * 12), nrow = n, ncol = 12)
  epsilon <- epsilon %*% diag(sqrt(tilde_Omega_diag))
  
  # Compute y_i^{[3:14]}
  y3_14 <- matrix(rep(tilde_alpha, each = n), nrow = n, byrow = FALSE) + y1_2 %*% tilde_beta + epsilon
  
  data <- cbind(y1_2, y3_14)
  colnames(data) <- paste0("V", 1:14)
  
  return(list(data=data, class=cluster_assignments))
}

simulate_sportisse_data <- function(
    n = 100,                  # number of observations
    K = 3,                    # number of clusters
    pik = c(0.5, 0.25, 0.25), # cluster probabilities
    d = 6,                    # number of variables
    tau = 2.31,              # signal strength
    probmiss_y = rep(c(1.45, 0.2, -3), 2),  # missing mechanism parameters for Y
    intercept_y = -1.38,     # intercept for missing mechanism
    missing_type = "MNARy"    # type of missing mechanism
) {
  if (length(pik) != K) stop("Length of pik must equal K")
  if (sum(pik) != 1) stop("pik probabilities must sum to 1")
  if (length(probmiss_y) != d) stop("Length of probmiss_y must equal d")
  
  # Define delta matrix (interaction signals)
  delta <- matrix(0, nrow = K, ncol = d)
  delta[1, 1] <- tau 
  delta[2, 2] <- tau
  delta[3, 3] <- tau
  delta[1, 4] <- tau
  delta[2, 5] <- tau
  delta[3, 6] <- tau
  
  # Missing mechanism parameters
  probmiss_z <- matrix(0, nrow = K, ncol = d, byrow = TRUE)
  
  # Simulate true cluster membership
  Z <- SimuZ(n = n, pik = pik)
  Partition_true <- apply(Z, 1, function(z) which(z == 1))
  
  # Generate data matrix Y
  Y <- matrix(NA, nrow = n, ncol = d)
  for (j in 1:d) {
    Y[, j] <- Z %*% delta[, j] + rnorm(n)
  }
  
  # Introduce missing data
  C <- SimuC(
    pik = pik, 
    Y = Y, 
    Z = Z, 
    mecha = missing_type,
    probmiss_z = probmiss_z, 
    probmiss_y = probmiss_y, 
    intercept_y = intercept_y
  )
  YNA <- Y
  YNA[C] <- NA
  
  return(list(
    Z = Z,                    # True cluster membership
    Partition_true = Partition_true,  # Cluster assignments
    Y = Y,                    # Complete data
    YNA = YNA,                # Data with missing values
    delta = delta,            # Signal matrix
    C = C                     # Missing data indicator
  ))
}

simulate_ESM_data <- function(
    n_samples = 450,
    relevant_means = list(c(-1, 2), c(2, -1)),
    relevant_sigma = list(matrix(c(1, 0.1, 0.1, 1), 2, 2), 
                          matrix(c(1, 0.1, 0.1, 1), 2, 2)),
    noise_features = list(
      f3 = list(mean = 1.5, sd = 1),
      f4 = list(mean = 3, sd = 0.5),
      f5 = list(mean = 1.8, sd = 0.9),
      f6 = list(mean = 0.3, sd = 0.5),
      f7 = list(mean = 2, sd = 1),
      f8 = list(mean = -2, sd = 2),
      f9 = list(mean = 4, sd = 3),
      f10 = list(mean = 1.5, sd = 0.1),
      f11 = list(mean = -4, sd = 2)
    ),
    correlated_pairs = list(
      f12_f13 = list(mean = c(2, 2), sigma = matrix(c(1, 0.3, 0.3, 1), 2, 2)),
      f14_f15 = list(mean = c(-3, 2), sigma = matrix(c(1, 0.1, 0.1, 1), 2, 2))
    ),
    add_correlations = TRUE,
    diagonal = FALSE
) {
  # Generate cluster assignments
  labels <- sample(1:2, n_samples, replace = TRUE)
  n1 <- sum(labels == 1)
  n2 <- sum(labels == 2)
  
  # Generate relevant features (f1, f2) from mixture of Gaussians
  relevant_features <- matrix(0, nrow = n_samples, ncol = 2)
  
  relevant_features[labels == 1, ] <- MASS::mvrnorm(
    n = n1,
    mu = relevant_means[[1]],
    Sigma = relevant_sigma[[1]]
  )
  
  relevant_features[labels == 2, ] <- MASS::mvrnorm(
    n = n2,
    mu = relevant_means[[2]],
    Sigma = relevant_sigma[[2]]
  )
  
  # Generate irrelevant features (f3-f11)
  irrelevant_features <- matrix(0, nrow = n_samples, ncol = length(noise_features))
  
  for (i in 1:length(noise_features)) {
    feature <- noise_features[[i]]
    irrelevant_features[, i] <- rnorm(n_samples, mean = feature$mean, sd = feature$sd)
  }
  
  # Generate correlated feature pairs (f12-f13, f14-f15)
  correlated_features <- matrix(0, nrow = n_samples, ncol = 4)
  
  if (add_correlations) {
    # Generate f12 and f13 (correlated)
    correlated_features[, 1:2] <- MASS::mvrnorm(
      n = n_samples,
      mu = correlated_pairs$f12_f13$mean,
      Sigma = correlated_pairs$f12_f13$sigma
    )
    
    # Generate f14 and f15 (correlated)
    correlated_features[, 3:4] <- MASS::mvrnorm(
      n = n_samples,
      mu = correlated_pairs$f14_f15$mean,
      Sigma = correlated_pairs$f14_f15$sigma
    )
  } else {
    # Generate uncorrelated features if add_correlations is FALSE
    correlated_features[, 1] <- rnorm(n_samples, mean = correlated_pairs$f12_f13$mean[1], sd = 1)
    correlated_features[, 2] <- rnorm(n_samples, mean = correlated_pairs$f12_f13$mean[2], sd = 1)
    correlated_features[, 3] <- rnorm(n_samples, mean = correlated_pairs$f14_f15$mean[1], sd = 1)
    correlated_features[, 4] <- rnorm(n_samples, mean = correlated_pairs$f14_f15$mean[2], sd = 1)
  }
  
  # Combine all features
  data <- cbind(relevant_features, irrelevant_features, correlated_features)
  colnames(data) <- c(paste0("f", 1:2), paste0("f", 3:15))
  
  # Calculate correlation matrix
  cor_matrix <- cor(data)
  
  # If diagonal=TRUE, decorrelate the features
  if (diagonal) {
    # Perform eigendecomposition of correlation matrix
    eigen_decomp <- eigen(cov(data))
    
    # Calculate the whitening transformation
    whitening_matrix <- eigen_decomp$vectors %*% diag(1/sqrt(eigen_decomp$values)) %*% t(eigen_decomp$vectors)
    
    # Apply the whitening transformation
    data_whitened <- scale(data, center = TRUE, scale = FALSE) %*% whitening_matrix
    
    # Restore original means and variances
    means_original <- colMeans(data)
    sds_original <- apply(data, 2, sd)
    
    # Scale back to original range
    data_decorrelated <- scale(data_whitened, center = FALSE, scale = 1/apply(data_whitened, 2, sd))
    data_decorrelated <- scale(data_decorrelated, center = -means_original/sds_original, scale = 1/sds_original)
    
    data <- data_decorrelated
    cor_matrix_after <- cor(data)
    
    cor_matrix <- list(original = cor_matrix, decorrelated = cor_matrix_after)
  }
  
  # Return results
  return(list(
    data = data,
    labels = labels,
    relevant_idx = 1:2,
    feature_info = list(
      relevant_means = relevant_means,
      relevant_sigma = relevant_sigma,
      noise_features = noise_features,
      correlated_pairs = correlated_pairs
    ),
    cor_matrix = cor_matrix
  ))
}
