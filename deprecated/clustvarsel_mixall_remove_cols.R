clvarselgrbkw_mixall_parallel <- function(X, G = 2:9, strategy = clusterStrategy(),
                                          samp = FALSE, sampsize = 2000,
                                          BIC.diff = 0, itermax = 100,
                                          verbose = interactive(),
                                          num_cores = parallel::detectCores() - 1)
{
  require(parallel)
  
  # Convert X to a data frame if it's not already
  X <- as.data.frame(X)
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(colnames(X))) 
    colnames(X) <- paste("X", 1:d, sep = "")
  G <- setdiff(G, 1)
  
  # If needed, sample the subset of observations
  if(samp) { sub <- sample(1:n, min(sampsize,n), replace = FALSE) }
  else     { sub <- seq.int(1,n) }
  
  # Function to determine data type
  get_data_type <- function(col) {
    if (is.numeric(col) && all(col %% 1 == 0, na.rm = TRUE)) {
      return("poisson")
    } else if (is.factor(col) || is.character(col)) {
      return("categorical")
    } else {
      return("gaussian")
    }
  }
  
  # Prepare data and determine clustering function
  prepare_data_and_cluster_function <- function(data) {
    data_types <- sapply(data, get_data_type)
    unique_types <- unique(data_types)
    
    if (length(unique_types) == 1) {
      cluster_fun <- switch(unique_types,
                            "poisson" = clusterPoisson,
                            "categorical" = clusterCategorical,
                            "gaussian" = clusterDiagGaussian)
      models <- lapply(unique_types, function(type) {
        switch(type,
               "poisson" = "poisson_pk_ljk",
               "categorical" = "categorical_pk_pjk",
               "gaussian" = "gaussian_pk_sjk")
      })
      return(list(data = data, models = models[[1]], 
                  cluster_fun = cluster_fun, is_mixed = FALSE))
    } else {
      data_list <- lapply(unique_types, function(type) {
        data[, data_types == type, drop = FALSE]
      })
      models <- lapply(unique_types, function(type) {
        switch(type,
               "poisson" = "poisson_pk_ljk",
               "categorical" = "categorical_pk_pjk",
               "gaussian" = "gaussian_pk_sjk")
      })
      return(list(data = data_list, models = models, 
                  cluster_fun = clusterMixedData, 
                  is_mixed = TRUE))
    }
  }
  
  # Function to fit clustering model
  fit_cluster_model <- function(data, G) {
    tryCatch({
      prepared <- prepare_data_and_cluster_function(data)
      if (prepared$is_mixed) {
        model <- prepared$cluster_fun(prepared$data, nbCluster = G, 
                                      strategy = strategy, criterion = "BIC", 
                                      models = prepared$models, nbCore = 0)
      } else {
        model <- prepared$cluster_fun(prepared$data, nbCluster = G, 
                                      strategy = strategy, criterion = "BIC",
                                      models = prepared$models, nbCore = 0)
      }
      return(list(model = model, BIC = -model@criterion, nbCluster = model@nbCluster))
    }, error = function(e) {
      return(list(model = NULL, BIC = -Inf, nbCluster = NA))
    })
  }
  
  # Initialize cluster for parallel processing
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))
  
  # Load required packages on all worker nodes
  clusterEvalQ(cl, {
    library(MixAll)
    library(BMA)
    library(nnet)
    library(mlogitBMA)
    NULL
  })
  
  # Export necessary functions and objects to the cluster
  clusterExport(cl, c("prepare_data_and_cluster_function", "get_data_type", 
                      "fit_cluster_model", "compute_bic_reg",
                      "strategy"), envir = environment())
  
  # Initialize variables
  S <- X
  NS <- data.frame(matrix(ncol = 0, nrow = nrow(X)))
  info <- data.frame(Var = character(), BIC = numeric(), BICdiff = numeric(), 
                     Step = character(), Decision = character(), 
                     Model = character(), G = integer(), 
                     stringsAsFactors = FALSE)
  
  # Initial clustering
  if (verbose) cat("Initialize model\n")
  init_model <- fit_cluster_model(S, G)
  BICS <- init_model$BIC
  
  criterion <- 1
  iter <- 0
  
  while((criterion == 1) & (iter < itermax) & (ncol(S) > 1)) {
    iter <- iter + 1
    if(verbose) cat(paste("iter", iter, "\n"))
    
    # Removing step - Parallel processing
    if(verbose) cat("- removing step\n")
    remove_results <- parLapply(cl, 1:ncol(S), function(i) {
      S_minus_i <- S[, -i, drop = FALSE]
      model_minus_i <- fit_cluster_model(S_minus_i, G)
      
      # Remove variables (columns) with missing values
      cols_with_missing <- colSums(is.na(S_minus_i)) > 0
      S_minus_i_obs <- S_minus_i[, !cols_with_missing, drop = FALSE]
      
      # Get the removed variable
      S_i_obs <- S[, i, drop = FALSE]
      # Remove observations where S_i_obs is missing
      obs_with_missing_S_i <- is.na(S_i_obs)
      S_minus_i_obs <- S_minus_i_obs[!obs_with_missing_S_i, , drop = FALSE]
      S_i_obs <- S_i_obs[!obs_with_missing_S_i, , drop = FALSE]
      
      BICreg <- compute_bic_reg(S_minus_i_obs, S_i_obs,
                                family = get_data_type(S_i_obs))
      
      list(BIC_total = model_minus_i$BIC + BICreg, 
           BICreg = BICreg,
           model = model_minus_i$model, 
           nbCluster = model_minus_i$nbCluster)
    })
    
    BIC_remove <- sapply(remove_results, function(x) x$BIC_total)
    BIC_reg <- sapply(remove_results, function(x) x$BICreg)
    cdiff_remove <- BICS - BIC_remove
    m_remove <- min(cdiff_remove)
    arg_remove <- which.min(cdiff_remove)
    
    if(m_remove < BIC.diff) {
      # Remove variable
      removed_var <- colnames(S)[arg_remove]
      BICS <- BICS - BIC_reg[arg_remove] - cdiff_remove[arg_remove]
      info <- rbind(info, data.frame(
        Var = removed_var, BIC = BICS, BICdiff = m_remove,
        Step = "Remove", Decision = "Accepted",
        Model = class(remove_results[[arg_remove]]$model)[1],
        G = remove_results[[arg_remove]]$nbCluster,
        stringsAsFactors = FALSE
      ))
      NS[[removed_var]] <- S[[removed_var]]
      S <- S[, -arg_remove, drop = FALSE]
    } else {
      info <- rbind(info, data.frame(
        Var = colnames(S)[arg_remove], BIC = BICS, BICdiff = m_remove,
        Step = "Remove", Decision = "Rejected",
        Model = class(remove_results[[arg_remove]]$model)[1],
        G = remove_results[[arg_remove]]$nbCluster,
        stringsAsFactors = FALSE
      ))
    }
    
    # Adding step - Parallel processing
    if(ncol(NS) > 2) {
      if(verbose) cat("+ adding step\n")
      add_results <- parLapply(cl, 1:ncol(NS), function(i) {
        S_plus_i <- cbind(S, NS[, i, drop = FALSE])
        model_plus_i <- fit_cluster_model(S_plus_i, G)
        
        # Remove variables (columns) with missing values
        cols_with_missing <- colSums(is.na(S)) > 0
        S_obs <- S[, !cols_with_missing, drop = FALSE]
        
        # Get the added variable
        NS_i_obs <- NS[, i, drop = FALSE]
        # Remove observations where NS_i_obs is missing
        obs_with_missing_NS_i <- is.na(NS_i_obs)
        S_obs <- S_obs[!obs_with_missing_NS_i, , drop = FALSE]
        NS_i_obs <- NS_i_obs[!obs_with_missing_NS_i, , drop = FALSE]
        
        BICreg <- compute_bic_reg(S_obs, NS_i_obs, 
                                  family = get_data_type(NS_i_obs))
        
        list(BIC_total = model_plus_i$BIC,
             BICreg = BICreg,
             model = model_plus_i$model,
             nbCluster = model_plus_i$nbCluster)
      })
      
      BIC_add <- sapply(add_results, function(x) x$BIC_total)
      BIC_reg <- sapply(add_results, function(x) x$BICreg)
      cdiff_add <- BIC_add - (BICS + BIC_reg)
      m_add <- max(cdiff_add)
      arg_add <- which.max(cdiff_add)
      
      if(m_add > BIC.diff) {
        # Add variable to S and update clustering BICS
        added_var <- colnames(NS)[arg_add]
        BICS <- BIC_add[arg_add]
        info <- rbind(info, data.frame(
          Var = added_var, BIC = BICS, BICdiff = m_add,
          Step = "Add", Decision = "Accepted",
          Model = class(add_results[[arg_add]]$model)[1],
          G = add_results[[arg_add]]$nbCluster,
          stringsAsFactors = FALSE
        ))
        S[[added_var]] <- NS[[added_var]]
        NS <- NS[, -arg_add, drop = FALSE]
      } else {
        info <- rbind(info, data.frame(
          Var = colnames(NS)[arg_add], BIC = BIC_add[arg_add], BICdiff = m_add,
          Step = "Add", Decision = "Rejected",
          Model = class(add_results[[arg_add]]$model)[1],
          G = add_results[[arg_add]]$nbCluster,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if(verbose) {
      print(info[nrow(info), c("Var", "BICdiff", "Step", "Decision")])
    }
    
    # Check if the variables in S have changed
    criterion <- if(ncol(S) == 0) 0 else 1
  }
  
  if(iter >= itermax) 
    warning("Algorithm stopped because maximum number of iterations was reached")
  
  # Prepare output
  varnames <- colnames(X)
  subset <- if(ncol(S) == 0) NULL else match(colnames(S), varnames)
  
  # Fit final model to get optimal number of clusters
  final_model <- fit_cluster_model(S, G)
  optimal_G <- final_model$nbCluster
  
  out <- list(
    variables = varnames,
    subset = subset,
    cluster_model = final_model,
    steps.info = info,
    optimal_G = optimal_G,
    search = "greedy",
    direction = "backward"
  )
  class(out) <- "clustvarsel_mixall"
  return(out)
}