library(MixAll)
library(MASS)
library(mclust)
library(clustvarsel)
library(stats)
library(VarSelLCM)
library(BMA)
library(mlogitBMA)
library(mlogit)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)
library(reshape2)
library(ggplot2)


compute_bic_reg <- function(x, y, family = "gaussian", full_res = FALSE) {
  # Ensure required packages are loaded
  required_packages <- c("BMA", "nnet", "mlogitBMA")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, 
                 "is required but not installed. Please install it."))
    }
  }
  
  x <- as.matrix(x)
  
  if (family != "multinomial") {
    y <- as.vector(unlist(y)) 
  } else {
    y <- as.factor(unlist(y)) 
  }
  
  # Check length
  if (nrow(x) != length(y)) {
    stop(paste("Mismatch in lengths: x has", nrow(x), 
               "rows but y has", length(y), "elements"))
  }
  
  n <- length(y)
  
  # Initialize variables
  selected_vars <- NULL
  bic_value <- NA
  fit <- NULL
  
  if (family != "multinomial") {
    bic_model <- bic.glm(x, y, glm.family = family, 
                         strict = FALSE, nbest = 1)
    
    # Select the best model
    best_model <- bic_model$which[1, ]
    selected_vars <- which(best_model)
    
    # Fit on selected vars
    if (family == "gaussian") {
      fit <- lm(y ~ ., data = data.frame(y = y, 
                                         x = x[, selected_vars, drop = FALSE]))
    } else if (family == "poisson") {
      fit <- glm(y ~ ., data = data.frame(y = y, 
                                          x = x[, selected_vars, drop = FALSE]),
                 family = poisson())
    }
    
    # Compute BIC
    bic_value <- -BIC(fit)
    
  } else {
    data_df <- as.data.frame(cbind(y, x))
    names(data_df)[1] <- "y"
    
    f_mnl <- as.formula(paste("y ~", paste(colnames(x), collapse = " + ")))
    
    # Perform BMA for Multinomial Logistic Regression
    bic_model <- tryCatch({
      bic.mlogit(f = f_mnl, data = data_df, 
                 base.choice = 1, varying = NULL, sep = ".", 
                 approx = TRUE, include.intercepts = TRUE, nbest = 1)
    }, error = function(e) {
      stop("Error in bic.mlogit: ", e$message)
    })
    
    # Extract the selected variables
    selected_models <- bic_model$bic.glm$which
    selected_vars <- which(selected_models)
    
    f_sel <- as.formula(paste("y ~", 
                              paste(colnames(x)[selected_vars], 
                                    collapse = " + ")))
    
    fit <- tryCatch({
      nnet::multinom(f_sel, data = data_df, trace = FALSE)
    }, error = function(e) {
      stop("Error in fitting multinomial model: ", e$message)
    })
    
    # Compute BIC for the multinomial model
    ll <- logLik(fit)
    # Number of parameters: (num preds + intercept) * (num classes - 1)
    p <- (length(selected_vars) + 1) * (length(levels(y)) - 1)
    bic_value <- -(-2 * as.numeric(ll) + log(n) * p)
  }
  
  if (full_res){
    return(list(BIC = bic_value, 
                selected_variables = selected_vars, 
                model = fit))
  }
  else {
    return(BIC = bic_value)
  }
}

clvarselgrbkw_mixall_parallel <- function(X, 
                                          G = 2:9, 
                                          strategy = clusterStrategy(),
                                          samp = FALSE, 
                                          sampsize = 2000,
                                          BIC.diff = 0, 
                                          itermax = 100,
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
      
      # Calculate BIC for regression of removed variable on remaining variables
      fully_observed <- complete.cases(S_minus_i)
      S_minus_i_obs <- S_minus_i[fully_observed, , drop = FALSE]
      S_i_obs <- S[fully_observed, i, drop = FALSE]
      BICreg <- compute_bic_reg(S_minus_i_obs, S_i_obs,
                                family = get_data_type(S[, i]))
      
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
    if(ncol(NS) > 0) {
      if(verbose) cat("+ adding step\n")
      add_results <- parLapply(cl, 1:ncol(NS), function(i) {
        S_plus_i <- cbind(S, NS[, i, drop = FALSE])
        model_plus_i <- fit_cluster_model(S_plus_i, G)
        
        # Calculate BIC for regression of added variable on current variables
        fully_observed <- complete.cases(S)
        S_obs <- S[fully_observed, , drop = FALSE]
        NS_i_obs <- NS[fully_observed, i, drop = FALSE]
        BICreg <- compute_bic_reg(S_obs, NS_i_obs, 
                                  family = get_data_type(NS[, i]))
        
        list(BIC_total = model_plus_i$BIC,
             BICreg = BICreg,
             model = model_plus_i$model,
             nbCluster = model_plus_i$nbCluster)
      })
      
      BIC_add <- sapply(add_results, function(x) x$BIC_total)
      BIC_reg <- sapply(add_results, function(x) x$BICreg)
      cdiff_add <- BIC_add - (BICS + BIC_reg)
      m_add <- min(cdiff_add)
      arg_add <- which.min(cdiff_add)
      
      if(m_add > BIC.diff) {
        # Add variable
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
    steps.info = info,
    optimal_G = optimal_G,
    search = "greedy",
    direction = "backward"
  )
  class(out) <- "clustvarsel"
  return(out)
}


