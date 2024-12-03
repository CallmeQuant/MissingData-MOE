source("clustvarsel_mixall.R")
source("exp_func.R")
source("amputation.R")
# run analysis for given missing rate
run_analysis <- function(data, missing_rate) {
  set.seed(123 + round(missing_rate * 100))
  
  data_missing_obj <- produce_NA(data, mechanism = "MAR",
                                 perc.missing = missing_rate)
  data_missing <- data_missing$data.incomp
  # Run variable selection methods
  result_mixall <- clvarselgrbkw_mixall_parallel(data_missing, G = 2:5, 
                                                 strategy = clusterStrategy(
                                                   nbTry = 5,
                                                   nbInit = 10,
                                                   initMethod = "class",
                                                   initAlgo = "SEM",
                                                   nbInitIteration = 25,
                                                   initEpsilon = 1e-3,
                                                   nbShortRun = 10,
                                                   shortRunAlgo = "EM",
                                                   nbShortIteration = 100,
                                                   shortEpsilon = 1e-05,
                                                   longRunAlgo = "EM",
                                                   nbLongIteration = 1000,
                                                   longEpsilon = 1e-05),
                                                 itermax = 10, 
                                                 verbose = TRUE)
  
  # Extract results
  selected_mixall <- if(inherits(result_mixall, "try-error")) NULL else result_mixall$subset
  optimal_G <- result_mixall$optimal_G
  # Return results
  return(list(
    missing_rate = missing_rate,
    mixall = selected_mixall,
    G = optimal_G
  ))
}


# Run analysis for different missing rates
# missing_rates <- c(0, 0.05, 0.10, 0.15, 0.20)
missing_rates <- c(0.05)
results <- lapply(missing_rates, function(rate) run_analysis(data, rate))

# Function to evaluate results
evaluate_results <- function(selected_vars, 
                             true_relevant = c(1, 2, 3),
                             true_redundant = c(4, 5), 
                             true_noise = c(6, 7)) {
  if(is.null(selected_vars)) return(list(
    correctly_identified = numeric(0),
    misclassified = numeric(0)
  ))
  correctly_identified <- intersect(selected_vars, true_relevant)
  misclassified <- setdiff(union(selected_vars, true_relevant),
                           intersect(selected_vars, true_relevant))
  
  return(list(
    correctly_identified = correctly_identified,
    misclassified = misclassified
  ))
}

# Print results for each missing rate
for(i in seq_along(results)) {
  cat(sprintf("\nResults for %d%% missing values:\n", results[[i]]$missing_rate * 100))
  
  # Evaluate each method
  methods <- c("mixall")
  for(method in methods) {
    cat(sprintf("\n%s:\n", method))
    eval_result <- evaluate_results(results[[i]][[method]])
    cat("Correctly identified num clusters:", 
        paste(results[[i]][['G']], collapse = ", "), "\n")
    cat("Correctly identified relevant variables:", 
        paste(eval_result$correctly_identified, collapse = ", "), "\n")
    cat("Misclassified variables:", 
        paste(eval_result$misclassified, collapse = ", "), "\n")
  }
}
