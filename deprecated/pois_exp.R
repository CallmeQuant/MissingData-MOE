source("clustvarsel_mixall.R")
source("exp_func.R")
source("amputation.R")
set.seed(123)

generate_count_data <- function(n = 1000, d = 10, k = 3, missing = FALSE) {
  relevant_vars <- paste0("X", 1:5)
  irrelevant_vars <- paste0("X", 6:10)
  
  data <- data.frame(matrix(0, nrow = n, ncol = d))
  names(data) <- c(relevant_vars, irrelevant_vars)
  
  cluster_sizes <- rep(n / k, k)
  start_idx <- 1
  for (i in 1:k) {
    end_idx <- start_idx + cluster_sizes[i] - 1
    # Relevant variables have cluster-specific means
    data[start_idx:end_idx, relevant_vars] <- matrix(rpois(length(relevant_vars) * cluster_sizes[i], lambda = 3 + i), ncol = length(relevant_vars))
    # Irrelevant variables have the same distribution across clusters
    data[start_idx:end_idx, irrelevant_vars] <- matrix(rpois(length(irrelevant_vars) * cluster_sizes[i], lambda = 3), ncol = length(irrelevant_vars))
    start_idx <- end_idx + 1
  }
  
  if (missing) {
    data <- produce_NA(data)$data.incomp
  }
  return(data)
}

# Intialize EM
EM_init <- clusterStrategy(nbTry = 10,
                           nbInit = 25,
                           initMethod = "class",
                           initAlgo = "SEM",
                           nbInitIteration = 50,
                           initEpsilon = 1e-3,
                           nbShortRun = 10,
                           shortRunAlgo = "EM",
                           nbShortIteration = 100,
                           shortEpsilon = 1e-05,
                           longRunAlgo = "EM",
                           nbLongIteration = 1000,
                           longEpsilon = 1e-05)

# Without missing 
count_data_wo_missing <- generate_count_data(n = 300)

count_result_wo_missing <- clvarselgrbkw_mixall_parallel(count_data_wo_missing, G = 3, 
                                                strategy = EM_init, 
                                                itermax = 8, verbose = TRUE)

print("Count Data w/o Missing Results:")
print(count_result_wo_missing$subset)
print(count_result_wo_missing$steps.info)

# With missing 
count_data_missing <- generate_count_data(n = 300, missing=TRUE)

count_result_missing <- clvarselgrbkw_mixall(count_data_missing, G = 3, 
                                             strategy = EM_init, 
                                             itermax = 8 , verbose = TRUE)

print("Count Data Missing Results:")
print(count_result_missing$subset)
print(count_result_missing$steps.info)