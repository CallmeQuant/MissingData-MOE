source("clustvarsel_mixall.R")
source("exp_func.R")
set.seed(123)

# Parameters
n <- 2000
p <- c(0.25, 0.25, 0.2, 0.3)
mu <- rbind(c(0, 0, 0),
            c(-6, 6, 0),
            c(0, 0, 6),
            c(-6, 6, 6))

# Create covariance matrix
A <- rotation_matrix_3d("z", pi/6) %*% rotation_matrix_3d("x", pi/3)
Sigma <- A %*% diag(c(6*sqrt(2), 1, 2)) %*% t(A)
diag_vals <- diag(Sigma)
Sigma <- diag(x = diag_vals)

# Generate mixture data
component <- sample(1:4, n, replace = TRUE, prob = p)
X <- matrix(0, nrow = n, ncol = 3)
for (k in 1:4) {
  idx <- which(component == k)
  X[idx,] <- mvrnorm(length(idx), mu = mu[k,], Sigma = Sigma)
}

# Generate fourth and fifth variables
epsilon <- mvrnorm(n, mu = c(0, 0), 
                   Sigma = rotation_matrix_2d(pi/6) %*% diag(c(1, 3)) %*% t(rotation_matrix_2d(pi/6)))
Y45 <- cbind(X[,1:2] %*% matrix(c(0.5, 1, 2, 0), nrow = 2) + epsilon + c(-1, 2))

# Generate two noisy variables
noise <- matrix(rnorm(n*2), ncol = 2)

# Combine all variables
data <- cbind(X, Y45, noise)
colnames(data) <- paste0("V", 1:7)


## Run variable selection
result <- clvarselgrbkw_mixall_parallel(data, G = 2:5, 
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
                                        itermax = 8, verbose = TRUE)

# Print results
print("Variable Selection Results:")
print(result$subset)
print(result$steps.info)

result2 <- clustvarsel(data, G=2:5, search="greedy",
                       direction = "backward", emModels1 = "V", 
                       emModels2 = "VVI", allow.EEE = FALSE, forcetwo = FALSE)

result3 <- VarSelCluster(data, gvals = 2:8,
                         vbleSelec = TRUE, crit.varsel = "BIC", nbcores = 4)

# Evaluate results
true_relevant <- c(1, 2, 3)  # True relevant variables
true_redundant <- c(4, 5)    # True redundant variables
true_noise <- c(6, 7)        # True noise variables

selected_mixall<- result$subset
selected_clustvarsel <- result2$subset
selected_varsellcm <- match(result3@model@names.relevant, colnames(data))

# Selected variables
print("Selected variables")
print("Clustvarsel: ")
print(selected_clustvarsel)

print("Mixall: ")
print(selected_mixall)

print("VarSelLCM: ")
print(selected_varsellcm)

# Correctly identified relevant variables
print("Correctly identified relevant variables")
print("Clustvarsel: ")
print(intersect(selected_clustvarsel, true_relevant))

print("Mixall: ")
print(intersect(selected_mixall, true_relevant))

print("VarSelLCM: ")
print(intersect(selected_varsellcm, true_relevant))

# Correctly identified irrelevant variables
print("Correctly identified irrelevant variables (redundant + noise)")
print("Clustvarsel: ")
print(setdiff(1:7, union(selected_clustvarsel, c(true_redundant, true_noise))))

print("Mixall: ")
print(setdiff(1:7, union(selected_mixall, c(true_redundant, true_noise))))

print("VarSelLCM: ")
print(setdiff(1:7, union(selected_varsellcm, c(true_redundant, true_noise))))

# Misclassify variables
print("Misclassified variables")
print("Clustvarsel: ")
print(setdiff(union(selected_clustvarsel, true_relevant), 
              intersect(selected_clustvarsel, true_relevant)))

print("Mixall: ")
print(setdiff(union(selected_mixall, true_relevant), 
              intersect(selected_mixall, true_relevant)))

print("VarSelLCM: ")
print(setdiff(union(selected_varsellcm, true_relevant),
              intersect(selected_varsellcm, true_relevant)))

