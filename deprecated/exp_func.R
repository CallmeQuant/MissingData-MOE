# Function to create 2D rotation matrix
rotation_matrix_2d <- function(angle) {
  matrix(c(cos(angle), -sin(angle),
           sin(angle), cos(angle)), nrow = 2)
}

# Function to create 3D rotation matrix
rotation_matrix_3d <- function(axis, angle) {
  if (axis == "z") {
    matrix(c(cos(angle), -sin(angle), 0,
             sin(angle), cos(angle), 0,
             0, 0, 1), nrow = 3)
  } else if (axis == "x") {
    matrix(c(1, 0, 0,
             0, cos(angle), -sin(angle),
             0, sin(angle), cos(angle)), nrow = 3)
  }
}

compute_error_metrics <- function(actual, predicted) {
  # Mean Squared Error
  mse <- mean((actual - predicted)^2)
  
  # Root Mean Square Error
  rmse <- sqrt(mse)
  
  # Normalized Root Mean Square Error
  nrmse <- rmse / sqrt(mean(actual^2))
  
  return(list(
    mse = mse,
    rmse = rmse,
    nrmse = nrmse
  ))
}

create_imputation_comparison <- function(data, imputed_values) {
  if (!is.data.frame(imputed_values)) {
    imputed_values <- as.data.frame(imputed_values)
  }
  
  # Initialize results dataframe
  comparison_df <- data.frame(
    row = imputed_values$row,
    col = imputed_values$col,
    actual_value = NA,
    imputed_value = imputed_values$value,
    abs_difference = NA
  )
  
  for (i in 1:nrow(comparison_df)) {
    row_idx <- comparison_df$row[i]
    col_idx <- comparison_df$col[i]
    actual_val <- data[row_idx, col_idx]
    comparison_df$actual_value[i] <- actual_val
    comparison_df$abs_difference[i] <- abs(actual_val - comparison_df$imputed_value[i])
  }
  
  error_metrics <- compute_error_metrics(comparison_df$actual_value, 
                                         comparison_df$imputed_value)
  
  attr(comparison_df, "mse") <- error_metrics$mse
  attr(comparison_df, "rmse") <- error_metrics$rmse
  attr(comparison_df, "nrmse") <- error_metrics$nrmse
  attr(comparison_df, "mean_abs_diff") <- mean(comparison_df$abs_difference)
  attr(comparison_df, "max_abs_diff") <- max(comparison_df$abs_difference)
  
  return(comparison_df)
}
