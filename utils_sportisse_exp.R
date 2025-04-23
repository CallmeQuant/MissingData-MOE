read_experiment_files <- function(base_dir, dimensions = c(6, 9), 
                                  missing_patterns = c("MNARz", "MNARy", "MNARyz"), 
                                  file_types = c("detailed", "summary")) {
  
  detailed_results <- list()
  summary_results <- list()
  
  for (d in dimensions) {
    for (pattern in missing_patterns) {
      
      # Construct file paths
      detailed_path <- file.path(base_dir, "mnar", 
                                 paste0("exp_", d, "_20_", pattern, 
                                        "_full_sportisse_detailed.csv"))
      summary_path <- file.path(base_dir, "mnar", 
                                paste0("exp_", d, "_20_", pattern, 
                                       "_full_sportisse_summary.csv"))
      
      if (file.exists(detailed_path)) {
        detailed_df <- read.csv(detailed_path)
        detailed_df$dimension <- d
        detailed_df$missing_pattern <- pattern
        list_name <- paste0("d", d, "_", pattern, "_detailed")
        detailed_results[[list_name]] <- detailed_df
      } else {
        warning(paste("File not found:", detailed_path))
      }
      
      if (file.exists(summary_path)) {
        summary_df <- read.csv(summary_path)
        summary_df$dimension <- d
        summary_df$missing_pattern <- pattern
        list_name <- paste0("d", d, "_", pattern, "_summary")
        summary_results[[list_name]] <- summary_df
      } else {
        warning(paste("File not found:", summary_path))
      }
    }
  }
  
  return(list(detailed = detailed_results, summary = summary_results))
}

prepare_plotting_data <- function(data, metric_type = "ARI") {
  metric_cols <- grep(paste0("^", metric_type, "_"), names(data), value = TRUE)
  
  plot_data <- data %>%
    select(Iteration, dimension, missing_pattern, all_of(metric_cols)) %>%
    pivot_longer(
      cols = all_of(metric_cols),
      names_to = "Mechanism",
      values_to = metric_type
    ) %>%
    mutate(Mechanism = str_replace(Mechanism, paste0(metric_type, "_"), ""))
  
  if (metric_type == "ARI" && !"True" %in% unique(plot_data$Mechanism)) {
    unique_combos <- plot_data %>%
      select(dimension, missing_pattern) %>%
      distinct()
    
    true_data <- unique_combos %>%
      crossing(Iteration = unique(plot_data$Iteration)) %>%
      mutate(Mechanism = "True", ARI = 0.9)  
    
    plot_data <- bind_rows(plot_data, true_data)
  }
  
  return(plot_data)
}

create_metric_plot <- function(plot_data, metric = "ARI", 
                               y_limit = NULL, 
                               h_line = NULL, 
                               h_line_color = "red", 
                               plot_title = NULL) {
  plot_data_filtered <- plot_data %>%
    filter(Mechanism != "True")
  
  plot_data_filtered$dimension <- factor(plot_data_filtered$dimension)
  plot_data_filtered$missing_pattern <- factor(plot_data_filtered$missing_pattern, 
                                               levels = c("MNARy", "MNARz", "MNARyz"))
  
  palette <- viridisLite::viridis(length(unique(plot_data_filtered$Mechanism)))
  
  if (is.null(plot_title)) {
    plot_title <- paste(metric, "by Missing Pattern and Dimension")
  }
  
  p <- ggplot(plot_data_filtered, aes(x = Mechanism, y = .data[[metric]], fill = Mechanism)) +
    geom_boxplot(outlier.size = 1, alpha = 0.8) +
    facet_grid(dimension ~ missing_pattern) +
    scale_fill_manual(values = palette) +
    theme_bw() + 
    labs(title = plot_title,
         y = metric,
         x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold", color = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  if (!is.null(h_line)) {
    p <- p + geom_hline(yintercept = h_line, linetype = "dashed", 
                        color = h_line_color, linewidth = 0.8)
  }
  
  if (!is.null(y_limit)) {
    p <- p + scale_y_continuous(limits = y_limit)
  }
  
  return(p)
}