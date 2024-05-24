# Additional Functions for the Clustering Procedure of Estimated Liquidity Bonds Time Series

# Data Exploration Functions ---------------------------------------------------

## Missing Values

# Function, which plots missing values in a data frame
plot_na <- function(df) {
  
  df_name <- (deparse(substitute(df)))
  
  plot <- vis_miss(df, warn_large_data = F) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) + 
    ylab("cusip_id") + 
    ggtitle(paste0("Missing Data in ", df_name))
  print(plot)
}

# Function which plots non missing values per column of data frame
plot_non_na_column <- function(df){
  # Calculate the number of missing values in each column
  missing_values <- df %>% 
    summarise_all(function(x) sum(!is.na(x))) %>% 
    gather(key = "Column", value = "Non_Missing_Values")
  
  # Plot the data
  ggplot(missing_values, aes(x = Column, y = Non_Missing_Values)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Date", y = "Number of Unique Bonds Measured", title = "Non Missing Values per Time Point")
  
}

# Function which plots non missing values per row of data frame
plot_non_na_row <- function(df){
  # Calculate the number of missing values in each row
  missing_values <- df %>% 
    mutate(Row = row_number(), Non_Missing_Values = rowSums(!is.na(.))) %>% 
    select(Row, Non_Missing_Values)
  
  ggplot(missing_values, aes(x = Row, y = Non_Missing_Values)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Bond", y = "Number of Observations per Bond", title = "Non Missing Values per Bond")
  
}

# Data Preprocessing -----------------------------------------------------------
# Function, which creates new data frames only with cusip_ids with number of 
# observations over threshold
over_x_observations <- function(list, x){
  
  list_over_x <- list()
  list_name <- deparse(substitute(list))
  
  for (i in seq_along(list)) {
    
    name <- colnames(ebl)[ei_size_restricted_cols][i]
    df <- list[[i]]
    df <- df[rowSums(!is.na(df)) >= x, ]
    list_over_x[[name]] <- df
  }
  
  list_over_x_name <- paste(list_name, "over", x, sep = "_")
  assign(list_over_x_name, list_over_x, envir = .GlobalEnv)
  print(list_over_x_name)
}

# Funtion which plots how many rows have number of observations over the threshold x
plot_over_x_observations_counts <- function(df){
  
  counts <- c()
  for (i in 1:ncol(df)) {
    # Get the data frame
    bonds_over_i <- df[rowSums(!is.na(df)) >= i, ]
    counts <- append(counts, nrow(bonds_over_i))
  }
  plot(counts, xlab = "Observations", ylab = "Number of Bonds with x Observations" )
}


# Function which selects the top y columns with the least missingness for a new df
top_y_rows <- function(list, y){
  list_top_y <- list()
  list_name <- deparse(substitute(list))
  
  for (i in seq_along(list)) {
    # Get the data frame
    df <- list[[i]]
    df_name <- names(list)[i]
    
    # Calculate the number of NA values in each column
    na_counts <- colSums(is.na(df))
    
    # Sort the columns by the number of NA values
    sorted_columns <- sort(na_counts, index.return = TRUE)$ix
    
    # Select the top y columns with the least NA values
    selected_columns <- sorted_columns[1:y]
    
    # Reorder the selected columns to match the original order
    selected_columns <- selected_columns[order(selected_columns)]
    
    # Subset the dataframe to the selected columns
    list_top_y[[df_name]] <- df[, selected_columns]
  }
  
  # Assign the modified list back to the global environment
  list_top_y_name <- paste(list_name, "top", y, sep = "_")
  assign(list_top_y_name, list_top_y, envir = .GlobalEnv)
  
  print(list_top_y_name)
}


# Function, which creates new data frames only with cusip_ids with a maximum of z NAs
remove_rows_with_na <- function(list, z) {
  
  list_z_na <- list()
  list_name <- deparse(substitute(list))
  
  for (i in seq_along(list)) {
    
    name <- colnames(ebl)[ei_size_restricted_cols][i]
    df <- list[[i]]
    na_counts <- rowSums(is.na(df))
    rows_to_remove <- which(na_counts > z)
    df <- df[-rows_to_remove, ]
    list_z_na[[name]] <- df
  }
  
  list_z_na_name <- paste(list_name, z, "NAs", sep = "_")
  assign(list_z_na_name, list_z_na, envir = .GlobalEnv)
  print(list_z_na_name)
  
}

# Clustering -------------------------------------------------------------------

# Hierarchical Clustering
hc_results_ebl <- function(df, raw_df, a, b){
  
  raw_df <- as.data.frame(raw_df)
  # global_y_range <- range(raw_df, na.rm = TRUE)
  df_name <- deparse(substitute(df))
  cut_results <- list()
  series_plots <- list()
  
  for (dist_method in c("euclidean", "DTW")){
    
    dist_matrix <- dist(df, dist_method)
    dist_matrix_name <- paste(df_name, dist_method, sep = "_")
    # assign(dist_matrix_name, dist_matrix, envir = .GlobalEnv)
    
    for(agglomeration_method in c("complete", "ward.D", "ward.D2")){
      clusters <- hierarch_cluster(dist_matrix = dist_matrix, 
                                   compute = F, agglomeration_method = agglomeration_method)
      clusters_name <- paste("hc", agglomeration_method, dist_matrix_name, sep = "_")
      
      # Cut Tree
      for (k in a:b){
        cut_clusters <- cutree(clusters, k = k)
        cut_name <- paste(clusters_name, k, sep = "_")
        cut_results[[cut_name]] <- cut_clusters
        # plots
        
        raw_df$cluster <- factor(cut_clusters)
        # Get the unique cluster ids
        cluster_ids <- unique(raw_df$cluster)
        par(mar = c(5, 4, 4, 2) + 0.1)
        for(i in seq_along(cluster_ids)) {
          # Subset the data frame for the current cluster
          raw_df_subset <- raw_df[raw_df$cluster == cluster_ids[i], ]
          raw_df_subset$cluster <- NULL
          xlab <- paste("Cluster", cluster_ids[i], sep = " ")
          
          
          # Plot all rows of the original data frame
          colors <- hcl.colors(nrow(raw_df_subset), palette = "Set 2", alpha = 1, rev = F) # "Geyser", "Zissou 1"
          matplot(t(raw_df_subset), type = "l", lty = 1, col = colors, main = "", xlab = xlab, ylab = "", cex = 0.6)
          axis(side = 1, at = seq(0, ncol(raw_df_subset), by = 10))
          # legend("topright", inset=c(-0.2,0), legend = rownames(raw_df_subset), col = colors, lty = 1, cex = 0.5, 
          #        xpd = TRUE, horiz = FALSE, x.intersp = 1, text.width = strwidth("17275R"),
          #        ncol = 1)
          
          # Add the plot to the list
          p <- recordPlot()
          series_plots[[paste(cut_name, i, sep = "_")]] <- p
          
        }
        par(mar = c(5, 4, 4, 2) + 0.1)
        raw_df_subset$cluster <- NULL
      }
    }
  }
  
  # save lists in global environment
  cut_results_name <- paste("cut_hc_results", df_name, sep = "_")
  print(cut_results_name)
  assign(cut_results_name, cut_results, envir = .GlobalEnv)
  
  series_plots_name <- paste("series_plots_hc", df_name, sep = "_")
  assign(series_plots_name, series_plots, envir = .GlobalEnv)
  print(series_plots_name)
}

# k-means
kmeans_results <- function (df, raw_df, a, b, initialization){
  df_name <- deparse(substitute(df))
  
  raw_df <- as.data.frame(raw_df)
  df_name <- deparse(substitute(df))
  results <- list()
  plots <- list()
  
  for (k in a:b){
    clusters <- kmeans_cluster(df, k, initialization = initialization, nstart = 50)
    clusters <- clusters$cluster
    clusters_name <- paste("kmeans", k, df_name, sep = "_")
    results[[clusters_name]] <- clusters
    
    raw_df$cluster <- factor(clusters)
    # Get the unique cluster ids
    cluster_ids <- unique(raw_df$cluster)
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(i in seq_along(cluster_ids)) {
      # Subset the data frame for the current cluster
      raw_df_subset <- raw_df[raw_df$cluster == cluster_ids[i], ]
      raw_df_subset$cluster <- NULL
      xlab <- paste("Cluster", cluster_ids[i], sep = " ")
      
      # Plot all rows of the original data frame
      colors <- hcl.colors(nrow(raw_df_subset), palette = "Set 2", alpha = 1, rev = F) # "Geyser", "Zissou 1"
      matplot(t(raw_df_subset), type = "l", lty = 1, col = colors, main = "", xlab = xlab, ylab = "", cex = 0.6)
      axis(side = 1, at = seq(0, ncol(raw_df_subset), by = 10))
      # legend("topright", inset=c(-0.2,0), legend = rownames(raw_df_subset), col = colors, lty = 1, cex = 0.5, 
      #        xpd = TRUE, horiz = FALSE, x.intersp = 1, text.width = strwidth("17275R"),
      #        ncol = 1)
      
      # Add the plot to the list
      p <- recordPlot()
      plots[[paste(clusters_name, i, sep = "_")]] <- p
      
    }
    par(mar = c(5, 4, 4, 2) + 0.1)
    raw_df_subset$cluster <- NULL
  }
  
  
  if(!is.null(initialization)){
    results_name <- paste("kmeans", initialization, "results", df_name, sep = "_")
    plots_name <- paste("kmeans", initialization, "plots", df_name, sep = "_")
  }
  else {
    results_name <- paste("kmeans_results", df_name, sep = "_")
    plots_name <- paste("kmeans_plots", df_name, sep = "_")
  }
  print(results_name)
  assign(results_name, results, envir = .GlobalEnv)
  print(plots_name)
  assign(plots_name, plots, envir = .GlobalEnv)
  
}


# Evaluation -------------------------------------------------------------------

internal_evaluation <- function(dist_matrix, clusters_object, toGlobEnv = T){
  
  clust_stats <- cqcluster.stats(dist_matrix, clusters_object, standardisation = "ave")
  
  int_eval_df <- data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int_eval_df) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
  
  int_eval_df["Calinski-Harabasz Index"] <- clust_stats$ch
  int_eval_df["Dunn Index"] <- clust_stats$dunn
  int_eval_df["Average Silhouette Width"] <- clust_stats$asw
  
  if (toGlobEnv == T){
    clusters_object_name <- deparse(substitute(clusters_object))
    df_name <- paste("int_eval", clusters_object_name, sep = "_")
    assign(df_name, int_eval_df, envir = .GlobalEnv)
    print(df_name)
  }
  
  return(int_eval_df)
  
}

# Function, which extracts a data frame consisting of the best-ranked partitions by each internal evaluation index
top_partitions <- function(ed_evaluation_df, dtw_evaluation_df, top){
  
  df_top <- data.frame()
  for (col in 1:ncol(ed_evaluation_df)){
    sorted_indices <- order(ed_evaluation_df[, col], decreasing = TRUE)
    temp <- ed_evaluation_df[sorted_indices, ]
    df_top <- rbind(df_top, head(temp, top))
  }
  for (col in 1:ncol(dtw_evaluation_df)){
    sorted_indices <- order(dtw_evaluation_df[, col], decreasing = TRUE)
    temp <- dtw_evaluation_df[sorted_indices, ]
    df_top <- rbind(df_top, head(temp, top))
  }
  
  row_names <- rownames(df_top)
  # Find the row names that end with "1" and exist in the data frame without the "1"
  duplicate_rows <- row_names %in% paste0(sub("1$", "", row_names), "1") & sub("1$", "", row_names) %in% row_names
  # Drop these rows
  df_top <- df_top[!duplicate_rows, ]
  
  return(df_top)
  
}

# Function, which saves the evaluation reults as a latex table
latex_table <- function(df, top = NULL, dir = "latex_tables"){
  
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  
  df_name <- deparse(substitute(df))
  
  if (!is.null(top)){
    df_top <- data.frame()
    
    for (col in 1:ncol(df)){
      sorted_indices <- order(df[, col], decreasing = TRUE)
      temp <- df[sorted_indices, ]
      df_top <- rbind(df_top, head(temp, top))
    }
    file_path <- file.path(dir, paste(df_name, top, ".tex", sep = ""))
  }
  else {
    df_top <- df
    file_path <- file.path(dir, paste(df_name, ".tex", sep = ""))
  }
  print(xtable::xtable(df_top), type = "latex", file = file_path)
}

