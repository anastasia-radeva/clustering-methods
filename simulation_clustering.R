# The following estimates clusters within the simulated data frame. 
# k-means and the agglomerative hierarchical clustering algorithms are initialized
# with several initialization functions and agglomeration criteria, respectively.
# The results are evaluated using internal evaluation criteria and compared using
# external evaluation criteria. The performance of the the algorithms is assessed
# based on the evaluation results.

# Dependencies -----------------------------------------------------------------
source("SARFIMA-GARCH_Simulation.R")
source("feature_extractions.R")
source("clustering_algorithms.R")

library("tidyverse")
library("parallel")
library("xtable")
library("graphics")
library("gridExtra")
library("reshape2")
library("ggplot2")
library("fpc")
library("ClusterR")
library("cluster")
library("ppclust")
library("dtwclust")
library("fclust")

# Extract Features -------------------------------------------------------------
extract_features(df_simulated, T) # returns df_simulated_features
df_plot(df_simulated_features)
latex_table(df_simulated_features)
rescale_features(df_simulated_features, tsa = T) # returns scaled_df_simulated_features
df_plot_features(scaled_df_simulated_features)
latex_table(scaled_df_simulated_features)

# Hierarchical Clustering Functions --------------------------------------------
# Function, which performs hierarchical clustering with different distances, agglomeration methods, tree depths
# Returns a list of the results
hc_results <- function(df, a, b){
    
  df_name <- deparse(substitute(df))
  
  results <- list()
  cut_results <- list()
  plots = list()
  cut_plots = list()
  
  for (dist_method in c("euclidean", "DTW")){
    
    dist_matrix <- dist(df, dist_method)
    dist_matrix_name <- paste("dist", df_name, dist_method, sep = "_")
    # assign(dist_matrix_name, dist_matrix, envir = .GlobalEnv)
    
    for(agglomeration_method in c("complete", "average", "ward.D", "ward.D2")){
      clusters <- hierarch_cluster(dist_matrix = dist_matrix, 
                                   compute = F, agglomeration_method = agglomeration_method)
      clusters_name <- paste("clusters_hc", agglomeration_method, dist_matrix_name, sep = "_")
      # Store the clusters in the results list
      results[[clusters_name]] <- clusters
      
      # Tree Plot
      par(mar = c(5, 4, 4, 2) + 0.1) 
      plot(clusters, main = "", 
           ylab = "Distance", xlab = "", cex.main = 0.8)
      tree <- recordPlot()
      plots[[clusters_name]] <- tree
      
      # Cut Tree
      for (k in a:b){
        cut_clusters <- cutree(clusters, k = k)
        cut_name <- paste(clusters_name, k, sep = "_")
        cut_results[[cut_name]] <- cut_clusters
        
        # Plots
        par(mar = c(120/(sqrt(2)*k), 4, 2, 2) + 0.1)
        # empty plot
        
        plot(1, 1, main = "", 
             xlim = c(1 - 0.5, length(cut_clusters) + 0.5), ylim = range(cut_clusters), 
             type = "n", xlab = "", ylab = "Cluster ID", xaxt = "n", cex.main = 0.8, ygap.axis = 0.3)
        # Add faint y-axis grid lines
        clip(1 - 2, length(cut_clusters) + 2, min(cut_clusters) - 1, max(cut_clusters) + 1)
        abline(h = seq(floor(min(cut_clusters)), ceiling(max(cut_clusters)), 
                       by = 1), col = "lightgray", lty = "longdash", xpd = F)
        # Add points
        points(cut_clusters, cex = 1) # col = colors, pch = 19,
        axis(1, at = seq(1,length(cut_clusters), by = 1), labels = rownames(df)[seq(1, nrow(df), by = 1)], las = 2, cex.axis = 0.7)
        axis(2, at = seq(floor(min(cut_clusters)), ceiling(max(cut_clusters)), by = 1))
        
        cut_tree <- recordPlot()
        cut_plots[[cut_name]] <- cut_tree
        
      }
    }
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  # save lists in global environment
  results_name <- paste("hc_results", df_name, sep = "_")
  print(results_name)
  assign(results_name, results, envir = .GlobalEnv)
  
  plots_name <- paste("tree_plots", df_name, sep = "_")
  print(plots_name)
  assign(plots_name, plots, envir = .GlobalEnv)
  
  cut_results_name <- paste("cut_hc_results", df_name, sep = "_")
  print(cut_results_name)
  assign(cut_results_name, cut_results, envir = .GlobalEnv)
  
  cut_plots_name <- paste("cut_tree_plots", df_name, sep = "_")
  print(cut_plots_name)
  assign(cut_plots_name, cut_plots, envir = .GlobalEnv)
  
}
# Function, which plots the raw data in clusters
plot_assigned_raw_ts <- function(assignments_list, raw_df){
  
  assignments_list_name <- deparse(substitute(assignments_list))
  assignments_raw <- list()
  
  for (name in names(assignments_list)){
    
    raw_df$cluster <- factor(unlist(assignments_list[[name]]))
    cluster_ids <- unique(raw_df$cluster)
    clusters_name <- name
    
    # For each cluster id
    for(i in seq_along(cluster_ids)) {
      # Subset the data frame for the current cluster
      raw_df_subset <- raw_df[raw_df$cluster == cluster_ids[i], ]
      raw_df_subset$cluster <- NULL
      
      xlab <- paste("Cluster", cluster_ids[i], sep = " ")
      # Plot all rows of the original data frame
      colors <- hcl.colors(nrow(raw_df_subset), palette = "Set 2", alpha = 1, rev = F) # "Geyser", "Zissou 1"
      matplot(t(raw_df_subset), type = "l", lty = 1, col = colors, main = "", xlab = xlab, ylab = "", ylim = c(-10, 10), cex = 0.8)
      legend("topright" , legend = rownames(raw_df_subset), col = colors, lty = 1, cex = 0.65)
      # Add the plot to the list
      plot <- recordPlot()
      plot_name <- paste(name, i, sep = "_")
      assignments_raw[[plot_name]] <- plot
    }
    
  }
  plots_name <- paste("raw_ts_plots", assignments_list_name, sep = "_")
  print(plots_name)
  assign(plots_name, assignments_raw, envir = .GlobalEnv)
  
}
# Function, which saves plots to specified directory
save_plots <- function(plots, dir = ".") {
  # Create the directory if it doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  
  # Iterate over each plot in the list
  for (i in seq_along(plots)) {
    # Open a PDF file
    pdf(file.path(dir, paste(names(plots)[i], ".pdf", sep = "")))
    # Replay the plot
    replayPlot(plots[[i]])
    # Close the PDF file
    dev.off()
  }
}
# Function, which evaluates hc clustering
evaluate_hc_clustering <- function(df, list){
  
  list_name <- deparse(substitute(list))
  
  euclidean <- list[grepl("euclidean", names(list))]
  dtw <- list[grepl("DTW", names(list))]
  
  internal_evaluation_df <- data.frame(matrix(nrow = length(list), ncol = 3))
  rownames(internal_evaluation_df) <- c(names(list))
  colnames(internal_evaluation_df) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
  
  for (i in 1:length(euclidean)){
    clusters_name <- names(euclidean)[i]
    internal_evaluation_df[clusters_name, ] <- unname(internal_evaluation(df, "euclidean", euclidean[[i]], toGlobEnv = F))
  }
  
  for (j in 1:length(dtw)){
    clusters_name <- names(dtw)[j]
    internal_evaluation_df[clusters_name, ] <- unname(internal_evaluation(df, "DTW", dtw[[j]], toGlobEnv = F))
  }
  
  int_df_name <- paste("internal_evaluation", list_name, sep = "_")
  print(int_df_name)
  assign(int_df_name, internal_evaluation_df, envir = .GlobalEnv)
  
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
# Function, which visualizes the results of the clustering evaluation
create_evaluation_plots <- function(df, top = nrow(df), dir = "evaluation_plots") {
  
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  
  for(i in 1:ncol(df)) {
    df_col <- data.frame(Column = rownames(df), Value = df[, i])
    df_col$Group <- colnames(df)[i]
    # top X
    df_col <- df_col[order(df_col$Value, decreasing = TRUE), ][1:top, ]
    
    # plot title
    if(top == nrow(df)){
      title = ""
    }
    else {
      title = paste("Top", top, "based on", df_col$Group, sep = " ")
    }
    
    p <- ggplot(df_col, aes(x = Value, y = reorder(Column, Value), fill = Group)) +
      geom_col() +
      ggtitle(title) +
      xlab(df_col$Group) +
      ylab("") +
      theme_minimal() +
      theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 8),
            legend.position = "none")
    
    print(p)
    # plot <- recordPlot()
    plot_name <- paste(colnames(df)[i], deparse(substitute(df)), top, sep = "")
    file_path <- file.path(dir, paste(plot_name, ".pdf", sep = ""))
    
    if(top == nrow(df)){
      ggsave(filename = file_path, plot = p)
    }
    else {
      ggsave(filename = file_path, plot = p, height = top * 0.3)
    }
    
  }
}

# HC Simulated Data Results ----------------------------------------------------
# scaled_df_simulated_features
hc_results(scaled_df_simulated_features, 4, 13)
save_plots(tree_plots_scaled_df_simulated_features, dir = "hc_tree_plots_scaled_df_simulated_features")
save_plots(cut_tree_plots_scaled_df_simulated_features, dir = "hc_cut_tree_plots_scaled_df_simulated_features")
evaluate_hc_clustering(scaled_df_simulated_features, cut_hc_results_scaled_df_simulated_features)
plot_assigned_raw_ts(cut_hc_results_scaled_df_simulated_features, as.data.frame(df_simulated))
save_plots(raw_ts_plots_cut_hc_results_scaled_df_simulated_features, "raw_ts_plots_cut_hc_results_scaled_df_simulated_features")
# evaluation outputs
hc_dtw_internal <- internal_evaluation_cut_hc_results_scaled_df_simulated_features[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_df_simulated_features)), ]
hc_euclidean_internal <- internal_evaluation_cut_hc_results_scaled_df_simulated_features[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_df_simulated_features)), ]

latex_table(hc_dtw_internal, top = 10)
create_evaluation_plots(internal_evaluation_cut_hc_results_scaled_df_simulated_features[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_df_simulated_features)), ], top = 5)
create_evaluation_plots(internal_evaluation_cut_hc_results_scaled_df_simulated_features[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_df_simulated_features)), ], top = 5)

# K-means Clustering Functions -------------------------------------------------

# Function, which performs k-means clustering with different number of clusters
# Returns a list of the results
kmeans_results <- function (df, a, b, initialization){
  df_name <- deparse(substitute(df))
  
  results <- list()
  plots = list()
  
  for (k in a:b){
    clusters <- kmeans_cluster(df, k, initialization = initialization, nstart = 50)
    clusters <- clusters$cluster
    clusters_name <- paste("clusters_kmeans", k, df_name, sep = "_")
    results[[clusters_name]] <- clusters
    
    par(mar = c(80/k, 4, 2, 2) + 0.1)
    # empty plot
    plot(1, 1, main = "", 
         xlim = c(1 - 0.05*length(clusters), length(clusters) + 0.05*length(clusters)),
         ylim = range(clusters), 
         type = "n", xlab = "", ylab = "Cluster ID", xaxt = "n", cex.main = 0.8, ygap.axis = 0.3)
    # Add faint y-axis grid lines
    clip(-1, length(clusters) + 2, min(clusters) - 1, max(clusters) + 1)
    abline(h = seq(min(clusters), max(clusters), by = 1), col = "lightgray", lty = "longdash", xpd = F)
    
    # Add points
    points(clusters, cex = 1) # col = colors, pch = 19,
    axis(1, at = 1:length(clusters), labels = rownames(df), las = 2, cex.axis = 0.7)
    # axis(2, at = seq(floor(min(clusters)), ceiling(max(clusters)), by = 1))
    axis(2, at = seq(1, length(unique(clusters)), by = 1))
    
    plot <- recordPlot()
    plots[[clusters_name]] <- plot
    par(mar = c(5, 4, 4, 2) + 0.1)
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
# Function, which evaluates kmeans clustering
evaluate_kmeans_clustering <- function(df, list){
  
  list_name <- deparse(substitute(list))
  
  internal_evaluation_df <- data.frame(matrix(nrow = length(list), ncol = 3))
  rownames(internal_evaluation_df) <- c(names(list))
  colnames(internal_evaluation_df) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
  
  for (j in 1:length(list)){
    clusters_name <- names(list)[j]
    internal_evaluation_df[clusters_name, ] <- unname(internal_evaluation(df, "euclidean", list[[j]], toGlobEnv = F))
  }
  
  int_df_name <- paste("internal_evaluation", list_name, sep = "_")
  print(int_df_name)
  assign(int_df_name, internal_evaluation_df, envir = .GlobalEnv)
}
# Function, which performs k-shape clustering with different number of clusters
# Returns a list of the results
kshape_results <- function (df, a, b){
  df_name <- deparse(substitute(df))
  
  results <- list()
  plots = list()
  
  for (k in a:b){
    clusters <- kshape_cluster(df, k, seed = 77)
    clusters <- clusters@cluster
    clusters_name <- paste("clusters_kshape", k, df_name, sep = "_")
    results[[clusters_name]] <- clusters
    
    par(mar = c(80/k, 4, 2, 2) + 0.1)
    # empty plot
    plot(1, 1, main = paste("Cluster Assignment:", clusters_name, sep = " "), 
         xlim = c(1 - 0.05*length(clusters), length(clusters) + 0.05*length(clusters)),
         ylim = range(clusters), 
         type = "n", xlab = "", ylab = "Cluster ID", xaxt = "n", cex.main = 0.8, ygap.axis = 0.3)
    # Add faint y-axis grid lines
    clip(-1, length(clusters) + 2, min(clusters) - 1, max(clusters) + 1)
    abline(h = seq(min(clusters), max(clusters), by = 1), col = "lightgray", lty = "longdash", xpd = F)
    
    # Add points
    points(clusters, cex = 1) # col = colors, pch = 19,
    axis(1, at = 1:length(clusters), labels = rownames(df), las = 2, cex.axis = 0.7)
    # axis(2, at = seq(floor(min(clusters)), ceiling(max(clusters)), by = 1))
    axis(2, at = seq(1, length(unique(clusters)), by = 1))
    
    plot <- recordPlot()
    plots[[clusters_name]] <- plot
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
  
  results_name <- paste("kshape_results", df_name, sep = "_")
  print(results_name)
  assign(results_name, results, envir = .GlobalEnv)
  
  plots_name <- paste("kshape_plots", df_name, sep = "_")
  print(plots_name)
  assign(plots_name, plots, envir = .GlobalEnv)
  
}
# Function, which evaluates k-shape clustering
evaluate_kshape_clustering <- function(df, list, true_labels = NULL){
  
  list_name <- deparse(substitute(list))
  
  internal_evaluation_df <- data.frame(matrix(nrow = length(list), ncol = 3))
  rownames(internal_evaluation_df) <- c(names(list))
  colnames(internal_evaluation_df) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
  
  for (j in 1:length(list)){
    clusters_name <- names(list)[j]
    internal_evaluation_df[clusters_name, ] <- unname(internal_evaluation(df, "sbd", list[[j]], toGlobEnv = F))
  }
  
  int_df_name <- paste("internal_evaluation", list_name, sep = "_")
  print(int_df_name)
  assign(int_df_name, internal_evaluation_df, envir = .GlobalEnv)
  
  # External evaluation
  if (!is.null(true_labels)){
    external_evaluation_df <- data.frame(matrix(nrow = length(list), ncol = 4))
    rownames(external_evaluation_df) <- c(names(list))
    colnames(external_evaluation_df) <- c("Adjusted Rand Index", "Jaccard Index", "Purity", "Normalized Mutual Information")
    for (l in 1:length(list)){
      clusters_name <- names(list)[l]
      external_evaluation_df[clusters_name, ] <- unname(external_evaluation(true_labels, list[[l]], toGlobEnv = F))
    }
    ext_df_name <- paste("external_evaluation", list_name, sep = "_")
    print(ext_df_name)
    assign(ext_df_name, external_evaluation_df, envir = .GlobalEnv)
  }
}

# K-means Simulated Data Results -----------------------------------------------
# scaled_df_simulated_features

# Determine k
elbow(scaled_df_simulated_features, 20)
save_plots(elbow_scaled_df_simulated_features)

# k-means
kmeans_results(scaled_df_simulated_features, 4, 13, initialization = NULL)
save_plots(kmeans_plots_scaled_df_simulated_features, "kmeans_plots_scaled_df_simulated_features")
# plot_assigned_raw_ts(kmeans_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(scaled_df_simulated_features, kmeans_results_scaled_df_simulated_features)
internal_evaluation_kmeans_results_scaled_df_simulated_features
# latex_table(internal_evaluation_kmeans_results_scaled_df_simulated_features)
# create_evaluation_plots(internal_evaluation_kmeans_results_scaled_df_simulated_features)

# k-means++
kmeans_results(scaled_df_simulated_features, 4, 13, initialization = "plus_plus")
save_plots(kmeans_plus_plus_plots_scaled_df_simulated_features, "kmeans_plus_plus_plots_scaled_df_simulated_features")
# plot_assigned_raw_ts(kmeans_plus_plus_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_plus_plus_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_plus_plus_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(scaled_df_simulated_features, kmeans_plus_plus_results_scaled_df_simulated_features)
rownames(internal_evaluation_kmeans_plus_plus_results_scaled_df_simulated_features) <- paste(rownames(internal_evaluation_kmeans_plus_plus_results_scaled_df_simulated_features), "_plus", sep = "")
# latex_table(internal_evaluation_kmeans_plus_plus_results_scaled_df_simulated_features)
# create_evaluation_plots(internal_evaluation_kmeans_plus_plus_results_scaled_df_simulated_features)

# k-means simfp
kmeans_results(scaled_df_simulated_features, 4, 13, initialization = "simfp_init")
save_plots(kmeans_simfp_init_plots_scaled_df_simulated_features, "kmeans_simfp_init_plots_scaled_df_simulated_features")
# plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_simfp_init_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_simfp_init_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(scaled_df_simulated_features, kmeans_simfp_init_results_scaled_df_simulated_features)
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_df_simulated_features) <- paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_df_simulated_features), "_simfp", sep = "")
# latex_table(internal_evaluation_kmeans_simfp_init_results_scaled_df_simulated_features)
# create_evaluation_plots(internal_evaluation_kmeans_simfp_init_results_scaled_df_simulated_features)

# k-means grep
kmeans_results(scaled_df_simulated_features, 4, 13, initialization = "grep_init")
save_plots(kmeans_grep_init_plots_scaled_df_simulated_features, "kmeans_grep_init_plots_scaled_df_simulated_features")
# plot_assigned_raw_ts(kmeans_grep_init_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_grep_init_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_grep_init_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(scaled_df_simulated_features, kmeans_grep_init_results_scaled_df_simulated_features)
rownames(internal_evaluation_kmeans_grep_init_results_scaled_df_simulated_features) <- paste(rownames(internal_evaluation_kmeans_grep_init_results_scaled_df_simulated_features), "_grep", sep = "")
# latex_table(internal_evaluation_kmeans_grep_init_results_scaled_df_simulated_features)
# create_evaluation_plots(internal_evaluation_kmeans_grep_init_results_scaled_df_simulated_features)

# k-shape
kshape_results(scaled_df_simulated_features, 4, 13)
save_plots(kshape_plots_scaled_df_simulated_features, "kshape_plots_scaled_df_simulated_features")
# plot_assigned_raw_ts(kshape_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kshape_results_scaled_df_simulated_features, "raw_ts_plots_kshape_results_scaled_df_simulated_features")
evaluate_kshape_clustering(scaled_df_simulated_features, kshape_results_scaled_df_simulated_features)
latex_table(internal_evaluation_kshape_results_scaled_df_simulated_features)

# Evaluate Euclidean distance partitions
simulated_euclidean_results <- rbind(hc_euclidean_internal, internal_evaluation_kmeans_results_scaled_df_simulated_features, 
                           internal_evaluation_kmeans_plus_plus_results_scaled_df_simulated_features, 
                           internal_evaluation_kmeans_simfp_init_results_scaled_df_simulated_features, 
                           internal_evaluation_kmeans_grep_init_results_scaled_df_simulated_features)

latex_table(simulated_euclidean_results, top = 15)

# robustness check using external evaluation criteria
best_rows_index <- c("DTW 13\\(CHI, DI)", "DTW 12\\(CHI, ASW)", "DTW 11\\(ASW, CHI)", "DTW 9\\(DI)", "DTW 10\\(DI)", "DTW 10\\(ASW)", "ED 13\\(CHI, DI)", "ED 12\\(CHI, ASW)", "ED 12\\(CHI)", "ED 9\\(DI)", "ED 10\\(DI)", "ED 11\\(ASW)", "ED 111\\(ASW)")

partitions <- rbind(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_13,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_12, 
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_DTW_11,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_9,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_10, 
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW_10,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_euclidean_13,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean_12, 
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_euclidean_12, 
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_euclidean_9,
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_average_dist_scaled_df_simulated_features_euclidean_10, 
                    cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean_11,
                    kmeans_grep_init_results_scaled_df_simulated_features$clusters_kmeans_11_scaled_df_simulated_features)

purity_matrix <- matrix(nrow = 13, ncol = 13)
for (row in 1:nrow(partitions)){
  for (row1 in 1:nrow(partitions)){
    purity_matrix[row, row1] <- external_validation(partitions[row, ], partitions[row1, ], method = "purity")
  }
}
purity_matrix <- as.data.frame(purity_matrix)
rownames(purity_matrix) <- best_rows_index
colnames(purity_matrix) <- best_rows_index
latex_table(purity_matrix)

jaccard_matrix <- matrix(nrow = 13, ncol = 13)
for (row in 1:nrow(partitions)){
  for (row1 in 1:nrow(partitions)){
    jaccard_matrix[row, row1] <- external_validation(partitions[row, ], partitions[row1, ], method = "jaccard_index")
  }
}
jaccard_matrix <- as.data.frame(jaccard_matrix)
rownames(jaccard_matrix) <- best_rows_index
colnames(jaccard_matrix) <- best_rows_index
latex_table(jaccard_matrix)

adjusted_rand_index_matrix <- matrix(nrow = 13, ncol = 13)
for (row in 1:nrow(partitions)){
  for (row1 in 1:nrow(partitions)){
    adjusted_rand_index_matrix[row, row1] <- external_validation(partitions[row, ], partitions[row1, ], method = "adjusted_rand_index")
  }
}
adjusted_rand_index_matrix <- as.data.frame(adjusted_rand_index_matrix)
rownames(adjusted_rand_index_matrix) <- best_rows_index
colnames(adjusted_rand_index_matrix) <- best_rows_index
latex_table(adjusted_rand_index_matrix)

# raw time series as clusters
# ch, dunn dtw 1
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_13, df_simulated, save = T)
# ch, dunn dtw 2
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_12, df_simulated, save = T)
# ch, dunn dtw 3
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_11, df_simulated, save = T)

# asw dtw 1
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_11, df_simulated, save = T)
# asw dtw 2
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_10, df_simulated, save = T)
# asw dtw 3
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_DTW_9, df_simulated, save = T)

# ch, dunn euclidean 1
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_euclidean_13, df_simulated, save = T)
# ch euclidean 2
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D2_dist_scaled_df_simulated_features_euclidean_12 , df_simulated, save = T)
# ch euclidean 3
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_complete_dist_scaled_df_simulated_features_euclidean_12, df_simulated, save = T)

#dunn euclidean 2
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_mcquitty_dist_scaled_df_simulated_features_euclidean_4, df_simulated, save = T)
#dunn euclidean 3
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_average_dist_scaled_df_simulated_features_euclidean_6, df_simulated, save = T)

# asw euclidean 1
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean_11, df_simulated, save = T)
# asw euclidean 2
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean_12, df_simulated, save = T)
# asw euclidean 3
clusters_plot(cut_hc_results_scaled_df_simulated_features$clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean_10, df_simulated, save = T)


# Alternative Methods ----------------------------------------------------------
# Raw-Based Clustering of df_simulated -----------------------------------------
hc_results(df_simulated, 4, 13)
save_plots(tree_plots_df_simulated, dir = "hc_tree_plots_df_simulated")
save_plots(cut_tree_plots_df_simulated, dir = "hc_cut_tree_plots_df_simulated")
evaluate_hc_clustering(df_simulated, cut_hc_results_df_simulated)
latex_table(internal_evaluation_cut_hc_results_df_simulated)
internal_evaluation_cut_hc_results_df_simulated
cut_tree_plots_df_simulated$clusters_hc_complete_dist_df_simulated_euclidean_13

# k-means
kmeans_results(df_simulated, 4, 13, initialization = NULL)
save_plots(kmeans_plots_df_simulated, "kmeans_plots_df_simulated")
# plot_assigned_raw_ts(kmeans_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(df_simulated, kmeans_results_df_simulated)
internal_evaluation_kmeans_results_df_simulated

# feature-based clustering of df_simulated with no tsa -------------------------
scaled_df_simulated_features_notsa <- scaled_df_simulated_features[ , !grepl("TSA", names(scaled_df_simulated_features))]
hc_results(scaled_df_simulated_features_notsa, 4, 13)
save_plots(tree_plots_scaled_df_simulated_features_notsa, dir = "hc_tree_plots_scaled_df_simulated_features_notsa")
save_plots(cut_tree_plots_scaled_df_simulated_features_notsa, dir = "hc_cut_tree_plots_scaled_df_simulated_features_notsa")
evaluate_hc_clustering(scaled_df_simulated_features_notsa, cut_hc_results_scaled_df_simulated_features_notsa)
latex_table(internal_evaluation_cut_hc_results_scaled_df_simulated_features_notsa)
internal_evaluation_cut_hc_results_scaled_df_simulated_features_notsa
cut_tree_plots_scaled_df_simulated_features_notsa$clusters_hc_complete_dist_scaled_df_simulated_features_notsa_DTW_11

# k-means
kmeans_results(scaled_df_simulated_features_notsa, 4, 13, initialization = NULL)
save_plots(kmeans_plots_scaled_df_simulated_features_notsa, "kmeans_plots_scaled_df_simulated_features_notsa")
# plot_assigned_raw_ts(kmeans_results_scaled_df_simulated_features, as.data.frame(df_simulated))
# save_plots(raw_ts_plots_kmeans_results_scaled_df_simulated_features, "raw_ts_plots_kmeans_results_scaled_df_simulated_features")
evaluate_kmeans_clustering(scaled_df_simulated_features_notsa, kmeans_results_scaled_df_simulated_features_notsa)
internal_evaluation_kmeans_results_scaled_df_simulated_features_notsa


## fuzzy c-means (not included in thesis) --------------------------------------
fuzzy_results <- function(df, a, b){
  
  df_name <- deparse(substitute(df))
  
  results <- list()
  plots = list()
  
  for(k in a:b){
    # cluster
    
    clusters <- fcm(df, k, nstart=50, iter.max=50, numseed = 7)
      # cmeans(df, centers = k, iter.max = 100)
      #Fclust(df, k, type = type, noise = T, stand = F, distance = dist)
    clusters_name <- paste("clusters_fuzzy", df_name, k, sep = "_")
    results[[clusters_name]] <- clusters
    
    # plot
    fuzzy_membership_matrix <- clusters$u
    membership_df <- as.data.frame(fuzzy_membership_matrix)
    membership_df$DataPoint <- rownames(membership_df)
    # Melt the data frame to a long format
    membership_long <- melt(membership_df, id.vars = "DataPoint")
    
    h_lines <- seq(0.5, length(unique(membership_long$DataPoint)) + 0.5, by = 1)
    v_lines <- seq(0.5, length(unique(membership_long$variable)) + 0.5, by = 1)
    # Plot the heatmap with a grid on top
    p <- ggplot(membership_long, aes(x = variable, y = DataPoint, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      labs(x = "Cluster", y = "", fill = "Membership Value", 
           title = "Heatmap of Membership Values")
    plots[[clusters_name]] <- p
    
    
  }
  
  results_name <- paste("fuzzy_results", df_name, sep = "_")
  plots_name <- paste("fuzzy_plots", df_name, sep = "_")
  
  assign(results_name, results, envir = .GlobalEnv)
  print(results_name)
  assign(plots_name, plots, envir = .GlobalEnv)
  print(plots_name)
}

fuzzy_results(scaled_df_simulated_features, 7, 13)
fuzzy_plots_scaled_df_simulated_features
fuzzy_results_scaled_df_simulated_features

fuzzy_cluster(scaled_df_simulated_features, 9)
clusters_fuzzy_standard_scaled_df_simulated_features_20$clus
fuzzy_membership_matrix <- clusters_fuzzy_standard_scaled_df_simulated_features_12$U
fclust::summary.fclust(clusters_fuzzy_standard_scaled_df_simulated_features_12)


average_membership <- colMeans(fuzzy_membership_matrix)
print(average_membership)
