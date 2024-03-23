# Dependencies -----------------------------------------------------------------
source("SARFIMA-GARCH_Simulation.R")
source("simulation_plots.R")
source("feature_extractions.R")
source("feature_df.R")
source("clustering_algorithms.R")

# Extract Features -------------------------------------------------------------
extract_features(df_simulated, T) # returns df_simulated_features
df_plot(df_simulated_features)

rescale_features(df_simulated_features, T) # returns scaled_df_simulated_features
df_plot(scaled_df_simulated_features)
# HIERARCHICAL -----------------------------------------------------------------

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
    
    for(agglomeration_method in c("complete", "average", "mcquitty", "ward.D", "ward.D2")){
      clusters <- hierarch_cluster(dist_matrix = dist_matrix, 
                       compute = F, agglomeration_method = agglomeration_method)
      clusters_name <- paste("clusters_hc", agglomeration_method, dist_matrix_name, sep = "_")
      # Store the clusters in the results list
      results[[clusters_name]] <- clusters
      
      # Tree Plot
      par(mar = c(5, 4, 4, 2) + 0.1) 
      plot(clusters, main = paste("Tree:", clusters_name, sep = " "), 
           ylab = "Distance", xlab = "", cex.main = 0.8)
      tree <- recordPlot()
      plots[[clusters_name]] <- tree
      
      # Cut Tree
      for (k in a:b){
        cut_clusters <- cutree(clusters, k = k)
        cut_name <- paste(clusters_name, k, sep = "_")
        cut_results[[cut_name]] <- cut_clusters
        
        # Plots
        par(mar = c(11, 4, 2, 2) + 0.1)
        # empty plot
        plot(1, 1, main = paste("Cut Tree:", cut_name, sep = " "), 
             xlim = c(1, length(cut_clusters)), ylim = range(cut_clusters), 
             type = "n", xlab = "", ylab = "Cluster ID", xaxt = "n", cex.main = 0.8)
        # Add faint y-axis grid lines
        abline(h = seq(floor(min(cut_clusters)), ceiling(max(cut_clusters)), 
                       by = 1), col = "lightgray", lty = "longdash")
        # Add points
        points(cut_clusters, cex = 1) # col = colors, pch = 19,
        axis(1, at = 1:length(cut_clusters), labels = rownames(df), las = 2)
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

# Function, which evaluates clustering
evaluate_clustering <- function(df, list, true_labels = NULL){
  
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

# Function, which visualizes the results of the clustering evaluation
# coming soon

# K-MEANS ----------------------------------------------------------------------
# scaled_df_simulated_features

## simple k-means
elbow(scaled_df_simulated_features, 20)
kmeans_cluster(scaled_df_simulated_features, 9, initialization = T, nstart = 50)

## k-shape
kshape_cluster(scaled_df_simulated_features, 9, seed = 77)

## k-means++
plus_plus(df, k, seed=8)
kmeans_cluster(df, num_clusters, initialization = F, nstart = 50)

## fuzzy c-means
fuzzy_cluster(scaled_df_simulated_features, 9, "euclidean")

## visualizations
cluster_assignment_plot(scaled_df_simulated_features, clusters_kshape_scaled_df_simulated_features_9@cluster)
clusters_plot(clusters_kmeans_scaled_df_simulated_features_9$cluster, as.data.frame(df_simulated))


# Simulated Data Results -------------------------------------------------------
# scaled_df_simulated_features
hc_results(scaled_df_simulated_features, 7, 12)
save_plots(tree_plots_scaled_df_simulated_features, dir = "hc_tree_plots_scaled_df_simulated_features")
save_plots(cut_tree_plots_scaled_df_simulated_features, dir = "hc_cut_tree_plots_scaled_df_simulated_features")
evaluate_clustering(scaled_df_simulated_features, cut_hc_results_scaled_df_simulated_features, true_labels = simulated_data_labels_numeric)

# df_simulated
hc_results(df_simulated, 7, 12)
save_plots(tree_plots_df_simulated, dir = "hc_tree_plots_df_simulated")
save_plots(cut_tree_plots_df_simulated, dir = "hc_cut_tree_plots_df_simulated")
