# Outline
# 1 READ ME
# 2 Dependencies
# 3 Clustering Algorithms
#   3.1 k-Means
#   3.2 Hierarchical
#   3.3 Model-based
# 4 Visualizations
#   4.1 
#   4.2 
#   4.3 
# 5 Evaluation Criteria
#   5.1 External Evaluation Criteria
#   5.2 Internal Evaluation Criteria 

# 1 READ ME --------------------------------------------------------------------
# The following are functions 
# - which perform clustering and save the result in the Global(soon to be virtual?) Environment
# - which visualize the clusters
# - which evaluate the clustering performance using internal/external methods

# 2 Dependencies ---------------------------------------------------------------
load_package <- function(pack){
  if (!is.element(pack, installed.packages()[,1])){
    install.packages(pack, dependencies = TRUE)
  }
  library(pack, character.only = TRUE)
}

# 3 Clustering Algorithms ------------------------------------------------------

## 3.1 Partitioning ------------------------------------------------------------
### k-Means --------------------------------------------------------------------
# Elbow method to determine number of clusters

elbow <- function(df, max_num_clusters){
  
  load_package("purrr")
  load_package("dplyr")
  load_package("ggplot2")
  
  wss <- map_dbl(1:max_num_clusters, ~{kmeans(df, ., nstart=50,iter.max = 15 )$tot.withinss})
  n_clust <- 1:max_num_clusters
  elbow_df <- as.data.frame(cbind("n_clust" = n_clust, "wss" = wss))
  
  plot <- ggplot(elbow_df) +
    geom_line(aes(y = wss, x = n_clust), colour = "#82518c") +
    scale_x_continuous(breaks = min(elbow_df$n_clust):max(elbow_df$n_clust)) +
    theme_minimal()
  
  df_name <- deparse(substitute(df))
  plot_name <- paste("elbow", df_name, sep = "_")
  assign(plot_name, plot, envir = .GlobalEnv)
  
  return(plot)
}
# Example:
# elbow(scaled_df_simulated_features, 20)

# k-means ++ initialization
plus_plus <- function(df, k, seed=8) {
  set.seed(seed)
  centroids <- list(df[1,])
  
  for (i in 2:k) {
    dist_sq <- apply(df, 1, function(x) min(sapply(centroids, function(c) sum((x - c)^2))))
    probs <- dist_sq / sum(dist_sq)
    cumulative_probs <- cumsum(probs)
    r <- runif(1)
    
    for (j in 1:length(cumulative_probs)) {
      if (r < cumulative_probs[j]) {
        centroids[[i]] <- df[j,]
        break
      }
    }
  }
  
  return(do.call(rbind, centroids))
}
# Example:
# test_initialization <- plus_plus(scaled_df_simulated_features, 11)

# k-means clustering algorithm with boolean "initialization" whether to use k-means ++ initialization method
kmeans_cluster <- function(df, num_clusters, initialization = F, nstart = 50){
  load_package("stats")
  
  df_name <- deparse(substitute(df))
  
  if (initialization == T) {
    centroids <- plus_plus(df, num_clusters)
    clusters <- kmeans(df, centers = centroids, nstart = nstart)
    
    clusters_name <- paste("clusters_kmeans", df_name, num_clusters, "plus",sep = "_")
  }
  else {
    clusters <- kmeans(df, centers = num_clusters, iter.max = 15, nstart = nstart)
    clusters_name <- paste("clusters_kmeans", df_name, num_clusters, sep = "_")
  }
  print(clusters_name)
  assign(clusters_name, clusters, envir = .GlobalEnv)
}
# Example:
# kmeans_cluster(scaled_df_simulated_features, 11, initialization = T) # clusters_kmeans_scaled_df_simulated_features_11_plus
# cluster_assignment_plot(scaled_df_simulated_features, clusters_kmeans_scaled_df_simulated_features_11_plus$cluster)
# kmeans_cluster(scaled_df_simulated_features, 11) # clusters_kmeans_scaled_df_simulated_features_10
# cluster_assignment_plot(scaled_df_simulated_features, clusters_kmeans_scaled_df_simulated_features_11$cluster)

### k-Shape --------------------------------------------------------------------
# tsclust function does not have an initialization parameter
kshape_cluster <- function(df, num_clusters, seed = 77){
  load_package("dtwclust")
  
  clusters <- tsclust(df, "partitional", k = num_clusters, distance = "sbd", centroid = "shape", seed = seed)
  
  df_name <- deparse(substitute(df))
  clusters_name <- paste("clusters_kshape", df_name, num_clusters,sep = "_")
  print(clusters_name)
  assign(clusters_name, clusters, envir = .GlobalEnv)
}
# Example:
# kshape_cluster(scaled_df_simulated_features, 11)
# cluster_assignment_plot(scaled_df_simulated_features, clusters_kshape_scaled_df_simulated_features_11@cluster)

## fuzzy c-means ---------------------------------------------------------------
fuzzy_cluster <- function(df, num_clusters, dist){
  load_package("dtwclust")
  
  clusters <- tsclust(
    series = scaled_df_simulated_features,
    type = "fuzzy",
    k = 11,
    preproc = NULL,
    distance = dist,
    seed = 7,
    trace = FALSE,
    error.check = TRUE
  )
  
  df_name <- deparse(substitute(df))
  clusters_name <- paste("clusters_fuzzy", df_name, num_clusters,sep = "_")
  print(clusters_name)
  assign(clusters_name, clusters, envir = .GlobalEnv)
}
# Example:
# fuzzy_cluster(scaled_df_simulated_features, 11, "dtw_basic")
# cluster_assignment_plot(scaled_df_simulated_features, clusters_fuzzy_scaled_df_simulated_features_11@cluster)
# membership_matrix <- clusters_fuzzy_scaled_df_simulated_features_11@fcluster


## 3.2 Hierarchical ------------------------------------------------------------

# distance method options: 
# stats::"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
# dtwclust:: "DTW"

# agglomerative method options: 
# stats::"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
# "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

# Function which computes distance matrix and stores it in the list
compute_dist_matrix <- function(df, dist_method) {
  load_package("stats")

  dist_matrix <- dist(df, method=dist_method)
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  # assign the new distance matrix to the global environment
  dist_matrix_name <- paste("dist", df_name, dist_method, sep = "_")
  # assign the distance matrix to the list in the global environment
  dist_matrix_list <- dist_matrix
  assign(dist_matrix_name, dist_matrix, envir = .GlobalEnv)
  print(dist_matrix_name)
  return(dist_matrix)
  
}
# Example:
# compute_dist_matrix(scaled_df_simulated_features, "DTW")

# Function which produces hierarchical clustering
# Name of list of clusters can be interpreted as follows:
# "clusters" + abbreviation of clustering algorithm 
# + agglomeration method for hierarchical clustering 
# + distance method + name of data frame
hierarch_cluster <- function(dist_matrix = NULL, 
                            compute = F, df = NULL, dist_method = NULL, 
                             agglomeration_method) {
  load_package("stats")
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  
  # compute distance matrix, if not given as input
  if (compute == T) {
    dist <- compute_dist_matrix(df, dist_method = dist_method)
    dist_name <- paste("dist", df_name, dist_method, sep = "_")
  }
  else{
    dist <- dist_matrix
    dist_name <- deparse(substitute(dist_matrix))
  }
  clusters <- hclust(dist, method=agglomeration_method)
  
  clusters_name <- paste("clusters_hc", agglomeration_method, dist_name, sep = "_")
  assign(clusters_name, clusters, envir = .GlobalEnv)
  print(clusters_name)
  return(clusters)
}
# Example:
# hierarch_cluster(compute = T, df = scaled_df_simulated_features, dist_method = "euclidean", agglomeration_method = "ward.D2")

# Function which cuts hierarchical tree 
cut_clusters <- function(clusters, num_clusters) {
  load_package("stats")
  
  cut_clusters <- cutree(clusters, k=num_clusters)

  # get the name of the input data frame
  clusters_name <- deparse(substitute(clusters))
  cut_name <- paste("cut", num_clusters, clusters_name, sep = "_")
  assign(cut_name, cut_clusters, envir = .GlobalEnv)
  print(cut_name)
  return(cut_clusters)
}
# Example:
# cut_clusters(clusters_hc_ward.D2_euclidean_scaled_df_simulated_features, 11)

## 3.4 Model-based -------------------------------------------------------------

# 4 Visualizations -------------------------------------------------------------

# Function, which plots a data frame
df_plot <- function(df) {
  load_package("ggplot2")
  
  # Determine the global y range
  global_y_range <- range(df, na.rm = TRUE)
  
  # For each row in the dataframe
  for(i in 1:nrow(df)) {
    # Create a new dataframe for this row
    df_row <- data.frame(Column = colnames(df), Value = unlist(df[i, ]))
    df_row$Group <- rownames(df)[i]
    p <- ggplot(df_row, aes(Column, Value, group = Group)) +
      geom_line() +  # Use geom_point() here
      ggtitle(paste("Row:", rownames(df)[i])) +
      xlab("") +
      ylab("") +
      ylim(global_y_range) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    print(p)
  }
}
# Example:
# df_plot(as.data.frame(df_simulated))

# Function which plots the cluster assignment
cluster_assignment_plot <- function(df, cut_clusters){
  
  load_package("RColorBrewer")
  load_package("graphics")
  
  par(mar = c(11, 4, 2, 2) + 0.1)
  
  # empty plot
  plot(1, 1, xlim = c(1, length(cut_clusters)), ylim = range(cut_clusters), type = "n", xlab = "", ylab = "Cluster ID", xaxt = "n")
  # Add faint y-axis grid lines
  abline(h = seq(floor(min(cut_clusters)), ceiling(max(cut_clusters)), by = 1), col = "lightgray", lty = "longdash")
  # Add points
  colors <- hcl.colors(length(cut_clusters), palette = "SunsetDark", alpha = 1, rev = T)
  points(cut_clusters, cex = 1) # col = colors, pch = 19,
  axis(1, at = 1:length(cut_clusters), labels = rownames(df), las = 2)
  axis(2, at = seq(floor(min(cut_clusters)), ceiling(max(cut_clusters)), by = 1))
  
  plot <- recordPlot()
  cut_clusters_name <- deparse(substitute(cut_clusters))
  plot_name <- paste("plot", cut_clusters_name, sep = "_")
  assign(plot_name, plot, envir = .GlobalEnv)
  par(mar = c(5, 4, 4, 2) + 0.1)  # Default
  return(plot)
}
# Example:
# cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_ward.D2_euclidean_scaled_df_simulated_features)

# Function which plots a hierarchical tree
hc_plot <- function(clusters){
  
  par(mar = c(5, 4, 4, 2) + 0.1) 
  plot(clusters, ylab = "Distance", xlab = "")
  
  plot <- recordPlot()
  clusters_name <- deparse(substitute(clusters))
  plot_name <- paste("plot", clusters_name, sep = "_")
  assign(plot_name, plot, envir = .GlobalEnv)
  par(mar = c(5, 4, 4, 2) + 0.1) # Default
  return(plot)
}
# Example:
# hc_plot(clusters_hc_ward.D2_euclidean_scaled_df_simulated_features)

# Function which plots clusters of the original time series
clusters_plot <- function(cluster_object, df) {
  load_package("ggplot2")
  load_package("gridExtra")
  
  # cluster assignment
  df$cluster <- factor(cluster_object)
  
  # Get the unique cluster ids
  cluster_ids <- unique(df$cluster)
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # For each cluster id
  for(i in seq_along(cluster_ids)) {
    # Subset the data frame for the current cluster
    df_subset <- df[df$cluster == cluster_ids[i], ]
    df_subset$cluster <- NULL
    
    xlab <- paste("Cluster", cluster_ids[i], sep = " ")
    # Plot all rows of the original data frame
    colors <- hcl.colors(nrow(df_subset), palette = "Set 2", alpha = 1, rev = F) # "Geyser", "Zissou 1"
    matplot(t(df_subset), type = "l", lty = 1, col = colors, xlab = xlab, ylab = "")
    legend("topright" , legend = rownames(df_subset), col = colors, lty = 1, cex = 0.65)
    # Add the plot to the list
    p <- recordPlot()
    plots[[i]] <- p
    
  }
  clusters_name <- deparse(substitute(cluster_object))
  plots_name <- paste("series_plots", clusters_name, sep = "_")
  assign(plots_name, plots, envir = .GlobalEnv)
  print(plots_name)
}
# Example:
# clusters_plot(cut_9_clusters_hc_average_dist_scaled_df_simulated_features_DTW, as.data.frame(df_simulated))

# 5 Evaluation Criteria --------------------------------------------------------

internal_evaluation <- function(df, distance_method, clusters_object, toGlobEnv = T){
  load_package("fpc")
  
  dist <- compute_dist_matrix(df, distance_method)
  clust_stats <- cluster.stats(dist, clusters_object)
  
  int_eval_df <- data.frame(matrix(ncol = 3, nrow = 1))
  colnames(int_eval_df) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
  
  int_eval_df["Calinski-Harabasz Index"] <- clust_stats$ch
  int_eval_df["Dunn Index"] <- clust_stats$dunn
  int_eval_df["Average Silhouette Width"] <- clust_stats$avg.silwidth
  
  if (toGlobEnv == T){
    clusters_object_name <- deparse(substitute(clusters_object))
    df_name <- paste("int_eval", clusters_object_name, sep = "_")
    assign(df_name, int_eval_df, envir = .GlobalEnv)
    print(df_name)
  }
  
  return(int_eval_df)
  
}
# Example:
# internal_evaluation(scaled_df_simulated_features, "euclidean", clusters_kshape_scaled_df_simulated_features_9@cluster)
# internal_evaluation(scaled_df_simulated_features, "DTW", cut_11_clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW)
# internal_evaluation(scaled_df_simulated_features, "euclidean", clusters_kmeans_scaled_df_simulated_features_11$cluster)
# internal_evaluation(scaled_df_simulated_features, "euclidean", clusters_kmeans_scaled_df_simulated_features_11_plus)

external_evaluation <- function(true_labels, clusters_object, toGlobEnv = T){
  load_package("ClusterR")
  
  ext_eval_df <- data.frame(matrix(ncol = 4, nrow = 1))
  colnames(ext_eval_df) <- c("Adjusted Rand Index", "Jaccard Index", "Purity", "Normalized Mutual Information")
  
  ext_eval_df["Adjusted Rand Index"] <- external_validation(true_labels, clusters_object, method = "adjusted_rand_index")
  ext_eval_df["Jaccard Index"] <- external_validation(true_labels, clusters_object, method = "jaccard_index")
  ext_eval_df["Purity"] <- external_validation(true_labels, clusters_object, method = "purity")
  ext_eval_df["Normalized Mutual Information"] <- external_validation(true_labels, clusters_object, method = "nmi")
  
  if (toGlobEnv == T){
    clusters_object_name <- deparse(substitute(clusters_object))
    df_name <- paste("ext_eval", clusters_object_name, sep = "_")
    assign(df_name, ext_eval_df, envir = .GlobalEnv)
    print(df_name)
  }
  
  return(ext_eval_df)
  
}

# Example:
# internal_evaluation(simulated_data_labels_numeric, clusters_kshape_scaled_df_simulated_features_9@cluster)
