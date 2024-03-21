# Dependencies -----------------------------------------------------------------
source("SARFIMA-GARCH_Simulation.R")
source("simulation_plots.R")
source("feature_extractions.R")
source("feature_df.R")
source("clustering_algorithms.R")

# Extract Features -------------------------------------------------------------
extract_features(df_simulated, T) # returns df_simulated_features
rescale_features(df_simulated_features, T)

# list to store all distance matrices
dist_matrix_list <- list()
# list to store all hc clusters
hc_clusters_list <- list()
# list to store all hc dendograms
hc_dendograms_list <- list()
# compute_dist_matrix <- function(df, dist_method)
# hierarch_cluster <- function(dist_matrix = NULL, 
#                                compute = F, df = NULL, dist_method = NULL, 
#                                agglomeration_method)
# hc_plot <- function(clusters)
# cut_clusters <- function(clusters, num_clusters)
# cluster_assignment_plot <- function(df, cut_clusters)

# distance method options: 
# stats::"euclidea", "maximum", "manhattan", "canberra", "binary" or "minkowski"
# dtwclust:: "DTW"

# agglomerative method options: 
# stats::"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
# "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

# df_simulated -----------------------------------------------------------------
dist_matrix_list <- compute_dist_matrix(df_simulated, "euclidean") # dist_df_simulated_euclidean
dist_matrix_list <- compute_dist_matrix(df_simulated, "DTW") # dist_df_simulated_DTW

## eucledian --------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_df_simulated_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_df_simulated_euclidean

hc_plot(clusters_hc_complete_dist_df_simulated_euclidean)
hc_plot(clusters_hc_average_dist_df_simulated_euclidean)
hc_plot(clusters_hc_mcquitty_dist_df_simulated_euclidean)
hc_plot(clusters_hc_ward.D_dist_df_simulated_euclidean)
hc_plot(clusters_hc_ward.D2_dist_df_simulated_euclidean)
hc_plot(clusters_hc_centroid_dist_df_simulated_euclidean)
hc_plot(clusters_hc_median_dist_df_simulated_euclidean)

cut_clusters(clusters_hc_complete_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_complete_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_average_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_average_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_mcquitty_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_mcquitty_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_ward.D_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_ward.D_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_ward.D2_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_ward.D2_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_centroid_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_centroid_dist_df_simulated_euclidean)
cut_clusters(clusters_hc_median_dist_df_simulated_euclidean, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_median_dist_df_simulated_euclidean)

## DTW --------------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_df_simulated_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_df_simulated_DTW

hc_plot(clusters_hc_complete_dist_df_simulated_DTW)
hc_plot(clusters_hc_average_dist_df_simulated_DTW)
hc_plot(clusters_hc_mcquitty_dist_df_simulated_DTW)
hc_plot(clusters_hc_ward.D_dist_df_simulated_DTW)
hc_plot(clusters_hc_ward.D2_dist_df_simulated_DTW)
hc_plot(clusters_hc_centroid_dist_df_simulated_DTW)
hc_plot(clusters_hc_median_dist_df_simulated_DTW)

cut_clusters(clusters_hc_complete_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_complete_dist_df_simulated_DTW)
cut_clusters(clusters_hc_average_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_average_dist_df_simulated_DTW)
cut_clusters(clusters_hc_mcquitty_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_mcquitty_dist_df_simulated_DTW)
cut_clusters(clusters_hc_ward.D_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_ward.D_dist_df_simulated_DTW)
cut_clusters(clusters_hc_ward.D2_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_ward.D2_dist_df_simulated_DTW)
cut_clusters(clusters_hc_centroid_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_centroid_dist_df_simulated_DTW)
cut_clusters(clusters_hc_median_dist_df_simulated_DTW, 11)
cluster_assignment_plot(df_simulated, cut_11_clusters_hc_median_dist_df_simulated_DTW)

# df_simulated_features --------------------------------------------------------

dist_matrix_list <- compute_dist_matrix(df_simulated_features, "euclidean") # dist_df_simulated_features_euclidean
dist_matrix_list <- compute_dist_matrix(df_simulated_features, "DTW") # dist_df_simulated_features_DTW

## eucledian --------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_df_simulated_features_euclidean

hc_plot(clusters_hc_complete_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_average_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_mcquitty_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_ward.D_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_ward.D2_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_centroid_dist_df_simulated_features_euclidean)
hc_plot(clusters_hc_median_dist_df_simulated_features_euclidean)

cut_clusters(clusters_hc_complete_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_complete_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_average_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_average_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_mcquitty_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_mcquitty_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_ward.D_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_ward.D_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_ward.D2_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_ward.D2_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_centroid_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_centroid_dist_df_simulated_features_euclidean)
cut_clusters(clusters_hc_median_dist_df_simulated_features_euclidean, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_median_dist_df_simulated_features_euclidean)

## DTW --------------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_df_simulated_features_DTW

hc_plot(clusters_hc_complete_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_average_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_mcquitty_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_ward.D_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_ward.D2_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_centroid_dist_df_simulated_features_DTW)
hc_plot(clusters_hc_median_dist_df_simulated_features_DTW)

cut_clusters(clusters_hc_complete_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_complete_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_average_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_average_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_mcquitty_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_mcquitty_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_ward.D_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_ward.D_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_ward.D2_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_ward.D2_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_centroid_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_centroid_dist_df_simulated_features_DTW)
cut_clusters(clusters_hc_median_dist_df_simulated_features_DTW, 11)
cluster_assignment_plot(df_simulated_features, cut_11_clusters_hc_median_dist_df_simulated_features_DTW)

# scaled_df_simulated_features -------------------------------------------------
dist_matrix_list <- compute_dist_matrix(scaled_df_simulated_features, "euclidean") # dist_scaled_df_simulated_features_euclidean
dist_matrix_list <- compute_dist_matrix(scaled_df_simulated_features, "DTW") # dist_scaled_df_simulated_features_DTW

# eucledian --------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_scaled_df_simulated_features_euclidean
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_euclidean, 
                 compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_scaled_df_simulated_features_euclidean

hc_plot(clusters_hc_complete_dist_scaled_df_simulated_features_euclidean) # meh
hc_plot(clusters_hc_average_dist_scaled_df_simulated_features_euclidean) # meh
hc_plot(clusters_hc_mcquitty_dist_scaled_df_simulated_features_euclidean) # bad
hc_plot(clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean) # meh
hc_plot(clusters_hc_ward.D2_dist_scaled_df_simulated_features_euclidean) # meh
hc_plot(clusters_hc_centroid_dist_scaled_df_simulated_features_euclidean) # bad
hc_plot(clusters_hc_median_dist_scaled_df_simulated_features_euclidean) # bad

cut_clusters(clusters_hc_complete_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_complete_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_average_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_average_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_mcquitty_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_mcquitty_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_ward.D_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_ward.D2_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_ward.D2_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_centroid_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_centroid_dist_scaled_df_simulated_features_euclidean)
cut_clusters(clusters_hc_median_dist_scaled_df_simulated_features_euclidean, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_median_dist_scaled_df_simulated_features_euclidean)

# DTW --------------------------------------------------------------------------
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "complete") # clusters_hc_complete_dist_scaled_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "average") # clusters_hc_average_dist_scaled_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "mcquitty") # clusters_hc_mcquitty_dist_scaled_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D") # clusters_hc_ward.D_dist_scaled_df_simulated_features_DTW
hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
                 compute = F, df = NULL, dist_method = NULL, "ward.D2") # clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW
# hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
#                  compute = F, df = NULL, dist_method = NULL, "centroid") # clusters_hc_centroid_dist_scaled_df_simulated_features_DTW
# hierarch_cluster(dist_matrix = dist_scaled_df_simulated_features_DTW, 
#                  compute = F, df = NULL, dist_method = NULL, "median") # clusters_hc_median_dist_scaled_df_simulated_features_DTW

hc_plot(clusters_hc_complete_dist_scaled_df_simulated_features_DTW)
hc_plot(clusters_hc_average_dist_scaled_df_simulated_features_DTW) # kind of
hc_plot(clusters_hc_mcquitty_dist_scaled_df_simulated_features_DTW)
hc_plot(clusters_hc_ward.D_dist_scaled_df_simulated_features_DTW) # this looks good i think
hc_plot(clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW)
# hc_plot(clusters_hc_centroid_dist_scaled_df_simulated_features_DTW)
# hc_plot(clusters_hc_median_dist_scaled_df_simulated_features_DTW)

cut_clusters(clusters_hc_complete_dist_scaled_df_simulated_features_DTW, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_complete_dist_scaled_df_simulated_features_DTW)
cut_clusters(clusters_hc_average_dist_scaled_df_simulated_features_DTW, 9)
cluster_assignment_plot(scaled_df_simulated_features, cut_9_clusters_hc_average_dist_scaled_df_simulated_features_DTW)
cut_clusters(clusters_hc_mcquitty_dist_scaled_df_simulated_features_DTW, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_mcquitty_dist_scaled_df_simulated_features_DTW)
cut_clusters(clusters_hc_ward.D_dist_scaled_df_simulated_features_DTW, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_ward.D_dist_scaled_df_simulated_features_DTW)
cut_clusters(clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW, 11)
cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_ward.D2_dist_scaled_df_simulated_features_DTW)
# cut_clusters(clusters_hc_centroid_dist_scaled_df_simulated_features_DTW, 11)
# cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_centroid_dist_scaled_df_simulated_features_DTW)
# cut_clusters(clusters_hc_median_dist_scaled_df_simulated_features_DTW, 11)
# cluster_assignment_plot(scaled_df_simulated_features, cut_11_clusters_hc_median_dist_scaled_df_simulated_features_DTW)


