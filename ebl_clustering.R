# The following estimates clusters of bonds. The provisioned bond liquidity forecasts 
# are prepared for comparison with the Diebold-Mariano Test. The forecast loss
# differenial of each pair of forecasts is extracted and clustered using the proposed
# feature-based method.

# k-means and the agglomerative hierarchical clustering algorithms are initialized
# with several initialization functions and agglomeration criteria, respectively.
# The results are evaluated using internal evaluation criteria and their similarity 
# is compared using external evaluation criteria. The more similar the produced 
# partitions are, the more confident we can be of the estimated clusters' presence
# in the data.

# Dependencies
source("feature_extractions.R")
source("clustering_alorithms.R")
source("ebl_clustering_functions.R")

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

# 1 DATA PREPARATION -----------------------------------------------------------
# Data Exploration
ebl <- read.csv("expected_bond_liquidity_all_models.csv") # replace with directory of data
head(ebl)
length(unique(ebl$cusip_id))
sort(names(ebl))

# Data Preprocessing of Forecasts
# Only consider restricted models (less missingness)
other_cols_to_keep <- c("date", "cusip_id")
ei_cols <- grepl("EI", colnames(ebl))
size_cols <- grepl("size", colnames(ebl))
restricted_cols <- grepl("restricted", colnames(ebl))
ei_size_restricted_cols <- ei_cols&size_cols&restricted_cols
all_cols <- c(other_cols_to_keep, colnames(ebl)[ei_size_restricted_cols])
restricted_ebl <- ebl[, all_cols]
names(restricted_ebl)

# Create a data frame out of every restricted model column
# rows: unique "cupis_id"
# cols: "date"
model_list <- list()

model_list <- mclapply(colnames(ebl)[ei_size_restricted_cols], function(col) {
  # data frame to long format
  long_df <- reshape2::melt(restricted_ebl, id.vars = c("date", "cusip_id"), measure.vars = col)
  # back to wide format
  wide_df <- reshape2::dcast(long_df, cusip_id ~ date, value.var = "value")
  
  # row names are cusip id + _model
  rownames(wide_df) <- wide_df$cusip_id
  wide_df$cusip_id <- NULL
  return(wide_df)
}, mc.cores = detectCores() - 1)
names(model_list) <- colnames(ebl)[ei_size_restricted_cols]

dim(model_list$EI_size_restricted_model)

# Missing Values
plot_na(model_list$EI_size_restricted_model) # only plot for one, missingness is the same for every df in list
plot_non_na_column(model_list$EI_size_restricted_model)
plot_non_na_row(model_list$EI_size_restricted_model)

## 1.1 Extract Balanced Data Panel ---------------------------------------------

# forecasts
over_x_observations(model_list, 30)
top_y_rows(model_list_over_30, 104)
remove_rows_with_na(model_list_over_30_top_104, 0)
dim(model_list_over_30_top_104_0_NAs$EI_size_restricted_model) # checks
# NAIVE one-step-ahead forecast
model_list_over_30_top_104_0_NAs[["EI_NAIVE_size_restricted_model"]] <- 
  t(apply(realized_bas_size_over_30_top_104_0_NAs, 1, function(x) rep(x[1], length(x))))

# realized_bas_size
# Only keep bond ids with more than x observations
realized_bas_size_over_30 <- realized_bas_size[rowSums(!is.na(realized_bas_size)) >= 30, ]
# Only keep the top y columns with the least missingness
na_counts <- colSums(is.na(realized_bas_size_over_30))
# Sort the columns by the number of NA values
sorted_columns <- sort(na_counts, index.return = TRUE)$ix
# Select the top y columns with the least NA values
selected_columns <- sorted_columns[1:104]
# Reorder the selected columns to match the original order
selected_columns <- selected_columns[order(selected_columns)]
# Subset the dataframe to the selected columns
realized_bas_size_over_30_top_104 <- realized_bas_size_over_30[, selected_columns]
# Remove rows with missing values
na_counts <- rowSums(is.na(realized_bas_size_over_30_top_104))
rows_to_remove <- which(na_counts > 0)
realized_bas_size_over_30_top_104_0_NAs <- realized_bas_size_over_30_top_104[-rows_to_remove, ]
# checks
any(is.na(realized_bas_size_over_30_top_104_0_NAs))
dim(realized_bas_size_over_30_top_104_0_NAs)
colnames(realized_bas_size_over_30_top_104_0_NAs)


## 1.2 Compute the Forecast Loss of Every Forecast -----------------------------
fc_loss_list_alt <- list()

for (model in names(model_list_over_30_top_104_0_NAs)){
  fc_loss_list_alt[[model]] <- (realized_bas_size_over_30_top_104_0_NAs - model_list_over_30_top_104_0_NAs[[model]])^2
}

## 1.3 Compute the Forecast Loss Differentials for All Model Pairs -------------
delta_list_alt <- list()

model_pairs <- combn(1:length(fc_loss_list_alt), 2)

for (i in seq(ncol(model_pairs))) {
  m1 <- model_pairs[1, i]
  m2 <- model_pairs[2, i]
  
  if (names(fc_loss_list_alt)[m1] == "EI_size_restricted_model") {
    m1_name <- "GBRT"
  } else {
    m1_name <- gsub("EI_|_size_restricted_model", "", names(fc_loss_list_alt)[m1])
  }
  m2_name <- gsub("EI_|_size_restricted_model", "", names(fc_loss_list_alt)[m2])
  delta_name <- paste(m1_name, m2_name, "delta", sep = "_")
  
  delta_list_alt[[delta_name]] <- fc_loss_list_alt[[m1]] - fc_loss_list_alt[[m2]]
}
# plot deltas
df_plot(delta_list_alt$MIDAS_NAIVE_delta)

# Remove Beginning of Corona From Data Panel
delta_list_alt_100 <- list()
for(i in seq_along(delta_list_alt)) {
  # Split the data frame into two based on the column index
  delta_list_alt_100[[i]] <- delta_list_alt[[i]][, 1:100]
}
names(delta_list_alt_100) <- names(delta_list_alt)

# peak
delta_list_alt_4 <- list()
for(i in seq_along(delta_list_alt)) {
  # Split the data frame into two based on the column index
  delta_list_alt_4[[i]] <- delta_list_alt[[i]][, 101:104]
}
names(delta_list_alt_4) <- names(delta_list_alt)

## 1.4 Extract Features --------------------------------------------------------
delta_list_alt_features <- list()
delta_list_alt_features[["GBRT_ENET_delta"]] <- extract_features(delta_list_alt_100[["GBRT_ENET_delta"]], tsa = T)
delta_list_alt_features[["GBRT_AR1_delta"]] <- extract_features(delta_list_alt_100[["GBRT_AR1_delta"]], tsa = T)
delta_list_alt_features[["GBRT_LIN_delta"]] <- extract_features(delta_list_alt_100[["GBRT_LIN_delta"]], tsa = T)
delta_list_alt_features[["GBRT_MIDAS_delta"]] <- extract_features(delta_list_alt_100[["GBRT_MIDAS_delta"]], tsa = T)
delta_list_alt_features[["GBRT_NAIVE_delta"]] <- extract_features(delta_list_alt_100[["GBRT_NAIVE_delta"]], tsa = T)
delta_list_alt_features[["ENET_AR1_delta"]] <- extract_features(delta_list_alt_100[["ENET_AR1_delta"]], tsa = T)
delta_list_alt_features[["ENET_LIN_delta"]] <- extract_features(delta_list_alt_100[["ENET_LIN_delta"]], tsa = T)
delta_list_alt_features[["ENET_MIDAS_delta"]] <- extract_features(delta_list_alt_100[["ENET_MIDAS_delta"]], tsa = T)
delta_list_alt_features[["ENET_NAIVE_delta"]] <- extract_features(delta_list_alt_100[["ENET_NAIVE_delta"]], tsa = T)
delta_list_alt_features[["AR1_LIN_delta"]] <- extract_features(delta_list_alt_100[["AR1_LIN_delta"]], tsa = T)
delta_list_alt_features[["AR1_MIDAS_delta"]] <- extract_features(delta_list_alt_100[["AR1_MIDAS_delta"]], tsa = T)
delta_list_alt_features[["AR1_NAIVE_delta"]] <- extract_features(delta_list_alt_100[["AR1_NAIVE_delta"]], tsa = T)
delta_list_alt_features[["LIN_MIDAS_delta"]] <- extract_features(delta_list_alt_100[["LIN_MIDAS_delta"]], tsa = T)
delta_list_alt_features[["LIN_NAIVE_delta"]] <- extract_features(delta_list_alt_100[["LIN_NAIVE_delta"]], tsa = T)
delta_list_alt_features[["MIDAS_NAIVE_delta"]] <- extract_features(delta_list_alt_100[["MIDAS_NAIVE_delta"]], tsa = T)

# rescale features
delta_list_alt_scaled_features <- lapply(delta_list_alt_features, rescale_features)
delta_list_alt_scaled_features$ # check
  
  # elbow method on all data frames
  for (delta_features in delta_list_alt_scaled_features){
    print(elbow(delta_features, nrow(delta_features)/2))
  }

# 2 CLUSTERING -----------------------------------------------------------------
## 2.1 GBRT_ENET_delta -------------------------------------------------------

### Feature-Based Clustering
elbow_scaled_GBRT_ENET_delta
# Extract Features as Objects
scaled_GBRT_ENET_delta_alt <- delta_list_alt_scaled_features$GBRT_ENET_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_GBRT_ENET_delta_alt <- dist(scaled_GBRT_ENET_delta_alt, "euclidean")
dtw_scaled_GBRT_ENET_delta_alt <- dist(scaled_GBRT_ENET_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_GBRT_ENET_delta_alt, delta_list_alt_100$GBRT_ENET_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_GBRT_ENET_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_GBRT_ENET_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_GBRT_ENET_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_GBRT_ENET_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_ENET_delta_alt, kmeans_results_scaled_GBRT_ENET_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_GBRT_ENET_delta_alt, delta_list_alt_100$GBRT_ENET_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_ENET_delta_alt, kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_GBRT_ENET_delta_alt, delta_list_alt_100$GBRT_ENET_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_GBRT_ENET_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_GBRT_ENET_delta_alt, cut_hc_results_scaled_GBRT_ENET_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_GBRT_ENET_delta_alt, cut_hc_results_scaled_GBRT_ENET_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt)), ]

#### Internal Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt_DTW, 5)

# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt), "_simfp", sep = "")

internal_evaluation_GBRT_ENET_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_GBRT_ENET_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt, internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt_ED)
latex_table(internal_evaluation_GBRT_ENET_delta_alt_ED, 5)

plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3, delta_list_alt_100$GBRT_ENET_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt, delta_list_alt_100$GBRT_ENET_delta, T, "kmeans_48_scaled_GBRT_ENET_delta_alt_simfp")

#### External Evaluation

GBRT_ENET_delta_top3 <- top_partitions(internal_evaluation_GBRT_ENET_delta_alt_ED, 
                                       internal_evaluation_cut_hc_results_scaled_GBRT_ENET_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

GBRT_ENET_delta_top3_partitions <- rbind(kmeans_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                         kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt,
                                         kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_49_scaled_GBRT_ENET_delta_alt,
                                         kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_41_scaled_GBRT_ENET_delta_alt,
                                         kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                         kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_4_scaled_GBRT_ENET_delta_alt,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_DTW_3,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_4,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_50,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_49,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_48,
                                         cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_3)


GBRT_ENET_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(GBRT_ENET_delta_top3_partitions)){
  for (row1 in 1:nrow(GBRT_ENET_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(GBRT_ENET_delta_top3_partitions[row, ], 
                                                              GBRT_ENET_delta_top3_partitions[row1, ])
  }
}

GBRT_ENET_purity_matrix <- as.data.frame(GBRT_ENET_purity_matrix)
rownames(GBRT_ENET_purity_matrix) <- best_rows_index
colnames(GBRT_ENET_purity_matrix) <- best_rows_index
latex_table(GBRT_ENET_purity_matrix)
# raw-based clustering ---------------------------------------------------------
GBRT_ENET_delta_alt_raw <- delta_list_alt_100[["GBRT_ENET_delta"]]
euclidean_GBRT_ENET_delta_alt_raw <- dist(GBRT_ENET_delta_alt_raw, "euclidean")
dtw_GBRT_ENET_delta_alt_raw <- dist(GBRT_ENET_delta_alt_raw, "DTW")
elbow(GBRT_ENET_delta_alt_raw, 70)

#### k-means with random centroids initialization
kmeans_results(GBRT_ENET_delta_alt_raw, delta_list_alt_100$GBRT_ENET_delta, 3, 20, NULL)
kmeans_plots_GBRT_ENET_delta_alt_raw$kmeans_4_GBRT_ENET_delta_alt_raw_1 # check
internal_evaluation_kmeans_results_GBRT_ENET_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_GBRT_ENET_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_results_GBRT_ENET_delta_alt_raw)) {
  
  internal_evaluation_kmeans_results_GBRT_ENET_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_GBRT_ENET_delta_alt_raw, kmeans_results_GBRT_ENET_delta_alt_raw[[cluster]], F)
  
}

#### k-means with simfp initialization
kmeans_results(GBRT_ENET_delta_alt_raw, delta_list_alt_100$GBRT_ENET_delta, 3, 20, "simfp_init")
kmeans_simfp_init_plots_GBRT_ENET_delta_alt_raw$kmeans_4_GBRT_ENET_delta_alt_raw_4 # check
internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw)) {
  
  internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_GBRT_ENET_delta_alt_raw, kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw[[cluster]], F)
  
}

kmeans_simfp_init_plots_GBRT_ENET_delta_alt_raw$kmeans_6_GBRT_ENET_delta_raw_6

#### Hierarchical
hc_results_ebl(GBRT_ENET_delta_alt_raw, delta_list_alt_100$GBRT_ENET_delta, 3, 20)
cut_hc_results_GBRT_ENET_delta_alt_raw$hc_ward.D_GBRT_ENET_delta_alt_raw_DTW_4

internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_GBRT_ENET_delta_alt_raw)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw[cluster, ] <- 
      internal_evaluation(euclidean_GBRT_ENET_delta_alt_raw, cut_hc_results_GBRT_ENET_delta_alt_raw[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw[cluster, ] <- 
      internal_evaluation(dtw_GBRT_ENET_delta_alt_raw, cut_hc_results_GBRT_ENET_delta_alt_raw[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw_ED <- internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw)), ]
internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw_DTW <- internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw[grepl("DTW", rownames(internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw_DTW, 5)

# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw), "_simfp", sep = "")

internal_evaluation_GBRT_ENET_delta_alt_raw_ED <- rbind(internal_evaluation_kmeans_results_GBRT_ENET_delta_alt_raw, 
                                                        internal_evaluation_kmeans_simfp_init_results_GBRT_ENET_delta_alt_raw, 
                                                        internal_evaluation_cut_hc_results_GBRT_ENET_delta_alt_raw_ED)
latex_table(internal_evaluation_GBRT_ENET_delta_alt_raw_ED, 5)

plot_assigned_raw_ts(kmeans_results_GBRT_ENET_delta_alt_raw$kmeans_3_GBRT_ENET_delta_alt_raw, delta_list_alt_100$GBRT_ENET_delta)

plot_assigned_raw_ts(cut_hc_results_GBRT_ENET_delta_alt_raw$hc_complete_GBRT_ENET_delta_alt_raw_DTW_8, delta_list_alt_100$GBRT_ENET_delta)
series_plots_hc_GBRT_ENET_delta_alt_raw$hc_complete_GBRT_ENET_delta_alt_raw_DTW_6_6

## 2.2 GBRT_AR1_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_GBRT_AR1_delta_alt <- delta_list_alt_scaled_features$GBRT_AR1_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_GBRT_AR1_delta_alt <- dist(scaled_GBRT_AR1_delta_alt, "euclidean")
dtw_scaled_GBRT_AR1_delta_alt <- dist(scaled_GBRT_AR1_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_GBRT_AR1_delta_alt, delta_list_alt_100$GBRT_AR1_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_GBRT_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_GBRT_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_GBRT_AR1_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_GBRT_AR1_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_AR1_delta_alt, kmeans_results_scaled_GBRT_AR1_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_GBRT_AR1_delta_alt, delta_list_alt_100$GBRT_AR1_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_AR1_delta_alt, kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_GBRT_AR1_delta_alt, delta_list_alt_100$GBRT_AR1_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_GBRT_AR1_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_GBRT_AR1_delta_alt, cut_hc_results_scaled_GBRT_AR1_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_GBRT_AR1_delta_alt, cut_hc_results_scaled_GBRT_AR1_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt)), ]

#### Internal Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt), "_simfp", sep = "")

internal_evaluation_GBRT_AR1_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_GBRT_AR1_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt, internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt_ED)
latex_table(internal_evaluation_GBRT_AR1_delta_alt_ED, 5)

plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_AR1_delta_alt$hc_complete_scaled_GBRT_AR1_delta_alt_DTW_33, delta_list_alt_100$GBRT_AR1_delta)
plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_AR1_delta_alt$hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19, delta_list_alt_100$GBRT_AR1_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt$kmeans_5_scaled_GBRT_AR1_delta_alt, delta_list_alt_100$GBRT_AR1_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_GBRT_AR1_delta_alt$kmeans_5_scaled_GBRT_AR1_delta_alt, delta_list_alt_100$GBRT_AR1_delta, 
                     grid_plot = T, output_dir = "kmeans_5_scaled_GBRT_AR1_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_AR1_delta_alt$hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19, delta_list_alt_100$GBRT_AR1_delta, grid_plot = T, output_dir = "hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19")
#### External Evaluation

GBRT_AR1_delta_top3 <- top_partitions(internal_evaluation_GBRT_AR1_delta_alt_ED, 
                                      internal_evaluation_cut_hc_results_scaled_GBRT_AR1_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

GBRT_AR1_delta_top3_partitions <- rbind(kmeans_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_49_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_41_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_4_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_4,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_50,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_49,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_48,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_3)


GBRT_AR1_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(GBRT_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(GBRT_AR1_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(GBRT_AR1_delta_top3_partitions[row, ], 
                                                              GBRT_AR1_delta_top3_partitions[row1, ])
  }
}

GBRT_AR1_purity_matrix <- as.data.frame(GBRT_AR1_purity_matrix)
rownames(GBRT_AR1_purity_matrix) <- best_rows_index
colnames(GBRT_AR1_purity_matrix) <- best_rows_index
latex_table(GBRT_AR1_purity_matrix)
# raw-based clustering ---------------------------------------------------------
GBRT_AR1_delta_alt_raw <- delta_list_alt_100[["GBRT_AR1_delta"]]
euclidean_GBRT_AR1_delta_alt_raw <- dist(GBRT_AR1_delta_alt_raw, "euclidean")
dtw_GBRT_AR1_delta_alt_raw <- dist(GBRT_AR1_delta_alt_raw, "DTW")
elbow(GBRT_AR1_delta_alt_raw, 70)

#### k-means with random centroids initialization
kmeans_results(GBRT_AR1_delta_alt_raw, delta_list_alt_100$GBRT_AR1_delta, 3, 20, NULL)
kmeans_plots_GBRT_AR1_delta_alt_raw$kmeans_4_GBRT_AR1_delta_alt_raw_1 # check
internal_evaluation_kmeans_results_GBRT_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_GBRT_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_results_GBRT_AR1_delta_alt_raw)) {
  
  internal_evaluation_kmeans_results_GBRT_AR1_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_GBRT_AR1_delta_alt_raw, kmeans_results_GBRT_AR1_delta_alt_raw[[cluster]], F)
  
}

#### k-means with simfp initialization
kmeans_results(GBRT_AR1_delta_alt_raw, delta_list_alt_100$GBRT_AR1_delta, 3, 20, "simfp_init")
kmeans_simfp_init_plots_GBRT_AR1_delta_alt_raw$kmeans_4_GBRT_AR1_delta_alt_raw_4 # check
internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw)) {
  
  internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_GBRT_AR1_delta_alt_raw, kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(GBRT_AR1_delta_alt_raw, delta_list_alt_100$GBRT_AR1_delta, 3, 20)
cut_hc_results_GBRT_AR1_delta_alt_raw$hc_ward.D_GBRT_AR1_delta_alt_raw_DTW_4

internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_GBRT_AR1_delta_alt_raw)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw[cluster, ] <- 
      internal_evaluation(euclidean_GBRT_AR1_delta_alt_raw, cut_hc_results_GBRT_AR1_delta_alt_raw[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw[cluster, ] <- 
      internal_evaluation(dtw_GBRT_AR1_delta_alt_raw, cut_hc_results_GBRT_AR1_delta_alt_raw[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw_ED <- internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw)), ]
internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw_DTW <- internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw[grepl("DTW", rownames(internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw_DTW, 5)

# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw), "_simfp", sep = "")

internal_evaluation_GBRT_AR1_delta_alt_raw_ED <- rbind(internal_evaluation_kmeans_results_GBRT_AR1_delta_alt_raw, 
                                                       internal_evaluation_kmeans_simfp_init_results_GBRT_AR1_delta_alt_raw, 
                                                       internal_evaluation_cut_hc_results_GBRT_AR1_delta_alt_raw_ED)
latex_table(internal_evaluation_GBRT_AR1_delta_alt_raw_ED, 15)

plot_assigned_raw_ts(kmeans_results_GBRT_AR1_delta_alt_raw$kmeans_3_GBRT_AR1_delta_alt_raw, delta_list_alt_100$GBRT_AR1_delta)

plot_assigned_raw_ts(cut_hc_results_GBRT_AR1_delta_alt_raw$hc_complete_GBRT_AR1_delta_alt_raw_euclidean_3, 
                     delta_list_alt_100$GBRT_AR1_delta)

## 2.3 GBRT_LIN_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_GBRT_LIN_delta_alt <- delta_list_alt_scaled_features$GBRT_LIN_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_GBRT_LIN_delta_alt <- dist(scaled_GBRT_LIN_delta_alt, "euclidean")
dtw_scaled_GBRT_LIN_delta_alt <- dist(scaled_GBRT_LIN_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_GBRT_LIN_delta_alt, delta_list_alt_100$GBRT_LIN_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_GBRT_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_GBRT_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_GBRT_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_GBRT_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_LIN_delta_alt, kmeans_results_scaled_GBRT_LIN_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_GBRT_LIN_delta_alt, delta_list_alt_100$GBRT_LIN_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_LIN_delta_alt, kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_GBRT_LIN_delta_alt, delta_list_alt_100$GBRT_LIN_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_GBRT_LIN_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_GBRT_LIN_delta_alt, cut_hc_results_scaled_GBRT_LIN_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_GBRT_LIN_delta_alt, cut_hc_results_scaled_GBRT_LIN_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt), "_simfp", sep = "")

internal_evaluation_GBRT_LIN_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_GBRT_LIN_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt, internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt_ED)
latex_table(internal_evaluation_GBRT_LIN_delta_alt_ED, 5)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt$kmeans_3_scaled_GBRT_LIN_delta_alt, delta_list_alt_100$GBRT_LIN_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_GBRT_LIN_delta_alt$kmeans_5_scaled_GBRT_LIN_delta_alt, delta_list_alt_100$GBRT_LIN_delta, 
                     grid_plot = T, output_dir = "kmeans_5_scaled_GBRT_AR1_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_LIN_delta_alt$hc_ward.D2_scaled_GBRT_LIN_delta_alt_DTW_19, delta_list_alt_100$GBRT_LIN_delta, grid_plot = T, output_dir = "hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19")
#### External Evaluation

GBRT_LIN_delta_top3 <- top_partitions(internal_evaluation_GBRT_LIN_delta_alt_ED, 
                                      internal_evaluation_cut_hc_results_scaled_GBRT_LIN_delta_alt_DTW, 3)

## 2.4 GBRT_MIDAS_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_GBRT_MIDAS_delta_alt <- delta_list_alt_scaled_features$GBRT_MIDAS_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_GBRT_MIDAS_delta_alt <- dist(scaled_GBRT_MIDAS_delta_alt, "euclidean")
dtw_scaled_GBRT_MIDAS_delta_alt <- dist(scaled_GBRT_MIDAS_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_GBRT_MIDAS_delta_alt, delta_list_alt_100$GBRT_MIDAS_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_GBRT_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_GBRT_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_GBRT_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_GBRT_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_MIDAS_delta_alt, kmeans_results_scaled_GBRT_MIDAS_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_GBRT_MIDAS_delta_alt, delta_list_alt_100$GBRT_MIDAS_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_MIDAS_delta_alt, kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_GBRT_MIDAS_delta_alt, delta_list_alt_100$GBRT_MIDAS_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_GBRT_MIDAS_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_GBRT_MIDAS_delta_alt, cut_hc_results_scaled_GBRT_MIDAS_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_GBRT_MIDAS_delta_alt, cut_hc_results_scaled_GBRT_MIDAS_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt), "_simfp", sep = "")

internal_evaluation_GBRT_MIDAS_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_GBRT_MIDAS_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_MIDAS_delta_alt, internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt_ED)
latex_table(internal_evaluation_GBRT_MIDAS_delta_alt_ED, 5)

plot_assigned_raw_ts(kmeans_results_scaled_GBRT_MIDAS_delta_alt$kmeans_4_scaled_GBRT_MIDAS_delta_alt, delta_list_alt_100$GBRT_MIDAS_delta)

plot_assigned_raw_ts(kmeans_results_scaled_GBRT_MIDAS_delta_alt$kmeans_39_scaled_GBRT_MIDAS_delta_alt, delta_list_alt_100$GBRT_MIDAS_delta, grid_plot = T, output_dir = "hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19")
#### External Evaluation

GBRT_MIDAS_delta_top3 <- top_partitions(internal_evaluation_GBRT_MIDAS_delta_alt_ED, 
                                        internal_evaluation_cut_hc_results_scaled_GBRT_MIDAS_delta_alt_DTW, 3)


## 2.5 GBRT_NAIVE_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_GBRT_NAIVE_delta_alt <- delta_list_alt_scaled_features$GBRT_NAIVE_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_GBRT_NAIVE_delta_alt <- dist(scaled_GBRT_NAIVE_delta_alt, "euclidean")
dtw_scaled_GBRT_NAIVE_delta_alt <- dist(scaled_GBRT_NAIVE_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_GBRT_NAIVE_delta_alt, delta_list_alt_100$GBRT_NAIVE_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_GBRT_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_GBRT_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_GBRT_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_GBRT_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_NAIVE_delta_alt, kmeans_results_scaled_GBRT_NAIVE_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_GBRT_NAIVE_delta_alt, delta_list_alt_100$GBRT_NAIVE_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_GBRT_NAIVE_delta_alt, kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_GBRT_NAIVE_delta_alt, delta_list_alt_100$GBRT_NAIVE_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_GBRT_NAIVE_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_GBRT_NAIVE_delta_alt, cut_hc_results_scaled_GBRT_NAIVE_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_GBRT_NAIVE_delta_alt, cut_hc_results_scaled_GBRT_NAIVE_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt_DTW, 5)
series_plots_hc_scaled_GBRT_NAIVE_delta_alt$hc_complete_scaled_GBRT_NAIVE_delta_alt_DTW_4_1
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt), "_simfp", sep = "")

internal_evaluation_GBRT_NAIVE_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_GBRT_NAIVE_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_GBRT_NAIVE_delta_alt, internal_evaluation_cut_hc_results_scaled_GBRT_NAIVE_delta_alt_ED)
latex_table(internal_evaluation_GBRT_NAIVE_delta_alt_ED, 5)

plot_assigned_raw_ts(cut_hc_results_scaled_GBRT_NAIVE_delta_alt$hc_ward.D2_scaled_GBRT_NAIVE_delta_alt_euclidean_3, 
                     delta_list_alt_100$GBRT_NAIVE_delta)

## 2.6 ENET_AR1_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_ENET_AR1_delta_alt <- delta_list_alt_scaled_features$ENET_AR1_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_ENET_AR1_delta_alt <- dist(scaled_ENET_AR1_delta_alt, "euclidean")
dtw_scaled_ENET_AR1_delta_alt <- dist(scaled_ENET_AR1_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_ENET_AR1_delta_alt, delta_list_alt_100$ENET_AR1_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_ENET_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_ENET_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_ENET_AR1_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_ENET_AR1_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_AR1_delta_alt, kmeans_results_scaled_ENET_AR1_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_ENET_AR1_delta_alt, delta_list_alt_100$ENET_AR1_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_AR1_delta_alt, kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_ENET_AR1_delta_alt, delta_list_alt_100$ENET_AR1_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_ENET_AR1_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_ENET_AR1_delta_alt, cut_hc_results_scaled_ENET_AR1_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_ENET_AR1_delta_alt, cut_hc_results_scaled_ENET_AR1_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt)), ]

#### Internal Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt), "_simfp", sep = "")

internal_evaluation_ENET_AR1_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_ENET_AR1_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt, internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt_ED)
latex_table(internal_evaluation_ENET_AR1_delta_alt_ED, 5)

series_plots_hc_scaled_ENET_AR1_delta_alt$hc_complete_scaled_ENET_AR1_delta_alt_euclidean_7_7
kmeans_simfp_init_plots_scaled_ENET_AR1_delta_alt$kmeans_10_scaled_ENET_AR1_delta_alt_10

plot_assigned_raw_ts(cut_hc_results_scaled_ENET_AR1_delta_alt$hc_complete_scaled_ENET_AR1_delta_alt_DTW_5, 
                     delta_list_alt_100$ENET_AR1_delta)

plot_assigned_raw_ts(kmeans_results_scaled_ENET_AR1_delta_alt$kmeans_3_scaled_ENET_AR1_delta_alt, delta_list_alt_100$ENET_AR1_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_ENET_AR1_delta_alt$kmeans_5_scaled_ENET_AR1_delta_alt, delta_list_alt_100$ENET_AR1_delta, 
                     grid_plot = T, output_dir = "kmeans_5_scaled_GBRT_AR1_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_ENET_AR1_delta_alt$hc_ward.D2_scaled_ENET_AR1_delta_alt_DTW_19, delta_list_alt_100$ENET_AR1_delta, grid_plot = T, output_dir = "hc_ward.D2_scaled_GBRT_AR1_delta_alt_DTW_19")

#### External Evaluation

ENET_AR1_delta_top3 <- top_partitions(internal_evaluation_ENET_AR1_delta_alt_ED, 
                                      internal_evaluation_cut_hc_results_scaled_ENET_AR1_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

ENET_AR1_delta_top3_partitions <- rbind(kmeans_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_49_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_41_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_4_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_4,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_50,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_49,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_48,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_3)


ENET_AR1_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(ENET_AR1_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(ENET_AR1_delta_top3_partitions[row, ], 
                                                              ENET_AR1_delta_top3_partitions[row1, ])
  }
}

ENET_AR1_purity_matrix <- as.data.frame(ENET_AR1_purity_matrix)
rownames(ENET_AR1_purity_matrix) <- best_rows_index
colnames(ENET_AR1_purity_matrix) <- best_rows_index
latex_table(ENET_AR1_purity_matrix)

# raw-based clustering ---------------------------------------------------------

ENET_AR1_delta_alt_raw <- delta_list_alt_100[["ENET_AR1_delta"]]
euclidean_ENET_AR1_delta_alt_raw <- dist(ENET_AR1_delta_alt_raw, "euclidean")
dtw_ENET_AR1_delta_alt_raw <- dist(ENET_AR1_delta_alt_raw, "DTW")
elbow(ENET_AR1_delta_alt_raw, 70)

# k-means with random centroids initialization
kmeans_results(ENET_AR1_delta_alt_raw, delta_list_alt_100$ENET_AR1_delta, 3, 20, NULL)
kmeans_plots_ENET_AR1_delta_alt_raw$kmeans_4_ENET_AR1_delta_alt_raw_1 # check
internal_evaluation_kmeans_results_ENET_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_ENET_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_results_ENET_AR1_delta_alt_raw)) {
  
  internal_evaluation_kmeans_results_ENET_AR1_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_ENET_AR1_delta_alt_raw, kmeans_results_ENET_AR1_delta_alt_raw[[cluster]], F)
  
}

# k-means with simfp initialization
kmeans_results(ENET_AR1_delta_alt_raw, delta_list_alt_100$ENET_AR1_delta, 3, 20, "simfp_init")
kmeans_simfp_init_plots_ENET_AR1_delta_alt_raw$kmeans_4_ENET_AR1_delta_alt_raw_4 # check
internal_evaluation_kmeans_simfp_init_results_ENET_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_ENET_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(kmeans_simfp_init_results_ENET_AR1_delta_alt_raw)) {
  
  internal_evaluation_kmeans_simfp_init_results_ENET_AR1_delta_alt_raw[cluster, ] <- 
    internal_evaluation(euclidean_ENET_AR1_delta_alt_raw, kmeans_simfp_init_results_ENET_AR1_delta_alt_raw[[cluster]], F)
  
}

kmeans_simfp_init_plots_ENET_AR1_delta_alt_raw$kmeans_6_ENET_AR1_delta_raw_6

# hierarchical
hc_results_ebl(ENET_AR1_delta_alt_raw, delta_list_alt_100$ENET_AR1_delta, 3, 20)
cut_hc_results_ENET_AR1_delta_alt_raw$hc_ward.D_ENET_AR1_delta_alt_raw_DTW_4

internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_ENET_AR1_delta_alt_raw)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw[cluster, ] <- 
      internal_evaluation(euclidean_ENET_AR1_delta_alt_raw, cut_hc_results_ENET_AR1_delta_alt_raw[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw[cluster, ] <- 
      internal_evaluation(dtw_ENET_AR1_delta_alt_raw, cut_hc_results_ENET_AR1_delta_alt_raw[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw_ED <- internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw)), ]
internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw_DTW <- internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw[grepl("DTW", rownames(internal_evaluation_cut_hc_results_ENET_AR1_delta_alt_raw)), ]


## 2.7 ENET_LIN_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_ENET_LIN_delta_alt <- delta_list_alt_scaled_features$ENET_LIN_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_ENET_LIN_delta_alt <- dist(scaled_ENET_LIN_delta_alt, "euclidean")
dtw_scaled_ENET_LIN_delta_alt <- dist(scaled_ENET_LIN_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_ENET_LIN_delta_alt, delta_list_alt_100$ENET_LIN_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_ENET_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_ENET_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_ENET_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_ENET_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_LIN_delta_alt, kmeans_results_scaled_ENET_LIN_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_ENET_LIN_delta_alt, delta_list_alt_100$ENET_LIN_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_LIN_delta_alt, kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_ENET_LIN_delta_alt, delta_list_alt_100$ENET_LIN_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_ENET_LIN_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_ENET_LIN_delta_alt, cut_hc_results_scaled_ENET_LIN_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_ENET_LIN_delta_alt, cut_hc_results_scaled_ENET_LIN_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt), "_simfp", sep = "")

internal_evaluation_ENET_LIN_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_ENET_LIN_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt, internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt_ED)
latex_table(internal_evaluation_ENET_LIN_delta_alt_ED, 5)

series_plots_hc_scaled_ENET_LIN_delta_alt$hc_ward.D_scaled_ENET_LIN_delta_alt_euclidean_3_2

plot_assigned_raw_ts(cut_hc_results_scaled_ENET_LIN_delta_alt$hc_complete_scaled_ENET_LIN_delta_alt_euclidean_50, 
                     delta_list_alt_100$ENET_LIN_delta)

plot_assigned_raw_ts(kmeans_results_scaled_ENET_LIN_delta_alt$kmeans_3_scaled_ENET_LIN_delta_alt, delta_list_alt_100$ENET_LIN_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_ENET_LIN_delta_alt$kmeans_5_scaled_ENET_LIN_delta_alt, delta_list_alt_100$ENET_LIN_delta, 
                     grid_plot = T, output_dir = "kmeans_5_scaled_GBRT_AR1_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_ENET_LIN_delta_alt$hc_complete_scaled_ENET_LIN_delta_alt_euclidean_50, delta_list_alt_100$ENET_LIN_delta, 
                     grid_plot = T, 
                     output_dir = "hc_complete_scaled_ENET_LIN_delta_alt_euclidean_50")

#### External Evaluation

ENET_LIN_delta_top3 <- top_partitions(internal_evaluation_ENET_LIN_delta_alt_ED, 
                                      internal_evaluation_cut_hc_results_scaled_ENET_LIN_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

ENET_LIN_delta_top3_partitions <- rbind(kmeans_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_49_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_41_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                        kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_4_scaled_GBRT_ENET_delta_alt,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_DTW_3,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_4,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_50,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_49,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_48,
                                        cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_3)


ENET_LIN_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(ENET_LIN_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(ENET_LIN_delta_top3_partitions[row, ], 
                                                              ENET_LIN_delta_top3_partitions[row1, ])
  }
}

ENET_LIN_purity_matrix <- as.data.frame(ENET_LIN_purity_matrix)
rownames(ENET_LIN_purity_matrix) <- best_rows_index
colnames(ENET_LIN_purity_matrix) <- best_rows_index
latex_table(ENET_LIN_purity_matrix)
## 2.8 ENET_MIDAS_delta --------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_ENET_MIDAS_delta_alt <- delta_list_alt_scaled_features$ENET_MIDAS_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_ENET_MIDAS_delta_alt <- dist(scaled_ENET_MIDAS_delta_alt, "euclidean")
dtw_scaled_ENET_MIDAS_delta_alt <- dist(scaled_ENET_MIDAS_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_ENET_MIDAS_delta_alt, delta_list_alt_100$ENET_MIDAS_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_ENET_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_ENET_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_ENET_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_ENET_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_MIDAS_delta_alt, kmeans_results_scaled_ENET_MIDAS_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_ENET_MIDAS_delta_alt, delta_list_alt_100$ENET_MIDAS_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_MIDAS_delta_alt, kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_ENET_MIDAS_delta_alt, delta_list_alt_100$ENET_MIDAS_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_ENET_MIDAS_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_ENET_MIDAS_delta_alt, cut_hc_results_scaled_ENET_MIDAS_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_ENET_MIDAS_delta_alt, cut_hc_results_scaled_ENET_MIDAS_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt), "_simfp", sep = "")

internal_evaluation_ENET_MIDAS_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_ENET_MIDAS_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt, internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt_ED)
latex_table(internal_evaluation_ENET_MIDAS_delta_alt_ED, 5)

series_plots_hc_scaled_ENET_MIDAS_delta_alt$hc_ward.D2_scaled_ENET_MIDAS_delta_alt_euclidean_12_12
kmeans_simfp_init_plots_scaled_ENET_MIDAS_delta_alt$kmeans_7_scaled_ENET_MIDAS_delta_alt_7


plot_assigned_raw_ts(cut_hc_results_scaled_ENET_MIDAS_delta_alt$hc_complete_scaled_ENET_MIDAS_delta_alt_euclidean_37, 
                     delta_list_alt_100$ENET_MIDAS_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt$kmeans_3_scaled_ENET_MIDAS_delta_alt, delta_list_alt_100$ENET_MIDAS_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_ENET_MIDAS_delta_alt$kmeans_3_scaled_ENET_MIDAS_delta_alt, delta_list_alt_100$ENET_MIDAS_delta, 
                     grid_plot = T, output_dir = "kmeans_3_scaled_ENET_MIDAS_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_ENET_MIDAS_delta_alt$hc_complete_scaled_ENET_MIDAS_delta_alt_euclidean_37, 
                     delta_list_alt_100$ENET_MIDAS_delta, grid_plot = T, 
                     output_dir = "hc_complete_scaled_ENET_MIDAS_delta_alt_euclidean_37")

#### External Evaluation

ENET_MIDAS_delta_top3 <- top_partitions(internal_evaluation_ENET_MIDAS_delta_alt_ED, 
                                        internal_evaluation_cut_hc_results_scaled_ENET_MIDAS_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

ENET_MIDAS_delta_top3_partitions <- rbind(kmeans_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_euclidean_3,
                                          kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_48_scaled_GBRT_ENET_delta_alt,
                                          kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_49_scaled_GBRT_ENET_delta_alt,
                                          kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_41_scaled_GBRT_ENET_delta_alt,
                                          kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_3_scaled_GBRT_ENET_delta_alt,
                                          kmeans_simfp_init_results_scaled_GBRT_ENET_delta_alt$kmeans_4_scaled_GBRT_ENET_delta_alt,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_3,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D_scaled_GBRT_ENET_delta_alt_DTW_3,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_ward.D2_scaled_GBRT_ENET_delta_alt_DTW_4,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_50,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_49,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_48,
                                          cut_hc_results_scaled_GBRT_ENET_delta_alt$hc_complete_scaled_GBRT_ENET_delta_alt_DTW_3)


ENET_MIDAS_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(ENET_MIDAS_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(ENET_MIDAS_delta_top3_partitions[row, ], 
                                                              ENET_MIDAS_delta_top3_partitions[row1, ])
  }
}

ENET_MIDAS_purity_matrix <- as.data.frame(ENET_MIDAS_purity_matrix)
rownames(ENET_MIDAS_purity_matrix) <- best_rows_index
colnames(ENET_MIDAS_purity_matrix) <- best_rows_index
latex_table(ENET_MIDAS_purity_matrix)

## 2.9 ENET_NAIVE_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_ENET_NAIVE_delta_alt <- delta_list_alt_scaled_features$ENET_NAIVE_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_ENET_NAIVE_delta_alt <- dist(scaled_ENET_NAIVE_delta_alt, "euclidean")
dtw_scaled_ENET_NAIVE_delta_alt <- dist(scaled_ENET_NAIVE_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_ENET_NAIVE_delta_alt, delta_list_alt_100$ENET_NAIVE_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_ENET_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_ENET_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_ENET_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_ENET_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_NAIVE_delta_alt, kmeans_results_scaled_ENET_NAIVE_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_ENET_NAIVE_delta_alt, delta_list_alt_100$ENET_NAIVE_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_ENET_NAIVE_delta_alt, kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_ENET_NAIVE_delta_alt, delta_list_alt_100$ENET_NAIVE_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_ENET_NAIVE_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_ENET_NAIVE_delta_alt, cut_hc_results_scaled_ENET_NAIVE_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_ENET_NAIVE_delta_alt, cut_hc_results_scaled_ENET_NAIVE_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt_DTW, 15)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt), "_simfp", sep = "")

internal_evaluation_ENET_NAIVE_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_ENET_NAIVE_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_ENET_NAIVE_delta_alt, internal_evaluation_cut_hc_results_scaled_ENET_NAIVE_delta_alt_ED)
latex_table(internal_evaluation_ENET_NAIVE_delta_alt_ED, 15)

## 2.10 AR1_LIN_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_AR1_LIN_delta_alt <- delta_list_alt_scaled_features$AR1_LIN_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_AR1_LIN_delta_alt <- dist(scaled_AR1_LIN_delta_alt, "euclidean")
dtw_scaled_AR1_LIN_delta_alt <- dist(scaled_AR1_LIN_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_AR1_LIN_delta_alt, delta_list_alt_100$AR1_LIN_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_AR1_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_AR1_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_AR1_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_AR1_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_LIN_delta_alt, kmeans_results_scaled_AR1_LIN_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_AR1_LIN_delta_alt, delta_list_alt_100$AR1_LIN_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_LIN_delta_alt, kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_AR1_LIN_delta_alt, delta_list_alt_100$AR1_LIN_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_AR1_LIN_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_AR1_LIN_delta_alt, cut_hc_results_scaled_AR1_LIN_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_AR1_LIN_delta_alt, cut_hc_results_scaled_AR1_LIN_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt), "_simfp", sep = "")

internal_evaluation_AR1_LIN_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_AR1_LIN_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt, internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt_ED)
latex_table(internal_evaluation_AR1_LIN_delta_alt_ED, 5)

series_plots_hc_scaled_AR1_LIN_delta_alt$hc_ward.D2_scaled_AR1_LIN_delta_alt_euclidean_7_1

plot_assigned_raw_ts(cut_hc_results_scaled_AR1_LIN_delta_alt$hc_ward.D2_scaled_AR1_LIN_delta_alt_DTW_48, 
                     delta_list_alt_100$AR1_LIN_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt$kmeans_49_scaled_AR1_LIN_delta_alt, delta_list_alt_100$AR1_LIN_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_AR1_LIN_delta_alt$kmeans_49_scaled_AR1_LIN_delta_alt, delta_list_alt_100$AR1_LIN_delta, 
                     grid_plot = T, output_dir = "kmeans_49_scaled_AR1_LIN_delta_alt_simfp")
plot_assigned_raw_ts(cut_hc_results_scaled_AR1_LIN_delta_alt$hc_ward.D2_scaled_AR1_LIN_delta_alt_DTW_48, 
                     delta_list_alt_100$AR1_LIN_delta, grid_plot = T, 
                     output_dir = "hc_ward.D2_scaled_AR1_LIN_delta_alt_DTW_48")

#### External Evaluation

AR1_LIN_delta_top3 <- top_partitions(internal_evaluation_AR1_LIN_delta_alt_ED, 
                                     internal_evaluation_cut_hc_results_scaled_AR1_LIN_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

AR1_LIN_delta_top3_partitions <- rbind()


AR1_LIN_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(AR1_LIN_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(AR1_LIN_delta_top3_partitions[row, ], 
                                                              AR1_LIN_delta_top3_partitions[row1, ])
  }
}

AR1_LIN_purity_matrix <- as.data.frame(AR1_LIN_purity_matrix)
rownames(AR1_LIN_purity_matrix) <- best_rows_index
colnames(AR1_LIN_purity_matrix) <- best_rows_index
latex_table(AR1_LIN_purity_matrix)


## 2.11 AR1_MIDAS_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_AR1_MIDAS_delta_alt <- delta_list_alt_scaled_features$AR1_MIDAS_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_AR1_MIDAS_delta_alt <- dist(scaled_AR1_MIDAS_delta_alt, "euclidean")
dtw_scaled_AR1_MIDAS_delta_alt <- dist(scaled_AR1_MIDAS_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_AR1_MIDAS_delta_alt, delta_list_alt_100$AR1_MIDAS_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_AR1_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_AR1_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_AR1_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_AR1_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_MIDAS_delta_alt, kmeans_results_scaled_AR1_MIDAS_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_AR1_MIDAS_delta_alt, delta_list_alt_100$AR1_MIDAS_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_MIDAS_delta_alt, kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_AR1_MIDAS_delta_alt, delta_list_alt_100$AR1_MIDAS_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_AR1_MIDAS_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_AR1_MIDAS_delta_alt, cut_hc_results_scaled_AR1_MIDAS_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_AR1_MIDAS_delta_alt, cut_hc_results_scaled_AR1_MIDAS_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt), "_simfp", sep = "")

internal_evaluation_AR1_MIDAS_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_AR1_MIDAS_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt, internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt_ED)
latex_table(internal_evaluation_AR1_MIDAS_delta_alt_ED, 5)

series_plots_hc_scaled_AR1_MIDAS_delta_alt$hc_ward.D_scaled_AR1_MIDAS_delta_alt_DTW_50_50

plot_assigned_raw_ts(cut_hc_results_scaled_AR1_MIDAS_delta_alt$hc_complete_scaled_AR1_MIDAS_delta_alt_euclidean_50, 
                     delta_list_alt_100$AR1_MIDAS_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt$kmeans_49_scaled_AR1_MIDAS_delta_alt, delta_list_alt_100$AR1_MIDAS_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_AR1_MIDAS_delta_alt$kmeans_3_scaled_AR1_MIDAS_delta_alt, delta_list_alt_100$AR1_MIDAS_delta, 
                     grid_plot = T, output_dir = "kmeans_3_scaled_AR1_MIDAS_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_AR1_MIDAS_delta_alt$hc_ward.D2_scaled_AR1_MIDAS_delta_alt_DTW_49, 
                     delta_list_alt_100$AR1_MIDAS_delta, grid_plot = T, 
                     output_dir = "hc_ward.D2_scaled_AR1_MIDAS_delta_alt_DTW_49")

#### External Evaluation

AR1_MIDAS_delta_top3 <- top_partitions(internal_evaluation_AR1_MIDAS_delta_alt_ED, 
                                       internal_evaluation_cut_hc_results_scaled_AR1_MIDAS_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

AR1_MIDAS_delta_top3_partitions <- rbind()


AR1_MIDAS_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(AR1_MIDAS_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(AR1_MIDAS_delta_top3_partitions[row, ], 
                                                              AR1_MIDAS_delta_top3_partitions[row1, ])
  }
}

AR1_MIDAS_purity_matrix <- as.data.frame(AR1_MIDAS_purity_matrix)
rownames(AR1_MIDAS_purity_matrix) <- best_rows_index
colnames(AR1_MIDAS_purity_matrix) <- best_rows_index
latex_table(AR1_MIDAS_purity_matrix)

## 2.12 AR1_NAIVE_delta ---------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_AR1_NAIVE_delta_alt <- delta_list_alt_scaled_features$AR1_NAIVE_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_AR1_NAIVE_delta_alt <- dist(scaled_AR1_NAIVE_delta_alt, "euclidean")
dtw_scaled_AR1_NAIVE_delta_alt <- dist(scaled_AR1_NAIVE_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_AR1_NAIVE_delta_alt, delta_list_alt_100$AR1_NAIVE_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_AR1_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_AR1_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_AR1_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_AR1_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_NAIVE_delta_alt, kmeans_results_scaled_AR1_NAIVE_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_AR1_NAIVE_delta_alt, delta_list_alt_100$AR1_NAIVE_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_AR1_NAIVE_delta_alt, kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_AR1_NAIVE_delta_alt, delta_list_alt_100$AR1_NAIVE_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_AR1_NAIVE_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_AR1_NAIVE_delta_alt, cut_hc_results_scaled_AR1_NAIVE_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_AR1_NAIVE_delta_alt, cut_hc_results_scaled_AR1_NAIVE_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt), "_simfp", sep = "")

internal_evaluation_AR1_NAIVE_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_AR1_NAIVE_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt, internal_evaluation_cut_hc_results_scaled_AR1_NAIVE_delta_alt_ED)
latex_table(internal_evaluation_AR1_NAIVE_delta_alt_ED, 5)

plot_assigned_raw_ts(cut_hc_results_scaled_AR1_NAIVE_delta_alt$hc_ward.D_scaled_AR1_NAIVE_delta_alt_DTW_20, 
                     delta_list_alt_100$AR1_NAIVE_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_AR1_NAIVE_delta_alt$kmeans_34_scaled_AR1_NAIVE_delta_alt, 
                     delta_list_alt_100$AR1_NAIVE_delta)

## 2.13 LIN_MIDAS_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_LIN_MIDAS_delta_alt <- delta_list_alt_scaled_features$LIN_MIDAS_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_LIN_MIDAS_delta_alt <- dist(scaled_LIN_MIDAS_delta_alt, "euclidean")
dtw_scaled_LIN_MIDAS_delta_alt <- dist(scaled_LIN_MIDAS_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_LIN_MIDAS_delta_alt, delta_list_alt_100$LIN_MIDAS_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_LIN_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_LIN_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_LIN_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_LIN_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_LIN_MIDAS_delta_alt, kmeans_results_scaled_LIN_MIDAS_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_LIN_MIDAS_delta_alt, delta_list_alt_100$LIN_MIDAS_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_LIN_MIDAS_delta_alt, kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_LIN_MIDAS_delta_alt, delta_list_alt_100$LIN_MIDAS_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_LIN_MIDAS_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_LIN_MIDAS_delta_alt, cut_hc_results_scaled_LIN_MIDAS_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_LIN_MIDAS_delta_alt, cut_hc_results_scaled_LIN_MIDAS_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt), "_simfp", sep = "")

internal_evaluation_LIN_MIDAS_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_LIN_MIDAS_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt, internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt_ED)
latex_table(internal_evaluation_LIN_MIDAS_delta_alt_ED, 5)

plot_assigned_raw_ts(cut_hc_results_scaled_LIN_MIDAS_delta_alt$hc_complete_scaled_LIN_MIDAS_delta_alt_euclidean_48, 
                     delta_list_alt_100$LIN_MIDAS_delta)
plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt$kmeans_3_scaled_LIN_MIDAS_delta_alt, delta_list_alt_100$LIN_MIDAS_delta)

plot_assigned_raw_ts(kmeans_simfp_init_results_scaled_LIN_MIDAS_delta_alt$kmeans_3_scaled_LIN_MIDAS_delta_alt, delta_list_alt_100$LIN_MIDAS_delta, 
                     grid_plot = T, output_dir = "kmeans_3_scaled_LIN_MIDAS_delta_alt")
plot_assigned_raw_ts(cut_hc_results_scaled_LIN_MIDAS_delta_alt$hc_complete_scaled_LIN_MIDAS_delta_alt_DTW_50, 
                     delta_list_alt_100$LIN_MIDAS_delta, grid_plot = T, 
                     output_dir = "hc_complete_scaled_LIN_MIDAS_delta_alt_DTW_50")

#### External Evaluation

LIN_MIDAS_delta_top3 <- top_partitions(internal_evaluation_LIN_MIDAS_delta_alt_ED, 
                                       internal_evaluation_cut_hc_results_scaled_LIN_MIDAS_delta_alt_DTW, 3)

best_rows_index <- c("ED 3\\(CHI, ASW)", "ED 3_1\\(CHI, ASW)", "ED 3\\(CHI)", 
                     "ED 48\\(DI)", "ED 49\\(DI)", "ED 41\\(ASW)", "ED 3\\(ASW)", "ED 4\\(ASW)",
                     "DTW 3\\(CHI, ASW)", "DTW 3_1\\(CHI, ASW)", "DTW 4\\(CHI)", 
                     "DTW 50\\(DI)", "DTW 49\\(DI)", "DTW 48\\(DI)", "DTW 3\\(ASW)")

LIN_MIDAS_delta_top3_partitions <- rbind()


LIN_MIDAS_purity_matrix <- matrix(nrow = 15, ncol = 15)
for (row in 1:nrow(ENET_AR_delta_top3_partitions)){
  for (row1 in 1:nrow(LIN_MIDAS_delta_top3_partitions)){
    GBRT_ENET_purity_matrix[row, row1] <- nonsymmetric_purity(LIN_MIDAS_delta_top3_partitions[row, ], 
                                                              LIN_MIDAS_delta_top3_partitions[row1, ])
  }
}

LIN_MIDAS_purity_matrix <- as.data.frame(LIN_MIDAS_purity_matrix)
rownames(LIN_MIDAS_purity_matrix) <- best_rows_index
colnames(LIN_MIDAS_purity_matrix) <- best_rows_index
latex_table(LIN_MIDAS_purity_matrix)

## 2.14 LIN_NAIVE_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_LIN_NAIVE_delta_alt <- delta_list_alt_scaled_features$LIN_NAIVE_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_LIN_NAIVE_delta_alt <- dist(scaled_LIN_NAIVE_delta_alt, "euclidean")
dtw_scaled_LIN_NAIVE_delta_alt <- dist(scaled_LIN_NAIVE_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_LIN_NAIVE_delta_alt, delta_list_alt_100$LIN_NAIVE_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_LIN_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_LIN_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_LIN_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_LIN_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_LIN_NAIVE_delta_alt, kmeans_results_scaled_LIN_NAIVE_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_LIN_NAIVE_delta_alt, delta_list_alt_100$LIN_NAIVE_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_LIN_NAIVE_delta_alt, kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_LIN_NAIVE_delta_alt, delta_list_alt_100$LIN_NAIVE_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_LIN_NAIVE_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_LIN_NAIVE_delta_alt, cut_hc_results_scaled_LIN_NAIVE_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_LIN_NAIVE_delta_alt, cut_hc_results_scaled_LIN_NAIVE_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt), "_simfp", sep = "")

internal_evaluation_LIN_NAIVE_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_LIN_NAIVE_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_LIN_NAIVE_delta_alt, internal_evaluation_cut_hc_results_scaled_LIN_NAIVE_delta_alt_ED)
latex_table(internal_evaluation_LIN_NAIVE_delta_alt_ED, 5)

## 2.15 MIDAS_NAIVE_delta ----------------------------------------------------------

### Feature-Based Clustering

# Extract Features as Objects
scaled_MIDAS_NAIVE_delta_alt <- delta_list_alt_scaled_features$MIDAS_NAIVE_delta
# Compute Distance Matrix for Evaluation
euclidean_scaled_MIDAS_NAIVE_delta_alt <- dist(scaled_MIDAS_NAIVE_delta_alt, "euclidean")
dtw_scaled_MIDAS_NAIVE_delta_alt <- dist(scaled_MIDAS_NAIVE_delta_alt, "DTW")

#### k-means with Random Centroids Initialization
kmeans_results(scaled_MIDAS_NAIVE_delta_alt, delta_list_alt_100$MIDAS_NAIVE_delta, 3, 50, NULL)
internal_evaluation_kmeans_results_scaled_MIDAS_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_results_scaled_MIDAS_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_results_scaled_MIDAS_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_results_scaled_MIDAS_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_MIDAS_NAIVE_delta_alt, kmeans_results_scaled_MIDAS_NAIVE_delta_alt[[cluster]], F)
  
}

#### k-means with Simfp Initialization
kmeans_results(scaled_MIDAS_NAIVE_delta_alt, delta_list_alt_100$MIDAS_NAIVE_delta, 3, 50, "simfp_init")
internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")
for(cluster in names(kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt)) {
  
  internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt[cluster, ] <- 
    internal_evaluation(euclidean_scaled_MIDAS_NAIVE_delta_alt, kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt[[cluster]], F)
  
}

#### Hierarchical
hc_results_ebl(scaled_MIDAS_NAIVE_delta_alt, delta_list_alt_100$MIDAS_NAIVE_delta, 3, 50)
internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt) <- c("Calinski-Harabasz Index", "Dunn Index", "Average Silhouette Width")

for(cluster in names(cut_hc_results_scaled_MIDAS_NAIVE_delta_alt)) {
  
  if(grepl("euclidean", cluster)) {
    internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(euclidean_scaled_MIDAS_NAIVE_delta_alt, cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[[cluster]], F)
  }
  
  if(grepl("DTW", cluster)) {
    internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[cluster, ] <- 
      internal_evaluation(dtw_scaled_MIDAS_NAIVE_delta_alt, cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[[cluster]], F)
  }
}

internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt_ED <- internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[grepl("euclidean", rownames(internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt)), ]
internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt_DTW <- internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt[grepl("DTW", rownames(internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt)), ]

#### Evaluation
# DTW
latex_table(internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt_DTW, 5)
# ED
# change rownames of simfp k-means df so they are different from k-means df
rownames(internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt) <- 
  paste(rownames(internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt), "_simfp", sep = "")

internal_evaluation_MIDAS_NAIVE_delta_alt_ED <- rbind(internal_evaluation_kmeans_results_scaled_MIDAS_NAIVE_delta_alt, internal_evaluation_kmeans_simfp_init_results_scaled_MIDAS_NAIVE_delta_alt, internal_evaluation_cut_hc_results_scaled_MIDAS_NAIVE_delta_alt_ED)
latex_table(internal_evaluation_MIDAS_NAIVE_delta_alt_ED, 5)