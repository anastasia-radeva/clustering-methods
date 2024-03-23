# Dependencies -----------------------------------------------------------------
load_package("parameters")
load_package("longmemo")
load_package("tseries")
load_package("rugarch")
load_package("e1071")
load_package("parallel")

# source("SARFIMA-GARCH_Simulation.R")

# One Time Definitions -------------------------------------------------------
# Volatility clustering
garch_model <- rugarch::ugarchspec(variance.model = list(model = "sGARCH", 
                                                         garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0),
                                                     include.mean = F), 
                                   distribution.model = "ged")
# Nonlinearity in Autoregressive Structure of Volatility
egarch_model <- rugarch::ugarchspec(variance.model = list(model = "eGARCH",
                                                          garchOrder = c(1, 1)), 
                                    mean.model = list(armaOrder = c(0, 0),
                                                      include.mean = F), 
                                    distribution.model = "norm")

# Function which extracts features using the defined functions -----------------
extract_features_old <- function(df, tsa) {
  
  # Dependencies
  source("feature_extractions.R")
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  
  # create data frame to store features in with the correct columns
  df_features <- data.frame()
  
  # Compute features for each time series (row) in the data frame --------------
  for(i in 1:nrow(df)) {
    
    # Mean
    df_features[i, "Mean"] <- mean(df[i, ])
    
    # Volatility
    df_features[i, "SD"] <- sd(df[i, ])
    
    # Skewness
    df_features[i, "Skew"] <- e1071::skewness(df[i, ], na.rm = F, type = "1")
    
    # Kurtosis
    df_features[i, "Kurtosis"] <- e1071::kurtosis(df[i, ], na.rm = F, type = "1")# returns excess kurtosis
    
    # Periodicity
    df_features[i, "Periodicity"] <- extract_frequency(df[i, ])
    
    # # Trend
    # df_features[i, "Trend"] <- 1 - sd(df[i, ])/sd(detrend(df[i, ]))
    # 
    # # Seasonality
    # df_features[i, "Seasonality"] <- 1 - sd(df[i, ])/sd(deseason(df[i, ]))
    
    # Autocorrelation
    df_features[i, "Autocorrelation"] <- extract_autocorrelation(df[i, ])
    
    # Long-range Dependence
    df_features[i, "Long_range_Dependence"] <- extract_long_range_dependence(df[i, ])
    
    # Non-linearity in Autoregressive Structure of Mean
    df_features[i, "Nonlinearity_Mean"] <- extract_nonlinearity_mean(df[i, ])
    
    # Volatility Clustering
    df_features[i, "Volatility_Clustering"] <- extract_vola_clustering(df[i, ], garch_model)
   
    # Non-linearity in Autoregressive Structure of Volatility: Leverage effect
    df_features[i, "Nonlinearity_Volatility"] <- extract_nonlinearity_vola(df[i, ], egarch_model)
    
    
    if (tsa == T){
      tsa_df <- ts_adjust(df)
      
      # Mean
      df_features[i, "TSA_Mean"] <- mean(tsa_df[i, ])
      
      # Volatility
      df_features[i, "TSA_SD"] <- sd(tsa_df[i, ])
      
      # Skewness
      df_features[i, "TSA_Skew"] <- e1071::skewness(tsa_df[i, ], na.rm = F, type = "1")
      
      # Kurtosis
      df_features[i, "TSA_Kurtosis"] <- e1071::kurtosis(tsa_df[i, ], na.rm = F, type = "1")# returns excess kurtosis
      
      # Trend
      df_features[i, "TSA_Trend"] <- 1 - sd(tsa_df[i, ])/sd(deseason(df[i, ]))
      
      # Seasonality
      df_features[i, "TSA_Seasonality"] <- 1 - sd(tsa_df[i, ])/sd(detrend(df[i, ]))
      
      # Autocorrelation
      df_features[i, "TSA_Autocorrelation"] <- extract_autocorrelation(tsa_df[i, ])
      
      # Long-range Dependence
      df_features[i, "TSA_Long_range_Dependence"] <- extract_long_range_dependence(tsa_df[i, ])
      
      # Non-linearity in Autoregressive Structure of Mean
      df_features[i, "TSA_Nonlinearity_Mean"] <- extract_nonlinearity_mean(tsa_df[i, ])
      
      # Volatility Clustering
      df_features[i, "TSA_Volatility_Clustering"] <- extract_vola_clustering(tsa_df[i, ], garch_model)
      
      # Non-linearity in Autoregressive Structure of Volatility: Leverage effect
      df_features[i, "TSA_Nonlinearity_Volatility"] <- extract_nonlinearity_vola(tsa_df[i, ], egarch_model)
      
    }
    
  }
  
  # set row names to name of series out of original df
  rownames(df_features) <- row.names(df)
  
  # assign the new data frame to the global environment
  df_features_name <- paste(df_name, "_features", sep = "")
  assign(df_features_name, df_features, envir = .GlobalEnv)
  print(df_features_name)
  return(df_features)
}

extract_features <- function(df, tsa) {
  
  # Dependencies
  source("feature_extractions.R")
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  
  # create data frame to store features in with the correct columns
  df_features <- data.frame()
  
  # Compute features for each time series (row) in the data frame --------------
  df_features <- do.call(rbind, mclapply(1:nrow(df), function(i) {
    features <- data.frame(
      Mean = mean(df[i, ]),
      SD = sd(df[i, ]),
      Skew = e1071::skewness(df[i, ], na.rm = F, type = "1"),
      Kurtosis = e1071::kurtosis(df[i, ], na.rm = F, type = "1"), # returns excess kurtosis
      Periodicity = extract_frequency(df[i, ]),
      Autocorrelation = extract_autocorrelation(df[i, ]),
      Long_range_Dependence = extract_long_range_dependence(df[i, ]),
      Nonlinearity_Mean = extract_nonlinearity_mean(df[i, ]),
      Volatility_Clustering = extract_vola_clustering(df[i, ], garch_model),
      Nonlinearity_Volatility = extract_nonlinearity_vola(df[i, ], egarch_model)
    )
    
    if (tsa == T){
      tsa_df <- ts_adjust(df)
      
      features <- cbind(features, data.frame(
        TSA_Mean = mean(tsa_df[i, ]),
        TSA_SD = sd(tsa_df[i, ]),
        TSA_Skew = e1071::skewness(tsa_df[i, ], na.rm = F, type = "1"),
        TSA_Kurtosis = e1071::kurtosis(tsa_df[i, ], na.rm = F, type = "1"), # returns excess kurtosis
        TSA_Trend = 1 - sd(tsa_df[i, ])/sd(deseason(df[i, ])),
        TSA_Seasonality = 1 - sd(tsa_df[i, ])/sd(detrend(df[i, ])),
        TSA_Autocorrelation = extract_autocorrelation(tsa_df[i, ]),
        TSA_Long_range_Dependence = extract_long_range_dependence(tsa_df[i, ]),
        TSA_Nonlinearity_Mean = extract_nonlinearity_mean(tsa_df[i, ]),
        TSA_Volatility_Clustering = extract_vola_clustering(tsa_df[i, ], garch_model),
        TSA_Nonlinearity_Volatility = extract_nonlinearity_vola(tsa_df[i, ], egarch_model)
      ))
    }
    
    return(features)
  }, mc.cores = detectCores()))
  
  # set row names to name of series out of original df
  rownames(df_features) <- row.names(df)
  
  # assign the new data frame to the global environment
  df_features_name <- paste(df_name, "_features", sep = "")
  assign(df_features_name, df_features, envir = .GlobalEnv)
  print(df_features_name)
  return(df_features)
}
# Example:
# extract_features(df_simulated, tsa = T)
# rescale_features(df_simulated_features, tsa = T)
