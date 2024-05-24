# In the following, functions that extract features from the raw time series data are defined.

# Dependencies -----------------------------------------------------------------
library("mgcv")
library("stats")
library("splines")
library("stats")
library("longmemo")
library("tseries")
library("stats")
library("TSA")
library("tsDyn")
library("rugarch")
library("scales")
library("parameters")
library("e1071")
library("parallel")
library("datawizard")

# Feature Extraction Functions -------------------------------------------------

# Function to adjust data for trend and seasonality
ts_adjust <- function(df) {
  
  # create data frame to store features in with the correct columns
  # Initialize tsa_df as a matrix with the same number of rows as df
  tsa_df <- matrix(nrow = nrow(df), ncol = ncol(df))
  
  for(i in 1:nrow(df)) {
    
    # Detrend
    detrended <- detrend(df[i, ])
    diff <- detrended
    
    # Remove Seasonal Component from detrended series
    # not with deseason function, because it will extract period from detrended series and not original
    period <- extract_frequency(df[i, ])
    
    if (period > 0 && length(df[i, ]) > period) {
      diff <- detrended[1:period]
      for (j in (period + 1):length(df[i, ])) {
        value <- detrended[j] - detrended[j - period]
        diff <- append(diff, value)
      }
    }
    deseasonalized <- diff
    
    tsa_df[i, ] = deseasonalized
  }
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  # set row names to name of series out of original df
  rownames(tsa_df) <- row.names(df)
  
  # assign the new data frame to the global environment
  tsa_df_name <- paste("tsa_", df_name, sep = "")
  assign(tsa_df_name, tsa_df, envir = .GlobalEnv)
  
  return(tsa_df)
}

# Detrend with regression spline
detrend <- function(ts){
  
  tt <- 1:length(ts)
  trend <- rep(NA,length(ts))
  trend[!is.na(ts)] <- fitted(gam(ts ~ s(tt)))
  detrended_ts <- ts-trend
  return(detrended_ts)
}



# Function to detrend time series and find frequency
extract_frequency <- function(ts, threshold = 0.2, acf_threshold = 0.6) {
  
  # Detrend time series using a regression spline
  detrended <- detrend(ts)
  
  # Autocorrelation function for all lags up to 2/3 of series length to look for two peaks
  acf_values <- acf(detrended, lag.max = length(ts)*2 / 3, plot = FALSE)$acf
  
  # Find peaks and troughs in autocorrelation function
  diff_acf <- diff(acf_values)
  peaks <- which(diff_acf[-length(diff_acf)] >= 0 & diff_acf[-1] <= 0) + 1
  troughs <- which(diff_acf[-length(diff_acf)] <= 0 & diff_acf[-1] >= 0) + 1
  
  # Frequency is the first peak with conditions
  for (peak in peaks) {
    if (peak > min(troughs) && (acf_values[peak] - min(acf_values[troughs])) >= threshold && abs(acf_values[peak]) > acf_threshold) {
      # Check if there is another peak around period*2
      period <- peak
      next_peak <- period * 2
      buffer <- length(ts) / 50
      if (any(peaks >= next_peak - buffer & peaks <= next_peak + buffer)) {
        return(period)
      }
    }
  }
  
  # If no such peak is found, frequency is set to 0 (equivalent to non-seasonal)
  return(0)
}

# Function which removes Seasonal Component
# Differencing the Periodical Series
deseason<- function(ts){
  
  period <- extract_frequency(ts)
  diff = ts[1:period]
  
  if (period > 0 && length(ts) > period) {
    for (i in (period + 1):length(ts)) {
      value = ts[i] - ts[i - period]
      diff = append(diff, value)
    }
  }
  else diff = ts
  
  return(diff)
}

# Function which extracts autocorrelation
# represented by Ljung-Box Test Statistic
# we use a test statistic as quantification because a higher test statistic 
# often indicates stronger evidence against the null hypothesis
extract_autocorrelation <- function(ts){
  
  lb_test <-  Box.test(ts, lag = 2*sqrt(length(ts)), type = "Ljung-Box")
  lb_statistic <- lb_test$statistic # rule of thumb: lag = 2sqrt(T)
  lb_significance <- as.numeric(lb_test$p.value < 0.05)
  # normalize the test statistic by the number of observations and the number of lags
  # more interpretable and comparable across different datasets or models 
  lb_statistic_normalized <- lb_statistic*lb_significance*10/(length(ts)*2*sqrt(length(ts)))
  return(lb_statistic_normalized)
  
}

# Function which extracts Long-range Dependence
# represented by an estimate of the Hurst parameter using the Whittle estimator (Stroe-Kunolds, 2009)
extract_long_range_dependence <- function(ts){
  
  #assume fGN for model
  hurst_par <- WhittleEst(ts, model = "fGn")$coefficients[1]
  return(hurst_par)
  
}

# Function which extracts non-linearity in Autoregressive Structure of Mean
# represented by Teraesvirta Neural Network Test for Nonlinearity - Test Statistic
# we use a test statistic as quantification because a higher test statistic 
# often indicates stronger evidence against the null hypothesis
extract_nonlinearity_mean <- function(ts){
  
  trv_test <- terasvirta.test(ts(ts), lag = 1, type = "F") # using conventional lag
  trv_statistic <- trv_test$statistic
  trv_significance <- as.numeric(trv_test$p.value < 0.05)
  # normalize the test statistic by the number of observations and the number of lags
  # more interpretable and comparable across different datasets or models
  trv_statistic_normalized <- trv_statistic*trv_significance*10/(length(ts))
  return(unname(trv_statistic_normalized))
  
}

# [-1, 1] Scaling --------------------------------------------------------------
rescale_features <- function(df, tsa = T){
  
  # Define custom scaling function
  # rescale_1 <- function(x) 
  rescale_2 <- function(x) {tanh(x/3)} # {2 * ((x - min(x)) / (max(x) - min(x))) - 1}
  rescale_3 <- function(x) {exp(-(ncol(df)/x)*0.01)}
  
  # Create a copy of the original data frame
  scaled_df <- df
  
  if (tsa == T){
    selected_cols_1 <- c("SD", "TSA_SD", "Nonlinearity_Mean", "TSA_Nonlinearity_Mean")
    selected_cols_2 <- c("Mean", "Skew", "Kurtosis", 
                         "TSA_Mean","TSA_Skew", "TSA_Kurtosis")
  }
  else{
    # Select only the "Mean", "SD", "Skew", and "Kurtosis" columns
    selected_cols_1 <- c("SD", "Nonlinearity_Mean")
    selected_cols_2 <- c("Mean","Skew", "Kurtosis")
  }
  
  selected_cols_3 <- c("Periodicity")
  
  # Apply the rescale function to the selected columns
  scaled_df[selected_cols_1] <- apply(df[selected_cols_1], 2, scales::rescale)
  scaled_df[selected_cols_2] <- apply(df[selected_cols_2], 2, rescale_2)
  scaled_df[selected_cols_3] <- apply(df[selected_cols_3], 2, rescale_3)
  
  # Get the name of the input data frame
  df_name <- deparse(substitute(df))
  
  # Set row names to name of series out of original df
  rownames(scaled_df) <- row.names(df)
  
  # Assign the new data frame to the global environment
  scaled_df_name <- paste("scaled_", df_name, sep = "")
  assign(scaled_df_name, scaled_df, envir = .GlobalEnv)
  
  return(scaled_df)
  
}
# Function which extracts features data frame using the defined functions ------
extract_features <- function(df, tsa) {
  
  # Dependencies
  #source("feature_extractions.R")
  
  # get the name of the input data frame
  df_name <- deparse(substitute(df))
  
  df <- as.matrix(df)
  # create data frame to store features in with the correct columns
  df_features <- data.frame()
  
  # Compute features for each time series (row) in the data frame
  df_features <- do.call(rbind, mclapply(1:nrow(df), function(i) {
    features <- data.frame(
      Mean = mean(df[i, ]),
      SD = sd(df[i, ]),
      Skew = as.numeric(timeDate::skewness(df[i, ], method = "moment")),
      Kurtosis = as.numeric(timeDate::kurtosis(df[i, ], method = "fisher")), # returns excess kurtosis
      Periodicity = extract_frequency(df[i, ]),
      Autocorrelation = extract_autocorrelation(df[i, ]),
      Long_range_Dependence = extract_long_range_dependence(df[i, ]),
      Nonlinearity_Mean = extract_nonlinearity_mean(df[i, ])
    )
    
    if (tsa == T){
      tsa_df <- ts_adjust(df)
      
      features <- cbind(features, data.frame(
        TSA_Mean = mean(tsa_df[i, ]),
        TSA_SD = sd(tsa_df[i, ]),
        TSA_Skew = as.numeric(timeDate::skewness(tsa_df[i, ], method = "moment")),
        TSA_Kurtosis = as.numeric(timeDate::kurtosis(tsa_df[i, ], method = "fisher")), # returns excess kurtosis
        TSA_Trend = 1 - sd(tsa_df[i, ])/sd(deseason(df[i, ])),
        TSA_Seasonality = 1 - sd(tsa_df[i, ])/sd(detrend(df[i, ])),
        TSA_Autocorrelation = extract_autocorrelation(tsa_df[i, ]),
        TSA_Long_range_Dependence = extract_long_range_dependence(tsa_df[i, ]),
        TSA_Nonlinearity_Mean = extract_nonlinearity_mean(tsa_df[i, ])
        
      ))
    }
    
    return(features)
  }, mc.cores = detectCores() - 1))
  
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