# In the following, functions to extract features from the raw time series data 
# are defined.

# ------------------------------------------------------------------------------
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

# Detrend with regression spline -----------------------------------------------
detrend <- function(ts){
    
  load_package("mgcv")
  load_package("stats")
    
  tt <- 1:length(ts)
  trend <- rep(NA,length(ts))
  trend[!is.na(ts)] <- fitted(gam(ts ~ s(tt)))
  detrended_ts <- ts-trend
    
  return(detrended_ts)
    
}
  
# Function to detrend time series and find frequency ---------------------------
extract_frequency <- function(ts, threshold = 0.2, acf_threshold = 0.6) {  # Adjust the thresholds here
 
  # Dependencies
  load_package("splines")
  
  # Detrend time series using a regression spline with 3 knots
  knots <- quantile(1:length(ts), probs = seq(0, 1, length.out = 4))[-c(1, 4)]
  detrended <- residuals(lm(ts ~ ns(1:length(ts), knots = knots)))
                
  # Autocorrelation function for all lags up to 1/3 of series length
  acf_values <- acf(detrended, lag.max = length(ts) / 3, plot = FALSE)$acf
                
  # Find peaks and troughs in autocorrelation function
  diff_acf <- diff(acf_values)
  peaks <- which(diff_acf[-length(diff_acf)] >= 0 & diff_acf[-1] <= 0) + 1
  troughs <- which(diff_acf[-length(diff_acf)] <= 0 & diff_acf[-1] >= 0) + 1
                
  # Frequency is the first peak with conditions
  for (peak in peaks) {
  if (peak > min(troughs) && (acf_values[peak] - min(acf_values[troughs])) >= threshold && abs(acf_values[peak]) > acf_threshold) {
  return(peak)
    }
  }
                
  # If no such peak is found, frequency is set to 0 (equivalent to non-seasonal)
  return(0)
}

# Function which removes Seasonal Component ------------------------------------
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

# Function which extracts autocorrelation --------------------------------------------------------------
# represented by Ljung-Box Test Statistic
# we use a test statistic as quantification because a higher test statistic 
# often indicates stronger evidence against the null hypothesis
extract_autocorrelation <- function(ts){
  
  load_package("stats")
  
  lb_test <-  Box.test(ts, lag = 2*sqrt(length(ts)), type = "Ljung-Box")
  lb_statistic <- lb_test$statistic # rule of thumb: lag = 2sqrt(T)
  lb_significance <- as.numeric(lb_test$p.value < 0.05)
  # normalize the test statistic by the number of observations and the number of lags
  # more interpretable and comparable across different datasets or models 
  lb_statistic_normalized <- lb_statistic*lb_significance/(length(ts)*2*sqrt(length(ts)))
  return(lb_statistic_normalized)
  
}

# Function which extracts Long-range Dependence --------------------------------
# represented by an estimate of the Hurst parameter using the Whittle estimator (Stroe-Kunolds, 2009)
extract_long_range_dependence <- function(ts){
 
  load_package("longmemo")
  
  #assume fGN for model
  hurst_par <- WhittleEst(ts, model = "fGn")$coefficients[1]
  return(1-hurst_par)
  
}

# Function which extracts non-linearity in Autoregressive Structure of Mean ----
# represented by Teraesvirta Neural Network Test for Nonlinearity - Test Statistic
# we use a test statistic as quantification because a higher test statistic 
# often indicates stronger evidence against the null hypothesis
extract_nonlinearity_mean <- function(ts){
  
  load_package("tseries")
  load_package("stats")
  load_package("TSA")
  load_package("tsDyn")
  
  coef_diff <- 0
  
  
  fit <- tryCatch({
    setar(ts(ts), m=1)
    
  }, warning = function(w) {
    # If a warning is raised, return NULL
    NULL
  }, error = function(e) {
    # If an error is raised, return NULL
    NULL
  })
  
  if (!is.null(fit)) {
    
    # fit_summary <- summary(fit)
    # std_errors <- fit_summary$coef[,2]
    # error <- unname(std_errors["phiL.1"])*unname(std_errors["phiH.1"])
    
    # absolute value difference
    if (coef(fit)["phiH.1"]*coef(fit)["phiL.1"] > 0) {
      coef_diff <- unname(abs(abs(coef(fit)["phiH.1"]) - abs(coef(fit)["phiL.1"])))
    }
    else {
      coef_diff <- unname(abs(coef(fit)["phiH.1"]) + abs(coef(fit)["phiL.1"]))
    }
   
  }
  
  return(coef_diff)
  
  # trv_test <- terasvirta.test(ts(ts), lag = 1, type = "F")
  # 
  # trv_statistic <- trv_test$statistic # using conventional lag
  # trv_significance <- as.numeric(trv_test$p.value < 0.05)
  # # normalize the test statistic by the number of observations and the number of lags
  # # more interpretable and comparable across different datasets or models
  # trv_statistic_normalized <- trv_statistic*trv_significance/(length(ts)*1)
  # return(trv_statistic_normalized)
  
}


# Function which extracts Volatility Clustering --------------------------------
# represented by the sum of fitted GARCH coefficients
extract_vola_clustering <- function(ts, garch_model){
  
  load_package("rugarch")
  load_package("stats")

  fit <- tryCatch({
    rugarch::ugarchfit(spec = garch_model, data = ts)

  }, warning = function(w) {
    # If a warning is raised, return NULL
    NULL
  }, error = function(e) {
    # If an error is raised, return NULL
    NULL
  })
  
  
  coef_sum <- 0
  
  # try to zero parameters with p-values greater than 0.05, and re-estimate the model
  red <- tryCatch({
    rugarch::reduce(fit, pvalue = 0.05)

  }, warning = function(w) {
    # If a warning is raised, return NULL
    NULL
  }, error = function(e) {
    # If an error is raised, return NULL
    NULL
  })
  
  # If red is NULL (i.e., the solver failed to converge), zero out the coefficients
  if (!is.null(red)) {
  coef_sum <- unname(coef(red)["alpha1"] + coef(red)["beta1"])
  }

  return(abs(coef_sum))
  
}

# for(i in 1:nrow(df_simulated)) {
#   print(extract_vola_clustering(df_simulated[i, ], garch_model))
# }
# print(extract_vola_clustering(tsa_df_simulated[7, ], garch_model))
# test <- tsa_df_simulated[7, ]
# plot(test, type = "l")
# fit_test <- rugarch::ugarchfit(spec = garch_model, data = test)
# extract_vola_clustering(test, garch_model)


# Function which extracts Non-linearity in the Volatility: Leverage effect -----
# represented by the gamma coefficient of a fitted exponential GARCH
extract_nonlinearity_vola <- function(ts, egarch_model){
  
  load_package("rugarch")
  
  efit <- tryCatch({
    rugarch::ugarchfit(spec = egarch_model, data = ts)
    
  }, warning = function(w) {
    # If a warning is raised, return NULL
    NULL
  }, error = function(e) {
    # If an error is raised, return NULL
    NULL
  })
  
  gamma_par <- 0
  # try to zero parameters with p-values greater than 0.05, and re-estimate the model
  
  ered <- tryCatch({
    rugarch::reduce(efit, pvalue = 0.05)
    
  }, warning = function(w) {
    # If a warning is raised, return NULL
    NULL
  }, error = function(e) {
    # If an error is raised, return NULL
    NULL
  })
  
  # If ered is NULL (i.e., the solver failed to converge), zero out the coefficients
  if (!is.null(ered)) {
    # If ered is not NULL, check if gamma1 is negative
    gamma1_negative <- as.numeric(coef(ered)["gamma1"] < 0)
    gamma_par <- coef(ered)["gamma1"]*gamma1_negative
  }
  
  return(unname(gamma_par))
}


# [-1, 1] Scaling --------------------------------------------------------------
rescale_features <- function(df, tsa){
  
  load_package("scales")
  
  # Define custom scaling function
  rescale_2 <- function(x) {2 * ((x - min(x)) / (max(x) - min(x))) - 1}
  rescale_3 <- function(x) {exp(-(n/x)*0.01)}
  
  # Create a copy of the original data frame
  scaled_df <- df
   
  if (tsa == T){
    selected_cols_1 <- c("SD", "TSA_SD", "Nonlinearity_Mean", "TSA_Nonlinearity_Mean")
    selected_cols_2 <- c("Mean","Skew", "Kurtosis", "TSA_Mean", "TSA_Skew", "TSA_Kurtosis")
  }
  else{
    # Select only the "Mean", "SD", "Skew", and "Kurtosis" columns
    selected_cols_1 <- c("SD", "Nonlinearity_Mean")
    selected_cols_2 <- c("Mean", "Skew", "Kurtosis")
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


