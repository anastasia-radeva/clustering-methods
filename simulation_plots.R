# Outline
# 1 READ ME
# 2 Dependencies
# 3 SARFIMA
#   3.2 Mean
#   3.3 Volatility
#   3.4 Skewness
#   3.5 Kurtosis
#   3.6 Drift
#   3.7 Seasonality
#   3.8 Periodicity
#   3.9 Autocorrelation
#   3.10 Long-range Dependence
#   3.11 Non-linearity in Mean
# 4 GARCH
#   4.1 Parameters
#   4.2 Non-linearity in Volatility
#   4.3 Volatility Clustering
# Function which plots data frames

# --------------------------------------------------------------------------------
# 1.READ ME
# Plotting the simple synthetic processes, each representing one time series feature 
# from either the clustering literature or the stylized facts about financial return series
# Simulated with ARFIMA and GARCH in SARFIMA-GARCH_Simulation.R
# Goal: simple process, which can be classified by a prominent feature
# *the simulations do not strive to represent realistic processes

# 2 Dependencies ---------------------------------------------------------------
load_package <- function(pack){
  if (!is.element(pack, installed.packages()[,1])){
    install.packages(pack, dependencies = TRUE)
  }
  library(pack, character.only = TRUE)
}

load_package("ggplot2")
load_package("e1071")

source("SARFIMA-GARCH_Simulation.R")

# 3 SARFIMA --------------------------------------------------------------------

# 3.1 Mean ---------------------------------------------------------------------

plot(mean_series_1, type = "l", main = "Mean 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(mean_series_2, type = "l", main = "Mean 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 3.2 Volatility ---------------------------------------------------------------

plot(sd_series_1, type = "l", main = "Standard Deviation 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(sd_series_2, type = "l", main = "Standard Deviation 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 3.3 Skewness -----------------------------------------------------------------

plot(skewed_series_1, type = "l", main = "Skewness 1", xlab = "", ylab = "", ylim = c(-10, 10))
# density plot
plot(density(skewed_series_1), type = "l")
# skewness value
# skewness(skewed_series_1) # around -1
# kurtosis(skewed_series_1) # around 3

plot(skewed_series_2, type = "l", main = "Skewness 2", xlab = "", ylab = "", ylim = c(-10, 10))
# density plot
plot(density(skewed_series_2), type = "l")
#skewness value
# skewness(skewed_series_2) # around -1
# kurtosis(skewed_series_2) # around 3

# 3.4 Kurtosis -----------------------------------------------------------------

plot(leptokurtic_series_1, type = "l", main = "Kurtosis 1", xlab = "", ylab = "", ylim = c(-10, 10))
# density plot
plot(density(leptokurtic_series_1))
# kurtosis value
# kurtosis(leptokurtic_series_1) # a lot higher than 3

plot(leptokurtic_series_2, type = "l", main = "Kurtosis 2", xlab = "", ylab = "", ylim = c(-10, 10))
# density plot
plot(density(leptokurtic_series_2))
# kurtosis value
# kurtosis(leptokurtic_series_2) # a lot higher than 3

# 3.5 Trend --------------------------------------------------------------------

plot(drift_series_1, type = "l", main = "Drift 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(drift_series_2, type = "l", main = "Drift 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 3.6 Seasonality --------------------------------------------------------------

plot(seasonal_series_1, type = "l", main = "Seasonality 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(seasonal_series_2, type = "l", main = "Seasonality 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 3.7 Periodicity --------------------------------------------------------------

plot(periodic_series_1, type = "l", main = "Periodicity 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(periodic_series_2, type = "l", main = "Periodicity 2", xlab = "", ylab = "", ylim = c(-10, 10))


# 3.8 Autocorrelation ----------------------------------------------------------

plot(autocorrelated_series_1, type = "l", main = "Autocorrelation 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(autocorrelated_series_2, type = "l", main = "Autocorrelation 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 3.9 Long-range Dependence ---------------------------------------------------
# On the plots, the long range effects is visible - through this amplification/dampening of the autocorrelation
#                                                 - better starting from the lag where the autocorrelation function 
#                                                   of the corresponding arma process is less prominent than the long memory effect, 
#                                                   so the latter prevails
# It is hard to detect the magnitude of the two effects, when they have the same sign. 
# ARFIMA(ar1_1, dfrac)

plot(long_range_series_1, type = "l", main = "Long-range Dependence 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(long_range_series_2, type = "l", main = "Long-range Dependence 2", xlab = "", ylab = "", ylim = c(-10, 10))

# slower decaying autocorrelation than ARMA(1, 0)
sum(abs(autocorrelations(long_range_series_2)))
sum(abs(autocorrelations(autocorrelated_series_1)))
autocorrelations(long_range_series_1, maxlag = 100)
autocorrelations(autocorrelated_series_1, maxlag = 100)

# 3.10 Non linearity in the Mean -----------------------------------------------

plot(nonlinear_mean_series_1, type = "l", main = "Nonlinearity in Mean 1", xlab = "", ylab = "", ylim = c(-10, 10))
plot(nonlinear_mean_series_2, type = "l", main = "Nonlinearity in Mean 2", xlab = "", ylab = "", ylim = c(-10, 10))

# 4 GARCH ----------------------------------------------------------------------

# 4.1 Non.linearity in the volatility ------------------------------------------

# plot conditional sd
plot(nonlinear_vola_series_1, type = "l", main = "Nonlinearity in Vola 1", xlab = "", ylab = "", ylim = c(-10, 10))
skewness(nonlinear_vola_series_1)
skewness(nonlinear_s4_1@path$sigmaSim)
# VERY DIFF, maybe consider a skewed distribution model, if you are going to keep the noise
# if not, find a way to scale the sigma
# plot sigma series
plot(nonlinear_s4_1, which = 1, type = "l")

# plot conditional sd
plot(nonlinear_vola_series_2, type = "l", main = "Nonlinearity in Vola 2", xlab = "", ylab = "", ylim = c(-10, 10))
# the conditional volatility series and the conditional+unconditional volatility series have different skewness  
skewness(nonlinear_vola_series_2)
skewness(nonlinear_s4_2@path$sigmaSim)
# plot series (in this case, same as the residuals series, no arma)
plot(nonlinear_s4_2, which = 1, type = "l")

# 4.2 VOLATILITY CLUSTERING ----------------------------------------------------

# plot series 
plot(vola_clustering_series_1, type = "l", main = "Volatility CLustering 1", xlab = "", ylab = "")
# plot conditional sd
plot(vola_clustering_s4_1, which = 1, type = "l")

# plot series
plot(vola_clustering_series_2, type = "l", main = "Volatility CLustering 2", xlab = "", ylab = "")
#plot conditional sd
plot(vola_clustering_s4_2, which = 1, type = "l")
