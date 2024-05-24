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

simulation_plots <- list()

# 3 SARFIMA --------------------------------------------------------------------

# 3.1 Unconditional Volatility -------------------------------------------------
# Create the initial plot with the first series

plot(wn_1, type = "l", col = "deepskyblue3", lwd = 1, main = "White Noise", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(wn_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("wn_1", "wn_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["sd"]] <- recordPlot()

# 3.2 Skewness -----------------------------------------------------------------
plot(skewed_wn_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Skew White Noise", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(skewed_wn_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("skewed_wn_1", "skewed_wn_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["skewness"]] <- recordPlot()

# density plot
plot(density(skewed_wn_1), type = "l")
# skewness value
# moments::skewness(skewed_series_1) # around -1
# moments::kurtosis(skewed_series_1) # around 3

# density plot
plot(density(skewed_wn_2), type = "l")
#skewness value
# moments::skewness(skewed_series_2) # around -1
# moments::kurtosis(skewed_series_2) # around 3

# 3.3 Kurtosis -----------------------------------------------------------------
plot(leptokurtic_wn_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Leptokurtic White Noise", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(leptokurtic_wn_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("leptokurtic_wn_1", "leptokurtic_wn_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["kurtosis"]] <- recordPlot()

# density plot
plot(density(leptokurtic_wn_1))
# kurtosis value
# moments::kurtosis(leptokurtic_series_1) # a lot higher than 3

# density plot
plot(density(leptokurtic_wn_2))
# kurtosis value
# moments::kurtosis(leptokurtic_series_2) # a lot higher than 3

# 3.4 Trend --------------------------------------------------------------------
plot(trend_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Upward Trend", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(trend_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("trend_1", "trend_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["trend"]] <- recordPlot()


# 3.5 Seasonality --------------------------------------------------------------
plot(seasonal_ar_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Seasonal AR(1) with Different Periods", xlab = "n", ylab = "")
lines(seasonal_ar_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("seasonal_ar_1", "seasonal_ar_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["seasonal"]] <- recordPlot()


# 3.6 Periodicity --------------------------------------------------------------
plot(seasonal_ar_3, type = "l", col = "deepskyblue3", lwd = 1, main = "Seasonal AR(1) with Identical Periods", xlab = "n", ylab = "")
lines(seasonal_ar_4, lwd = 1, col = "lightpink3")
legend("topright", legend = c("seasonal_ar_3", "seasonal_ar_4"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["period"]] <- recordPlot()

# 3.7 Autocorrelation ----------------------------------------------------------
plot(ar_1, type = "l", col = "deepskyblue3", lwd = 1, main = "AR(1)", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(ar_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("ar_1", "ar_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["autocorrelation"]] <- recordPlot()


# 3.8 Long-range Dependence ---------------------------------------------------
# On the plots, the long range effects is visible - through this amplification/dampening of the autocorrelation
#                                                 - better starting from the lag where the autocorrelation function 
#                                                   of the corresponding arma process is less prominent than the long memory effect, 
#                                                   so the latter prevails
# It is hard to detect the magnitude of the two effects, when they have the same sign. 
# ARFIMA(ar1_1, dfrac)
plot(fractional_diff_ar_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Fractionally Integrated AR(1)", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(fractional_diff_ar_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("fractional_diff_ar_1", "fractional_diff_ar_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["long_range_dependence"]] <- recordPlot()


# slower decaying autocorrelation than ARMA(1, 0)
sum(abs(autocorrelations(long_range_series_2)))
sum(abs(autocorrelations(autocorrelated_series_1)))
autocorrelations(long_range_series_1, maxlag = 100)
autocorrelations(autocorrelated_series_1, maxlag = 100)

# 3.9 Non linearity in the Mean -----------------------------------------------
plot(tar_1, type = "l", col = "deepskyblue3", lwd = 1, main = "Threshold Regime Switching AR(1)", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(tar_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("tar_1", "tar_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["nonlinear_mean"]] <- recordPlot()

# 4 GARCH ----------------------------------------------------------------------

# 4.1 Nonlinearity in the volatility ------------------------------------------

plot(egarch_2, type = "l", col = "deepskyblue3", lwd = 1, main = "EGARCH(1, 1)", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(egarch_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("egarch_1", "egarch_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["nonlinear_vola"]] <- recordPlot()

# the conditional volatility series and the conditional+unconditional volatility series have different skewness  
skewness(nonlinear_s4_2@path$sigmaSim)

# 4.2 Volatility Clustering ----------------------------------------------------
plot(garch_1, type = "l", col = "deepskyblue3", lwd = 1, main = "GARCH(1, 1)", xlab = "n", ylab = "", ylim = c(-10, 10))
lines(garch_2, lwd = 1, col = "lightpink3")
legend("topright", legend = c("garch_1", "garch_2"), col = c("deepskyblue3", "lightpink3"), lty = 1, lwd = 2, cex = 0.8)
simulation_plots[["vola_clustering"]] <- recordPlot()

# save plots -------------------------------------------------------------------
save_plots(simulation_plots, "simulation_plots")
df_plot(as.data.frame(df_simulated))

# separate plots ---------------------------------------------------------------
global_y_range <- range(df_simulated, na.rm = TRUE)
# Create a data frame
df <- data.frame(Column = 1:length(trend_2), Value = trend_2)

df_row <- data.frame(Column = colnames(df), Value = unlist(df[i, ]))
df_row$Group <- rownames(df)[i]
p <- ggplot(df, aes(Column, Value)) +
  geom_line() + 
  ggtitle("") +
  xlab("") +
  ylab("") +
  ylim(global_y_range) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"))
print(p)
