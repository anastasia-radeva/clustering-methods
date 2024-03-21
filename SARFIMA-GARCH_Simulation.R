# Outline
# 1 READ ME
# 2 Dependencies
# 3 SARFIMA
#   3.1 Parameters
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
# 
# Note: 3.x SARFIMA and 4.y GARCH are both structured in the following way.
#     - Data generation for first time series
#     - Data generation for second time series

# 1 READ ME --------------------------------------------------------------------
# Simple synthetic processes, each representing one time series feature 
# from either the clustering literature or the stylized facts about financial return series
# Simulated with ARFIMA and GARCH 
# Goal: simple process, which can be classified by a prominent feature
# *the simulations do not strive to represent realistic processes

# 2 Dependencies ---------------------------------------------------------------
load_package("sarima")
load_package("astsa")
load_package("arfima")
load_package("rugarch")
load_package("TSA")
load_package("fGarch")

# 3 SARFIMA --------------------------------------------------------------------

# 3.1 SARFIMA Parameters

# time series have the same length
n <- 500
# constants
c_1 <- 0
c_2 <- 2
# drift
drift <- 0.01
# standard deviation
sd <- 0.1
# skew
skew_asym = 0.2
# degrees of freedom in Skew Student-t distribution
shape_lepto <- 3
# seasonal autocorrelation coefficient
# An SAR signature usually occurs when the autocorrelation at the seasonal period is positive,
# therefore both coefficients are positive.
# whereas an SMA signature usually occurs when the seasonal autocorrelation is negative.
sar1_1 <- 0.97
sar1_2 <- 0.95
# Autocorregressive Coefficient
ar1_1 <- -0.8
ar1_2 <- 0.9
# Fractional Differencing Parameter
d <- 0.3
# Threshold Autoregressive Coefficients
tar1_1 <- 0.85
tar1_2 <- 0.25
# Threshold
r <- 1
# Delay
del <- 1

# 3.2 Mean ---------------------------------------------------------------------
# two base constant processes with different means 
# ARFIMA(-), e ~ N(c_x, 0)
# mean and sd are arguments of the rand.gen 
# so they don't determine the whole series' mean and sd, just the noise's
# we use the arfima.sim function for consistency and comparability later
set.seed(1)
mean_series_1 <- arfima.sim(n = n, muHat = c_1, sd = 0)
set.seed(2)
mean_series_2 <- arfima.sim(n = n, muHat = c_2, sd = 0)

# 3.3 Volatility ---------------------------------------------------------------
# similar standard deviation, but different means
# ARFIMA(-), e ~ N(c_x, sd)
set.seed(1)
sd_series_1 <- arfima.sim(n = n, muHat = c_1)
set.seed(2)
sd_series_2 <- arfima.sim(n = n, muHat = c_1)

# 3.4 Skewness -----------------------------------------------------------------
# same skew parameter, but different means
# ARFIMA(-), skewed N(c_x sd, xi)
# skewed normal distribution
# The skewness Xi in fGarch is based on Fernandez & Steel 2000
# it is an inverse scale factor in the positive and the negative directions
# for Xi = 1, the distribution is symmetrical. When 0 < Xi < 1, negative skew. When Xi > 1, positive skew.
# However, it is important that the probability has still zero mean and unit variance. 
# Otherwise, it would be impossible, or at least difficult, to separate the fluctuations in the mean and variance 
# from the fluctuations in the shape of the density. 
# we use our existing mean and sd for comparability

set.seed(1)
skewed_series_1 <- arfima.sim(n = n, rand.gen = rsnorm, muHat = c_1, xi = skew_asym)
set.seed(2)
skewed_series_2 <- arfima.sim(n = n, rand.gen = rsnorm, muHat = c_1, xi = skew_asym)

# 3.5 Kurtosis -----------------------------------------------------------------
# same shape parameter, but different means
# ARFIMA (-), T(c_x, sd_x, nu)
# student t-distribution
# shape parameter nu > 2, the higher nu, the lighter the tails, approaches normal distribution

set.seed(3)
leptokurtic_series_1 <- arfima.sim(n = n, rand.gen = rstd, muHat = c_1, nu = shape_lepto)
set.seed(9)
leptokurtic_series_2 <- arfima.sim(n = n, rand.gen = rstd, muHat = c_1, nu = shape_lepto)

# 3.6 Trend --------------------------------------------------------------------
# same drift, but different intercept
# ARFIMA(d = 1)
# mean parameter determines the trend, when d > 0

set.seed(1)
drift_series_1 <- arfima.sim(n = n, model = list(dint = 1), muHat = drift, sd = 0)
set.seed(2)
drift_series_2 <- c_2 + arfima.sim(n = n, model = list(dint = 1), muHat = drift, sd = 0)

# 3.7 Seasonality --------------------------------------------------------------
# same seasonal autocorregressive coefficient, but different periodicity
# SARFIMA(ar1_1, sar1_1, period = n/x)

# the first 1/x*n components of the innovation play a role similar to the intercept, 
# as arfima.sim does not have a parameter for this
set.seed(1)
seasonal_series_1 <- arfima.sim(model = list(seasonal = list(phi = sar1_1, period = n/4)), 
                                n = n, innov = c(rep(c_2, 1/4*n), rep(0, (3/4)*n)) + rnorm(n))
set.seed(2)
seasonal_series_2 <- arfima.sim(model = list(seasonal = list(phi = sar1_1, period = n/5)), 
                                n = n, innov = c(rep(c_2, 1/5*n), rep(0, 4/5*n)) + rnorm(n))

# 3.8 Periodicity --------------------------------------------------------------
# same periodicity, but different seasonal autoregressive coefficients and shape
# SARFIMA(sar_x)

set.seed(1)
periodic_series_1 <- arfima.sim(model = list(seasonal = list(phi = sar1_1, period = n/10)), 
                                n = n, innov = c(rep(c_2, 1/20*n), rep(0, 19/20*n)) + rnorm(n))
set.seed(2)
periodic_series_2 <- arfima.sim(model = list(seasonal = list(phi = sar1_2, period = n/10)), 
                                n = n, innov = c(rep(c_2, 1/20*n), rep(0, 19/20*n)) + rnorm(n))

# 3.9 Autocorrelation ----------------------------------------------------------
# same AR coefficients, but different means
# ARFIMA(a1_1), first component of innov is used as intercept, as arfima.sim does not have a parameter for this

# autocorrelated_series_1 <- arfima.sim(model = list(phi = ar1_1), n = n, innov = c(c_2, rep(0, n-1)))
# autocorrelated_series_2 <- arfima.sim(model = list(phi = ar1_1), n = n, innov = c(c_2*2, rep(0, n-1)))
set.seed(1)
autocorrelated_series_1 <- arfima.sim(model = list(phi = ar1_1), n = n, sigma2 = 1)
set.seed(2)
autocorrelated_series_2 <- arfima.sim(model = list(phi = ar1_1), n = n, sigma2 = 1)

# 3.10 Long-range Dependence ---------------------------------------------------
# same fractional differencing parameter 
# but different intercept
# Definition in Vera-Valdes (2021): for -0.5 < d < 0.5, the series exhibits long-range dependence
# 0 < d < 0.5 long memory
# autocorrelation function decays at a hyperbolic rate and the series is said to be persistent
# -0.5 < d < 0 negative long memory or anti-persistence
# large values are likely to be followed by small values and small values by large values
# Interaction b/w phi, theta and d: If the signs are the same, the arma effects are amplified by d. 
#                                   Otherwise, if the signs are different, the arma effects are dampened. 
# It is hard to detect the magnitude of the two effects, when they have the same sign. 
# ARFIMA(ar1_1, dfrac)

# long_range_series_1 <- arfima.sim(model = list(phi = ar1_1, dfrac = -d), n = n, innov = c(c_2, rep(0, n-1)))
# long_range_series_2 <- arfima.sim(model = list(phi = ar1_1, dfrac = -d), n = n, innov = c(c_2*2, rep(0, n-1)))
set.seed(1)
long_range_series_1 <- arfima.sim(model = list(phi = ar1_1, dfrac = -d), n = n, sigma2 = 1)
set.seed(2)
long_range_series_2 <- arfima.sim(model = list(phi = ar1_1, dfrac = -d), n = n, sigma2 = 1)

# 3.11 Non linearity in the Mean -----------------------------------------------
# Threshold Regime Switching Autoregressive Process

set.seed(1)
nonlinear_mean_series_1 <- tar.sim(ntransient = 0, n = n, Phi1 = c(0, tar1_1), Phi2 = c(0, tar1_2), 
                                  thd = r, p = 1, d = del, sigma1 = 1, sigma2 = 1)$y
set.seed(2)
nonlinear_mean_series_2 <- tar.sim(ntransient = 0, n = n, Phi1 = c(0, tar1_1), Phi2 = c(0, tar1_2), 
                                   thd = r, p = 1, d = del, sigma1 = 1, sigma2 = 1)$y

# 4 GARCH ----------------------------------------------------------------------

# 4.1 ARCH-GARCH Parameters
# constants
omega_1 <- 0.1
omega_2 <- 0.15
# ARCH parameters
alpha1_1 <- 0.1
alpha1_2 <- 0.5
# GARCH parameters
beta1_1 <- 0.2 
beta1_2 <- 0.49 # sum of alpha1_2 and beta1_2 close to 1
# asymmetry parameter for eGARCH 
gamma_1 <- -0.6
# fractional differencing parameter

# 4.2 Nonlinearity in the volatility -------------------------------------------
# same gamma and alpha
# but different residuals
# leverage effect - negative shocks have a stronger impact on volatility
# EGARCH(alpha1_1, beta1_1), ARMA(-), e ~ N(0, 1)

nonlinear_spec = ugarchspec(variance.model = list(model = "eGARCH", 
                                                  garchOrder = c(1, 1), 
                                                  variance.targeting = FALSE), 
                            mean.model = list(armaOrder = c(0, 0),
                                              include.mean = F), 
                            distribution.model = "norm")

# parameter assignment for first series
setfixed(nonlinear_spec) <- list(omega = omega_1, alpha1 = alpha1_1, beta1 = beta1_1, gamma1 = gamma_1)
set.seed(2)
nonlinear_s4_1 = ugarchpath(nonlinear_spec, n.sim = n)
nonlinear_vola_series_1 <- nonlinear_s4_1@path$seriesSim

set.seed(3)
nonlinear_s4_2 = ugarchpath(nonlinear_spec, n.sim = n)
nonlinear_vola_series_2 <- nonlinear_s4_2@path$seriesSim

# 4.3 Volatility Clustering ----------------------------------------------------
# same volatility clustering coefficients (alpha1 and beta1)
# but different constant terms omega
# GARCH(alpha1_2, beta1_2), ARMA(-), e ~ N(0, 1)
vola_clustering_spec = ugarchspec(variance.model = list(model = "sGARCH", 
                                        garchOrder = c(1, 1), 
                                        variance.targeting = FALSE), 
                  mean.model = list(armaOrder = c(0, 0),
                                    include.mean = F), 
                  distribution.model = "norm")
# parameter assignment for first series
# alpha1 + beta1 close to 1 to simulate vola clustering
setfixed(vola_clustering_spec) <- list(omega = omega_1, alpha1 = alpha1_2, beta1 = beta1_2)
set.seed(2)
vola_clustering_s4_1 = ugarchpath(vola_clustering_spec, n.sim = n)
vola_clustering_series_1 <- vola_clustering_s4_1@path$seriesSim

set.seed(3)
vola_clustering_s4_2 = ugarchpath(vola_clustering_spec, n.sim = n)
vola_clustering_series_2 <- vola_clustering_s4_2@path$seriesSim

# Data Frame -------------------------------------------------------------------

simulated_data <- data.frame(sd_series_1, sd_series_2, 
                             skewed_series_1, skewed_series_2, 
                             leptokurtic_series_1, leptokurtic_series_2, 
                             drift_series_1, drift_series_2, 
                             seasonal_series_1, seasonal_series_2, 
                             periodic_series_1, periodic_series_2, 
                             autocorrelated_series_1, autocorrelated_series_2, 
                             long_range_series_1, long_range_series_2, 
                             nonlinear_mean_series_1, nonlinear_mean_series_2, 
                             nonlinear_vola_series_1, nonlinear_vola_series_2, 
                             vola_clustering_series_1, vola_clustering_series_2)



df_simulated <- t(simulated_data)

# create data labels for later
simulated_data_labels <- c("sd", "sd", 
                          "skewed", "skewed", 
                          "leptokurtic", "leptokurtic", 
                          "drift", "drift", 
                          "seasonal", "seasonal", 
                          "periodic", "periodic", 
                          "autocorrelated", "autocorrelated", 
                          "long_range", "long_range", 
                          "nonlinear_mean", "nonlinear_mean", 
                          "nonlinear_vola", "nonlinear_vola", 
                          "vola_clustering", "vola_clustering")
simulated_data_labels_factor <- as.factor(simulated_data_labels)

# Convert factor to numeric
simulated_data_labels_numeric <- as.numeric(simulated_data_labels_factor)

