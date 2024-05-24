# clustering-time-series-data

- feature-based clustering
- simulation study
- estimating clusters in forecast loss differential panels of bond liquidity

This repository contains the implementation accompanying my bachelor's thesis on the topic of 'Clustering Time-Series Data'.
The study aims to estimate clusters within panels of time-series data. The particular time series are forecast loss differentials of models of bond liquidity. The cluster estimation is needed to enable a consistent estimation of the panels’ variance, used in the Diebold-Mariano test for comparison of one-step forecasts. We propose a method that groups data by a set of extracted features. This approach addresses several challenges posed by time-series data during the clustering procedure, including high dimensionality, noise, and temporal and amplitude translations of patterns. To assess the efficacy of the proposed method in combination with various clustering algorithms, we conduct a simulation study. Following this, we evaluate the performance of these procedures and apply the most effective methods on the real time-series data panels.
