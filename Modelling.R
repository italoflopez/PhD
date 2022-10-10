library(ROracle)
library(keyring)
library(forecast)
library(zoo)
library(tseries)
library(vars)
library(portes)
library(plotly)
library(quantmod)
library(Metrics)
library(fGarch)
library(rugarch)
library(zoo)
library(tis)
library(forecast)
library(astsa)
library(plotly)
library(quantmod)
library(sarima)
library(pracma)

# Time Series ML
library(tidymodels)
library(modeltime)
library(modeltime.ensemble)

# Timing & Parallel Processing
library(tictoc)
library(future)
library(doFuture)

# Core 
library(tidyquant)
library(tidyverse)
library(timetk)

library(bvartools)

final_data

pca <- prcomp(final_data)

summary(pca)
str(pca)
