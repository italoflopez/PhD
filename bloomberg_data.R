library(forecast)
library(zoo)
library(tseries)
library(vars)
library(plotly)
library(quantmod)
library(Metrics)
library(fGarch)
library(rugarch)
library(zoo)
library(tis)
library(forecast)
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

library(data.table)

setwd("C://Users//Italo//OneDrive//Documents//PhD//PhD//First")

bloomberg_data <- read.csv("data_bloomberg.csv",header = T)

bloomberg_data%>%select(VARIABLE)%>%unique()

variables_bloomberg <- read.csv("variables_bloomberg.csv", skip = 1, header = T)

slow_moving_bloomberg_variables <- variables_bloomberg$Slow.moving.variables

fast_moving_variables <- variables_bloomberg$Fast.moving.variables

pce_core <- bloomberg_data[bloomberg_data$VARIABLE==slow_moving_bloomberg_variables[1],c("FECHA","PX_LAST")]

pce_core$FECHA <- ifelse(as.numeric(substr(pce_core$FECHA, 8, 9)) > 24, gsub("-(\\d{2})$", "-19\\1", pce_core$FECHA), gsub("-(\\d{2})$", "-20\\1", pce_core$FECHA))


# Create a zoo object
pce_core <- zoo(pce_core$PX_LAST, order.by = as.yearmon(pce_core$FECHA, format = "%d-%b-%Y"))

# Convert the index to Date format with the desired day (01)
#index(pce_core) <- as.Date(format(index(zoo_object), "%Y-%m-01"))

manu_pmi<- bloomberg_data[bloomberg_data$VARIABLE==slow_moving_bloomberg_variables[15],c("FECHA","PX_LAST")]

manu_pmi <- zoo_bloomberg_monthly(manu_pmi)



us_high_yield<- bloomberg_data[bloomberg_data$VARIABLE==fast_moving_variables[8],c("FECHA","PX_LAST")]

us_high_yield <- zoo_bloomberg_monthly(us_high_yield)



commodities_price <- bloomberg_data[bloomberg_data$VARIABLE=="Commodity Research Bureau BLS/US Spot All Commodities",c("FECHA","PX_LAST")]
