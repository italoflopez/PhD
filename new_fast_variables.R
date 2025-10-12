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

library(readxl)

library(lubridate)

setwd("C://Users//Italo//OneDrive//Documents//PhD//PhD//First")

sector_sp500 <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//sectoriales_sp500.xlsx", sheet = 1)

# ordering the data data frame
sector_sp500 <- sector_sp500 %>%
  arrange(VARIABLE, FECHA)

# selecting the columns we want
sector_sp500 <- sector_sp500 %>%
  select(FECHA,VARIABLE,PX_LAST)


# Pivoting it wider
sector_sp500 <- sector_sp500 %>%
  pivot_wider(names_from = VARIABLE, values_from = PX_LAST)


# Convert daily data to monthly averages
sector_sp500 <- sector_sp500 %>%
  mutate(Month = floor_date(FECHA, "month")) %>% # Add a Month column
  group_by(Month) %>% # Group by Month
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Calculating log differences
sector_sp500 <- sector_sp500 %>%
  # Apply log difference to all columns except FECHA
  mutate(across(-Month, ~ c(NA, diff(log(.))*100)))


# Converting to a zoo object
sector_sp500 <- zoo(sector_sp500, order.by = as.Date(sector_sp500$Month))

# Specify the columns to drop by name
sector_sp500 <- sector_sp500[, !(colnames(sector_sp500) %in% c("Month", "SP500 Real State"))]

#Droppping NAs
sector_sp500 <- sector_sp500%>%na.omit()


vix <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//vix.xlsx", sheet = 1)

# ordering the data data frame
vix <- vix %>%
  arrange(VARIABLE, FECHA)

# selecting the columns we want
vix <- vix %>%
  select(FECHA,VARIABLE,PX_LAST)


# Pivoting it wider
vix <- vix %>%
  pivot_wider(names_from = VARIABLE, values_from = PX_LAST)


# Convert daily data to monthly averages
vix <- vix %>%
  mutate(Month = floor_date(FECHA, "month")) %>% # Add a Month column
  group_by(Month) %>% # Group by Month
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Converting to a zoo object
vix <- zoo(vix, order.by = as.Date(vix$Month))

# Specify the columns to drop by name
vix <- vix[, !(colnames(vix) %in% c("Month"))]


#All new fast variables together
new_fast_variables <- merge(sector_sp500,vix)


#Droppping NAs
new_fast_variables <- new_fast_variables%>%na.omit()

write.csv(as.data.frame(new_fast_variables),"new_fast_variables.csv")
