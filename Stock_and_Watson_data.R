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

source("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//plotly_simple.R")

urate_job_losers <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//66_U Job Losers.xls", sheet = 1, skip = 11)

urate_job_losers <- zoo(urate_job_losers$`1.6000000000000001`, order.by = urate_job_losers$`24473`)

data_stock_and_watson <- urate_job_losers

ulevel_reentrants <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//67_U LF Reenty.xls", sheet = 1, skip = 11)

ulevel_reentrants <- zoo(ulevel_reentrants$`1098`, order.by = ulevel_reentrants$`24473`)

data_stock_and_watson <- merge(data_stock_and_watson,ulevel_reentrants)

colnames(data_stock_and_watson) <- c("Unemployment Rate - Job Losers","Unemployment Level - Reentrants to Labor Force, Thousands of Persons, Monthly, Seasonally Adjusted")

u_job_leavers <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//68_U Job Leavers.xls", sheet = 1, skip = 11)

u_job_leavers <- zoo(u_job_leavers$`404`, order.by = u_job_leavers$`24473`)

data_stock_and_watson <- merge(data_stock_and_watson,u_job_leavers)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Unemployment Level - Job Leavers, Thousands of Persons, Monthly, Seasonally Adjusted"

u_new_entrants <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//69_U New Entrants.xls", sheet = 1, skip = 11)

u_new_entrants <- zoo(u_new_entrants$`409`, order.by = u_new_entrants$`24473`)

data_stock_and_watson <- merge(data_stock_and_watson,u_new_entrants)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Unemployment Level - New Entrants, Thousands of Persons, Monthly, Seasonally Adjusted"

emp_part_time <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//70_Emp_SlackWk.xls", sheet = 1, skip = 11)

emp_part_time <- zoo(emp_part_time$`2063`, order.by = emp_part_time$`20210`)

data_stock_and_watson <- merge(data_stock_and_watson,emp_part_time)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Employment Level - Part-Time for Economic Reasons, All Industries, Thousands of Persons, Monthly, Seasonally Adjusted"

ppi_finished_consumer_goods <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//124_PPI_FinConsGds.xls", sheet = 1, skip = 11)

ppi_finished_consumer_goods%>%str()

ppi_finished_consumer_goods <- zoo(ppi_finished_consumer_goods$`28.300000000000001`, order.by = ppi_finished_consumer_goods$`17258`)

data_stock_and_watson <- merge(data_stock_and_watson,ppi_finished_consumer_goods)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Producer Price Index by Commodity: Final Demand: Personal Consumption Goods (Finished Consumer Goods), Index 1982=100, Monthly, Seasonally Adjusted"



ppi_industrial_commodities <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//126_PPI_IndCom.xls", sheet = 1, skip = 11)

ppi_industrial_commodities%>%str()

ppi_industrial_commodities <- zoo(ppi_industrial_commodities$`12.300000000000001`, order.by = ppi_industrial_commodities$`4750`)

data_stock_and_watson <- merge(data_stock_and_watson,ppi_industrial_commodities)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Producer Price Index by Commodity: Industrial Commodities, Index 1982=100, Monthly, Not Seasonally Adjusted"



ppi_natural_gas <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//131_Real_Price_NatGas_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 11)

ppi_natural_gas%>%str()

ppi_natural_gas <- zoo(ppi_natural_gas$`7.4000000000000004`, order.by = ppi_natural_gas$`24473`)



cpi<-getSymbols('CPIAUCSL',src='FRED')

cpi<-CPIAUCSL

cpi%>%str()

cpi <- zoo(cpi%>%coredata(), order.by = cpi%>%index()%>%as.POSIXct())

ppi_natural_gas_real <- ppi_natural_gas/cpi

data_stock_and_watson <- merge(data_stock_and_watson,ppi_natural_gas_real)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Producer Price Index by Commodity: Fuels and Related Products and Power: Natural Gas, Index 1982=100, Monthly, Not Seasonally Adjusted, Deflated"


gs1_tb3m <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//154_GS1_Tb3m.xls", sheet = 1, skip = 5)

gs1_tb3m%>%str()

gs1_tb3m <- zoo(gs1_tb3m$`...2`, order.by = as.Date(paste0(gs1_tb3m$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

data_stock_and_watson <- merge(data_stock_and_watson,gs1_tb3m)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "GS1_Tb3m"



cp_tb3m <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//156_CP_Tbill Spread.xls", sheet = 1, skip = 5)

cp_tb3m%>%str()

cp_tb3m <- zoo(cp_tb3m$`...2`, order.by = as.Date(paste0(cp_tb3m$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

data_stock_and_watson <- merge(data_stock_and_watson,cp_tb3m)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "cp_tb3m"



commercial_and_industrial_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//165_Real_C&Lloand_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 11)

commercial_and_industrial_loans%>%str()

commercial_and_industrial_loans <- zoo(commercial_and_industrial_loans$`11.289400000000001`, order.by = commercial_and_industrial_loans$`17168`)

commercial_and_industrial_loans_real <- commercial_and_industrial_loans/cpi

data_stock_and_watson <- merge(data_stock_and_watson,commercial_and_industrial_loans_real)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Commercial and Industrial Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



consumer_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//166_Real_ConsLoans_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

consumer_loans%>%str()

consumer_loans <- zoo(consumer_loans$`CONSUMER`, order.by = consumer_loans$`observation_date`)

consumer_loans_real <- consumer_loans/cpi

data_stock_and_watson <- merge(data_stock_and_watson,consumer_loans_real)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Consumer Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



pct_change_non_revolving_consumer_credit <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//167_Real_NonRevCredit_ESTÁ EN TASAS DE VARIACIÓN_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

pct_change_non_revolving_consumer_credit%>%str()

pct_change_non_revolving_consumer_credit <- zoo(pct_change_non_revolving_consumer_credit$`NONREVSLAR`, order.by = pct_change_non_revolving_consumer_credit$`observation_date`)

data_stock_and_watson <- merge(data_stock_and_watson,pct_change_non_revolving_consumer_credit)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Percent Change of Total Nonrevolving Consumer Credit, Percent Change at Annual Rate, Monthly, Seasonally Adjusted Annual Rate"


real_state_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//168_Real_LoansRealEst_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

real_state_loans%>%str()

real_state_loans <- zoo(real_state_loans$`REALLN`, order.by = real_state_loans$`observation_date`)

real_state_loans_real <- real_state_loans/cpi

data_stock_and_watson <- merge(data_stock_and_watson,real_state_loans_real)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Real Estate Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



revolving_credit_outstanding <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//169_Real_RevolvCredit_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 5)

revolving_credit_outstanding%>%str()

revolving_credit_outstanding <- zoo(revolving_credit_outstanding$`DTCTLR_N.M`, order.by = as.Date(paste0(revolving_credit_outstanding$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

revolving_credit_outstanding_real <- revolving_credit_outstanding/cpi

data_stock_and_watson <- merge(data_stock_and_watson,revolving_credit_outstanding_real)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "Total Revolving credit outstanding, not seasonally adjusted level, Real"



gs10 <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//GS10.xls", sheet = 1, skip = 11)

gs10%>%str()

df_monthly <- gs10 %>%
  group_by(Month = format(observation_date, "%Y-%m")) %>%
  summarise(Value = mean(DGS10_20240531, na.rm = TRUE))

# Convert Month column to Date format
df_monthly$Month <- as.Date(paste0(df_monthly$Month, "-01"))

# Create a zoo object
gs10_monthly <- zoo(df_monthly$Value, df_monthly$Month%>%as.POSIXct())



baa_bond_yield<-getSymbols('DBAA',src='FRED')

baa_bond_yield<-DBAA

baa_bond_yield<-aggregate(baa_bond_yield,as.yearmon,mean,
                          na.rm=T)
index(baa_bond_yield)<-as.yearmon(index(baa_bond_yield))

baa_bond_yield%>%str()

baa_bond_yield <- zoo(baa_bond_yield%>%coredata(), order.by = baa_bond_yield%>%index()%>%as.POSIXct())

baa_minus_gs10 <- baa_bond_yield-gs10_monthly

data_stock_and_watson <- merge(data_stock_and_watson,baa_minus_gs10)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "BAA Corporate bond minus GS10"



mortgage_30 <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//MORT.xls", sheet = 1, skip = 10)

mortgage_30%>%str()

mortgage_30$Month <- floor_date(mortgage_30$observation_date, "month")

df_monthly <- mortgage_30 %>%
  group_by(Month) %>%
  summarise(Value = mean(MORTGAGE30US, na.rm = TRUE))

# Create a zoo object
mortgage_30 <- zoo(df_monthly$Value, df_monthly$Month)


data_stock_and_watson <- merge(data_stock_and_watson,mortgage_30)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "30-Year Fixed Rate Mortgage Average in the United States, Percent, Weekly, Not Seasonally Adjusted"



mortgage_30_minus_gs10 <- mortgage_30-gs10_monthly

data_stock_and_watson <- merge(data_stock_and_watson,mortgage_30_minus_gs10)

colnames(data_stock_and_watson)[ncol(data_stock_and_watson)] <- "30-year US Mortgage Rate bond minus GS10"

