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

library(stringr)

setwd("C://Users//Italo//OneDrive//Documents//PhD//PhD//First")

#First run Stock_and_Watson_data.R

urate_job_losers <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//66_U Job Losers.xls", sheet = 1, skip = 11)

urate_job_losers <- zoo(urate_job_losers$`1.6000000000000001`, order.by = urate_job_losers$`24473`)

data_stock_and_watson_stationary <- urate_job_losers%>%diff()


ulevel_reentrants <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//67_U LF Reenty.xls", sheet = 1, skip = 11)

ulevel_reentrants <- zoo(ulevel_reentrants$`1098`, order.by = ulevel_reentrants$`24473`)

ulevel_reentrants <- ulevel_reentrants%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,ulevel_reentrants)

colnames(data_stock_and_watson_stationary) <- c("Unemployment Rate - Job Losers","Unemployment Level - Reentrants to Labor Force, Thousands of Persons, Monthly, Seasonally Adjusted")



u_job_leavers <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//68_U Job Leavers.xls", sheet = 1, skip = 11)

u_job_leavers <- zoo(u_job_leavers$`404`, order.by = u_job_leavers$`24473`)

u_job_leavers <- u_job_leavers%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,u_job_leavers)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Unemployment Level - Job Leavers, Thousands of Persons, Monthly, Seasonally Adjusted"



u_new_entrants <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//69_U New Entrants.xls", sheet = 1, skip = 11)

u_new_entrants <- zoo(u_new_entrants$`409`, order.by = u_new_entrants$`24473`)

u_new_entrants <- u_new_entrants%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,u_new_entrants)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Unemployment Level - New Entrants, Thousands of Persons, Monthly, Seasonally Adjusted"




emp_part_time <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//70_Emp_SlackWk.xls", sheet = 1, skip = 11)

emp_part_time <- zoo(emp_part_time$`2063`, order.by = emp_part_time$`20210`)

emp_part_time <- emp_part_time%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,emp_part_time)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Employment Level - Part-Time for Economic Reasons, All Industries, Thousands of Persons, Monthly, Seasonally Adjusted"



ppi_finished_consumer_goods <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//124_PPI_FinConsGds.xls", sheet = 1, skip = 11)

ppi_finished_consumer_goods%>%str()

ppi_finished_consumer_goods <- zoo(ppi_finished_consumer_goods$`28.300000000000001`, order.by = ppi_finished_consumer_goods$`17258`)

ppi_finished_consumer_goods <- ppi_finished_consumer_goods%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,ppi_finished_consumer_goods)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Producer Price Index by Commodity: Final Demand: Personal Consumption Goods (Finished Consumer Goods), Index 1982=100, Monthly, Seasonally Adjusted"


ppi_industrial_commodities <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//126_PPI_IndCom.xls", sheet = 1, skip = 11)

ppi_industrial_commodities%>%str()

ppi_industrial_commodities <- zoo(ppi_industrial_commodities$`12.300000000000001`, order.by = ppi_industrial_commodities$`4750`)

ppi_industrial_commodities <- ppi_industrial_commodities%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,ppi_industrial_commodities)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Producer Price Index by Commodity: Industrial Commodities, Index 1982=100, Monthly, Not Seasonally Adjusted"



ppi_natural_gas <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//131_Real_Price_NatGas_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 11)

ppi_natural_gas%>%str()

ppi_natural_gas <- zoo(ppi_natural_gas$`7.4000000000000004`, order.by = ppi_natural_gas$`24473`)



cpi<-getSymbols('CPIAUCSL',src='FRED')

cpi<-CPIAUCSL

cpi%>%str()

cpi <- zoo(cpi%>%coredata(), order.by = cpi%>%index()%>%as.POSIXct())

ppi_natural_gas_real <- ppi_natural_gas/cpi

ppi_natural_gas_real <- ppi_natural_gas_real%>%log()%>%diff()*100


data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,ppi_natural_gas_real)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Producer Price Index by Commodity: Fuels and Related Products and Power: Natural Gas, Index 1982=100, Monthly, Not Seasonally Adjusted, Deflated"



gs1_tb3m <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//154_GS1_Tb3m.xls", sheet = 1, skip = 5)

gs1_tb3m%>%str()

gs1_tb3m <- zoo(gs1_tb3m$`...2`, order.by = as.Date(paste0(gs1_tb3m$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

gs1_tb3m <- gs1_tb3m#%>%diff()

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,gs1_tb3m)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "GS1_Tb3m"


cp_tb3m <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//156_CP_Tbill Spread.xls", sheet = 1, skip = 5)

cp_tb3m%>%str()

cp_tb3m <- zoo(cp_tb3m$`...2`, order.by = as.Date(paste0(cp_tb3m$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

cp_tb3m <- cp_tb3m#%>%diff()

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,cp_tb3m)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "cp_tb3m"



commercial_and_industrial_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//165_Real_C&Lloand_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 11)

commercial_and_industrial_loans%>%str()

commercial_and_industrial_loans <- zoo(commercial_and_industrial_loans$`11.289400000000001`, order.by = commercial_and_industrial_loans$`17168`)

commercial_and_industrial_loans_real <- commercial_and_industrial_loans/cpi

commercial_and_industrial_loans_real <- commercial_and_industrial_loans_real%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,commercial_and_industrial_loans_real)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Commercial and Industrial Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



consumer_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//166_Real_ConsLoans_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

consumer_loans%>%str()

consumer_loans <- zoo(consumer_loans$`CONSUMER`, order.by = consumer_loans$`observation_date`)

consumer_loans_real <- consumer_loans/cpi

consumer_loans_real <- consumer_loans_real%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,consumer_loans_real)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Consumer Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



pct_change_non_revolving_consumer_credit <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//167_Real_NonRevCredit_ESTÁ EN TASAS DE VARIACIÓN_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

pct_change_non_revolving_consumer_credit%>%str()

pct_change_non_revolving_consumer_credit <- zoo(pct_change_non_revolving_consumer_credit$`NONREVSLAR`, order.by = pct_change_non_revolving_consumer_credit$`observation_date`)

pct_change_non_revolving_consumer_credit <- pct_change_non_revolving_consumer_credit#%>%diff()

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,pct_change_non_revolving_consumer_credit)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Percentage Change of Total Nonrevolving Consumer Credit,  at Annual Rate, Monthly, Seasonally Adjusted Annual Rate"



real_state_loans <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//168_Real_LoansRealEst_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 10)

real_state_loans%>%str()

real_state_loans <- zoo(real_state_loans$`REALLN`, order.by = real_state_loans$`observation_date`)

real_state_loans_real <- real_state_loans/cpi

real_state_loans_real <- real_state_loans_real%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,real_state_loans_real)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Real Estate Loans, All Commercial Banks, Billions of U.S. Dollars, Monthly, Seasonally Adjusted, Real"



revolving_credit_outstanding <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//Datos Stock y Watson//169_Real_RevolvCredit_HAY QUE DEFLACTARLO.xls", sheet = 1, skip = 5)

revolving_credit_outstanding%>%str()

revolving_credit_outstanding <- zoo(revolving_credit_outstanding$`DTCTLR_N.M`, order.by = as.Date(paste0(revolving_credit_outstanding$`Time Period`, "-01"), format = "%Y-%m-%d")%>%as.POSIXct())

revolving_credit_outstanding_real <- revolving_credit_outstanding/cpi

revolving_credit_outstanding_real <- revolving_credit_outstanding_real%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,revolving_credit_outstanding_real)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Total Revolving credit outstanding, not seasonally adjusted level, Real"


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

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,baa_minus_gs10)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "BAA Corporate bond minus GS10"


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

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,mortgage_30_minus_gs10)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "30-year US Mortgage Rate bond minus GS10"



bloomberg_data <- read.csv("data_bloomberg.csv",header = T)

bloomberg_variables <- bloomberg_data%>%pull(VARIABLE)%>%unique()%>%sort()

bloomberg_variables[grep("PMI", bloomberg_variables)]

bloomberg_data%>%filter(VARIABLE=="ISM Manufacturing PMI")

pmi <- bloomberg_data%>%filter(VARIABLE=="ISM Manufacturing PMI")


# Function to convert dates
convert_dates <- function(date_str) {
  date <- dmy(date_str)
  year <- year(date)
  if (year >= 2050 && year < 2099) {
    year <- year - 100 # Assume 20th century for years >= 50
  } else if (year < 2049) {
    year <- year  # Assume 21st century for years < 50
  }
  return(format(floor_date(update(date, year = year), "month"), "%Y-%m-%d"))
}

# Apply the function to the dates column
pmi <- pmi %>%
  mutate(dates = sapply(FECHA, convert_dates))



# Ensure dates column is in character format
pmi$dates <- as.character(pmi$dates)

# Convert dates to Date class with the correct format
pmi$dates <- as.Date(pmi$dates, format = "%Y-%m-%d")

# Create the zoo object
pmi <- zoo(pmi$PX_LAST, order.by = pmi$dates%>%as.POSIXct())


data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,pmi)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "US PMI Manufacturing"



av_api_key("QXJ8WVBH1O3U46U6")

#Get S&P 500 data using Alpha Vantage (SPY is commonly used for S&P 500)
getSymbols("SPY", src = "av", api.key = "QXJ8WVBH1O3U46U6", output.size = "full")

sp_500 <- (SPY$SPY.High+SPY$SPY.Low)/2

sp_500<-aggregate(sp_500,as.yearmon,mean,
                           na.rm=T)
index(sp_500)<-as.yearmon(index(sp_500))

sp_500 <- zoo(sp_500%>%coredata(), order.by = sp_500%>%index()%>%as.POSIXct())

sp_500 <- sp_500%>%log()%>%diff()*100

sp_500 <- read_excel("C://Users//Italo//OneDrive//Documents//PhD//PhD//First//SP500.xls", sheet = 1)

sp_500 <- sp_500$Close

sp_500 <- sp_500%>%ts(start = c(1985,1),frequency = 12)

sp_500 <- zoo(sp_500)

index(sp_500)<-as.yearmon(index(sp_500))

sp_500 <- zoo(sp_500%>%coredata(), order.by = sp_500%>%index()%>%as.POSIXct())

sp_500 <- sp_500%>%log()%>%diff()*100
  
data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,sp_500)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "S&P 500"

high_yield <-getSymbols('BAMLH0A0HYM2EY',src='FRED')
high_yield<-BAMLH0A0HYM2EY

high_yield<-aggregate(high_yield,as.yearmon,mean,
                  na.rm=T)
index(high_yield)<-as.yearmon(index(high_yield))

high_yield <- zoo(high_yield%>%coredata(), order.by = high_yield%>%index()%>%as.POSIXct())

high_yield <- high_yield%>%diff()

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,high_yield)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Yield on High-Yield Bonds"

bank_credit <-getSymbols('TOTBKCR',src='FRED')
bank_credit<-TOTBKCR

bank_credit<-aggregate(bank_credit,as.yearmon,sum,
                      na.rm=T)
index(bank_credit)<-as.yearmon(index(bank_credit))

bank_credit <- zoo(bank_credit%>%coredata(), order.by = bank_credit%>%index()%>%as.POSIXct())

bank_credit <- bank_credit%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,bank_credit)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Bank Credit"

#New

price_commodities <- read.csv("price_commodities.csv",  header = T)

price_commodities <- zoo(price_commodities$PX_LAST, order.by = price_commodities$FECHA%>%as.Date("%d/%m/%y"))

price_commodities<-aggregate(price_commodities,as.yearmon,mean,
                  na.rm=T)
index(price_commodities)<-as.yearmon(index(price_commodities))

price_commodities <- zoo(price_commodities%>%coredata(), order.by = price_commodities%>%index()%>%as.POSIXct())

price_commodities <- price_commodities%>%log()%>%diff()*100

data_stock_and_watson_stationary <- merge(data_stock_and_watson_stationary,price_commodities)

colnames(data_stock_and_watson_stationary)[ncol(data_stock_and_watson_stationary)] <- "Commodities Index"


#End New
data_stock_and_watson_stationary <- data_stock_and_watson_stationary[complete.cases(data_stock_and_watson_stationary),]

#write.csv(as.data.frame(data_stock_and_watson_stationary),"data_stock_and_watson_stationary_for_macro_factors_4.csv")
