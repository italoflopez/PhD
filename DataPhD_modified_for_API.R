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
library(fredr)
#bla

setwd("C://Users//Italo//OneDrive//Documents//PhD//PhD//First")

fredr_set_key("233ade3cce49bb63a9fd83f4fb88bec7")

ip_total <- fredr("INDPRO")

ip_total <- ip_total%>%select(date,value)

ip_final_nonind<-fredr('IPFPNSS')

ip_final_nonind <- ip_final_nonind%>%select(date,value)

data <- full_join(ip_total, ip_final_nonind, by = "date") %>%
  arrange(date)

colnames(data) <-c("date","Index of IP: total","Index of IP: final products and nonindustrial
supplies")

ip_final<-fredr('IPFINAL')

ip_final <- ip_final%>%select(date,value)

data <- full_join(data, ip_final, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)]<- "Index of IP: final products"

ip_cons<-fredr('IPCONGD')

ip_cons <- ip_cons%>%select(date,value)


data<-full_join(data,ip_cons, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)]<- "Index of IP: consumer goods"

ip_dur_cons <- fredr('IPDCONGD')

ip_dur_cons <- ip_dur_cons %>% select(date, value)

data <- full_join(data, ip_dur_cons, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: durable consumer goods"

ip_nondur_cons <- fredr('IPNCONGD')

ip_nondur_cons <- ip_nondur_cons %>% select(date, value)

data <- full_join(data, ip_nondur_cons, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: nondurable consumer
goods"

ip_bus_equip <- fredr('IPBUSEQ')

ip_bus_equip <- ip_bus_equip %>% select(date, value)

data <- full_join(data, ip_bus_equip, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: business equipment"

ip_mat <- fredr('IPMAT')

ip_mat <- ip_mat %>% select(date, value)

data <- full_join(data, ip_mat, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: materials"

ip_manu <- fredr('IPMAN')

ip_manu <- ip_manu %>% select(date, value)

data <- full_join(data, ip_manu, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: mfg"

ip_manu_dur <- fredr('IPDMAN')

ip_manu_dur <- ip_manu_dur %>% select(date, value)

data <- full_join(data, ip_manu_dur, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: mfg, durables"

ip_manu_nondur <- fredr('IPNMAN')

ip_manu_nondur <- ip_manu_nondur %>% select(date, value)

data <- full_join(data, ip_manu_nondur, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: mfg, nondurables"

ip_energy <- fredr('IPB50089S')

ip_energy <- ip_energy %>% select(date, value)

data <- full_join(data, ip_energy, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: energy, total"

ip_nonenergy <- fredr('IPX5001ES')

ip_nonenergy <- ip_nonenergy %>% select(date, value)

data <- full_join(data, ip_nonenergy, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: nonenergy, total"

ip_comp_equip_semi <- fredr('IPHITEK2S')

ip_comp_equip_semi <- ip_comp_equip_semi %>% select(date, value)

data <- full_join(data, ip_comp_equip_semi, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: computers, comm. equip.,
and semiconductors"

ip_nonenergy_excl <- fredr('IPX50EHVS')

ip_nonenergy_excl <- ip_nonenergy_excl %>% select(date, value)

data <- full_join(data, ip_nonenergy_excl, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Index of IP: nonenergy excl CCS and
MVP"

cap_util <- fredr('TCU')

cap_util <- cap_util %>% select(date, value)

data <- full_join(data, cap_util, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Capacity utilization: total"

cap_util_manu <- fredr('MCUMFN')

cap_util_manu <- cap_util_manu %>% select(date, value)

data <- full_join(data, cap_util_manu, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Capacity utilization: mfg"

cap_util_manu_dur <- fredr('CAPUTLGMFDS')

cap_util_manu_dur <- cap_util_manu_dur %>% select(date, value)

data <- full_join(data, cap_util_manu_dur, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Capacity utilization: mfg, durables"

cap_util_manu_nondur <- fredr('CAPUTLGMFNS')

cap_util_manu_nondur <- cap_util_manu_nondur %>% select(date, value)

data <- full_join(data, cap_util_manu_nondur, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Capacity utilization: mfg, nondurables"

cap_util_exc_ccs <- fredr('CAPUTLX50HTKS')

cap_util_exc_ccs <- cap_util_exc_ccs %>% select(date, value)

data <- full_join(data, cap_util_exc_ccs, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Capacity utilization: computers, comm. equip., and semiconductors"

employ_total <- fredr('PAYEMS')

employ_total <- employ_total %>% select(date, value)

data <- full_join(data, employ_total, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "CLF employed: total"

mean_duration_unemp <- fredr('UEMPMEAN')

mean_duration_unemp <- mean_duration_unemp %>% select(date, value)

data <- full_join(data, mean_duration_unemp, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Mean duration of unemployment"

unemp_5_weeks <- fredr('UEMPLT5')

unemp_5_weeks <- unemp_5_weeks %>% select(date, value)

data <- full_join(data, unemp_5_weeks, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Persons unemployed less than 5 weeks"

unemp_5_14_weeks <- fredr('UEMP5TO14')

unemp_5_14_weeks <- unemp_5_14_weeks %>% select(date, value)

data <- full_join(data, unemp_5_14_weeks, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Persons unemployed 5 to 14 weeks"

unemp_15_26_weeks <- fredr('UEMP15T26')

unemp_15_26_weeks <- unemp_15_26_weeks %>% select(date, value)

data <- full_join(data, unemp_15_26_weeks, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Persons unemployed 15 to 26 weeks"

unemp_plus_15_weeks <- fredr('UEMP15OV')

unemp_plus_15_weeks <- unemp_plus_15_weeks %>% select(date, value)

data <- full_join(data, unemp_plus_15_weeks, by = "date") %>%
  arrange(date)

colnames(data)[ncol(data)] <- "Persons unemployed 15+ weeks"

init_claims<-fredr('ICSA')

init_claims<-init_claims%>%select(date,value)
# Aggregate to monthly average
init_claims <- init_claims %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month) %>%
  summarise(value = sum(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)


data <- full_join(data, init_claims, by = "date") %>%
  arrange(date)
colnames(data)[ncol(data)]<- "Avg weekly initial claims"

employ_priv <- fredr('ADPWNUSNERSA')

employ_priv <- employ_priv %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
employ_priv <- employ_priv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_priv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: total private"
data%>%tail(20)

employ_mining <- fredr('CES1021000001')

employ_mining <- employ_mining %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_mining <- employ_mining %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_mining, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: mining"
data%>%tail(20)

employ_construction <- fredr('USCONS')

employ_construction <- employ_construction %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_construction <- employ_construction %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_construction, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: construction"
data%>%tail(20)

employ_manu <- fredr('MANEMP')

employ_manu <- employ_manu %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_manu <- employ_manu %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_manu, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: manufacturing"
data%>%tail(20)

employ_manu_dur <- fredr('CES3133900001')

employ_manu_dur <- employ_manu_dur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_manu_dur <- employ_manu_dur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_manu_dur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: manufacturing, durables"
data%>%tail(20)

employ_manu_nondur <- fredr('CES3232900001')

employ_manu_nondur <- employ_manu_nondur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_manu_nondur <- employ_manu_nondur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_manu_nondur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: manufacturing, nondurables"
data%>%tail(20)

employ_util <- fredr('CES4422000001')

employ_util <- employ_util %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_util <- employ_util %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_util, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: utilities"
data%>%tail(20)

employ_retail_trade <- fredr('USTRADE')

employ_retail_trade <- employ_retail_trade %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_retail_trade <- employ_retail_trade %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_retail_trade, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: retail trade"
data%>%tail(20)

employ_wholesale_trade <- fredr('USWTRADE')

employ_wholesale_trade <- employ_wholesale_trade %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_wholesale_trade <- employ_wholesale_trade %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_wholesale_trade, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: wholesale trade"
data%>%tail(20)

employ_fin <- fredr('USFIRE')

employ_fin <- employ_fin %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_fin <- employ_fin %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_fin, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: financial activities"
data%>%tail(20)

employ_pro_bus_serv <- fredr('USPBS')

employ_pro_bus_serv <- employ_pro_bus_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_pro_bus_serv <- employ_pro_bus_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_pro_bus_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: professional and business services"
data%>%tail(20)

employ_educ_health_serv <- fredr('USEHS')

employ_educ_health_serv <- employ_educ_health_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_educ_health_serv <- employ_educ_health_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_educ_health_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: education and health services"
data%>%tail(20)

employ_leisure <- fredr('USLAH')

employ_leisure <- employ_leisure %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_leisure <- employ_leisure %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_leisure, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: lesiure and hospitality"
data%>%tail(20)

employ_other_serv <- fredr('USSERV')

employ_other_serv <- employ_other_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
employ_other_serv <- employ_other_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, employ_other_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Employment on nonag payrolls: other services"
data%>%tail(20)

sales_manu_trade <- fredr('CMRMTSPL')

sales_manu_trade <- sales_manu_trade %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
sales_manu_trade <- sales_manu_trade %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, sales_manu_trade, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Manufacturing and Trade Industries Sales"
data%>%tail(20)

sales_manu <- fredr('MNFCTRSMSA')

sales_manu <- sales_manu %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
sales_manu <- sales_manu %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, sales_manu, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Manufacturers Sales"
data%>%tail(20)

sales_whole <- fredr('WHLSLRSMSA')

sales_whole <- sales_whole %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
sales_whole <- sales_whole %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, sales_whole, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Merchant Wholesalers Sales"
data%>%tail(20)

sales_retail <- fredr('MRTSSM44000USS')

sales_retail <- sales_retail %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
sales_retail <- sales_retail %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, sales_retail, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Retail Sales: Retail Trade"
data%>%tail(20)

real_pce <- fredr('PCEC96')

real_pce <- real_pce %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
real_pce <- real_pce %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, real_pce, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Personal Consumption Expenditures"
data%>%tail(20)

real_pce_dur <- fredr('PCEDGC96')

real_pce_dur <- real_pce_dur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
real_pce_dur <- real_pce_dur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, real_pce_dur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Personal Consumption Expenditures: Durable Goods"
data%>%tail(20)

real_pce_nondur <- fredr('PCENDC96')

real_pce_nondur <- real_pce_nondur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
real_pce_nondur <- real_pce_nondur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, real_pce_nondur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Personal Consumption Expenditures: Nondurable Goods"
data%>%tail(20)

real_pce_serv <- fredr('PCESC96')

real_pce_serv <- real_pce_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
real_pce_serv <- real_pce_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, real_pce_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Personal Consumption Expenditures: Services"
data%>%tail(20)

house_started <- fredr('HOUST')

house_started <- house_started %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
house_started <- house_started %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, house_started, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "New Privately-Owned Housing Units Started: Total Units"
data%>%tail(20)

house_permit <- fredr('PERMIT')

house_permit <- house_permit %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
house_permit <- house_permit %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, house_permit, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- " New Privately-Owned Housing Units Authorized in Permit-Issuing Places: Total Units"
data%>%tail(20)

house_sold <- fredr('HSN1F')

house_sold <- house_sold %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
house_sold <- house_sold %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, house_sold, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "New One Family Houses Sold: United States"
data%>%tail(20)

real_manu_trade_invent <- fredr('INVCMRMT')

real_manu_trade_invent <- real_manu_trade_invent %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
real_manu_trade_invent <- real_manu_trade_invent %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, real_manu_trade_invent, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Real Manufacturing and Trade Inventories"
data%>%tail(20)

dow_jones <- fredr('DJIA')

dow_jones <- dow_jones %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
dow_jones <- dow_jones %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, dow_jones, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Dow Jones Industrial Average"
data%>%tail(20)

s_p_500 <- fredr('SP500')

s_p_500 <- s_p_500 %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
s_p_500 <- s_p_500 %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, s_p_500, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "S&P 500"
data%>%tail(20)

dollar_euro_spot <- fredr('DEXUSEU')

dollar_euro_spot <- dollar_euro_spot %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
dollar_euro_spot <- dollar_euro_spot %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, dollar_euro_spot, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "U.S. Dollars to Euro Spot Exchange Rate"
data%>%tail(20)

swiss_franc_spot <- fredr('DEXSZUS')

swiss_franc_spot <- swiss_franc_spot %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
swiss_franc_spot <- swiss_franc_spot %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, swiss_franc_spot, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Swiss Francs to U.S. Dollar Spot Exchange Rate"
data%>%tail(20)

yen_spot <- fredr('DEXJPUS')

yen_spot <- yen_spot %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
yen_spot <- yen_spot %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, yen_spot, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Japanese Yen to U.S. Dollar Spot Exchange Rate"
data%>%tail(20)

pound_spot <- fredr('DEXUSUK')

pound_spot <- pound_spot %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
pound_spot <- pound_spot %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pound_spot, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "U.S. Dollars to U.K. Pound Sterling Spot Exchange Rate"
data%>%tail(20)

commercial_paper <- fredr('COMPOUT')

commercial_paper <- commercial_paper %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
commercial_paper <- commercial_paper %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, commercial_paper, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Commercial Paper Outstanding"
data%>%tail(20)

aaa_bond_yield <- fredr('DAAA')

aaa_bond_yield <- aaa_bond_yield %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
aaa_bond_yield <- aaa_bond_yield %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, aaa_bond_yield, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Moody's Seasoned Aaa Corporate Bond Yield"
data%>%tail(20)

baa_bond_yield <- fredr('DBAA')

baa_bond_yield <- baa_bond_yield %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
baa_bond_yield <- baa_bond_yield %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, baa_bond_yield, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Moody's Seasoned Baa Corporate Bond Yield"
data%>%tail(20)

m1 <- fredr('WM1NS')

m1 <- m1 %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
m1 <- m1 %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, m1, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "M1"
data%>%tail(20)

m2 <- fredr('WM2NS')

m2 <- m2 %>% select(date, value)
# Aggregate to monthly mean or sum depending on original function
m2 <- m2 %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, m2, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "M2"
data%>%tail(20)

m3 <- fredr('MABMM301USM189S')

m3 <- m3 %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
m3 <- m3 %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, m3, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "M3"
data%>%tail(20)

mb <- fredr('BOGMBASE')

mb <- mb %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
mb <- mb %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, mb, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Monetary Base; Total"
data%>%tail(20)

reser_depo_ins_total <- fredr('TOTRESNS')

reser_depo_ins_total <- reser_depo_ins_total %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
reser_depo_ins_total <- reser_depo_ins_total %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, reser_depo_ins_total, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Reserves of Depository Institutions: Total"
data%>%tail(20)

reser_depo_ins_nonborrowed <- fredr('NONBORRES')

reser_depo_ins_nonborrowed <- reser_depo_ins_nonborrowed %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
reser_depo_ins_nonborrowed <- reser_depo_ins_nonborrowed %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, reser_depo_ins_nonborrowed, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Reserves of Depository Institutions, Nonborrowed"
data%>%tail(20)

ppi_comm_finished <- fredr('WPSFD49207')

ppi_comm_finished <- ppi_comm_finished %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
ppi_comm_finished <- ppi_comm_finished %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, ppi_comm_finished, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Producer Price Index by Commodity: Final Demand: Finished Goods"
data%>%tail(20)

ppi_finished_less_food_energy <- fredr('WPSFD4131')

ppi_finished_less_food_energy <- ppi_finished_less_food_energy %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
ppi_finished_less_food_energy <- ppi_finished_less_food_energy %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, ppi_finished_less_food_energy, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Producer Price Index by Commodity: Final Demand: Finished Goods Less Foods and Energy"
data%>%tail(20)

cpi <- fredr('CPIAUCSL')

cpi <- cpi %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi <- cpi %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: All Items in U.S. City Average"
data%>%tail(20)

cpi_food_beverages <- fredr('CPIFABSL')

cpi_food_beverages <- cpi_food_beverages %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_food_beverages <- cpi_food_beverages %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_food_beverages, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Food and Beverages in U.S. City Average"
data%>%tail(20)

cpi_housing <- fredr('CPIHOSSL')

cpi_housing <- cpi_housing %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_housing <- cpi_housing %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_housing, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Housing in U.S. City Average"
data%>%tail(20)

cpi_apparel <- fredr('CPIAPPSL')

cpi_apparel <- cpi_apparel %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_apparel <- cpi_apparel %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_apparel, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Apparel in U.S. City Average"
data%>%tail(20)

cpi_transportation <- fredr('CPITRNSL')

cpi_transportation <- cpi_transportation %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_transportation <- cpi_transportation %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_transportation, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Transportation in U.S. City Average"
data%>%tail(20)

cpi_med_care <- fredr('CPIMEDSL')

cpi_med_care <- cpi_med_care %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_med_care <- cpi_med_care %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_med_care, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Medical Care in U.S. City Average"
data%>%tail(20)

cpi_comm <- fredr('CUSR0000SAC')

cpi_comm <- cpi_comm %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_comm <- cpi_comm %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_comm, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Commodities in U.S. City Average"
data%>%tail(20)

cpi_serv <- fredr('CUSR0000SAS')

cpi_serv <- cpi_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_serv <- cpi_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: Services in U.S. City Average"
data%>%tail(20)

cpi_all_less_food <- fredr('CPIULFSL')

cpi_all_less_food <- cpi_all_less_food %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_all_less_food <- cpi_all_less_food %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_all_less_food, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: All Items Less Food in U.S. City Average"
data%>%tail(20)

cpi_all_less_shelter <- fredr('CUSR0000SA0L2')

cpi_all_less_shelter <- cpi_all_less_shelter %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_all_less_shelter <- cpi_all_less_shelter %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_all_less_shelter, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: All Items Less Shelter in U.S. City Average"
data%>%tail(20)

cpi_all_less_med <- fredr('CUSR0000SA0L5')

cpi_all_less_med <- cpi_all_less_med %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_all_less_med <- cpi_all_less_med %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_all_less_med, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: All Items Less Medical Care in U.S. City Average"
data%>%tail(20)

cpi_all_less_food_and_energy <- fredr('CPILFESL')

cpi_all_less_food_and_energy <- cpi_all_less_food_and_energy %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
cpi_all_less_food_and_energy <- cpi_all_less_food_and_energy %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, cpi_all_less_food_and_energy, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Consumer Price Index for All Urban Consumers: All Items Less Food and Energy in U.S. City Average"
data%>%tail(20)

pce <- fredr('PCEPI')

pce <- pce %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
pce <- pce %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pce, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Personal Consumption Expenditures: Chain-type Price Index"
data%>%tail(20)

pce_less_food_and_energy <- fredr('PCEPILFE')

pce_less_food_and_energy <- pce_less_food_and_energy %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
pce_less_food_and_energy <- pce_less_food_and_energy %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pce_less_food_and_energy, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Personal Consumption Expenditures Excluding Food and Energy (Chain-Type Price Index)"
data%>%tail(20)

pce_dur <- fredr('DDURRG3M086SBEA')

pce_dur <- pce_dur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
pce_dur <- pce_dur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pce_dur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Personal consumption expenditures: Durable goods (chain-type price index)"
data%>%tail(20)

pce_nondur <- fredr('DNDGRG3M086SBEA')

pce_nondur <- pce_nondur %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
pce_nondur <- pce_nondur %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pce_nondur, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Personal consumption expenditures: Nondurable goods (chain-type price index)"
data%>%tail(20)

pce_serv <- fredr('DSERRG3M086SBEA')

pce_serv <- pce_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
pce_serv <- pce_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, pce_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Personal consumption expenditures: Services (chain-type price index)"
data%>%tail(20)

avg_hour_earn_total_priv <- fredr('CES0500000003')

avg_hour_earn_total_priv <- avg_hour_earn_total_priv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_total_priv <- avg_hour_earn_total_priv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_total_priv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Total Private"
data%>%tail(20)

avg_hour_earn_constr <- fredr('CES2000000003')

avg_hour_earn_constr <- avg_hour_earn_constr %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_constr <- avg_hour_earn_constr %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_constr, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Construction"
data%>%tail(20)

avg_hour_earn_manu <- fredr('CES3000000003')

avg_hour_earn_manu <- avg_hour_earn_manu %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_manu <- avg_hour_earn_manu %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_manu, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Manufacturing"
data%>%tail(20)

avg_hour_earn_prof_and_bus_serv <- fredr('CES6000000003')

avg_hour_earn_prof_and_bus_serv <- avg_hour_earn_prof_and_bus_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_prof_and_bus_serv <- avg_hour_earn_prof_and_bus_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_prof_and_bus_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Professional and Business Services"
data%>%tail(20)

avg_hour_earn_educ_and_health_serv <- fredr('CES6500000003')

avg_hour_earn_educ_and_health_serv <- avg_hour_earn_educ_and_health_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_educ_and_health_serv <- avg_hour_earn_educ_and_health_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_educ_and_health_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Education and Health Services"
data%>%tail(20)

avg_hour_earn_other_serv <- fredr('CES8000000003')

avg_hour_earn_other_serv <- avg_hour_earn_other_serv %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_other_serv <- avg_hour_earn_other_serv %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_other_serv, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of All Employees, Other Services"
data%>%tail(20)

federal_surplus <- fredr('MTSDS133FMS')

federal_surplus <- federal_surplus %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
federal_surplus <- federal_surplus %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, federal_surplus, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Federal Surplus or Deficit [-]"
data%>%tail(20)

#New
avg_hour_earn_production_and_nonsupervisory <- fredr('AHETPI')

avg_hour_earn_production_and_nonsupervisory <- avg_hour_earn_production_and_nonsupervisory %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_production_and_nonsupervisory <- avg_hour_earn_production_and_nonsupervisory %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_production_and_nonsupervisory, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of Production and Nonsupervisory Employees, Total Private"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_manu <- fredr('CES3000000008')

avg_hour_earn_production_and_nonsupervisory_manu <- avg_hour_earn_production_and_nonsupervisory_manu %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_production_and_nonsupervisory_manu <- avg_hour_earn_production_and_nonsupervisory_manu %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_production_and_nonsupervisory_manu, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of Production and Nonsupervisory Employees, Manufacturing"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_durables <- fredr('CES3100000008')

avg_hour_earn_production_and_nonsupervisory_durables <- avg_hour_earn_production_and_nonsupervisory_durables %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_production_and_nonsupervisory_durables <- avg_hour_earn_production_and_nonsupervisory_durables %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_production_and_nonsupervisory_durables, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of Production and Nonsupervisory Employees, Durable Goods"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_construction <- fredr('CES2000000008')

avg_hour_earn_production_and_nonsupervisory_construction <- avg_hour_earn_production_and_nonsupervisory_construction %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
avg_hour_earn_production_and_nonsupervisory_construction <- avg_hour_earn_production_and_nonsupervisory_construction %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, avg_hour_earn_production_and_nonsupervisory_construction, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Average Hourly Earnings of Production and Nonsupervisory Employees, Construction"
data%>%tail(20)


# consumer_credit <- fredr('TOTALSL')

# consumer_credit <- consumer_credit %>% select(date, value)
# # Convert to monthly frequency (interpolated if necessary)
# consumer_credit <- consumer_credit %>%
#   mutate(year_month = floor_date(date, 'month')) %>%
#   group_by(year_month) %>%
#   summarise(value = mean(value, na.rm = FALSE)) %>%
#   mutate(date = as.Date(year_month)) %>%
#   select(date, value)
# 
# data <- full_join(data, consumer_credit, by = 'date') %>%
#   arrange(date)
# colnames(data)[ncol(data)] <- "Total Consumer Credit Owned and Securitized"
# data%>%tail(20)




# rate_credit_card <- fredr('TERMCBCCALLNS')

# rate_credit_card <- rate_credit_card %>% select(date, value)
# # Convert to monthly frequency (interpolated if necessary)
# rate_credit_card <- rate_credit_card %>%
#   mutate(year_month = floor_date(date, 'month')) %>%
#   group_by(year_month) %>%
#   summarise(value = mean(value, na.rm = FALSE)) %>%
#   mutate(date = as.Date(year_month)) %>%
#   select(date, value)
# 
# data <- full_join(data, rate_credit_card, by = 'date') %>%
#   arrange(date)
# colnames(data)[ncol(data)] <- "Commercial Bank Interest Rate on Credit Card Plans, All Accounts"
# data%>%tail(20)
# 
# 
rate_3_month_cd <- fredr('IR3TCD01USM156N')

rate_3_month_cd <- rate_3_month_cd %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
rate_3_month_cd <- rate_3_month_cd %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, rate_3_month_cd, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Interest Rates: 3-Month or 90-Day Rates and Yields: Certificates of Deposit: Total for United States"
data%>%tail(20)


rate_personal_loans_24_months <- fredr('TERMCBPER24NS')

rate_personal_loans_24_months <- rate_personal_loans_24_months %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
rate_personal_loans_24_months <- rate_personal_loans_24_months %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, rate_personal_loans_24_months, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Finance Rate on Personal Loans at Commercial Banks, 24 Month Loan"
data%>%tail(20)


unemployment_rate <- fredr('UNRATE')

unemployment_rate <- unemployment_rate %>% select(date, value)
# Convert to monthly frequency (interpolated if necessary)
unemployment_rate <- unemployment_rate %>%
  mutate(year_month = floor_date(date, 'month')) %>%
  group_by(year_month) %>%
  summarise(value = mean(value, na.rm = FALSE)) %>%
  mutate(date = as.Date(year_month)) %>%
  select(date, value)

data <- full_join(data, unemployment_rate, by = 'date') %>%
  arrange(date)
colnames(data)[ncol(data)] <- "Unemployment Rate"
data%>%tail(20)


# write.csv(as.data.frame(data),"data_level_FRED_25-10-2025.csv")

library(ggforce)

data_long <- data %>%
  pivot_longer(
    cols = -date, 
    names_to = "variable", 
    values_to = "value"
  )

data_long <-data_long[complete.cases(data_long),]

# p <- ggplot(data_long, aes(x = date, y = value)) +
#   geom_line(color = "steelblue", linewidth = 0.4) +
#   facet_wrap_paginate(~ variable, ncol = 6, nrow = 6, scales = "free_y") +
#   scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
#   labs(
#     title = "Multivariate Time Series (Monthly)",
#     subtitle = "Each panel shows one variable",
#     x = NULL, y = NULL
#   ) +
#   theme_minimal(base_size = 10) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# ggplotly(p)
# 
# facet <- p + facet_wrap_paginate(~ variable, ncol = 6, nrow = 6, scales = "free_y", page = 1)
# 
# facet
# 
# ggplotly(facet)
# 
# 
# vars <- unique(data_long$variable)
# n_per_page <- 36
# n_pages <- ceiling(length(vars) / n_per_page)
# n_pages
# 
# 
# page <- 1
# vars_page <- vars[((page - 1) * n_per_page + 1):min(page * n_per_page, length(vars))]
# data_page <- dplyr::filter(data_long, variable %in% vars_page)
# 
# 
# p_page <- ggplot(data_page, aes(x = date, y = value)) +
#   geom_line(color = "steelblue", linewidth = 0.4) +
#   facet_wrap(~ variable, ncol = 6, nrow = 6, scales = "free_y") +
#   theme_minimal(base_size = 10) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# ggplotly(p_page)