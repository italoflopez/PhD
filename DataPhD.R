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

options(download.file.method = "libcurl")
ip_total <- getSymbols("INDPRO", src = "FRED", auto.assign = FALSE)

ip_total<-INDPRO

ip_final_nonind<-getSymbols('IPFPNSS',src='FRED')
ip_final_nonind<-IPFPNSS
data <- merge(ip_total,ip_final_nonind)
colnames(data) <-c("Index of IP: total","Index of IP: final products and nonindustrial
supplies")

ip_final<-getSymbols('IPFINAL',src='FRED')
ip_final<-IPFINAL
data <- merge(data,ip_final)
colnames(data)[ncol(data)]<- "Index of IP: final products"

ip_cons<-getSymbols('IPCONGD',src='FRED')
ip_cons<-IPCONGD
data<-merge(data,ip_cons)
colnames(data)[ncol(data)]<- "Index of IP: consumer goods"

ip_dur_cons<-getSymbols('IPDCONGD',src='FRED')
ip_dur_cons<-IPDCONGD
data<-merge(data,ip_dur_cons)
colnames(data)[ncol(data)]<- "Index of IP: durable consumer goods"

ip_nondur_cons<-getSymbols('IPNCONGD',src='FRED')
ip_nondur_cons<-IPNCONGD
data<-merge(data,ip_nondur_cons)
colnames(data)[ncol(data)]<- "Index of IP: nondurable consumer
goods"

ip_bus_equip<-getSymbols('IPBUSEQ',src='FRED')
ip_bus_equip<-IPBUSEQ
data<-merge(data,ip_bus_equip)
colnames(data)[ncol(data)]<- "Index of IP: business equipment"

ip_mat<-getSymbols('IPMAT',src='FRED')
ip_mat<-IPMAT
data<-merge(data,ip_mat)
colnames(data)[ncol(data)]<- "Index of IP: materials"

ip_manu<-getSymbols('IPMAN',src='FRED')
ip_manu<-IPMAN
data<-merge(data,ip_manu)
colnames(data)[ncol(data)]<- "Index of IP: mfg"

ip_manu_dur<-getSymbols('IPDMAN',src='FRED')
ip_manu_dur<-IPDMAN
data<-merge(data,ip_manu_dur)
colnames(data)[ncol(data)]<- "Index of IP: mfg, durables"

ip_manu_nondur<-getSymbols('IPNMAN',src='FRED')
ip_manu_nondur<-IPNMAN
data<-merge(data,ip_manu_nondur)
colnames(data)[ncol(data)]<- "Index of IP: mfg, nondurables"

ip_energy<-getSymbols('IPB50089S',src='FRED')
ip_energy<-IPB50089S
data<-merge(data,ip_energy)
colnames(data)[ncol(data)]<- "Index of IP: energy, total"

ip_nonenergy<-getSymbols('IPX5001ES',src='FRED')
ip_nonenergy<-IPX5001ES
data<-merge(data,ip_nonenergy)
colnames(data)[ncol(data)]<- "Index of IP: nonenergy, total"

ip_comp_equip_semi<-getSymbols('IPHITEK2S',src='FRED')
ip_comp_equip_semi<-IPHITEK2S
data<-merge(data,ip_comp_equip_semi)
colnames(data)[ncol(data)]<- "Index of IP: computers, comm. equip.,
and semiconductors"

ip_nonenergy_excl<-getSymbols('IPX50EHVS',src='FRED')
ip_nonenergy_excl<-IPX50EHVS
data<-merge(data,ip_nonenergy_excl)
colnames(data)[ncol(data)]<- "Index of IP: nonenergy excl CCS and
MVP"

cap_util<-getSymbols('TCU',src='FRED')
cap_util<-TCU
data<-merge(data,cap_util)
colnames(data)[ncol(data)]<- "Capacity utilization: total"

cap_util_manu<-getSymbols('MCUMFN',src='FRED')
cap_util_manu<-MCUMFN
data<-merge(data,cap_util_manu)
colnames(data)[ncol(data)]<- "Capacity utilization: mfg"

cap_util_manu_dur<-getSymbols('CAPUTLGMFDS',src='FRED')
cap_util_manu_dur<-CAPUTLGMFDS
data<-merge(data,cap_util_manu_dur)
colnames(data)[ncol(data)]<- "Capacity utilization: mfg, durables"

cap_util_manu_nondur<-getSymbols('CAPUTLGMFNS',src='FRED')
cap_util_manu_nondur<-CAPUTLGMFNS
data<-merge(data,cap_util_manu_nondur)
colnames(data)[ncol(data)]<- "Capacity utilization: mfg, nondurables"

cap_util_exc_ccs<-getSymbols('CAPUTLX50HTKS',src='FRED')
cap_util_exc_ccs<-CAPUTLX50HTKS
data<-merge(data,cap_util_exc_ccs)
colnames(data)[ncol(data)]<- "Capacity utilization: computers, comm. equip., and semiconductors"

employ_total<-getSymbols('PAYEMS',src='FRED')
employ_total<-PAYEMS
data<-merge(data,employ_total)
colnames(data)[ncol(data)]<- "CLF employed: total"

mean_duration_unemp<-getSymbols('UEMPMEAN',src='FRED')
mean_duration_unemp<-UEMPMEAN
data<-merge(data,mean_duration_unemp)
colnames(data)[ncol(data)]<- "Mean duration of unemployment"

unemp_5_weeks<-getSymbols('UEMPLT5',src='FRED')
unemp_5_weeks<-UEMPLT5
data<-merge(data,unemp_5_weeks)
colnames(data)[ncol(data)]<- "Persons unemployed less than 5 weeks"

unemp_5_14_weeks<-getSymbols('UEMP5TO14',src='FRED')
unemp_5_14_weeks<-UEMP5TO14
data<-merge(data,unemp_5_14_weeks)
colnames(data)[ncol(data)]<- "Persons unemployed 5 to 14 weeks"

unemp_15_26_weeks<-getSymbols('UEMP15T26',src='FRED')
unemp_15_26_weeks<-UEMP15T26
data<-merge(data,unemp_15_26_weeks)
colnames(data)[ncol(data)]<- "Persons unemployed 15 to 26 weeks"

unemp_plus_15_weeks<-getSymbols('UEMP15OV',src='FRED')
unemp_plus_15_weeks<-UEMP15OV
data<-merge(data,unemp_plus_15_weeks)
colnames(data)[ncol(data)]<- "Persons unemployed 15+ weeks"

init_claims<-getSymbols('ICSA',src='FRED')
init_claims<-ICSA
init_claims<-aggregate(init_claims,as.yearmon,mean,
          na.rm=T)
index(data)<-as.yearmon(index(data))
data<-merge(data,init_claims)
colnames(data)[ncol(data)]<- "Avg weekly initial claims"

employ_priv<-getSymbols('ADPWNUSNERSA',src='FRED')
employ_priv<-ADPWNUSNERSA
employ_priv<-aggregate(employ_priv,as.yearmon,mean,
                       na.rm=T)
data<-merge(data,employ_priv)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: total private"
data%>%tail(20)

employ_mining<-getSymbols('CES1021000001',src='FRED')
employ_mining<-CES1021000001
index(employ_mining)<-as.yearmon(index(employ_mining))
data<-merge(data,employ_mining)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: mining"
data%>%tail(20)

employ_construction<-getSymbols('USCONS',src='FRED')
employ_construction<-USCONS
index(employ_construction)<-as.yearmon(index(employ_construction))
data<-merge(data,employ_construction)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: construction"
data%>%tail(20)

employ_manu<-getSymbols('MANEMP',src='FRED')
employ_manu<-MANEMP
index(employ_manu)<-as.yearmon(index(employ_manu))
data<-merge(data,employ_manu)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: manufacturing"
data%>%tail(20)

employ_manu_dur<-getSymbols('CES3133900001',src='FRED')
employ_manu_dur<-CES3133900001
index(employ_manu_dur)<-as.yearmon(index(employ_manu_dur))
data<-merge(data,employ_manu_dur)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: manufacturing, durables"
data%>%tail(20)

employ_manu_nondur<-getSymbols('CES3232900001',src='FRED')
employ_manu_nondur<-CES3232900001
index(employ_manu_nondur)<-as.yearmon(index(employ_manu_nondur))
data<-merge(data,employ_manu_nondur)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: manufacturing, nondurables"
data%>%tail(20)

employ_util<-getSymbols('CES4422000001',src='FRED')
employ_util<-CES4422000001
index(employ_util)<-as.yearmon(index(employ_util))
data<-merge(data,employ_util)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: utilities"
data%>%tail(20)

employ_retail_trade<-getSymbols('USTRADE',src='FRED')
employ_retail_trade<-USTRADE
index(employ_retail_trade)<-as.yearmon(index(employ_retail_trade))
data<-merge(data,employ_retail_trade)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: retail trade"
data%>%tail(20)

employ_wholesale_trade<-getSymbols('USWTRADE',src='FRED')
employ_wholesale_trade<-USWTRADE
index(employ_wholesale_trade)<-as.yearmon(index(employ_wholesale_trade))
data<-merge(data,employ_wholesale_trade)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: wholesale trade"
data%>%tail(20)

employ_fin<-getSymbols('USFIRE',src='FRED')
employ_fin<-USFIRE
index(employ_fin)<-as.yearmon(index(employ_fin))
data<-merge(data,employ_fin)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: financial activities"
data%>%tail(20)

employ_pro_bus_serv<-getSymbols('USPBS',src='FRED')
employ_pro_bus_serv<-USPBS
index(employ_pro_bus_serv)<-as.yearmon(index(employ_pro_bus_serv))
data<-merge(data,employ_pro_bus_serv)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: professional and business services"
data%>%tail(20)

employ_educ_health_serv<-getSymbols('USEHS',src='FRED')
employ_educ_health_serv<-USEHS
index(employ_educ_health_serv)<-as.yearmon(index(employ_educ_health_serv))
data<-merge(data,employ_educ_health_serv)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: education and health services"
data%>%tail(20)

employ_leisure<-getSymbols('USLAH',src='FRED')
employ_leisure<-USLAH
index(employ_leisure)<-as.yearmon(index(employ_leisure))
data<-merge(data,employ_leisure)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: lesiure and hospitality"
data%>%tail(20)

employ_other_serv<-getSymbols('USSERV',src='FRED')
employ_other_serv<-USSERV
index(employ_other_serv)<-as.yearmon(index(employ_other_serv))
data<-merge(data,employ_other_serv)
colnames(data)[ncol(data)]<- "Employment on nonag payrolls: other services"
data%>%tail(20)

sales_manu_trade<-getSymbols('CMRMTSPL',src='FRED')
sales_manu_trade<-CMRMTSPL
index(sales_manu_trade)<-as.yearmon(index(sales_manu_trade))
data<-merge(data,sales_manu_trade)
colnames(data)[ncol(data)]<- "Real Manufacturing and Trade Industries Sales"
data%>%tail(20)

sales_manu<-getSymbols('MNFCTRSMSA',src='FRED')
sales_manu<-MNFCTRSMSA
index(sales_manu)<-as.yearmon(index(sales_manu))
data<-merge(data,sales_manu)
colnames(data)[ncol(data)]<- "Manufacturers Sales"
data%>%tail(20)

sales_whole<-getSymbols('WHLSLRSMSA',src='FRED')
sales_whole<-WHLSLRSMSA
index(sales_whole)<-as.yearmon(index(sales_whole))
data<-merge(data,sales_whole)
colnames(data)[ncol(data)]<- "Merchant Wholesalers Sales"
data%>%tail(20)

sales_retail<-getSymbols('MRTSSM44000USS',src='FRED')
sales_retail<-MRTSSM44000USS
index(sales_retail)<-as.yearmon(index(sales_retail))
data<-merge(data,sales_retail)
colnames(data)[ncol(data)]<- "Retail Sales: Retail Trade"
data%>%tail(20)

real_pce<-getSymbols('PCEC96',src='FRED')
real_pce<-PCEC96
index(real_pce)<-as.yearmon(index(real_pce))
data<-merge(data,real_pce)
colnames(data)[ncol(data)]<- "Real Personal Consumption Expenditures"
data%>%tail(20)

real_pce_dur<-getSymbols('PCEDGC96',src='FRED')
real_pce_dur<-PCEDGC96
index(real_pce_dur)<-as.yearmon(index(real_pce_dur))
data<-merge(data,real_pce_dur)
colnames(data)[ncol(data)]<- "Real Personal Consumption Expenditures: Durable Goods"
data%>%tail(20)

real_pce_nondur<-getSymbols('PCENDC96',src='FRED')
real_pce_nondur<-PCENDC96
index(real_pce_nondur)<-as.yearmon(index(real_pce_nondur))
data<-merge(data,real_pce_nondur)
colnames(data)[ncol(data)]<- "Real Personal Consumption Expenditures: Nondurable Goods"
data%>%tail(20)

real_pce_serv<-getSymbols('PCESC96',src='FRED')
real_pce_serv<-PCESC96
index(real_pce_serv)<-as.yearmon(index(real_pce_serv))
data<-merge(data,real_pce_serv)
colnames(data)[ncol(data)]<- "Real Personal Consumption Expenditures: Services"
data%>%tail(20)

house_started<-getSymbols('HOUST',src='FRED')
house_started<-HOUST
index(house_started)<-as.yearmon(index(house_started))
data<-merge(data,house_started)
colnames(data)[ncol(data)]<- "New Privately-Owned Housing Units Started: Total Units"
data%>%tail(20)

house_permit<-getSymbols('PERMIT',src='FRED')
house_permit<-PERMIT
index(house_permit)<-as.yearmon(index(house_permit))
data<-merge(data,house_permit)
colnames(data)[ncol(data)]<- " New Privately-Owned Housing Units Authorized in Permit-Issuing Places: Total Units"
data%>%tail(20)

house_sold<-getSymbols('HSN1F',src='FRED')
house_sold<-HSN1F
index(house_sold)<-as.yearmon(index(house_sold))
data<-merge(data,house_sold)
colnames(data)[ncol(data)]<- "New One Family Houses Sold: United States"
data%>%tail(20)

real_manu_trade_invent<-getSymbols('INVCMRMT',src='FRED')
real_manu_trade_invent<-INVCMRMT
index(real_manu_trade_invent)<-as.yearmon(index(real_manu_trade_invent))
data<-merge(data,real_manu_trade_invent)
colnames(data)[ncol(data)]<- "Real Manufacturing and Trade Inventories"
data%>%tail(20)

dow_jones<-getSymbols('DJIA',src='FRED')
dow_jones<-DJIA
dow_jones<-aggregate(dow_jones,as.yearmon,mean,
                       na.rm=T)
index(dow_jones)<-as.yearmon(index(dow_jones))
data<-merge(data,dow_jones)
colnames(data)[ncol(data)]<- "Dow Jones Industrial Average"
data%>%tail(20)

s_p_500<-getSymbols('SP500',src='FRED')
s_p_500<-SP500
s_p_500<-aggregate(s_p_500,as.yearmon,mean,
                     na.rm=T)
index(s_p_500)<-as.yearmon(index(s_p_500))
data<-merge(data,s_p_500)
colnames(data)[ncol(data)]<- "S&P 500"
data%>%tail(20)

dollar_euro_spot<-getSymbols('DEXUSEU',src='FRED')
dollar_euro_spot<-DEXUSEU
dollar_euro_spot<-aggregate(dollar_euro_spot,as.yearmon,mean,
                   na.rm=T)
index(dollar_euro_spot)<-as.yearmon(index(dollar_euro_spot))
data<-merge(data,dollar_euro_spot)
colnames(data)[ncol(data)]<- "U.S. Dollars to Euro Spot Exchange Rate"
data%>%tail(20)

swiss_franc_spot<-getSymbols('DEXSZUS',src='FRED')
swiss_franc_spot<-DEXSZUS
swiss_franc_spot<-aggregate(swiss_franc_spot,as.yearmon,mean,
                           na.rm=T)
index(swiss_franc_spot)<-as.yearmon(index(swiss_franc_spot))
data<-merge(data,swiss_franc_spot)
colnames(data)[ncol(data)]<- "Swiss Francs to U.S. Dollar Spot Exchange Rate"
data%>%tail(20)

yen_spot<-getSymbols('DEXJPUS',src='FRED')
yen_spot<-DEXJPUS
yen_spot<-aggregate(yen_spot,as.yearmon,mean,
                            na.rm=T)
index(yen_spot)<-as.yearmon(index(yen_spot))
data<-merge(data,yen_spot)
colnames(data)[ncol(data)]<- "Japanese Yen to U.S. Dollar Spot Exchange Rate"
data%>%tail(20)

pound_spot<-getSymbols('DEXUSUK',src='FRED')
pound_spot<-DEXUSUK
pound_spot<-aggregate(pound_spot,as.yearmon,mean,
                    na.rm=T)
index(pound_spot)<-as.yearmon(index(pound_spot))
data<-merge(data,pound_spot)
colnames(data)[ncol(data)]<- "U.S. Dollars to U.K. Pound Sterling Spot Exchange Rate"
data%>%tail(20)

commercial_paper<-getSymbols('COMPOUT',src='FRED')
commercial_paper<-COMPOUT
commercial_paper<-aggregate(commercial_paper,as.yearmon,mean,
                      na.rm=T)
index(commercial_paper)<-as.yearmon(index(commercial_paper))
data<-merge(data,commercial_paper)
colnames(data)[ncol(data)]<- "Commercial Paper Outstanding"
data%>%tail(20)

aaa_bond_yield<-getSymbols('DAAA',src='FRED')
aaa_bond_yield<-DAAA
aaa_bond_yield<-aggregate(aaa_bond_yield,as.yearmon,mean,
                            na.rm=T)
index(aaa_bond_yield)<-as.yearmon(index(aaa_bond_yield))
data<-merge(data,aaa_bond_yield)
colnames(data)[ncol(data)]<- "Moody's Seasoned Aaa Corporate Bond Yield"
data%>%tail(20)

baa_bond_yield<-getSymbols('DBAA',src='FRED')
baa_bond_yield<-DBAA
baa_bond_yield<-aggregate(baa_bond_yield,as.yearmon,mean,
                          na.rm=T)
index(baa_bond_yield)<-as.yearmon(index(baa_bond_yield))
data<-merge(data,baa_bond_yield)
colnames(data)[ncol(data)]<- "Moody's Seasoned Baa Corporate Bond Yield"
data%>%tail(20)

m1<-getSymbols('WM1NS',src='FRED')
m1<-WM1NS
m1<-aggregate(m1,as.yearmon,mean,
                          na.rm=T)
index(m1)<-as.yearmon(index(m1))
data<-merge(data,m1)
colnames(data)[ncol(data)]<- "M1"
data%>%tail(20)

m2<-getSymbols('WM2NS',src='FRED')
m2<-WM2NS
m2<-aggregate(m2,as.yearmon,mean,
              na.rm=T)
index(m2)<-as.yearmon(index(m2))
data<-merge(data,m2)
colnames(data)[ncol(data)]<- "M2"
data%>%tail(20)

m3<-getSymbols('MABMM301USM189S',src='FRED')
m3<-MABMM301USM189S
index(m3)<-as.yearmon(index(m3))
data<-merge(data,m3)
colnames(data)[ncol(data)]<- "M3"
data%>%tail(20)

mb<-getSymbols('BOGMBASE',src='FRED')
mb<-BOGMBASE
index(mb)<-as.yearmon(index(mb))
data<-merge(data,mb)
colnames(data)[ncol(data)]<- "Monetary Base; Total"
data%>%tail(20)

reser_depo_ins_total<-getSymbols('TOTRESNS',src='FRED')
reser_depo_ins_total<-TOTRESNS
index(reser_depo_ins_total)<-as.yearmon(index(reser_depo_ins_total))
data<-merge(data,reser_depo_ins_total)
colnames(data)[ncol(data)]<- "Reserves of Depository Institutions: Total"
data%>%tail(20)

reser_depo_ins_nonborrowed<-getSymbols('NONBORRES',src='FRED')
reser_depo_ins_nonborrowed<-NONBORRES
index(reser_depo_ins_nonborrowed)<-as.yearmon(index(reser_depo_ins_nonborrowed))
data<-merge(data,reser_depo_ins_nonborrowed)
colnames(data)[ncol(data)]<- "Reserves of Depository Institutions, Nonborrowed"
data%>%tail(20)

ppi_comm_finished<-getSymbols('WPSFD49207',src='FRED')
ppi_comm_finished<-WPSFD49207
index(ppi_comm_finished)<-as.yearmon(index(ppi_comm_finished))
data<-merge(data,ppi_comm_finished)
colnames(data)[ncol(data)]<- "Producer Price Index by Commodity: Final Demand: Finished Goods"
data%>%tail(20)

ppi_finished_less_food_energy<-getSymbols('WPSFD4131',src='FRED')
ppi_finished_less_food_energy<-WPSFD4131
index(ppi_finished_less_food_energy)<-as.yearmon(index(ppi_finished_less_food_energy))
data<-merge(data,ppi_finished_less_food_energy)
colnames(data)[ncol(data)]<- "Producer Price Index by Commodity: Final Demand: Finished Goods Less Foods and Energy"
data%>%tail(20)

cpi<-getSymbols('CPIAUCSL',src='FRED')
cpi<-CPIAUCSL
index(cpi)<-as.yearmon(index(cpi))
data<-merge(data,cpi)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: All Items in U.S. City Average"
data%>%tail(20)

cpi_food_beverages<-getSymbols('CPIFABSL',src='FRED')
cpi_food_beverages<-CPIFABSL
index(cpi_food_beverages)<-as.yearmon(index(cpi_food_beverages))
data<-merge(data,cpi_food_beverages)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Food and Beverages in U.S. City Average"
data%>%tail(20)

cpi_housing<-getSymbols('CPIHOSSL',src='FRED')
cpi_housing<-CPIHOSSL
index(cpi_housing)<-as.yearmon(index(cpi_housing))
data<-merge(data,cpi_housing)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Housing in U.S. City Average"
data%>%tail(20)

cpi_apparel<-getSymbols('CPIAPPSL',src='FRED')
cpi_apparel<-CPIAPPSL
index(cpi_apparel)<-as.yearmon(index(cpi_apparel))
data<-merge(data,cpi_apparel)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Apparel in U.S. City Average"
data%>%tail(20)

cpi_transportation<-getSymbols('CPITRNSL',src='FRED')
cpi_transportation<-CPITRNSL
index(cpi_transportation)<-as.yearmon(index(cpi_transportation))
data<-merge(data,cpi_transportation)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Transportation in U.S. City Average"
data%>%tail(20)

cpi_med_care<-getSymbols('CPIMEDSL',src='FRED')
cpi_med_care<-CPIMEDSL
index(cpi_med_care)<-as.yearmon(index(cpi_med_care))
data<-merge(data,cpi_med_care)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Medical Care in U.S. City Average"
data%>%tail(20)

cpi_comm<-getSymbols('CUSR0000SAC',src='FRED')
cpi_comm<-CUSR0000SAC
index(cpi_comm)<-as.yearmon(index(cpi_comm))
data<-merge(data,cpi_comm)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Commodities in U.S. City Average"
data%>%tail(20)

cpi_serv<-getSymbols('CUSR0000SAS',src='FRED')
cpi_serv<-CUSR0000SAS
index(cpi_serv)<-as.yearmon(index(cpi_serv))
data<-merge(data,cpi_serv)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: Services in U.S. City Average"
data%>%tail(20)

cpi_all_less_food<-getSymbols('CPIULFSL',src='FRED')
cpi_all_less_food<-CPIULFSL
index(cpi_all_less_food)<-as.yearmon(index(cpi_all_less_food))
data<-merge(data,cpi_all_less_food)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: All Items Less Food in U.S. City Average"
data%>%tail(20)

cpi_all_less_shelter<-getSymbols('CUSR0000SA0L2',src='FRED')
cpi_all_less_shelter<-CUSR0000SA0L2
index(cpi_all_less_shelter)<-as.yearmon(index(cpi_all_less_shelter))
data<-merge(data,cpi_all_less_shelter)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: All Items Less Shelter in U.S. City Average"
data%>%tail(20)

cpi_all_less_med<-getSymbols('CUSR0000SA0L5',src='FRED')
cpi_all_less_med<-CUSR0000SA0L5
index(cpi_all_less_med)<-as.yearmon(index(cpi_all_less_med))
data<-merge(data,cpi_all_less_med)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: All Items Less Medical Care in U.S. City Average"
data%>%tail(20)

cpi_all_less_food_and_energy<-getSymbols('CPILFESL',src='FRED')
cpi_all_less_food_and_energy<-CPILFESL
index(cpi_all_less_food_and_energy)<-as.yearmon(index(cpi_all_less_food_and_energy))
data<-merge(data,cpi_all_less_food_and_energy)
colnames(data)[ncol(data)]<- "Consumer Price Index for All Urban Consumers: All Items Less Food and Energy in U.S. City Average"
data%>%tail(20)

pce<-getSymbols('PCEPI',src='FRED')
pce<-PCEPI
index(pce)<-as.yearmon(index(pce))
data<-merge(data,pce)
colnames(data)[ncol(data)]<- "Personal Consumption Expenditures: Chain-type Price Index"
data%>%tail(20)

pce_less_food_and_energy<-getSymbols('PCEPILFE',src='FRED')
pce_less_food_and_energy<-PCEPILFE
index(pce_less_food_and_energy)<-as.yearmon(index(pce_less_food_and_energy))
data<-merge(data,pce_less_food_and_energy)
colnames(data)[ncol(data)]<- "Personal Consumption Expenditures Excluding Food and Energy (Chain-Type Price Index)"
data%>%tail(20)

pce_dur<-getSymbols('DDURRG3M086SBEA',src='FRED')
pce_dur<-DDURRG3M086SBEA
index(pce_dur)<-as.yearmon(index(pce_dur))
data<-merge(data,pce_dur)
colnames(data)[ncol(data)]<- "Personal consumption expenditures: Durable goods (chain-type price index)"
data%>%tail(20)

pce_nondur<-getSymbols('DNDGRG3M086SBEA',src='FRED')
pce_nondur<-DNDGRG3M086SBEA
index(pce_nondur)<-as.yearmon(index(pce_nondur))
data<-merge(data,pce_nondur)
colnames(data)[ncol(data)]<- "Personal consumption expenditures: Nondurable goods (chain-type price index)"
data%>%tail(20)

pce_serv<-getSymbols('DSERRG3M086SBEA',src='FRED')
pce_serv<-DSERRG3M086SBEA
index(pce_serv)<-as.yearmon(index(pce_serv))
data<-merge(data,pce_serv)
colnames(data)[ncol(data)]<- "Personal consumption expenditures: Services (chain-type price index)"
data%>%tail(20)

avg_hour_earn_total_priv<-getSymbols('CES0500000003',src='FRED')
avg_hour_earn_total_priv<-CES0500000003
index(avg_hour_earn_total_priv)<-as.yearmon(index(avg_hour_earn_total_priv))
data<-merge(data,avg_hour_earn_total_priv)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Total Private"
data%>%tail(20)

avg_hour_earn_constr<-getSymbols('CES2000000003',src='FRED')
avg_hour_earn_constr<-CES2000000003
index(avg_hour_earn_constr)<-as.yearmon(index(avg_hour_earn_constr))
data<-merge(data,avg_hour_earn_constr)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Construction"
data%>%tail(20)

avg_hour_earn_manu<-getSymbols('CES3000000003',src='FRED')
avg_hour_earn_manu<-CES3000000003
index(avg_hour_earn_manu)<-as.yearmon(index(avg_hour_earn_manu))
data<-merge(data,avg_hour_earn_manu)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Manufacturing"
data%>%tail(20)

avg_hour_earn_prof_and_bus_serv<-getSymbols('CES6000000003',src='FRED')
avg_hour_earn_prof_and_bus_serv<-CES6000000003
index(avg_hour_earn_prof_and_bus_serv)<-as.yearmon(index(avg_hour_earn_prof_and_bus_serv))
data<-merge(data,avg_hour_earn_prof_and_bus_serv)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Professional and Business Services"
data%>%tail(20)

avg_hour_earn_educ_and_health_serv<-getSymbols('CES6500000003',src='FRED')
avg_hour_earn_educ_and_health_serv<-CES6500000003
index(avg_hour_earn_educ_and_health_serv)<-as.yearmon(index(avg_hour_earn_educ_and_health_serv))
data<-merge(data,avg_hour_earn_educ_and_health_serv)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Education and Health Services"
data%>%tail(20)

avg_hour_earn_other_serv<-getSymbols('CES8000000003',src='FRED')
avg_hour_earn_other_serv<-CES8000000003
index(avg_hour_earn_other_serv)<-as.yearmon(index(avg_hour_earn_other_serv))
data<-merge(data,avg_hour_earn_other_serv)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of All Employees, Other Services"
data%>%tail(20)

federal_surplus<-getSymbols('MTSDS133FMS',src='FRED')
federal_surplus<-MTSDS133FMS
index(federal_surplus)<-as.yearmon(index(federal_surplus))
data<-merge(data,federal_surplus)
colnames(data)[ncol(data)]<- "Federal Surplus or Deficit [-]"
data%>%tail(20)

#New
avg_hour_earn_production_and_nonsupervisory <- getSymbols('AHETPI',src='FRED')
avg_hour_earn_production_and_nonsupervisory<-AHETPI
index(avg_hour_earn_production_and_nonsupervisory)<-as.yearmon(index(avg_hour_earn_production_and_nonsupervisory))
data<-merge(data,avg_hour_earn_production_and_nonsupervisory)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of Production and Nonsupervisory Employees, Total Private"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_manu <- getSymbols('CES3000000008',src='FRED')
avg_hour_earn_production_and_nonsupervisory_manu<-CES3000000008
index(avg_hour_earn_production_and_nonsupervisory_manu)<-as.yearmon(index(avg_hour_earn_production_and_nonsupervisory_manu))
data<-merge(data,avg_hour_earn_production_and_nonsupervisory_manu)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of Production and Nonsupervisory Employees, Manufacturing"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_durables <- getSymbols('CES3100000008',src='FRED')
avg_hour_earn_production_and_nonsupervisory_durables<-CES3100000008
index(avg_hour_earn_production_and_nonsupervisory_durables)<-as.yearmon(index(avg_hour_earn_production_and_nonsupervisory_durables))
data<-merge(data,avg_hour_earn_production_and_nonsupervisory_durables)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of Production and Nonsupervisory Employees, Durable Goods"
data%>%tail(20)


avg_hour_earn_production_and_nonsupervisory_construction <- getSymbols('CES2000000008',src='FRED')
avg_hour_earn_production_and_nonsupervisory_construction<-CES2000000008
index(avg_hour_earn_production_and_nonsupervisory_construction)<-as.yearmon(index(avg_hour_earn_production_and_nonsupervisory_construction))
data<-merge(data,avg_hour_earn_production_and_nonsupervisory_construction)
colnames(data)[ncol(data)]<- "Average Hourly Earnings of Production and Nonsupervisory Employees, Construction"
data%>%tail(20)


# consumer_credit <- getSymbols('TOTALSL',src='FRED')
# consumer_credit<-TOTALSL
# index(consumer_credit)<-as.yearmon(index(consumer_credit))
# data<-merge(data,consumer_credit)
# colnames(data)[ncol(data)]<- "Total Consumer Credit Owned and Securitized"
# data%>%tail(20)




# rate_credit_card <- getSymbols('TERMCBCCALLNS',src='FRED')
# rate_credit_card<-TERMCBCCALLNS
# index(rate_credit_card)<-as.yearmon(index(rate_credit_card))
# data<-merge(data,rate_credit_card)
# colnames(data)[ncol(data)]<- "Commercial Bank Interest Rate on Credit Card Plans, All Accounts"
# data%>%tail(20)
# 
# 
rate_3_month_cd <- getSymbols('IR3TCD01USM156N',src='FRED')
rate_3_month_cd<-IR3TCD01USM156N
index(rate_3_month_cd)<-as.yearmon(index(rate_3_month_cd))
rate_3_month_cd <- na.approx(rate_3_month_cd)
data<-merge(data,rate_3_month_cd)
colnames(data)[ncol(data)]<- "Interest Rates: 3-Month or 90-Day Rates and Yields: Certificates of Deposit: Total for United States"
data%>%tail(20)


rate_personal_loans_24_months <- getSymbols('TERMCBPER24NS',src='FRED')
rate_personal_loans_24_months<-TERMCBPER24NS
index(rate_personal_loans_24_months)<-as.yearmon(index(rate_personal_loans_24_months))
rate_personal_loans_24_months <- na.approx(rate_personal_loans_24_months)
data<-merge(data,rate_personal_loans_24_months)
colnames(data)[ncol(data)]<- "Finance Rate on Personal Loans at Commercial Banks, 24 Month Loan"
data%>%tail(20)


unemployment_rate <- getSymbols('UNRATE',src='FRED')
unemployment_rate<-UNRATE
index(unemployment_rate)<-as.yearmon(index(unemployment_rate))
unemployment_rate <- na.approx(unemployment_rate)
data<-merge(data,unemployment_rate)
colnames(data)[ncol(data)]<- "Unemployment Rate"
data%>%tail(20)
