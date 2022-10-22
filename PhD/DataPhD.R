library(ROracle)
library(keyring)
library(forecast)
library(zoo)
library(tseries)
library(vars)
library(portes)
library(plotly)
library(quantmod)
library(matrixStats)

#Data for the PhD
ip_total <- getSymbols("INDPRO",src='FRED')
ip_total <- INDPRO
ip_total <- as.zoo(ip_total)
plotly_simple(ip_total,"Industrial Production: Total Index")
colnames(ip_total) <- "Industrial Production: Total Index"

data <- ip_total 

ip_fpnss <- getSymbols("IPFPNSS",src='FRED')
ip_fpnss <- IPFPNSS
ip_fpnss <- as.zoo(ip_fpnss)
plotly_simple(ip_fpnss,"Industrial Production: Final Products and Nonindustrial Supplies")
colnames(ip_fpnss) <- "Industrial Production: Final Products and Nonindustrial Supplies"

data <- merge(data,ip_fpnss)

ip_final <- getSymbols("IPFINAL",src='FRED')
ip_final <- IPFINAL
ip_final <- as.zoo(ip_final)
plotly_simple(ip_final,"Industrial Production: Final Products")
colnames(ip_final) <- "Industrial Production: Final Products"

data <- merge(data,ip_final)

ip_cons <- getSymbols("IPCONGD",src='FRED')
ip_cons <- IPCONGD
ip_cons <- as.zoo(ip_cons)
plotly_simple(ip_cons,"Industrial Production: Consumer Goods")
colnames(ip_cons) <- "Industrial Production: Consumer Goods"

data <- merge(data,ip_cons)

ip_dur_cons <- getSymbols("IPDCONGD",src='FRED')
ip_dur_cons <- IPDCONGD
ip_dur_cons <- as.zoo(ip_dur_cons)
plotly_simple(ip_dur_cons,"Industrial Production: Durable Consumer Goods")
colnames(ip_dur_cons) <- "Industrial Production: Durable Consumer Goods"

data <- merge(data,ip_dur_cons)

ip_non_dur_cons <- getSymbols("IPNCONGD",src='FRED')
ip_non_dur_cons <- IPNCONGD
ip_non_dur_cons <- as.zoo(ip_non_dur_cons)
plotly_simple(ip_non_dur_cons,"Industrial Production: Non-Durable Consumer Goods")
colnames(ip_non_dur_cons) <- "Industrial Production: Non-Durable Consumer Goods"

data <- merge(data,ip_non_dur_cons)

ip_bus_eq <- getSymbols("IPBUSEQ",src='FRED')
ip_bus_eq <- IPBUSEQ
ip_bus_eq <- as.zoo(ip_bus_eq)
plotly_simple(ip_bus_eq,"Industrial Production: Equipment: Business Equipment")
colnames(ip_bus_eq) <- "Industrial Production: Equipment: Business Equipment"

data <- merge(data,ip_bus_eq)

ip_mat <- getSymbols("IPMAT",src='FRED')
ip_mat <- IPMAT
ip_mat <- as.zoo(ip_mat)
plotly_simple(ip_mat,"Industrial Production: Materials")
colnames(ip_mat) <- "Industrial Production: Materials"

data <- merge(data,ip_mat)

ip_manu <- getSymbols("IPMAN",src='FRED')
ip_manu <- IPMAN
ip_manu <- as.zoo(ip_manu)
plotly_simple(ip_manu,"Industrial Production: Manufacturing")
colnames(ip_manu) <- "Industrial Production: Manufacturing"

data <- merge(data,ip_manu)

ip_ener <- getSymbols("IPB50089S",src='FRED')
ip_ener <- IPB50089S
ip_ener <- as.zoo(ip_ener)
plotly_simple(ip_ener,"Industrial Production: Energy, Total")
colnames(ip_ener) <- "Industrial Production: Energy, Total"

data <- merge(data,ip_ener)

ip_non_ener <- getSymbols("IPX5001ES",src='FRED')
ip_non_ener <- IPX5001ES
ip_non_ener <- as.zoo(ip_non_ener)
plotly_simple(ip_non_ener,"Industrial Production: Non-Energy, Total")
colnames(ip_non_ener) <- "Industrial Production: Non-Energy, Total"

data <- merge(data,ip_non_ener)

ip_mvp <- getSymbols("IPG3361T3S",src='FRED')
ip_mvp <- IPG3361T3S
ip_mvp <- as.zoo(ip_mvp)
plotly_simple(ip_mvp,"Industrial Production: Manufacturing: Durable Goods: Motor Vehicles and Parts")
colnames(ip_mvp) <- "Industrial Production: Manufacturing: Durable Goods: Motor Vehicles and Parts"

data <- merge(data,ip_mvp)

ip_pccoms <- getSymbols("IPHITEK2S",src='FRED')
ip_pccoms <- IPHITEK2S
ip_pccoms <- as.zoo(ip_pccoms)
plotly_simple(ip_pccoms,"Industrial Production: Manufacturing: Durable Goods: Computers, Communications Equipment, and Semiconductors")
colnames(ip_pccoms) <- "Industrial Production: Manufacturing: Durable Goods: Computers, Communications Equipment, and Semiconductors"

data <- merge(data,ip_pccoms)

cap_total <- getSymbols("TCU",src='FRED')
cap_total <- TCU
cap_total <- as.zoo(cap_total)
plotly_simple(cap_total,"Capacity Utilization: Total Index")
colnames(cap_total) <- "Capacity Utilization: Total Index"

data <- merge(data,cap_total)

cap_manu <- getSymbols("MCUMFN",src='FRED')
cap_manu <- MCUMFN
cap_manu <- as.zoo(cap_manu)
plotly_simple(cap_manu,"Capacity Utilization: Manufacturing")
colnames(cap_manu) <- "Capacity Utilization: Manufacturing"

data <- merge(data,cap_manu)

cap_manu_dur <- getSymbols("CAPUTLGMFDS",src='FRED')
cap_manu_dur <- CAPUTLGMFDS
cap_manu_dur <- as.zoo(cap_manu_dur)
plotly_simple(cap_manu_dur,"Capacity Utilization: Durable Manufacturing")
colnames(cap_manu_dur) <- "Capacity Utilization: Durable Manufacturing"

data <- merge(data,cap_manu_dur)

cap_manu_non_dur <- getSymbols("CAPUTLGMFNS",src='FRED')
cap_manu_non_dur <- CAPUTLGMFNS
cap_manu_non_dur <- as.zoo(cap_manu_non_dur)
plotly_simple(cap_manu_non_dur,"Capacity Utilization: Non-Durable Manufacturing")
colnames(cap_manu_non_dur) <- "Capacity Utilization: Non-Durable Manufacturing"

data <- merge(data,cap_manu_non_dur)

cap_comsemi <- getSymbols("CAPUTLHITEK2S",src='FRED')
cap_comsemi <- CAPUTLHITEK2S
cap_comsemi <- as.zoo(cap_comsemi)
plotly_simple(cap_comsemi,"Capacity Utilization: : Computers, Communications Equipment, and Semiconductors")
colnames(cap_comsemi) <- "Capacity Utilization: : Computers, Communications Equipment, and Semiconductors"

data <- merge(data,cap_comsemi)

help_wanted <- getSymbols("M0882BUSM350NNBR",src='FRED')
help_wanted <- M0882BUSM350NNBR
help_wanted <- as.zoo(help_wanted)
plotly_simple(help_wanted,"Help-Wanted Advertising in Newspapers for United States")
colnames(help_wanted) <- "Help-Wanted Advertising in Newspapers for United States"

data <- merge(data,help_wanted)

unemployed <- getSymbols("UNEMPLOY",src='FRED')
unemployed  <- UNEMPLOY
unemployed  <- as.zoo(unemployed)
plotly_simple(unemployed ,"Unemployment Level")
colnames(unemployed) <- "Unemployment Level"

data <- merge(data,unemployed )

employees <- getSymbols("PAYEMS",src='FRED')
employees  <- PAYEMS
employees  <- as.zoo(employees)
plotly_simple(employees,"All Employees, Total Nonfarm")
plot(employees)
colnames(employees) <- "All Employees, Total Nonfarm"

data <- merge(data,employees)

dur_un <- getSymbols("UEMPMEAN",src='FRED')
dur_un  <- UEMPMEAN
dur_un  <- as.zoo(dur_un)
plotly_simple(dur_un,"Average Weeks Unemployed")
plot(dur_un)
colnames(dur_un) <- "Average Weeks Unemployed"

data <- merge(data,dur_un)

un_5 <- getSymbols("UEMPLT5",src='FRED')
un_5  <- UEMPLT5
un_5  <- as.zoo(un_5)
plotly_simple(un_5,"Number Unemployed for Less Than 5 Weeks")
plot(un_5)
colnames(un_5) <- "Number Unemployed for Less Than 5 Weeks"

data <- merge(data,un_5)

un_5_14 <- getSymbols("UEMP5TO14",src='FRED')
un_5_14  <- UEMP5TO14
un_5_14  <- as.zoo(un_5_14)
plotly_simple(un_5_14,"Number Unemployed for 5-14 Weeks")
plot(un_5_14)
colnames(un_5_14) <- "Number Unemployed for 5-14 Weeks"

data <- merge(data,un_5_14)

un_15_26 <- getSymbols("UEMP15T26",src='FRED')
un_15_26  <- UEMP15T26
un_15_26  <- as.zoo(un_15_26)
plotly_simple(un_15_26,"Number Unemployed for 15-26 Weeks")
plot(un_15_26)
colnames(un_15_26) <- "Number Unemployed for 15-26 Weeks"

data <- merge(data,un_15_26)

un_15_plus <- getSymbols("UEMP15OV",src='FRED')
un_15_plus <- UEMP15OV
un_15_plus  <- as.zoo(un_15_plus)
plotly_simple(un_15_plus,"Number Unemployed for 15 Weeks & over")
plot(un_15_plus)
colnames(un_15_plus) <- "Number Unemployed for 15 Weeks & over"

data <- merge(data,un_15_plus)

ini_claim <- getSymbols("ICSA",src='FRED')
ini_claim <- ICSA
ini_claim  <- as.zoo(ini_claim)
plotly_simple(ini_claim,"Initial Claims")
plot(ini_claim)
colnames(ini_claim) <- "Initial Claims"

ini_claim <- aggregate(ini_claim,as.yearmon ,sum)
plot(ini_claim)

index(data) <- as.yearmon(index(data))
data <- merge(data,ini_claim)

private_payroll <- getSymbols("NPPTTL",src='FRED')
private_payroll<- NPPTTL
private_payroll  <- as.zoo(private_payroll)
plotly_simple(private_payroll,"Total Nonfarm Private Payroll Employment")
plot(private_payroll)
colnames(private_payroll) <- "Total Nonfarm Private Payroll Employment"
index(private_payroll) <- as.yearmon(index(private_payroll))

data <- merge(data,private_payroll)

sales_m_t <- getSymbols("CMRMTSPL",src='FRED')
sales_m_t <- CMRMTSPL
sales_m_t  <- as.zoo(sales_m_t)
plotly_simple(sales_m_t,"Real Manufacturing and Trade Industries Sales")
plot(sales_m_t)
colnames(sales_m_t) <- "Real Manufacturing and Trade Industries Sales"
index(sales_m_t) <- as.yearmon(index(sales_m_t))

data <- merge(data,sales_m_t)

real_cons_exp <- getSymbols("PCEC96",src='FRED')
real_cons_exp <- PCEC96
real_cons_exp  <- as.zoo(real_cons_exp)
plotly_simple(real_cons_exp,"Real Personal Consumption Expenditures")
plot(real_cons_exp)
colnames(real_cons_exp) <- "Real Personal Consumption Expenditures"
index(real_cons_exp) <- as.yearmon(index(real_cons_exp))

data <- merge(data,real_cons_exp)

real_cons_exp_dur <- getSymbols("PCEDGC96",src='FRED')
real_cons_exp_dur <- PCEDGC96
real_cons_exp_dur  <- as.zoo(real_cons_exp_dur)
plotly_simple(real_cons_exp_dur,"Real Personal Consumption Expenditures: Durable Goods")
plot(real_cons_exp_dur)
colnames(real_cons_exp_dur) <- "Real Personal Consumption Expenditures: Durable Goods"
index(real_cons_exp_dur) <- as.yearmon(index(real_cons_exp_dur))

data <- merge(data,real_cons_exp_dur)

real_cons_exp_non_dur <- getSymbols("PCENDC96",src='FRED')
real_cons_exp_non_dur <- PCENDC96
real_cons_exp_non_dur  <- as.zoo(real_cons_exp_non_dur)
plotly_simple(real_cons_exp_non_dur,"Real Personal Consumption Expenditures: Nondurable Goods")
plot(real_cons_exp_non_dur)
colnames(real_cons_exp_non_dur) <- "Real Personal Consumption Expenditures: Nondurable Goods"
index(real_cons_exp_non_dur) <- as.yearmon(index(real_cons_exp_non_dur))

data <- merge(data,real_cons_exp_non_dur)

real_cons_exp_ser <- getSymbols("PCESC96",src='FRED')
real_cons_exp_ser <- PCESC96
real_cons_exp_ser  <- as.zoo(real_cons_exp_ser)
plotly_simple(real_cons_exp_ser,"Real Personal Consumption Expenditures: Services")
plot(real_cons_exp_ser)
colnames(real_cons_exp_ser) <- "Real Personal Consumption Expenditures: Services"
index(real_cons_exp_ser) <- as.yearmon(index(real_cons_exp_ser))

data <- merge(data,real_cons_exp_ser)

real_cons_exp_mvp <- getSymbols("DMOTRX1Q020SBEA",src='FRED')
real_cons_exp_mvp  <- DMOTRX1Q020SBEA
real_cons_exp_mvp   <- as.zoo(real_cons_exp_mvp)
plotly_simple(real_cons_exp_mvp," Real personal consumption expenditures: Durable goods: Motor vehicles and parts")
plot(real_cons_exp_mvp)
colnames(real_cons_exp_mvp) <- " Real personal consumption expenditures: Durable goods: Motor vehicles and parts"
index(real_cons_exp_mvp) <- as.yearmon(index(real_cons_exp_mvp))

data <- merge(data,real_cons_exp_mvp)

house_starts <- getSymbols("HOUST",src='FRED')
house_starts  <- HOUST
house_starts   <- as.zoo(house_starts)
plotly_simple(house_starts," New Privately-Owned Housing Units Started: Total Units")
plot(house_starts)
colnames(house_starts) <- " New Privately-Owned Housing Units Started: Total Units"
index(house_starts) <- as.yearmon(index(house_starts))

data <- merge(data,house_starts)

house_permit <- getSymbols("PERMIT",src='FRED')
house_permit  <- PERMIT
house_permit   <- as.zoo(house_permit)
plotly_simple(house_permit," New Privately-Owned Housing Units Authorized in Permit-Issuing Places: Total Units")
plot(house_permit)
colnames(house_permit) <- " New Privately-Owned Housing Units Authorized in Permit-Issuing Places: Total Units"
index(house_permit) <- as.yearmon(index(house_permit))

data <- merge(data,house_permit)

invent_ma_tr <- getSymbols("INVCMRMT",src='FRED')
invent_ma_tr  <- INVCMRMT
invent_ma_tr   <- as.zoo(invent_ma_tr)
plotly_simple(invent_ma_tr,"Real Manufacturing and Trade Inventories")
plot(invent_ma_tr)
colnames(invent_ma_tr) <- "Real Manufacturing and Trade Inventories"
index(invent_ma_tr) <- as.yearmon(index(invent_ma_tr))

data <- merge(data,invent_ma_tr)

manu_new_all <- getSymbols("M06070USM144NNBR",src='FRED')
manu_new_all  <- M06070USM144NNBR
manu_new_all   <- as.zoo(manu_new_all)
plotly_simple(manu_new_all,"Manufacturers' New Orders, All Industries for United States")
plot(manu_new_all)
colnames(manu_new_all) <- "Manufacturers' New Orders, All Industries for United States"
index(manu_new_all) <- as.yearmon(index(manu_new_all))

data <- merge(data,manu_new_all)

manu_new_dur <- getSymbols("DGORDER",src='FRED')
manu_new_dur  <- DGORDER
manu_new_dur   <- as.zoo(manu_new_dur)
plotly_simple(manu_new_dur,"Manufacturers' New Orders: Durable Goods")
plot(manu_new_dur)
colnames(manu_new_dur) <- "Manufacturers' New Orders: Durable Goods"
index(manu_new_dur) <- as.yearmon(index(manu_new_dur))

data <- merge(data,manu_new_dur)

manu_new_non_dur <- getSymbols("AMNMNO",src='FRED')
manu_new_non_dur  <- AMNMNO
manu_new_non_dur   <- as.zoo(manu_new_non_dur)
plotly_simple(manu_new_non_dur,"Manufacturers' New Orders: Nondurable Goods")
plot(manu_new_non_dur)
colnames(manu_new_non_dur) <- "Manufacturers' New Orders: Nondurable Goods"
index(manu_new_non_dur) <- as.yearmon(index(manu_new_non_dur))

data <- merge(data,manu_new_non_dur)

getSymbols('^GSPC',src='yahoo',from="1980-01-02")
s_p_500 <- GSPC
s_p_500   <- as.zoo(s_p_500)
# s_p_500 <- zoo(s_p_500[,2],order.by = as.Date(s_p_500[,1]))
s_p_500 <- aggregate(s_p_500,as.yearmon,mean,
                     na.rm=T)
s_p_500 <- s_p_500[,"GSPC.Close"]

plotly_simple(s_p_500,"S&P 500")
plot(s_p_500)
colnames(s_p_500) <- "S&P 500"
index(s_p_500) <- as.yearmon(index(s_p_500))

data <- merge(data,s_p_500)
colnames(data)[ncol(data)] <- "S&P 500"

euro <- getSymbols("DEXUSEU",src='FRED')
euro  <- DEXUSEU
euro   <- as.zoo(euro)
euro <- aggregate(euro,as.yearmon,mean,
                     na.rm=T)

plotly_simple(euro,"U.S. Dollars to Euro Spot Exchange Rate")
plot(euro)
colnames(euro) <- "U.S. Dollars to Euro Spot Exchange Rate"
index(euro) <- as.yearmon(index(euro))

data <- merge(data,euro)

swiss <- getSymbols("DEXSZUS",src='FRED')
swiss  <- DEXSZUS
swiss   <- as.zoo(swiss)
swiss <- aggregate(swiss,as.yearmon,mean,
                  na.rm=T)

plotly_simple(swiss,"Swiss Francs to U.S. Dollar Spot Exchange Rate")
plot(swiss)
colnames(swiss) <- "Swiss Francs to U.S. Dollar Spot Exchange Rate"
index(swiss) <- as.yearmon(index(swiss))

data <- merge(data,swiss)

yen <- getSymbols("DEXJPUS",src='FRED')
yen  <- DEXJPUS
yen   <- as.zoo(yen)
yen <- aggregate(yen,as.yearmon,mean,
                   na.rm=T)

plotly_simple(yen,"Japanese Yen to U.S. Dollar Spot Exchange Rate")
plot(yen)
colnames(yen) <- "Japanese Yen to U.S. Dollar Spot Exchange Rate"
index(yen) <- as.yearmon(index(yen))

data <- merge(data,yen)

pound <- getSymbols("DEXUSUK",src='FRED')
pound  <- DEXUSUK
pound   <- as.zoo(pound)
pound <- aggregate(pound,as.yearmon,mean,
                 na.rm=T)

plotly_simple(pound,"U.S. Dollars to U.K. Pound Sterling Spot Exchange Rate")
plot(pound)
colnames(pound) <- "U.S. Dollars to U.K. Pound Sterling Spot Exchange Rate"
index(pound) <- as.yearmon(index(pound))

data <- merge(data,pound)
