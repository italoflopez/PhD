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

data%>%tail(30)

#First remove variables that have NAs at the end of the sample
data_updated <- data[, colSums(is.na(data%>%tail(30))) < 12]

data_updated%>%tail(30)

ncol(data_updated)

data_updated <- data_updated[complete.cases(data_updated),]

data_updated%>%tail(30)
ncol(data_updated)
nrow(data_updated)

colnames(data_updated)[1]

plotly_simple(data_updated[,1],colnames(data_updated)[1])
data_updated[,1]%>%adf.test()
data_updated[,1]%>%kpss.test()
data_updated[,1]%>%pp.test()

data_updated[,1]%>%log()%>%diff()
plotly_simple(data_updated[,1]%>%log()%>%diff(),colnames(data_updated)[1])

data_updated[,1]%>%log()%>%diff()%>%adf.test()
data_updated[,1]%>%log()%>%diff()%>%kpss.test()
data_updated[,1]%>%log()%>%diff()%>%pp.test()

ip_stationary <- data_updated[,1]%>%log()%>%diff()*100
colnames(ip_stationary)<-colnames(data_updated)[1]

data_stationary <- ip_stationary
colnames(data_stationary)<-colnames(data_updated)[1]

colnames(data_updated)[2]

plotly_simple(data_updated[,2],colnames(data_updated)[2])
data_updated[,2]%>%adf.test()
data_updated[,2]%>%kpss.test()
data_updated[,2]%>%pp.test()

data_updated[,2]%>%log()%>%diff()
plotly_simple(data_updated[,2]%>%log()%>%diff(),colnames(data_updated)[2])

data_updated[,2]%>%log()%>%diff()%>%adf.test()
data_updated[,2]%>%log()%>%diff()%>%kpss.test()
data_updated[,2]%>%log()%>%diff()%>%pp.test()

ip_final_non_stationary <- data_updated[,2]%>%log()%>%diff()*100
colnames(ip_final_non_stationary)<-colnames(data_updated)[2]

data_stationary <- merge(data_stationary,ip_final_non_stationary)
colnames(data_stationary)<-colnames(data_updated)[1:2]

colnames(data_updated)[3]

plotly_simple(data_updated[,3],colnames(data_updated)[3])
data_updated[,3]%>%adf.test()
data_updated[,3]%>%kpss.test()
data_updated[,3]%>%pp.test()

data_updated[,3]%>%log()%>%diff()
plotly_simple(data_updated[,3]%>%log()%>%diff(),colnames(data_updated)[3])

data_updated[,3]%>%log()%>%diff()%>%adf.test()
data_updated[,3]%>%log()%>%diff()%>%kpss.test()
data_updated[,3]%>%log()%>%diff()%>%pp.test()

ip_final_stationary <- data_updated[,3]%>%log()%>%diff()*100
colnames(ip_final_stationary)<-colnames(data_updated)[3]

data_stationary <- merge(data_stationary,ip_final_stationary)
colnames(data_stationary)[3]<-colnames(data_updated)[3]

colnames(data_updated)[4]

plotly_simple(data_updated[,4],colnames(data_updated)[4])
data_updated[,4]%>%adf.test()
data_updated[,4]%>%kpss.test()
data_updated[,4]%>%pp.test()

data_updated[,4]%>%log()%>%diff()
plotly_simple(data_updated[,4]%>%log()%>%diff(),colnames(data_updated)[4])

data_updated[,4]%>%log()%>%diff()%>%adf.test()
data_updated[,4]%>%log()%>%diff()%>%kpss.test()
data_updated[,4]%>%log()%>%diff()%>%pp.test()

ip_cons_stationary <- data_updated[,4]%>%log()%>%diff()*100
colnames(ip_cons_stationary)<-colnames(data_updated)[4]

data_stationary <- merge(data_stationary,ip_cons_stationary)
colnames(data_stationary)[4]<-colnames(data_updated)[4]

colnames(data_updated)[5]

plotly_simple(data_updated[,5],colnames(data_updated)[5])
data_updated[,5]%>%adf.test()
data_updated[,5]%>%kpss.test()
data_updated[,5]%>%pp.test()

data_updated[,5]%>%log()%>%diff()
plotly_simple(data_updated[,5]%>%log()%>%diff(),colnames(data_updated)[5])

data_updated[,5]%>%log()%>%diff()%>%adf.test()
data_updated[,5]%>%log()%>%diff()%>%kpss.test()
data_updated[,5]%>%log()%>%diff()%>%pp.test()

ip_cons_dur_stationary <- data_updated[,5]%>%log()%>%diff()*100
colnames(ip_cons_dur_stationary )<-colnames(data_updated)[5]

data_stationary <- merge(data_stationary,ip_cons_dur_stationary )
colnames(data_stationary)[5]<-colnames(data_updated)[5]

colnames(data_updated)[6]

plotly_simple(data_updated[,6],colnames(data_updated)[6])
data_updated[,6]%>%adf.test()
data_updated[,6]%>%kpss.test()
data_updated[,6]%>%pp.test()

data_updated[,6]%>%log()%>%diff()
plotly_simple(data_updated[,6]%>%log()%>%diff(),colnames(data_updated)[6])

data_updated[,6]%>%log()%>%diff()%>%adf.test()
data_updated[,6]%>%log()%>%diff()%>%kpss.test()
data_updated[,6]%>%log()%>%diff()%>%pp.test()

ip_cons_non_dur_stationary <- data_updated[,6]%>%log()%>%diff()*100
colnames(ip_cons_non_dur_stationary)<-colnames(data_updated)[6]

data_stationary <- merge(data_stationary,ip_cons_non_dur_stationary)
colnames(data_stationary)[6]<-colnames(data_updated)[6]

colnames(data_updated)[7]

plotly_simple(data_updated[,7],colnames(data_updated)[7])
data_updated[,7]%>%adf.test()
data_updated[,7]%>%kpss.test()
data_updated[,7]%>%pp.test()

data_updated[,7]%>%log()%>%diff()
plotly_simple(data_updated[,7]%>%log()%>%diff(),colnames(data_updated)[7])

data_updated[,7]%>%log()%>%diff()%>%adf.test()
data_updated[,7]%>%log()%>%diff()%>%kpss.test()
data_updated[,7]%>%log()%>%diff()%>%pp.test()

ip_equip_bus_stationary <- data_updated[,7]%>%log()%>%diff()*100
colnames(ip_equip_bus_stationary)<-colnames(data_updated)[7]

data_stationary <- merge(data_stationary,ip_equip_bus_stationary)
colnames(data_stationary)[7]<-colnames(data_updated)[7]

colnames(data_updated)[8]

plotly_simple(data_updated[,8],colnames(data_updated)[8])
data_updated[,8]%>%adf.test()
data_updated[,8]%>%kpss.test()
data_updated[,8]%>%pp.test()

data_updated[,8]%>%log()%>%diff()
plotly_simple(data_updated[,8]%>%log()%>%diff(),colnames(data_updated)[8])

data_updated[,8]%>%log()%>%diff()%>%adf.test()
data_updated[,8]%>%log()%>%diff()%>%kpss.test()
data_updated[,8]%>%log()%>%diff()%>%pp.test()

ip_mat_stationary <- data_updated[,8]%>%log()%>%diff()*100
colnames(ip_mat_stationary)<-colnames(data_updated)[8]

data_stationary <- merge(data_stationary,ip_mat_stationary)
colnames(data_stationary)[8]<-colnames(data_updated)[8]

colnames(data_updated)[9]

plotly_simple(data_updated[,9],colnames(data_updated)[9])
data_updated[,9]%>%adf.test()
data_updated[,9]%>%kpss.test()
data_updated[,9]%>%pp.test()

data_updated[,9]%>%log()%>%diff()
plotly_simple(data_updated[,9]%>%log()%>%diff(),colnames(data_updated)[9])

data_updated[,9]%>%log()%>%diff()%>%adf.test()
data_updated[,9]%>%log()%>%diff()%>%kpss.test()
data_updated[,9]%>%log()%>%diff()%>%pp.test()

ip_manu_stationary <- data_updated[,9]%>%log()%>%diff()*100
colnames(ip_manu_stationary)<-colnames(data_updated)[9]

data_stationary <- merge(data_stationary,ip_manu_stationary)
colnames(data_stationary)[9]<-colnames(data_updated)[9]

colnames(data_updated)[10]

plotly_simple(data_updated[,10],colnames(data_updated)[10])
data_updated[,10]%>%adf.test()
data_updated[,10]%>%kpss.test()
data_updated[,10]%>%pp.test()

data_updated[,10]%>%log()%>%diff()
plotly_simple(data_updated[,10]%>%log()%>%diff(),colnames(data_updated)[10])

data_updated[,10]%>%log()%>%diff()%>%adf.test()
data_updated[,10]%>%log()%>%diff()%>%kpss.test()
data_updated[,10]%>%log()%>%diff()%>%pp.test()

ip_ener_stationary <- data_updated[,10]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[10]

data_stationary <- merge(data_stationary,ip_ener_stationary)
colnames(data_stationary)[10]<-colnames(data_updated)[10]

colnames(data_updated)[11]

plotly_simple(data_updated[,11],colnames(data_updated)[11])
data_updated[,11]%>%adf.test()
data_updated[,11]%>%kpss.test()
data_updated[,11]%>%pp.test()

data_updated[,11]%>%log()%>%diff()
plotly_simple(data_updated[,11]%>%log()%>%diff(),colnames(data_updated)[11])

data_updated[,11]%>%log()%>%diff()%>%adf.test()
data_updated[,11]%>%log()%>%diff()%>%kpss.test()
data_updated[,11]%>%log()%>%diff()%>%pp.test()

ip_non_ener_stationary <- data_updated[,11]%>%log()%>%diff()*100
colnames(ip_non_ener_stationary)<-colnames(data_updated)[11]

data_stationary <- merge(data_stationary,ip_non_ener_stationary)
colnames(data_stationary)[11]<-colnames(data_updated)[11]

colnames(data_updated)[12]

plotly_simple(data_updated[,12],colnames(data_updated)[12])
data_updated[,12]%>%adf.test()
data_updated[,12]%>%kpss.test()
data_updated[,12]%>%pp.test()

data_updated[,12]%>%log()%>%diff()
plotly_simple(data_updated[,12]%>%log()%>%diff(),colnames(data_updated)[12])

data_updated[,12]%>%log()%>%diff()%>%adf.test()
data_updated[,12]%>%log()%>%diff()%>%kpss.test()
data_updated[,12]%>%log()%>%diff()%>%pp.test()

ip_manu_dur_mvp_stationary <- data_updated[,12]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[12]

data_stationary <- merge(data_stationary,ip_manu_dur_mvp_stationary)
colnames(data_stationary)[12]<-colnames(data_updated)[12]

colnames(data_updated)[13]

plotly_simple(data_updated[,13],colnames(data_updated)[13])
data_updated[,13]%>%adf.test()
data_updated[,13]%>%kpss.test()
data_updated[,13]%>%pp.test()

data_updated[,13]%>%log()%>%diff()
plotly_simple(data_updated[,13]%>%log()%>%diff(),colnames(data_updated)[13])

data_updated[,13]%>%log()%>%diff()%>%adf.test()
data_updated[,13]%>%log()%>%diff()%>%kpss.test()
data_updated[,13]%>%log()%>%diff()%>%pp.test()

ip_manu_dur_comp <- data_updated[,13]%>%log()%>%diff()*100
colnames(ip_manu_dur_comp)<-colnames(data_updated)[13]

data_stationary <- merge(data_stationary,ip_manu_dur_comp)
colnames(data_stationary)[13]<-colnames(data_updated)[13]

colnames(data_updated)[14]

plotly_simple(data_updated[,14],colnames(data_updated)[14])
data_updated[,14]%>%adf.test()
data_updated[,14]%>%kpss.test()
data_updated[,14]%>%pp.test()

data_updated[,14]%>%log()%>%diff()
plotly_simple(data_updated[,14]%>%log()%>%diff(),colnames(data_updated)[14])

data_updated[,14]%>%log()%>%diff()%>%adf.test()
data_updated[,14]%>%log()%>%diff()%>%kpss.test()
data_updated[,14]%>%log()%>%diff()%>%pp.test()

cap_util_stationary <- data_updated[,14]%>%log()%>%diff()*100
colnames(cap_util_stationary)<-colnames(data_updated)[14]

data_stationary <- merge(data_stationary,cap_util_stationary)
colnames(data_stationary)[14]<-colnames(data_updated)[14]

colnames(data_updated)[15]

plotly_simple(data_updated[,15],colnames(data_updated)[15])
data_updated[,15]%>%adf.test()
data_updated[,15]%>%kpss.test()
data_updated[,15]%>%pp.test()

data_updated[,15]%>%log()%>%diff()
plotly_simple(data_updated[,15]%>%log()%>%diff(),colnames(data_updated)[15])

data_updated[,15]%>%log()%>%diff()%>%adf.test()
data_updated[,15]%>%log()%>%diff()%>%kpss.test()
data_updated[,15]%>%log()%>%diff()%>%pp.test()

cap_util_manu_stationary <- data_updated[,15]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[15]

data_stationary <- merge(data_stationary,cap_util_manu_stationary)
colnames(data_stationary)[15]<-colnames(data_updated)[15]

colnames(data_updated)[16]

plotly_simple(data_updated[,16],colnames(data_updated)[16])
data_updated[,16]%>%adf.test()
data_updated[,16]%>%kpss.test()
data_updated[,16]%>%pp.test()

data_updated[,16]%>%log()%>%diff()
plotly_simple(data_updated[,16]%>%log()%>%diff(),colnames(data_updated)[16])

data_updated[,16]%>%log()%>%diff()%>%adf.test()
data_updated[,16]%>%log()%>%diff()%>%kpss.test()
data_updated[,16]%>%log()%>%diff()%>%pp.test()

cap_util_manu_2_stationary <- data_updated[,16]%>%log()%>%diff()*100
colnames(cap_util_manu_2_stationary)<-colnames(data_updated)[16]

data_stationary <- merge(data_stationary,cap_util_manu_2_stationary)
colnames(data_stationary)[16]<-colnames(data_updated)[16]

colnames(data_updated)[17]

plotly_simple(data_updated[,17],colnames(data_updated)[17])
data_updated[,17]%>%adf.test()
data_updated[,17]%>%kpss.test()
data_updated[,17]%>%pp.test()

data_updated[,17]%>%log()%>%diff()
plotly_simple(data_updated[,17]%>%log()%>%diff(),colnames(data_updated)[17])

data_updated[,17]%>%log()%>%diff()%>%adf.test()
data_updated[,17]%>%log()%>%diff()%>%kpss.test()
data_updated[,17]%>%log()%>%diff()%>%pp.test()

cap_util_manu_dur_stationary <- data_updated[,17]%>%log()%>%diff()*100
colnames(cap_util_manu_dur_stationary)<-colnames(data_updated)[17]

data_stationary <- merge(data_stationary,cap_util_manu_dur_stationary)
colnames(data_stationary)[17]<-colnames(data_updated)[17]

colnames(data_updated)[18]

plotly_simple(data_updated[,18],colnames(data_updated)[18])
data_updated[,18]%>%adf.test()
data_updated[,18]%>%kpss.test()
data_updated[,18]%>%pp.test()

data_updated[,18]%>%log()%>%diff()
plotly_simple(data_updated[,18]%>%log()%>%diff(),colnames(data_updated)[18])

data_updated[,18]%>%log()%>%diff()%>%adf.test()
data_updated[,18]%>%log()%>%diff()%>%kpss.test()
data_updated[,18]%>%log()%>%diff()%>%pp.test()

cap_util_manu_non_dur_stationary <- data_updated[,18]%>%log()%>%diff()*100
colnames(cap_util_manu_non_dur_stationary)<-colnames(data_updated)[18]

data_stationary <- merge(data_stationary,cap_util_manu_non_dur_stationary)
colnames(data_stationary)[18]<-colnames(data_updated)[18]

colnames(data_updated)[19]

plotly_simple(data_updated[,19],colnames(data_updated)[19])
data_updated[,19]%>%adf.test()
data_updated[,19]%>%kpss.test()
data_updated[,19]%>%pp.test()

data_updated[,19]%>%log()%>%diff()
plotly_simple(data_updated[,19]%>%log()%>%diff(),colnames(data_updated)[19])

data_updated[,19]%>%log()%>%diff()%>%adf.test()
data_updated[,19]%>%log()%>%diff()%>%kpss.test()
data_updated[,19]%>%log()%>%diff()%>%pp.test()

cap_util_comp_stationary <- data_updated[,19]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[19]

data_stationary <- merge(data_stationary,cap_util_comp_stationary)
colnames(data_stationary)[19]<-colnames(data_updated)[19]

colnames(data_updated)[20]

plotly_simple(data_updated[,20],colnames(data_updated)[20])
data_updated[,20]%>%adf.test()
data_updated[,20]%>%kpss.test()
data_updated[,20]%>%pp.test()

data_updated[,20]%>%log()%>%diff()
plotly_simple(data_updated[,20]%>%log()%>%diff(),colnames(data_updated)[20])

data_updated[,20]%>%log()%>%diff()%>%adf.test()
data_updated[,20]%>%log()%>%diff()%>%kpss.test()
data_updated[,20]%>%log()%>%diff()%>%pp.test()

unemploy_level_stationary <- data_updated[,20]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[20]

data_stationary <- merge(data_stationary,unemploy_level_stationary)
colnames(data_stationary)[20]<-colnames(data_updated)[20]

colnames(data_updated)[21]

plotly_simple(data_updated[,21],colnames(data_updated)[21])
data_updated[,21]%>%adf.test()
data_updated[,21]%>%kpss.test()
data_updated[,21]%>%pp.test()

data_updated[,21]%>%log()%>%diff()
plotly_simple(data_updated[,21]%>%log()%>%diff(),colnames(data_updated)[21])

data_updated[,21]%>%log()%>%diff()%>%adf.test()
data_updated[,21]%>%log()%>%diff()%>%kpss.test()
data_updated[,21]%>%log()%>%diff()%>%pp.test()

employ_total_stationary <- data_updated[,21]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[21]

data_stationary <- merge(data_stationary,employ_total_stationary)
colnames(data_stationary)[21]<-colnames(data_updated)[21]

colnames(data_updated)[22]

plotly_simple(data_updated[,22],colnames(data_updated)[22])
data_updated[,22]%>%adf.test()
data_updated[,22]%>%kpss.test()
data_updated[,22]%>%pp.test()

data_updated[,22]%>%log()%>%diff()
plotly_simple(data_updated[,22]%>%log()%>%diff(),colnames(data_updated)[22])

data_updated[,22]%>%log()%>%diff()%>%adf.test()
data_updated[,22]%>%log()%>%diff()%>%kpss.test()
data_updated[,22]%>%log()%>%diff()%>%pp.test()

avg_week_unemploy_stationary <- data_updated[,22]%>%log()%>%diff()*100
colnames(avg_week_unemploy_stationary)<-colnames(data_updated)[22]

data_stationary <- merge(data_stationary,avg_week_unemploy_stationary)
colnames(data_stationary)[22]<-colnames(data_updated)[22]

colnames(data_updated)[23]

plotly_simple(data_updated[,23],colnames(data_updated)[23])
data_updated[,23]%>%adf.test()
data_updated[,23]%>%kpss.test()
data_updated[,23]%>%pp.test()

data_updated[,23]%>%log()%>%diff()
plotly_simple(data_updated[,23]%>%log()%>%diff(),colnames(data_updated)[23])

data_updated[,23]%>%log()%>%diff()%>%adf.test()
data_updated[,23]%>%log()%>%diff()%>%kpss.test()
data_updated[,23]%>%log()%>%diff()%>%pp.test()

unemploy_less_5_stationary <- data_updated[,23]
colnames(unemploy_less_5_stationary)<-colnames(data_updated)[23]

data_stationary <- merge(data_stationary,unemploy_less_5_stationary)
colnames(data_stationary)[23]<-colnames(data_updated)[23]

colnames(data_updated)[24]

plotly_simple(data_updated[,24],colnames(data_updated)[24])
data_updated[,24]%>%adf.test()
data_updated[,24]%>%kpss.test()
data_updated[,24]%>%pp.test()

data_updated[,24]%>%log()%>%diff()
plotly_simple(data_updated[,24]%>%log()%>%diff(),colnames(data_updated)[24])

data_updated[,24]%>%log()%>%diff()%>%adf.test()
data_updated[,24]%>%log()%>%diff()%>%kpss.test()
data_updated[,24]%>%log()%>%diff()%>%pp.test()

unemploy_5_14_stationary <- data_updated[,24]
colnames(unemploy_5_14_stationary)<-colnames(data_updated)[24]

data_stationary <- merge(data_stationary,unemploy_5_14_stationary)
colnames(data_stationary)[24]<-colnames(data_updated)[24]

colnames(data_updated)[25]

plotly_simple(data_updated[,25],colnames(data_updated)[25])
data_updated[,25]%>%adf.test()
data_updated[,25]%>%kpss.test()
data_updated[,25]%>%pp.test()

data_updated[,25]%>%log()%>%diff()
plotly_simple(data_updated[,25]%>%log()%>%diff(),colnames(data_updated)[25])

data_updated[,25]%>%log()%>%diff()%>%adf.test()
data_updated[,25]%>%log()%>%diff()%>%kpss.test()
data_updated[,25]%>%log()%>%diff()%>%pp.test()

unemploy_15_26_stationary <- data_updated[,25]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,unemploy_15_26_stationary)
colnames(data_stationary)[25]<-colnames(data_updated)[25]

colnames(data_updated)[26]

plotly_simple(data_updated[,26],colnames(data_updated)[26])
data_updated[,26]%>%adf.test()
data_updated[,26]%>%kpss.test()
data_updated[,26]%>%pp.test()

data_updated[,26]%>%log()%>%diff()
plotly_simple(data_updated[,26]%>%log()%>%diff(),colnames(data_updated)[26])

data_updated[,26]%>%log()%>%diff()%>%adf.test()
data_updated[,26]%>%log()%>%diff()%>%kpss.test()
data_updated[,26]%>%log()%>%diff()%>%pp.test()

unemploy_15_plus_stationary <- data_updated[,26]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,unemploy_15_plus_stationary)
colnames(data_stationary)[26]<-colnames(data_updated)[26]

colnames(data_updated)[27]

plotly_simple(data_updated[,27],colnames(data_updated)[27])
data_updated[,27]%>%adf.test()
data_updated[,27]%>%kpss.test()
data_updated[,27]%>%pp.test()

data_updated[,27]%>%log()%>%diff()
plotly_simple(data_updated[,27]%>%log()%>%diff(),colnames(data_updated)[27])

data_updated[,27]%>%log()%>%diff()%>%adf.test()
data_updated[,27]%>%log()%>%diff()%>%kpss.test()
data_updated[,27]%>%log()%>%diff()%>%pp.test()

initial_stationary <- data_updated[,27]

data_stationary <- merge(data_stationary,initial_stationary)
colnames(data_stationary)[27]<-colnames(data_updated)[27]

colnames(data_updated)[28]

plotly_simple(data_updated[,28],colnames(data_updated)[28])
data_updated[,28]%>%adf.test()
data_updated[,28]%>%kpss.test()
data_updated[,28]%>%pp.test()

data_updated[,28]%>%log()%>%diff()
plotly_simple(data_updated[,28]%>%log()%>%diff(),colnames(data_updated)[28])

data_updated[,28]%>%log()%>%diff()%>%adf.test()
data_updated[,28]%>%log()%>%diff()%>%kpss.test()
data_updated[,28]%>%log()%>%diff()%>%pp.test()

employ_private_stationary <- data_updated[,28]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_private_stationary)
colnames(data_stationary)[28]<-colnames(data_updated)[28]

colnames(data_updated)[29]

plotly_simple(data_updated[,29],colnames(data_updated)[29])
data_updated[,29]%>%adf.test()
data_updated[,29]%>%kpss.test()
data_updated[,29]%>%pp.test()

data_updated[,29]%>%log()%>%diff()
plotly_simple(data_updated[,29]%>%log()%>%diff(),colnames(data_updated)[29])

data_updated[,29]%>%log()%>%diff()%>%adf.test()
data_updated[,29]%>%log()%>%diff()%>%kpss.test()
data_updated[,29]%>%log()%>%diff()%>%pp.test()

sales_manu_stationary <- data_updated[,29]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,sales_manu_stationary)
colnames(data_stationary)[29]<-colnames(data_updated)[29]

colnames(data_updated)[30]

plotly_simple(data_updated[,30],colnames(data_updated)[30])
data_updated[,30]%>%adf.test()
data_updated[,30]%>%kpss.test()
data_updated[,30]%>%pp.test()

data_updated[,30]%>%log()%>%diff()
plotly_simple(data_updated[,30]%>%log()%>%diff(),colnames(data_updated)[30])

data_updated[,30]%>%log()%>%diff()%>%adf.test()
data_updated[,30]%>%log()%>%diff()%>%kpss.test()
data_updated[,30]%>%log()%>%diff()%>%pp.test()

pers_comp_stationary <- data_updated[,30]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pers_comp_stationary)
colnames(data_stationary)[30]<-colnames(data_updated)[30]

colnames(data_updated)[31]

plotly_simple(data_updated[,31],colnames(data_updated)[31])
data_updated[,31]%>%adf.test()
data_updated[,31]%>%kpss.test()
data_updated[,31]%>%pp.test()

data_updated[,31]%>%log()%>%diff()
plotly_simple(data_updated[,31]%>%log()%>%diff(),colnames(data_updated)[31])

data_updated[,31]%>%log()%>%diff()%>%adf.test()
data_updated[,31]%>%log()%>%diff()%>%kpss.test()
data_updated[,31]%>%log()%>%diff()%>%pp.test()

pers_cons_dur_stationary <- data_updated[,31]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pers_cons_dur_stationary)
colnames(data_stationary)[31]<-colnames(data_updated)[31]

colnames(data_updated)[32]

plotly_simple(data_updated[,32],colnames(data_updated)[32])
data_updated[,32]%>%adf.test()
data_updated[,32]%>%kpss.test()
data_updated[,32]%>%pp.test()

data_updated[,32]%>%log()%>%diff()
plotly_simple(data_updated[,32]%>%log()%>%diff(),colnames(data_updated)[32])

data_updated[,32]%>%log()%>%diff()%>%adf.test()
data_updated[,32]%>%log()%>%diff()%>%kpss.test()
data_updated[,32]%>%log()%>%diff()%>%pp.test()

pers_cons_non_dur_stationary <- data_updated[,32]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pers_cons_non_dur_stationary)
colnames(data_stationary)[32]<-colnames(data_updated)[32]

colnames(data_updated)[33]

plotly_simple(data_updated[,33],colnames(data_updated)[33])
data_updated[,33]%>%adf.test()
data_updated[,33]%>%kpss.test()
data_updated[,33]%>%pp.test()

data_updated[,33]%>%log()%>%diff()
plotly_simple(data_updated[,33]%>%log()%>%diff(),colnames(data_updated)[33])

data_updated[,33]%>%log()%>%diff()%>%adf.test()
data_updated[,33]%>%log()%>%diff()%>%kpss.test()
data_updated[,33]%>%log()%>%diff()%>%pp.test()

pers_cons_ser_stationary <- data_updated[,33]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pers_cons_ser_stationary)
colnames(data_stationary)[33]<-colnames(data_updated)[33]

colnames(data_updated)[34]

plotly_simple(data_updated[,34],colnames(data_updated)[34])
data_updated[,34]%>%adf.test()
data_updated[,34]%>%kpss.test()
data_updated[,34]%>%pp.test()

data_updated[,34]%>%log()%>%diff()
plotly_simple(data_updated[,34]%>%log()%>%diff(),colnames(data_updated)[34])

data_updated[,34]%>%log()%>%diff()%>%adf.test()
data_updated[,34]%>%log()%>%diff()%>%kpss.test()
data_updated[,34]%>%log()%>%diff()%>%pp.test()

house_started_stationary <- data_updated[,34]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,house_started_stationary)
colnames(data_stationary)[34]<-colnames(data_updated)[34]

colnames(data_updated)[35]

plotly_simple(data_updated[,35],colnames(data_updated)[35])
data_updated[,35]%>%adf.test()
data_updated[,35]%>%kpss.test()
data_updated[,35]%>%pp.test()

data_updated[,35]%>%log()%>%diff()
plotly_simple(data_updated[,35]%>%log()%>%diff(),colnames(data_updated)[35])

data_updated[,35]%>%log()%>%diff()%>%adf.test()
data_updated[,35]%>%log()%>%diff()%>%kpss.test()
data_updated[,35]%>%log()%>%diff()%>%pp.test()

house_authorized_stationary <- data_updated[,35]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,house_authorized_stationary)
colnames(data_stationary)[35]<-colnames(data_updated)[35]

colnames(data_updated)[36]

plotly_simple(data_updated[,36],colnames(data_updated)[36])
data_updated[,36]%>%adf.test()
data_updated[,36]%>%kpss.test()
data_updated[,36]%>%pp.test()

data_updated[,36]%>%log()%>%diff()
plotly_simple(data_updated[,36]%>%log()%>%diff(),colnames(data_updated)[36])

data_updated[,36]%>%log()%>%diff()%>%adf.test()
data_updated[,36]%>%log()%>%diff()%>%kpss.test()
data_updated[,36]%>%log()%>%diff()%>%pp.test()

manu_invent_stationary <- data_updated[,36]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,manu_invent_stationary)
colnames(data_stationary)[36]<-colnames(data_updated)[36]

colnames(data_updated)[37]

plotly_simple(data_updated[,37],colnames(data_updated)[37])
data_updated[,37]%>%adf.test()
data_updated[,37]%>%kpss.test()
data_updated[,37]%>%pp.test()

data_updated[,37]%>%log()%>%diff()
plotly_simple(data_updated[,37]%>%log()%>%diff(),colnames(data_updated)[37])

data_updated[,37]%>%log()%>%diff()%>%adf.test()
data_updated[,37]%>%log()%>%diff()%>%kpss.test()
data_updated[,37]%>%log()%>%diff()%>%pp.test()

manu_new_dur_stationary <- data_updated[,37]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,manu_new_dur_stationary)
colnames(data_stationary)[37]<-colnames(data_updated)[37]

colnames(data_updated)[38]

plotly_simple(data_updated[,38],colnames(data_updated)[38])
data_updated[,38]%>%adf.test()
data_updated[,38]%>%kpss.test()
data_updated[,38]%>%pp.test()

data_updated[,38]%>%log()%>%diff()
plotly_simple(data_updated[,38]%>%log()%>%diff(),colnames(data_updated)[38])

data_updated[,38]%>%log()%>%diff()%>%adf.test()
data_updated[,38]%>%log()%>%diff()%>%kpss.test()
data_updated[,38]%>%log()%>%diff()%>%pp.test()

manu_new_non_dur_stationary <- data_updated[,38]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,manu_new_non_dur_stationary)
colnames(data_stationary)[38]<-colnames(data_updated)[38]

colnames(data_updated)[39]

plotly_simple(data_updated[,39],colnames(data_updated)[39])
data_updated[,39]%>%adf.test()
data_updated[,39]%>%kpss.test()
data_updated[,39]%>%pp.test()

data_updated[,39]%>%log()%>%diff()
plotly_simple(data_updated[,39]%>%log()%>%diff(),colnames(data_updated)[39])

data_updated[,39]%>%log()%>%diff()%>%adf.test()
data_updated[,39]%>%log()%>%diff()%>%kpss.test()
data_updated[,39]%>%log()%>%diff()%>%pp.test()

s_p_500_stationary <- data_updated[,39]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,s_p_500_stationary)
colnames(data_stationary)[39]<-colnames(data_updated)[39]

colnames(data_updated)[40]

plotly_simple(data_updated[,40],colnames(data_updated)[40])
data_updated[,40]%>%adf.test()
data_updated[,40]%>%kpss.test()
data_updated[,40]%>%pp.test()

data_updated[,40]%>%log()%>%diff()
plotly_simple(data_updated[,40]%>%log()%>%diff(),colnames(data_updated)[40])

data_updated[,40]%>%log()%>%diff()%>%adf.test()
data_updated[,40]%>%log()%>%diff()%>%kpss.test()
data_updated[,40]%>%log()%>%diff()%>%pp.test()

dollar_euro_stationary <- data_updated[,40]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,dollar_euro_stationary)
colnames(data_stationary)[40]<-colnames(data_updated)[40]

colnames(data_updated)[41]

plotly_simple(data_updated[,41],colnames(data_updated)[41])
data_updated[,41]%>%adf.test()
data_updated[,41]%>%kpss.test()
data_updated[,41]%>%pp.test()

data_updated[,41]%>%log()%>%diff()
plotly_simple(data_updated[,41]%>%log()%>%diff(),colnames(data_updated)[41])

data_updated[,41]%>%log()%>%diff()%>%adf.test()
data_updated[,41]%>%log()%>%diff()%>%kpss.test()
data_updated[,41]%>%log()%>%diff()%>%pp.test()

swiss_dollar_stationary <- data_updated[,41]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,swiss_dollar_stationary)
colnames(data_stationary)[41]<-colnames(data_updated)[41]

colnames(data_updated)[42]

plotly_simple(data_updated[,42],colnames(data_updated)[42])
data_updated[,42]%>%adf.test()
data_updated[,42]%>%kpss.test()
data_updated[,42]%>%pp.test()

data_updated[,42]%>%log()%>%diff()
plotly_simple(data_updated[,42]%>%log()%>%diff(),colnames(data_updated)[42])

data_updated[,42]%>%log()%>%diff()%>%adf.test()
data_updated[,42]%>%log()%>%diff()%>%kpss.test()
data_updated[,42]%>%log()%>%diff()%>%pp.test()

yen_dollar_stationary <- data_updated[,42]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,yen_dollar_stationary)
colnames(data_stationary)[42]<-colnames(data_updated)[42]

colnames(data_updated)[43]

plotly_simple(data_updated[,43],colnames(data_updated)[43])
data_updated[,43]%>%adf.test()
data_updated[,43]%>%kpss.test()
data_updated[,43]%>%pp.test()

data_updated[,43]%>%log()%>%diff()
plotly_simple(data_updated[,43]%>%log()%>%diff(),colnames(data_updated)[43])

data_updated[,43]%>%log()%>%diff()%>%adf.test()
data_updated[,43]%>%log()%>%diff()%>%kpss.test()
data_updated[,43]%>%log()%>%diff()%>%pp.test()

dollar_pound_stationary <- data_updated[,43]%>%log()%>%diff()*100


data_stationary <- merge(data_stationary,dollar_pound_stationary)
colnames(data_stationary)[43]<-colnames(data_updated)[43]

colnames(data_updated)[44]

plotly_simple(data_updated[,44],colnames(data_updated)[44])
data_updated[,44]%>%adf.test()
data_updated[,44]%>%kpss.test()
data_updated[,44]%>%pp.test()

colnames(data_updated)[45]

plotly_simple(data_updated[,45],colnames(data_updated)[45])
data_updated[,45]%>%adf.test()
data_updated[,45]%>%kpss.test()
data_updated[,45]%>%pp.test()

data_updated[,45]%>%diff()
plotly_simple(data_updated[,45]%>%diff(),colnames(data_updated)[45])

data_updated[,45]%>%diff()%>%adf.test()
data_updated[,45]%>%diff()%>%kpss.test()
data_updated[,45]%>%diff()%>%pp.test()

aaa_corporate_stationary <- data_updated[,45]%>%diff()


data_stationary <- merge(data_stationary,aaa_corporate_stationary)
colnames(data_stationary)[44]<-colnames(data_updated)[45]

colnames(data_updated)[46]

plotly_simple(data_updated[,46],colnames(data_updated)[46])
data_updated[,46]%>%adf.test()
data_updated[,46]%>%kpss.test()
data_updated[,46]%>%pp.test()

data_updated[,46]%>%diff()
plotly_simple(data_updated[,46]%>%diff(),colnames(data_updated)[46])

data_updated[,46]%>%diff()%>%adf.test()
data_updated[,46]%>%diff()%>%kpss.test()
data_updated[,46]%>%diff()%>%pp.test()

baa_corporate_stationary <- data_updated[,46]%>%diff()


data_stationary <- merge(data_stationary,baa_corporate_stationary)
colnames(data_stationary)[45]<-colnames(data_updated)[46]

colnames(data_updated)[47]

plotly_simple(data_updated[,47],colnames(data_updated)[47])
data_updated[,47]%>%adf.test()
data_updated[,47]%>%kpss.test()
data_updated[,47]%>%pp.test()

data_updated[,47]%>%log()%>%diff()
plotly_simple(data_updated[,47]%>%log()%>%diff(),colnames(data_updated)[47])

data_updated[,47]%>%log()%>%diff()%>%adf.test()
data_updated[,47]%>%log()%>%diff()%>%kpss.test()
data_updated[,47]%>%log()%>%diff()%>%pp.test()

m1_stationary <- data_updated[,47]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,m1_stationary)
colnames(data_stationary)[46]<-colnames(data_updated)[47]

colnames(data_updated)[48]

plotly_simple(data_updated[,48],colnames(data_updated)[48])
data_updated[,48]%>%adf.test()
data_updated[,48]%>%kpss.test()
data_updated[,48]%>%pp.test()

data_updated[,48]%>%log()%>%diff()
plotly_simple(data_updated[,48]%>%log()%>%diff(),colnames(data_updated)[48])

data_updated[,48]%>%log()%>%diff()%>%adf.test()
data_updated[,48]%>%log()%>%diff()%>%kpss.test()
data_updated[,48]%>%log()%>%diff()%>%pp.test()

cpi_stationary <- data_updated[,48]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_stationary)
colnames(data_stationary)[47]<-colnames(data_updated)[48]

colnames(data_updated)[49]

plotly_simple(data_updated[,49],colnames(data_updated)[49])
data_updated[,49]%>%adf.test()
data_updated[,49]%>%kpss.test()
data_updated[,49]%>%pp.test()

data_updated[,49]%>%log()%>%diff()
plotly_simple(data_updated[,49]%>%log()%>%diff(),colnames(data_updated)[49])

data_updated[,49]%>%log()%>%diff()%>%adf.test()
data_updated[,49]%>%log()%>%diff()%>%kpss.test()
data_updated[,49]%>%log()%>%diff()%>%pp.test()

m2_stationary <- data_updated[,49]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,m2_stationary)
colnames(data_stationary)[48]<-colnames(data_updated)[49]

colnames(data_updated)[50]

plotly_simple(data_updated[,50],colnames(data_updated)[50])
data_updated[,50]%>%adf.test()
data_updated[,50]%>%kpss.test()
data_updated[,50]%>%pp.test()

data_updated[,50]%>%log()%>%diff()
plotly_simple(data_updated[,50]%>%log()%>%diff(),colnames(data_updated)[50])

data_updated[,50]%>%log()%>%diff()%>%adf.test()
data_updated[,50]%>%log()%>%diff()%>%kpss.test()
data_updated[,50]%>%log()%>%diff()%>%pp.test()

reser_depo_ins_stationary <- data_updated[,50]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,reser_depo_ins_stationary)
colnames(data_stationary)[49]<-colnames(data_updated)[50]

data_stationary <- data_stationary[complete.cases(data_stationary),]

data_transformed <- scale(data_stationary)
colMeans(data_transformed)
sapply(data_transformed, sd)

final_data <- sapply(data_transformed, ts_clean_vec,period=12)
final_data <- as.zoo(final_data)
index(final_data) <- index(data_transformed)
final_data
final_data%>%tail(20)

final_data%>%colMeans()

sapply(final_data, sd)

write.csv(data_transformed,"standardize_macro_data.csv")

write.csv(ffr,"ffr.csv")

