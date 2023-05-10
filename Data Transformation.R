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

data%>%tail(30)

#First remove variables that have NAs at the end of the sample
data_updated <- data[, colSums(is.na(data%>%tail(30))) < 12]

data_updated%>%tail(30)

ncol(data_updated)

which(!is.na(data_updated[,50]))[1]<900

bla<-lapply(data_updated, function(x) which(!is.na(x))[1]<900)
eliminated_variables <- colnames(data)[which(bla!=TRUE)%>%as.vector()]
data_updated <- data_updated[,which(bla==TRUE)%>%as.vector()]
ncol(data_updated)

data_updated <- data_updated[complete.cases(data_updated),]

data_updated <-data_updated%>%as.zoo()

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

ip_manu_dur_stationary <- data_updated[,10]%>%log()%>%diff()*100
colnames(ip_manu_dur_stationary)<-colnames(data_updated)[10]

data_stationary <- merge(data_stationary,ip_manu_dur_stationary)
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

ip_manu_nondur_stationary <- data_updated[,11]%>%log()%>%diff()*100
colnames(ip_manu_nondur_stationary)<-colnames(data_updated)[11]

data_stationary <- merge(data_stationary,ip_manu_nondur_stationary)
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

ip_ener_stationary <- data_updated[,12]%>%log()%>%diff()*100
colnames(ip_ener_stationary)<-colnames(data_updated)[12]

data_stationary <- merge(data_stationary,ip_ener_stationary)
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

ip_non_ener_stationary <- data_updated[,13]%>%log()%>%diff()*100
colnames(ip_non_ener_stationary)<-colnames(data_updated)[13]

data_stationary <- merge(data_stationary,ip_non_ener_stationary)
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

ip_comp_comm_semi_stationary <- data_updated[,14]%>%log()%>%diff()*100
colnames(ip_comp_comm_semi_stationary)<-colnames(data_updated)[14]

data_stationary <- merge(data_stationary,ip_comp_comm_semi_stationary)
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

ip_non_ener_excl_ccs_and_mvp_stationary <- data_updated[,15]%>%log()%>%diff()*100
colnames(ip_non_ener_excl_ccs_and_mvp_stationary)<-colnames(data_updated)[15]

data_stationary <- merge(data_stationary,ip_non_ener_excl_ccs_and_mvp_stationary)
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

cap_util_total_stationary <- data_updated[,16]%>%log()%>%diff()*100
colnames(cap_util_total_stationary)<-colnames(data_updated)[16]

data_stationary <- merge(data_stationary,cap_util_total_stationary)
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

cap_util_manu_stationary <- data_updated[,17]%>%log()%>%diff()*100
colnames(cap_util_manu_stationary)<-colnames(data_updated)[17]

data_stationary <- merge(data_stationary,cap_util_manu_stationary)
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

cap_util_manu_dur_stationary <- data_updated[,18]%>%log()%>%diff()*100
colnames(cap_util_manu_dur_stationary)<-colnames(data_updated)[18]

data_stationary <- merge(data_stationary,cap_util_manu_dur_stationary)
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

cap_util_manu_nondur_stationary <- data_updated[,19]%>%log()%>%diff()*100
colnames(cap_util_manu_nondur_stationary)<-colnames(data_updated)[19]

data_stationary <- merge(data_stationary,cap_util_manu_nondur_stationary)
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

cap_util_comp_comm_semi_stationary <- data_updated[,20]%>%log()%>%diff()*100
colnames(cap_util_comp_comm_semi_stationary)<-colnames(data_updated)[20]

data_stationary <- merge(data_stationary,cap_util_comp_comm_semi_stationary)
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
colnames(employ_total_stationary)<-colnames(data_updated)[21]

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

plotly_simple(data_updated[,28]%>%as.numeric(),colnames(data_updated)[28])
data_updated[,28]%>%adf.test()
data_updated[,28]%>%kpss.test()
data_updated[,28]%>%pp.test()

data_updated[,28]%>%log()%>%diff()
plotly_simple(data_updated[,28]%>%as.ts()%>%log()%>%diff(),colnames(data_updated)[28])

data_updated[,28]%>%log()%>%diff()%>%na.omit()%>%adf.test()
data_updated[,28]%>%log()%>%diff()%>%na.omit()%>%kpss.test()
data_updated[,28]%>%log()%>%diff()%>%na.omit()%>%pp.test()

employ_mining_stationary <- data_updated[,28]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_mining_stationary)
colnames(data_stationary)[28]<-colnames(data_updated)[28]

colnames(data_updated)[29]

plotly_simple(data_updated[,29]%>%as.numeric(),colnames(data_updated)[29])
data_updated[,29]%>%adf.test()
data_updated[,29]%>%kpss.test()
data_updated[,29]%>%pp.test()

data_updated[,29]%>%log()%>%diff()
plotly_simple(data_updated[,29]%>%log()%>%diff()%>%as.numeric(),colnames(data_updated)[29])

data_updated[,29]%>%log()%>%diff()%>%na.omit()%>%adf.test()
data_updated[,29]%>%log()%>%diff()%>%na.omit()%>%kpss.test()
data_updated[,29]%>%log()%>%diff()%>%na.omit()%>%pp.test()

employ_construc_stationary <- data_updated[,29]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_construc_stationary)
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

employ_manu_stationary <- data_updated[,30]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_manu_stationary)
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

employ_manu_dur_stationary <- data_updated[,31]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_manu_dur_stationary)
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

employ_manu_non_dur_stationary <- data_updated[,32]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_manu_non_dur_stationary)
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

employ_util_stationary <- data_updated[,33]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_util_stationary)
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

employ_retail_stationary <- data_updated[,34]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_retail_stationary)
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

employ_wholesale_stationary <- data_updated[,35]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_wholesale_stationary)
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

employ_finan_stationary <- data_updated[,36]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_finan_stationary)
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

employ_prof_busi_stationary <- data_updated[,37]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_prof_busi_stationary)
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

employ_educ_health_stationary <- data_updated[,38]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_educ_health_stationary)
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

employ_leisure_stationary <- data_updated[,39]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_leisure_stationary)
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

employ_other_serv_stationary <- data_updated[,40]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,employ_other_serv_stationary)
colnames(data_stationary)[40]<-colnames(data_updated)[40]

colnames(data_updated)[41]

plotly_simple(data_updated[,41]%>%as.numeric(),colnames(data_updated)[41])
data_updated[,41]%>%adf.test()
data_updated[,41]%>%kpss.test()
data_updated[,41]%>%pp.test()

data_updated[,41]%>%log()%>%diff()
plotly_simple(data_updated[,41]%>%log()%>%diff(),colnames(data_updated)[41])

data_updated[,41]%>%log()%>%diff()%>%adf.test()
data_updated[,41]%>%log()%>%diff()%>%kpss.test()
data_updated[,41]%>%log()%>%diff()%>%pp.test()

real_manu_trade_sales_stationary <- data_updated[,41]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,real_manu_trade_sales_stationary)
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

manu_sales_stationary <- data_updated[,42]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,manu_sales_stationary)
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

merch_whole_sales_stationary <- data_updated[,43]%>%log()%>%diff()*100


data_stationary <- merge(data_stationary,merch_whole_sales_stationary)
colnames(data_stationary)[43]<-colnames(data_updated)[43]

colnames(data_updated)[44]

plotly_simple(data_updated[,44],colnames(data_updated)[44])
data_updated[,44]%>%adf.test()
data_updated[,44]%>%kpss.test()
data_updated[,44]%>%pp.test()

retail_sales_stationary <- data_updated[,44]%>%log()%>%diff()*100


data_stationary <- merge(data_stationary,retail_sales_stationary)
colnames(data_stationary)[44]<-colnames(data_updated)[44]

colnames(data_updated)[45]

plotly_simple(data_updated[,45]%>%as.zoo(),colnames(data_updated)[45])
data_updated[,45]%>%adf.test()
data_updated[,45]%>%kpss.test()
data_updated[,45]%>%pp.test()

data_updated[,45]%>%diff()
plotly_simple(data_updated[,45]%>%diff(),colnames(data_updated)[45])

data_updated[,45]%>%diff()%>%adf.test()
data_updated[,45]%>%diff()%>%kpss.test()
data_updated[,45]%>%diff()%>%pp.test()

housing_units_started_stationary <- data_updated[,45]%>%log()%>%diff()*100


data_stationary <- merge(data_stationary,housing_units_started_stationary)
colnames(data_stationary)[45]<-colnames(data_updated)[45]

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

house_permit_stationary <- data_updated[,46]%>%log()%>%diff()*100


data_stationary <- merge(data_stationary,house_permit_stationary)
colnames(data_stationary)[46]<-colnames(data_updated)[46]

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

house_sold_stationary <- data_updated[,47]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,house_sold_stationary)
colnames(data_stationary)[47]<-colnames(data_updated)[47]

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

swiss_franc_stationary <- data_updated[,48]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,swiss_franc_stationary)
colnames(data_stationary)[48]<-colnames(data_updated)[48]

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

yen_stationary <- data_updated[,49]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,yen_stationary)
colnames(data_stationary)[49]<-colnames(data_updated)[49]

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

sterling_pound_stationary <- data_updated[,50]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,sterling_pound_stationary)
colnames(data_stationary)[50]<-colnames(data_updated)[50]

colnames(data_updated)[51]

plotly_simple(data_updated[,51],colnames(data_updated)[51])
data_updated[,51]%>%adf.test()
data_updated[,51]%>%kpss.test()
data_updated[,51]%>%pp.test()

data_updated[,51]%>%diff()
plotly_simple(data_updated[,51]%>%diff(),colnames(data_updated)[51])

data_updated[,51]%>%log()%>%diff()%>%adf.test()
data_updated[,51]%>%log()%>%diff()%>%kpss.test()
data_updated[,51]%>%log()%>%diff()%>%pp.test()

aaa_bond_stationary <- data_updated[,51]%>%diff()*100

data_stationary <- merge(data_stationary,aaa_bond_stationary)
colnames(data_stationary)[51]<-colnames(data_updated)[51]

colnames(data_updated)[52]

plotly_simple(data_updated[,52],colnames(data_updated)[52])
data_updated[,52]%>%adf.test()
data_updated[,52]%>%kpss.test()
data_updated[,52]%>%pp.test()

data_updated[,52]%>%diff()
plotly_simple(data_updated[,52]%>%diff(),colnames(data_updated)[52])

data_updated[,52]%>%log()%>%diff()%>%adf.test()
data_updated[,52]%>%log()%>%diff()%>%kpss.test()
data_updated[,52]%>%log()%>%diff()%>%pp.test()

baa_bond_stationary <- data_updated[,52]%>%diff()*100

data_stationary <- merge(data_stationary,baa_bond_stationary)
colnames(data_stationary)[52]<-colnames(data_updated)[52]

colnames(data_updated)[53]

plotly_simple(data_updated[,53],colnames(data_updated)[53])
data_updated[,53]%>%adf.test()
data_updated[,53]%>%kpss.test()
data_updated[,53]%>%pp.test()

data_updated[,53]%>%log()%>%diff()
plotly_simple(data_updated[,53]%>%log()%>%diff(),colnames(data_updated)[53])

data_updated[,53]%>%log()%>%diff()%>%adf.test()
data_updated[,53]%>%log()%>%diff()%>%kpss.test()
data_updated[,53]%>%log()%>%diff()%>%pp.test()

m1_stationary <- data_updated[,53]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,m1_stationary)
colnames(data_stationary)[53]<-colnames(data_updated)[53]

colnames(data_updated)[54]

plotly_simple(data_updated[,54],colnames(data_updated)[54])
data_updated[,54]%>%adf.test()
data_updated[,54]%>%kpss.test()
data_updated[,54]%>%pp.test()

data_updated[,54]%>%log()%>%diff()
plotly_simple(data_updated[,54]%>%log()%>%diff(),colnames(data_updated)[54])

data_updated[,54]%>%log()%>%diff()%>%adf.test()
data_updated[,54]%>%log()%>%diff()%>%kpss.test()
data_updated[,54]%>%log()%>%diff()%>%pp.test()

m2_stationary <- data_updated[,54]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,m2_stationary)
colnames(data_stationary)[54]<-colnames(data_updated)[54]

colnames(data_updated)[55]

plotly_simple(data_updated[,55],colnames(data_updated)[55])
data_updated[,55]%>%adf.test()
data_updated[,55]%>%kpss.test()
data_updated[,55]%>%pp.test()

data_updated[,55]%>%log()%>%diff()
plotly_simple(data_updated[,55]%>%log()%>%diff(),colnames(data_updated)[55])

data_updated[,55]%>%log()%>%diff()%>%adf.test()
data_updated[,55]%>%log()%>%diff()%>%kpss.test()
data_updated[,55]%>%log()%>%diff()%>%pp.test()

m3_stationary <- data_updated[,55]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,m3_stationary)
colnames(data_stationary)[55]<-colnames(data_updated)[55]

colnames(data_updated)[56]

plotly_simple(data_updated[,56],colnames(data_updated)[56])
data_updated[,56]%>%adf.test()
data_updated[,56]%>%kpss.test()
data_updated[,56]%>%pp.test()

data_updated[,56]%>%log()%>%diff()
plotly_simple(data_updated[,56]%>%log()%>%diff(),colnames(data_updated)[56])

data_updated[,56]%>%log()%>%diff()%>%adf.test()
data_updated[,56]%>%log()%>%diff()%>%kpss.test()
data_updated[,56]%>%log()%>%diff()%>%pp.test()

mb_stationary <- data_updated[,56]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,mb_stationary)
colnames(data_stationary)[56]<-colnames(data_updated)[56]

colnames(data_updated)[57]

plotly_simple(data_updated[,57],colnames(data_updated)[57])
data_updated[,57]%>%adf.test()
data_updated[,57]%>%kpss.test()
data_updated[,57]%>%pp.test()

data_updated[,57]%>%log()%>%diff()
plotly_simple(data_updated[,57]%>%log()%>%diff(),colnames(data_updated)[57])

data_updated[,57]%>%log()%>%diff()%>%adf.test()
data_updated[,57]%>%log()%>%diff()%>%kpss.test()
data_updated[,57]%>%log()%>%diff()%>%pp.test()

reser_depo_inst_stationary <- data_updated[,57]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,reser_depo_inst_stationary)
colnames(data_stationary)[57]<-colnames(data_updated)[57]

colnames(data_updated)[58]

plotly_simple(data_updated[,58],colnames(data_updated)[58])
data_updated[,58]%>%adf.test()
data_updated[,58]%>%kpss.test()
data_updated[,58]%>%pp.test()

data_updated[,58]%>%log()%>%diff()
plotly_simple(data_updated[,58]%>%log()%>%diff(),colnames(data_updated)[58])

data_updated[,58]%>%log()%>%diff()%>%adf.test()
data_updated[,58]%>%log()%>%diff()%>%kpss.test()
data_updated[,58]%>%log()%>%diff()%>%pp.test()

reser_depo_inst_nonborrowed_stationary <- data_updated[,58]%>%growth.rate(lag = 1,simple = T)
reser_depo_inst_nonborrowed_stationary <-reser_depo_inst_nonborrowed_stationary%>%as.zoo()
index(reser_depo_inst_nonborrowed_stationary) <- as.yearmon(index(reser_depo_inst_nonborrowed_stationary))  

data_stationary <- merge(data_stationary,reser_depo_inst_nonborrowed_stationary)
colnames(data_stationary)[58]<-colnames(data_updated)[58]

colnames(data_updated)[59]

plotly_simple(data_updated[,59],colnames(data_updated)[59])
data_updated[,59]%>%adf.test()
data_updated[,59]%>%kpss.test()
data_updated[,59]%>%pp.test()

data_updated[,59]%>%log()%>%diff()
plotly_simple(data_updated[,59]%>%log()%>%diff(),colnames(data_updated)[59])

data_updated[,59]%>%log()%>%diff()%>%adf.test()
data_updated[,59]%>%log()%>%diff()%>%kpss.test()
data_updated[,59]%>%log()%>%diff()%>%pp.test()

ppi_commodities_finished_goods_stationary <- data_updated[,59]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,ppi_commodities_finished_goods_stationary)
colnames(data_stationary)[59]<-colnames(data_updated)[59]

colnames(data_updated)[60]

plotly_simple(data_updated[,60],colnames(data_updated)[60])
data_updated[,60]%>%adf.test()
data_updated[,60]%>%kpss.test()
data_updated[,60]%>%pp.test()

data_updated[,60]%>%log()%>%diff()
plotly_simple(data_updated[,60]%>%log()%>%diff(),colnames(data_updated)[60])

data_updated[,60]%>%log()%>%diff()%>%adf.test()
data_updated[,60]%>%log()%>%diff()%>%kpss.test()
data_updated[,60]%>%log()%>%diff()%>%pp.test()

ppi_less_food_and_energy_stationary <- data_updated[,60]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,ppi_less_food_and_energy_stationary)
colnames(data_stationary)[60]<-colnames(data_updated)[60]

colnames(data_updated)[61]

plotly_simple(data_updated[,61],colnames(data_updated)[61])
data_updated[,61]%>%adf.test()
data_updated[,61]%>%kpss.test()
data_updated[,61]%>%pp.test()

data_updated[,61]%>%log()%>%diff()
plotly_simple(data_updated[,61]%>%log()%>%diff(),colnames(data_updated)[61])

data_updated[,61]%>%log()%>%diff()%>%adf.test()
data_updated[,61]%>%log()%>%diff()%>%kpss.test()
data_updated[,61]%>%log()%>%diff()%>%pp.test()

cpi_stationary <- data_updated[,61]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_stationary)
colnames(data_stationary)[61]<-colnames(data_updated)[61]

colnames(data_updated)[62]

plotly_simple(data_updated[,62],colnames(data_updated)[62])
data_updated[,62]%>%adf.test()
data_updated[,62]%>%kpss.test()
data_updated[,62]%>%pp.test()

data_updated[,62]%>%log()%>%diff()
plotly_simple(data_updated[,62]%>%log()%>%diff(),colnames(data_updated)[62])

data_updated[,62]%>%log()%>%diff()%>%adf.test()
data_updated[,62]%>%log()%>%diff()%>%kpss.test()
data_updated[,62]%>%log()%>%diff()%>%pp.test()

cpi_food_and_beverages_stationary <- data_updated[,62]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_food_and_beverages_stationary)
colnames(data_stationary)[62]<-colnames(data_updated)[62]

colnames(data_updated)[63]

plotly_simple(data_updated[,63],colnames(data_updated)[63])
data_updated[,63]%>%adf.test()
data_updated[,63]%>%kpss.test()
data_updated[,63]%>%pp.test()

data_updated[,63]%>%log()%>%diff()
plotly_simple(data_updated[,63]%>%log()%>%diff(),colnames(data_updated)[63])

data_updated[,63]%>%log()%>%diff()%>%adf.test()
data_updated[,63]%>%log()%>%diff()%>%kpss.test()
data_updated[,63]%>%log()%>%diff()%>%pp.test()

cpi_housing_stationary <- data_updated[,63]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_housing_stationary)
colnames(data_stationary)[63]<-colnames(data_updated)[63]

colnames(data_updated)[64]

plotly_simple(data_updated[,64],colnames(data_updated)[64])
data_updated[,64]%>%adf.test()
data_updated[,64]%>%kpss.test()
data_updated[,64]%>%pp.test()

data_updated[,64]%>%log()%>%diff()
plotly_simple(data_updated[,64]%>%log()%>%diff(),colnames(data_updated)[64])

data_updated[,64]%>%log()%>%diff()%>%adf.test()
data_updated[,64]%>%log()%>%diff()%>%kpss.test()
data_updated[,64]%>%log()%>%diff()%>%pp.test()

cpi_apparel_stationary <- data_updated[,64]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_apparel_stationary)
colnames(data_stationary)[64]<-colnames(data_updated)[64]

colnames(data_updated)[65]

plotly_simple(data_updated[,65],colnames(data_updated)[65])
data_updated[,65]%>%adf.test()
data_updated[,65]%>%kpss.test()
data_updated[,65]%>%pp.test()

data_updated[,65]%>%log()%>%diff()
plotly_simple(data_updated[,65]%>%log()%>%diff(),colnames(data_updated)[65])

data_updated[,65]%>%log()%>%diff()%>%adf.test()
data_updated[,65]%>%log()%>%diff()%>%kpss.test()
data_updated[,65]%>%log()%>%diff()%>%pp.test()

cpi_transportation_stationary <- data_updated[,65]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_transportation_stationary)
colnames(data_stationary)[65]<-colnames(data_updated)[65]

colnames(data_updated)[66]

plotly_simple(data_updated[,66],colnames(data_updated)[66])
data_updated[,66]%>%adf.test()
data_updated[,66]%>%kpss.test()
data_updated[,66]%>%pp.test()

data_updated[,66]%>%log()%>%diff()
plotly_simple(data_updated[,66]%>%log()%>%diff(),colnames(data_updated)[66])

data_updated[,66]%>%log()%>%diff()%>%adf.test()
data_updated[,66]%>%log()%>%diff()%>%kpss.test()
data_updated[,66]%>%log()%>%diff()%>%pp.test()

cpi_med_care_stationary <- data_updated[,66]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_med_care_stationary)
colnames(data_stationary)[66]<-colnames(data_updated)[66]

colnames(data_updated)[67]

plotly_simple(data_updated[,67],colnames(data_updated)[67])
data_updated[,67]%>%adf.test()
data_updated[,67]%>%kpss.test()
data_updated[,67]%>%pp.test()

data_updated[,67]%>%log()%>%diff()
plotly_simple(data_updated[,67]%>%log()%>%diff(),colnames(data_updated)[67])

data_updated[,67]%>%log()%>%diff()%>%adf.test()
data_updated[,67]%>%log()%>%diff()%>%kpss.test()
data_updated[,67]%>%log()%>%diff()%>%pp.test()

cpi_commodities_stationary <- data_updated[,67]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_commodities_stationary)
colnames(data_stationary)[67]<-colnames(data_updated)[67]

colnames(data_updated)[68]

plotly_simple(data_updated[,68],colnames(data_updated)[68])
data_updated[,68]%>%adf.test()
data_updated[,68]%>%kpss.test()
data_updated[,68]%>%pp.test()

data_updated[,68]%>%log()%>%diff()
plotly_simple(data_updated[,68]%>%log()%>%diff(),colnames(data_updated)[68])

data_updated[,68]%>%log()%>%diff()%>%adf.test()
data_updated[,68]%>%log()%>%diff()%>%kpss.test()
data_updated[,68]%>%log()%>%diff()%>%pp.test()

cpi_services_stationary <- data_updated[,68]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_services_stationary)
colnames(data_stationary)[68]<-colnames(data_updated)[68]

colnames(data_updated)[69]

plotly_simple(data_updated[,69],colnames(data_updated)[69])
data_updated[,69]%>%adf.test()
data_updated[,69]%>%kpss.test()
data_updated[,69]%>%pp.test()

data_updated[,69]%>%log()%>%diff()
plotly_simple(data_updated[,69]%>%log()%>%diff(),colnames(data_updated)[69])

data_updated[,69]%>%log()%>%diff()%>%adf.test()
data_updated[,69]%>%log()%>%diff()%>%kpss.test()
data_updated[,69]%>%log()%>%diff()%>%pp.test()

cpi_all_less_food_stationary <- data_updated[,69]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_all_less_food_stationary)
colnames(data_stationary)[69]<-colnames(data_updated)[69]

colnames(data_updated)[70]

plotly_simple(data_updated[,70],colnames(data_updated)[70])
data_updated[,70]%>%adf.test()
data_updated[,70]%>%kpss.test()
data_updated[,70]%>%pp.test()

data_updated[,70]%>%log()%>%diff()
plotly_simple(data_updated[,70]%>%log()%>%diff(),colnames(data_updated)[70])

data_updated[,70]%>%log()%>%diff()%>%adf.test()
data_updated[,70]%>%log()%>%diff()%>%kpss.test()
data_updated[,70]%>%log()%>%diff()%>%pp.test()

cpi_all_less_shelter_stationary <- data_updated[,70]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_all_less_shelter_stationary)
colnames(data_stationary)[70]<-colnames(data_updated)[70]

colnames(data_updated)[71]

plotly_simple(data_updated[,71],colnames(data_updated)[71])
data_updated[,71]%>%adf.test()
data_updated[,71]%>%kpss.test()
data_updated[,71]%>%pp.test()

data_updated[,71]%>%log()%>%diff()
plotly_simple(data_updated[,71]%>%log()%>%diff(),colnames(data_updated)[71])

data_updated[,71]%>%log()%>%diff()%>%adf.test()
data_updated[,71]%>%log()%>%diff()%>%kpss.test()
data_updated[,71]%>%log()%>%diff()%>%pp.test()

cpi_all_less_med_care_stationary <- data_updated[,71]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_all_less_med_care_stationary)
colnames(data_stationary)[71]<-colnames(data_updated)[71]

colnames(data_updated)[72]

plotly_simple(data_updated[,72],colnames(data_updated)[72])
data_updated[,72]%>%adf.test()
data_updated[,72]%>%kpss.test()
data_updated[,72]%>%pp.test()

data_updated[,72]%>%log()%>%diff()
plotly_simple(data_updated[,72]%>%log()%>%diff(),colnames(data_updated)[72])

data_updated[,72]%>%log()%>%diff()%>%adf.test()
data_updated[,72]%>%log()%>%diff()%>%kpss.test()
data_updated[,72]%>%log()%>%diff()%>%pp.test()

cpi_all_less_food_and_energy_stationary <- data_updated[,72]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,cpi_all_less_food_and_energy_stationary)
colnames(data_stationary)[72]<-colnames(data_updated)[72]

colnames(data_updated)[73]

plotly_simple(data_updated[,73],colnames(data_updated)[73])
data_updated[,73]%>%adf.test()
data_updated[,73]%>%kpss.test()
data_updated[,73]%>%pp.test()

data_updated[,73]%>%log()%>%diff()
plotly_simple(data_updated[,73]%>%log()%>%diff(),colnames(data_updated)[73])

data_updated[,73]%>%log()%>%diff()%>%adf.test()
data_updated[,73]%>%log()%>%diff()%>%kpss.test()
data_updated[,73]%>%log()%>%diff()%>%pp.test()

pce_stationary <- data_updated[,73]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pce_stationary)
colnames(data_stationary)[73]<-colnames(data_updated)[73]

colnames(data_updated)[74]

plotly_simple(data_updated[,74],colnames(data_updated)[74])
data_updated[,74]%>%adf.test()
data_updated[,74]%>%kpss.test()
data_updated[,74]%>%pp.test()

data_updated[,74]%>%log()%>%diff()
plotly_simple(data_updated[,74]%>%log()%>%diff(),colnames(data_updated)[74])

data_updated[,74]%>%log()%>%diff()%>%adf.test()
data_updated[,74]%>%log()%>%diff()%>%kpss.test()
data_updated[,74]%>%log()%>%diff()%>%pp.test()

pce_less_food_and_energy_stationary <- data_updated[,74]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pce_less_food_and_energy_stationary)
colnames(data_stationary)[74]<-colnames(data_updated)[74]

colnames(data_updated)[75]

plotly_simple(data_updated[,75],colnames(data_updated)[75])
data_updated[,75]%>%adf.test()
data_updated[,75]%>%kpss.test()
data_updated[,75]%>%pp.test()

data_updated[,75]%>%log()%>%diff()
plotly_simple(data_updated[,75]%>%log()%>%diff(),colnames(data_updated)[75])

data_updated[,75]%>%log()%>%diff()%>%adf.test()
data_updated[,75]%>%log()%>%diff()%>%kpss.test()
data_updated[,75]%>%log()%>%diff()%>%pp.test()

pce_durables_stationary <- data_updated[,75]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pce_durables_stationary)
colnames(data_stationary)[75]<-colnames(data_updated)[75]

colnames(data_updated)[76]

plotly_simple(data_updated[,76],colnames(data_updated)[76])
data_updated[,76]%>%adf.test()
data_updated[,76]%>%kpss.test()
data_updated[,76]%>%pp.test()

data_updated[,76]%>%log()%>%diff()
plotly_simple(data_updated[,76]%>%log()%>%diff(),colnames(data_updated)[76])

data_updated[,76]%>%log()%>%diff()%>%adf.test()
data_updated[,76]%>%log()%>%diff()%>%kpss.test()
data_updated[,76]%>%log()%>%diff()%>%pp.test()

pce_nondurables_stationary <- data_updated[,76]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pce_nondurables_stationary)
colnames(data_stationary)[76]<-colnames(data_updated)[76]

colnames(data_updated)[77]

plotly_simple(data_updated[,77],colnames(data_updated)[77])
data_updated[,77]%>%adf.test()
data_updated[,77]%>%kpss.test()
data_updated[,77]%>%pp.test()

data_updated[,77]%>%log()%>%diff()
plotly_simple(data_updated[,77]%>%log()%>%diff(),colnames(data_updated)[77])

data_updated[,77]%>%log()%>%diff()%>%adf.test()
data_updated[,77]%>%log()%>%diff()%>%kpss.test()
data_updated[,77]%>%log()%>%diff()%>%pp.test()

pce_services_stationary <- data_updated[,77]%>%log()%>%diff()*100

data_stationary <- merge(data_stationary,pce_services_stationary)
colnames(data_stationary)[77]<-colnames(data_updated)[77]

colnames(data_updated)[78]

plotly_simple(data_updated[,78],colnames(data_updated)[78])
data_updated[,78]%>%adf.test()
data_updated[,78]%>%kpss.test()
data_updated[,78]%>%pp.test()

data_updated[,78]%>%log()%>%diff()
plotly_simple(data_updated[,78]%>%growth.rate(lag = 1,simple = T),colnames(data_updated)[78])

data_updated[,78]%>%growth.rate(lag = 1,simple = T)%>%adf.test()
data_updated[,78]%>%growth.rate(lag = 1,simple = T)%>%kpss.test()
data_updated[,78]%>%growth.rate(lag = 1,simple = T)%>%pp.test()

federal_surplus_stationary <- data_updated[,78]%>%growth.rate(lag = 1,simple = T)
federal_surplus_stationary <- federal_surplus_stationary%>%as.zoo()
index(federal_surplus_stationary) <- as.yearmon(index(federal_surplus_stationary))

data_stationary <- merge(data_stationary,federal_surplus_stationary)
colnames(data_stationary)[78]<-colnames(data_updated)[78]

data_stationary <- data_stationary[complete.cases(data_stationary),]
index(data_stationary)<-as.Date(index(data_stationary))
# data_transformed <- scale(data_stationary)
# index(data_transformed)<-as.Date(index(data_transformed))
# colMeans(data_transformed)
# sapply(data_transformed, sd)
# 
# final_data <- sapply(data_transformed, ts_clean_vec,period=12)
# final_data <- as.zoo(final_data)
# index(final_data) <- index(data_transformed)
# final_data
# final_data%>%tail(20)
# 
# final_data%>%colMeans()
# 
# sapply(final_data, sd)
# 
setwd("C://Users//Italo//Documents//PhD//PhD//First")
# 
write.csv(as.data.frame(data_stationary),"stationary_data_for_macro_factors.csv")
# index(ffr)<-as.Date(index(ffr))

write.csv(as.data.frame(ffr),"ffr.csv")

#data_transformed_attributes<-data_transformed%>%attributes()
