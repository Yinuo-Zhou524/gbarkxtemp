setwd("E://LSU//research//gbark_temperature")
library(MASS)
library(leaps)
library(dplyr)
latempspring1<-read.csv("temperature x gbark dry down 202403.csv")
#colnames(latempspring1)<-c("tag","date","sp",	"treatment"	,"segment_location",	"segment"	,"mass_w_parafilm",	"hour",	"minute")

#latempspring1$tag<-as.numeric(latempspring1$tag)
#latempspring1$sp[is.na(latempspring1$sp)] = "NA"
#latempspring1$ob<-seq(1,length(latempspring1$tag),1)


#make date and time readable
latempspring1$time<-paste(latempspring1$hour,latempspring1$minute,sep=":")
latempspring1$date_and_time1<-paste(latempspring1$measurement_date,latempspring1$time,sep=" ")
latempspring1$rtime<-as.POSIXct(latempspring1$date_and_time1,format="%Y.%m.%d %H:%M")

###calculate hours of dry down 
###What is end included time?
latempspring_dry1<-data.frame(matrix(nrow=0,ncol=3,data=NA))
names(latempspring_dry1)<-c("segment","start_time","end_time")


for(i in unique(latempspring1$segment[!is.na(latempspring1$mass_w_parafilm)]))
{latempspring_dry1[i,1]<-paste(i)
latempspring_dry1[i,2]<-paste(min(latempspring1$rtime[latempspring1$segment==i]))
#latempspring_dry1[i,3]<-paste(max(latempspring1$rtime[latempspring1$exclude==0 & latempspring1$tag==i]))
latempspring_dry1[i,3]<-paste(max(latempspring1$rtime[latempspring1$segment==i]))}

latempspring_dry1$start_time<-as.POSIXct(latempspring_dry1$start_time,format="%Y-%m-%d %H:%M")
#latempspring_dry1$end_included_time<-as.POSIXct(latempspring_dry1$end_included_time,format="%Y-%m-%d %H:%M")
latempspring_dry1$end_time<-as.POSIXct(latempspring_dry1$end_time,format="%Y-%m-%d %H:%M")
latempspring_dry1$total_time<-as.numeric(latempspring_dry1$end_time-latempspring_dry1$start_time)


#
latempspring<-merge(latempspring1,latempspring_dry1,by="segment")
latempspring$dry_time<-(as.numeric(latempspring$rtime-latempspring$start_time))/3600


#get the R^2 for each segment
x<-latempspring$mass_w_parafilm[latempspring$segment=="021A"]
y<-as.numeric(latempspring$dry_time[latempspring$segment=="021A"])

for(i in 1:length(x))
{xi<-x[x<= c(sort(x)[length(x)+1-i])]
yi<-y[y>= c(sort(y,decreasing=T)[length(y)+1-i])]
lm<-lm(xi~yi)
print(summary(lm)$r.squared)
}

##Make it for each segment, a function?
r2_exclude=function(a){
  x<-latempspring$mass_w_parafilm[latempspring$segment==a]
  y<-as.numeric(latempspring$dry_time[latempspring$segment==a])
  
  for(i in 1:length(x))
  {xi<-x[x<= c(sort(x)[length(x)+1-i])]
  yi<-y[y>= c(sort(y,decreasing=T)[length(y)+1-i])]
  lm<-lm(xi~yi)
  print(summary(lm)$r.squared)
  }
}

#R2>0.99

###Input the excluded data
latempspring2<-read.csv("temperature x gbark dry down 202403_exclude.csv")
colnames(latempspring2)<-c("tag","collection_date","round","date","sp",	"treatment"	,"treatment_ID",	"segment"	,	"exclude", "mass_w_parafilm","hour","minute")


#make date and time readable
latempspring2$time<-paste(latempspring2$hour,latempspring2$minute,sep=":")
latempspring2$date_and_time1<-paste(latempspring2$date,latempspring2$time,sep=" ")
latempspring2$rtime<-as.POSIXct(latempspring2$date_and_time1,format="%Y.%m.%d %H:%M")

###calculate hours of dry down 
###What is end included time?
latempspring_dry2<-data.frame(matrix(nrow=0,ncol=4,data=NA))
names(latempspring_dry2)<-c("segment","start_time","end_included_time","end_time")


for(i in unique(latempspring2$segment[!is.na(latempspring2$mass_w_parafilm)]))
{latempspring_dry2[i,1]<-paste(i)
latempspring_dry2[i,2]<-paste(min(latempspring2$rtime[latempspring2$segment==i]))
latempspring_dry2[i,3]<-paste(max(latempspring2$rtime[latempspring2$exclude==0 & latempspring2$segment==i]))
latempspring_dry2[i,4]<-paste(max(latempspring2$rtime[latempspring2$segment==i]))}

latempspring_dry2$start_time<-as.POSIXct(latempspring_dry2$start_time,format="%Y-%m-%d %H:%M")
latempspring_dry2$end_included_time<-as.POSIXct(latempspring_dry2$end_included_time,format="%Y-%m-%d %H:%M")
latempspring_dry2$end_time<-as.POSIXct(latempspring_dry2$end_time,format="%Y-%m-%d %H:%M")
latempspring_dry2$total_time<-as.numeric(latempspring_dry2$end_time-latempspring_dry2$start_time)


#
latempspring<-merge(latempspring2,latempspring_dry2,by="segment")
latempspring$dry_time<-(as.numeric(latempspring$rtime-latempspring$start_time))/3600

#exlude outlier points within stems that are probably due to wet stems or lack of water in the segment
latempspring<-latempspring[latempspring$exclude!=1,]



#input air temp and RH
latempspring_tRH10C<-read.csv("temp x gbark tRH 10C spring.csv")
latempspring_tRH25C<-read.csv("temp x gbark tRH 25C spring.csv")
latempspring_tRH35C<-read.csv("temp x gbark tRH 35C spring.csv")
latempspring_tRH45C<-read.csv("temp x gbark tRH 45C spring.csv")
latempspring_tRH55C<-read.csv("temp x gbark tRH 55C spring.csv")

#separate segment by temperature id 
latempspring_A<-latempspring[latempspring$treatment_ID=="A",]
latempspring_B<-latempspring[latempspring$treatment_ID=="B",]
latempspring_C<-latempspring[latempspring$treatment_ID=="C",]
latempspring_D<-latempspring[latempspring$treatment_ID=="D",]
latempspring_E<-latempspring[latempspring$treatment_ID=="E",]

#Sys.setenv(TZ="CDT")
#Sys.getenv("TZ")
#options(tz="CDT")

###10C---A----
latempspring_tRH10C$rtime<-as.POSIXct(latempspring_tRH10C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference

#Log10 ew =  10.79574 (1-273.16/T)                                                          [3]
#- 5.02800 Log10(T/273.16)
#+ 1.50475 10-4 (1 - 10^(-8.2969*(T/273.16-1)))
#+ 0.42873 10-3 (10^(+4.76955*(1-273.16/T)) - 1)
#+ 0.78614
#with T in [K] and ew in [hPa] (Goff, 1957)

#latempspring_tRH10C$svp_1<-0.6135*exp(17.502*latempspring_tRH10C$Temp_Avg/(240.97+latempspring_tRH10C$Temp_Avg)) #from LI610 manual
latempspring_tRH10C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latempspring_tRH10C$Temp_Avg+273.15)) - 5.02800*log10((latempspring_tRH10C$Temp_Avg+273.15)/273.16)+
                                      1.50475*10^(-4)*(1-10^(-8.2969*((latempspring_tRH10C$Temp_Avg+273.15)/273.16-1)))+
                                      0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latempspring_tRH10C$Temp_Avg+273.15))))+
                                      0.78614)) # Golf 1957 unit:kPa
#latempspring_tRH10C$svp_1<-0.1*(10^(-7.90298*(373.16/(latempspring_tRH10C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latempspring_tRH10C$Temp_Avg+273.15))-
#                           1.3816*10^(-7)*(10^(11.344*(1-(latempspring_tRH10C$Temp_Avg+273.15)/373.16))-1)+
#                           8.1328*10^(-3)*(10^(-3.49149*(373.16/(latempspring_tRH10C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latempspring_tRH10C$vpd_1<-((100-latempspring_tRH10C$RH_Avg)/100)*latempspring_tRH10C$svp_1
latempspring_tRH10C$D1<-(latempspring_tRH10C$vpd_1/101.3)
#mass flow correction
latempspring_tRH10C$ws_1<-(latempspring_tRH10C$svp_1/101.3)
latempspring_tRH10C$wa_1<-(latempspring_tRH10C$RH_Avg/100)*latempspring_tRH10C$ws_1
latempspring_tRH10C$wmean_1<- (latempspring_tRH10C$ws_1+latempspring_tRH10C$wa_1)/2
latempspring_tRH10C$massflow_factor_1<- 1 - latempspring_tRH10C$wmean_1


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latempspring_regs2A<-data.frame(matrix(nrow=0,ncol=11,data=NA))
names(latempspring_regs2A)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd")

for (i in unique(latempspring_A$segment))
{
  lm1<-lm(latempspring_A$mass_w_parafilm[latempspring_A$segment==i]~as.numeric(latempspring_A$dry_time[latempspring_A$segment==i]))
  
  #print(lm1$coefficients[2])
  latempspring_regs2A[i,1]<-paste(i)
  latempspring_regs2A[i,2]<-paste(lm1$coefficients[1])
  latempspring_regs2A[i,3]<-paste(lm1$coefficients[2])
  latempspring_regs2A[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latempspring_tRH10C$rtime<-as.numeric(latempspring_tRH10C$rtime)
for (i in unique(latempspring_A$segment))
{
  latempspring_regs2A[i,5]<-paste(mean(latempspring_tRH10C$D1[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,6]<-paste(sd(latempspring_tRH10C$D1[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,7]<-paste(mean(latempspring_tRH10C$Temp_Avg[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,8]<-paste(sd(latempspring_tRH10C$Temp_Avg[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,9]<-paste(mean(latempspring_tRH10C$massflow_factor_1[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,10]<-paste(sd(latempspring_tRH10C$massflow_factor_1[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
  latempspring_regs2A[i,11]<-paste(mean(latempspring_tRH10C$vpd_1[latempspring_tRH10C$rtime>=latempspring_A$start_time[latempspring_A$segment==i]&latempspring_tRH10C$rtime<=latempspring_A$end_time[latempspring_A$segment==i]],na.rm=T))
}


###25C---B----
latempspring_tRH25C$rtime<-as.POSIXct(latempspring_tRH25C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latempspring_tRH25C$svp_1<-0.6135*exp(17.502*latempspring_tRH25C$Temp_Avg/(240.97+latempspring_tRH25C$Temp_Avg)) #from LI610 manual
latempspring_tRH25C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latempspring_tRH25C$Temp_Avg+273.15)) - 5.02800*log10((latempspring_tRH25C$Temp_Avg+273.15)/273.16)+
                                      1.50475*10^(-4)*(1-10^(-8.2969*((latempspring_tRH25C$Temp_Avg+273.15)/273.16-1)))+
                                      0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latempspring_tRH25C$Temp_Avg+273.15))))+
                                      0.78614))#Golf 1957
#latempspring_tRH25C$svp_1<-0.1*(10^(-7.90298*(373.16/(latempspring_tRH25C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latempspring_tRH25C$Temp_Avg+273.15))-
#                          1.3816*10^(-7)*(10^(11.344*(1-(latempspring_tRH25C$Temp_Avg+273.15)/373.16))-1)+
#                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latempspring_tRH25C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latempspring_tRH25C$vpd_1<-((100-latempspring_tRH25C$RH_Avg)/100)*latempspring_tRH25C$svp_1
latempspring_tRH25C$D1<-(latempspring_tRH25C$vpd_1/101.3)

#mass flow correction
latempspring_tRH25C$ws_1<-(latempspring_tRH25C$svp_1/101.3)
latempspring_tRH25C$wa_1<-(latempspring_tRH25C$RH_Avg/100)*latempspring_tRH25C$ws_1
latempspring_tRH25C$wmean_1<- (latempspring_tRH25C$ws_1+latempspring_tRH25C$wa_1)/2
latempspring_tRH25C$massflow_factor_1<- 1 - latempspring_tRH25C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latempspring_regs2B<-data.frame(matrix(nrow=0,ncol=11,data=NA))
names(latempspring_regs2B)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd")

for (i in unique(latempspring_B$segment))
{
  lm1<-lm(latempspring_B$mass_w_parafilm[latempspring_B$segment==i]~as.numeric(latempspring_B$dry_time[latempspring_B$segment==i]))
  
  #print(lm1$coefficients[2])
  latempspring_regs2B[i,1]<-paste(i)
  latempspring_regs2B[i,2]<-paste(lm1$coefficients[1])
  latempspring_regs2B[i,3]<-paste(lm1$coefficients[2])
  latempspring_regs2B[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latempspring_tRH25C$rtime<-as.numeric(latempspring_tRH25C$rtime)
for (i in unique(latempspring_B$segment))
{
  latempspring_regs2B[i,5]<-paste(mean(latempspring_tRH25C$D1[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,6]<-paste(sd(latempspring_tRH25C$D1[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,7]<-paste(mean(latempspring_tRH25C$Temp_Avg[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,8]<-paste(sd(latempspring_tRH25C$Temp_Avg[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,9]<-paste(mean(latempspring_tRH25C$massflow_factor_1[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,10]<-paste(sd(latempspring_tRH25C$massflow_factor_1[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
  latempspring_regs2B[i,11]<-paste(mean(latempspring_tRH25C$vpd_1[latempspring_tRH25C$rtime>=latempspring_B$start_time[latempspring_B$segment==i]&latempspring_tRH25C$rtime<=latempspring_B$end_time[latempspring_B$segment==i]],na.rm=T))
}

###35C---C----
latempspring_tRH35C$rtime<-as.POSIXct(latempspring_tRH35C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latempspring_tRH35C$svp_1<-0.6135*exp(17.502*latempspring_tRH35C$Temp_Avg/(240.97+latempspring_tRH35C$Temp_Avg)) #from LI610 manual
latempspring_tRH35C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latempspring_tRH35C$Temp_Avg+273.15)) - 5.02800*log10((latempspring_tRH35C$Temp_Avg+273.15)/273.16)+
                                      1.50475*10^(-4)*(1-10^(-8.2969*((latempspring_tRH35C$Temp_Avg+273.15)/273.16-1)))+
                                      0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latempspring_tRH35C$Temp_Avg+273.15))))+
                                      0.78614))#Golf 1957
#latempspring_tRH35C$svp_1<-0.1*(10^(-7.90298*(373.16/(latempspring_tRH35C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latempspring_tRH35C$Temp_Avg+273.15))-
#                          1.3816*10^(-7)*(10^(11.344*(1-(latempspring_tRH35C$Temp_Avg+273.15)/373.16))-1)+
#                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latempspring_tRH35C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latempspring_tRH35C$vpd_1<-((100-latempspring_tRH35C$RH_Avg)/100)*latempspring_tRH35C$svp_1
latempspring_tRH35C$D1<-(latempspring_tRH35C$vpd_1/101.3)
#mass flow correction
latempspring_tRH35C$ws_1<-(latempspring_tRH35C$svp_1/101.3)
latempspring_tRH35C$wa_1<-(latempspring_tRH35C$RH_Avg/100)*latempspring_tRH35C$ws_1
latempspring_tRH35C$wmean_1<- (latempspring_tRH35C$ws_1+latempspring_tRH35C$wa_1)/2
latempspring_tRH35C$massflow_factor_1<- 1 - latempspring_tRH35C$wmean_1


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latempspring_regs2C<-data.frame(matrix(nrow=0,ncol=11,data=NA))
names(latempspring_regs2C)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd")

for (i in unique(latempspring_C$segment))
{
  lm1<-lm(latempspring_C$mass_w_parafilm[latempspring_C$segment==i]~as.numeric(latempspring_C$dry_time[latempspring_C$segment==i]))
  
  #print(lm1$coefficients[2])
  latempspring_regs2C[i,1]<-paste(i)
  latempspring_regs2C[i,2]<-paste(lm1$coefficients[1])
  latempspring_regs2C[i,3]<-paste(lm1$coefficients[2])
  latempspring_regs2C[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latempspring_tRH35C$rtime<-as.numeric(latempspring_tRH35C$rtime)
for (i in unique(latempspring_C$segment))
{
  latempspring_regs2C[i,5]<-paste(mean(latempspring_tRH35C$D1[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,6]<-paste(sd(latempspring_tRH35C$D1[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,7]<-paste(mean(latempspring_tRH35C$Temp_Avg[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,8]<-paste(sd(latempspring_tRH35C$Temp_Avg[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,9]<-paste(mean(latempspring_tRH35C$massflow_factor_1[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,10]<-paste(sd(latempspring_tRH35C$massflow_factor_1[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
  latempspring_regs2C[i,11]<-paste(mean(latempspring_tRH35C$vpd_1[latempspring_tRH35C$rtime>=latempspring_C$start_time[latempspring_C$segment==i]&latempspring_tRH35C$rtime<=latempspring_C$end_time[latempspring_C$segment==i]],na.rm=T))
}

###45C---D----
latempspring_tRH45C$rtime<-as.POSIXct(latempspring_tRH45C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latempspring_tRH45C$svp_1<-0.6135*exp(17.502*latempspring_tRH45C$Temp_Avg/(240.97+latempspring_tRH45C$Temp_Avg)) #from LI610 manual
latempspring_tRH45C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latempspring_tRH45C$Temp_Avg+273.15)) - 5.02800*log10((latempspring_tRH45C$Temp_Avg+273.15)/273.16)+
                                      1.50475*10^(-4)*(1-10^(-8.2969*((latempspring_tRH45C$Temp_Avg+273.15)/273.16-1)))+
                                      0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latempspring_tRH45C$Temp_Avg+273.15))))+
                                      0.78614))
#latempspring_tRH45C$svp_1<-0.1*(10^(-7.90298*(373.16/(latempspring_tRH45C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latempspring_tRH45C$Temp_Avg+273.15))-
#                          1.3816*10^(-7)*(10^(11.344*(1-(latempspring_tRH45C$Temp_Avg+273.15)/373.16))-1)+
#                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latempspring_tRH45C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latempspring_tRH45C$vpd_1<-((100-latempspring_tRH45C$RH_Avg)/100)*latempspring_tRH45C$svp_1
latempspring_tRH45C$D1<-(latempspring_tRH45C$vpd_1/101.3)
#mass flow correction
latempspring_tRH45C$ws_1<-(latempspring_tRH45C$svp_1/101.3)
latempspring_tRH45C$wa_1<-(latempspring_tRH45C$RH_Avg/100)*latempspring_tRH45C$ws_1
latempspring_tRH45C$wmean_1<- (latempspring_tRH45C$ws_1+latempspring_tRH45C$wa_1)/2
latempspring_tRH45C$massflow_factor_1<- 1 - latempspring_tRH45C$wmean_1


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latempspring_regs2D<-data.frame(matrix(nrow=0,ncol=11,data=NA))
names(latempspring_regs2D)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd")

for (i in unique(latempspring_D$segment))
{
  lm1<-lm(latempspring_D$mass_w_parafilm[latempspring_D$segment==i]~as.numeric(latempspring_D$dry_time[latempspring_D$segment==i]))
  
  #print(lm1$coefficients[2])
  latempspring_regs2D[i,1]<-paste(i)
  latempspring_regs2D[i,2]<-paste(lm1$coefficients[1])
  latempspring_regs2D[i,3]<-paste(lm1$coefficients[2])
  latempspring_regs2D[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latempspring_tRH45C$rtime<-as.numeric(latempspring_tRH45C$rtime)
for (i in unique(latempspring_D$segment))
{
  latempspring_regs2D[i,5]<-paste(mean(latempspring_tRH45C$D1[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,6]<-paste(sd(latempspring_tRH45C$D1[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,7]<-paste(mean(latempspring_tRH45C$Temp_Avg[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,8]<-paste(sd(latempspring_tRH45C$Temp_Avg[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,9]<-paste(mean(latempspring_tRH45C$massflow_factor_1[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,10]<-paste(sd(latempspring_tRH45C$massflow_factor_1[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
  latempspring_regs2D[i,11]<-paste(mean(latempspring_tRH45C$vpd_1[latempspring_tRH45C$rtime>=latempspring_D$start_time[latempspring_D$segment==i]&latempspring_tRH45C$rtime<=latempspring_D$end_time[latempspring_D$segment==i]],na.rm=T))
}


###55C---E----
latempspring_tRH55C$rtime<-as.POSIXct(latempspring_tRH55C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latempspring_tRH55C$svp_1<-0.6135*exp(17.502*latempspring_tRH55C$Temp_Avg/(240.97+latempspring_tRH55C$Temp_Avg)) #from LI610 manual
latempspring_tRH55C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latempspring_tRH55C$Temp_Avg+273.15)) - 5.02800*log10((latempspring_tRH55C$Temp_Avg+273.15)/273.16)+
                                      1.50475*10^(-4)*(1-10^(-8.2969*((latempspring_tRH55C$Temp_Avg+273.15)/273.16-1)))+
                                      0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latempspring_tRH55C$Temp_Avg+273.15))))+
                                      0.78614))
#latempspring_tRH55C$svp_1<-0.1*(10^(-7.90298*(373.16/(latempspring_tRH55C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latempspring_tRH55C$Temp_Avg+273.15))-
#                          1.3816*10^(-7)*(10^(11.344*(1-(latempspring_tRH55C$Temp_Avg+273.15)/373.16))-1)+
#                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latempspring_tRH55C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latempspring_tRH55C$vpd_1<-((100-latempspring_tRH55C$RH_Avg)/100)*latempspring_tRH55C$svp_1
latempspring_tRH55C$D1<-(latempspring_tRH55C$vpd_1/101.3)
#mass flow correction
latempspring_tRH55C$ws_1<-(latempspring_tRH55C$svp_1/101.3)
latempspring_tRH55C$wa_1<-(latempspring_tRH55C$RH_Avg/100)*latempspring_tRH55C$ws_1
latempspring_tRH55C$wmean_1<- (latempspring_tRH55C$ws_1+latempspring_tRH55C$wa_1)/2
latempspring_tRH55C$massflow_factor_1<- 1 - latempspring_tRH55C$wmean_1


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latempspring_regs2E<-data.frame(matrix(nrow=0,ncol=11,data=NA))
names(latempspring_regs2E)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd")

for (i in unique(latempspring_E$segment))
{
  lm1<-lm(latempspring_E$mass_w_parafilm[latempspring_E$segment==i]~as.numeric(latempspring_E$dry_time[latempspring_E$segment==i]))
  
  #print(lm1$coefficients[2])
  latempspring_regs2E[i,1]<-paste(i)
  latempspring_regs2E[i,2]<-paste(lm1$coefficients[1])
  latempspring_regs2E[i,3]<-paste(lm1$coefficients[2])
  latempspring_regs2E[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latempspring_tRH55C$rtime<-as.numeric(latempspring_tRH55C$rtime)
for (i in unique(latempspring_E$segment))
{
  latempspring_regs2E[i,5]<-paste(mean(latempspring_tRH55C$D1[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,6]<-paste(sd(latempspring_tRH55C$D1[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,7]<-paste(mean(latempspring_tRH55C$Temp_Avg[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,8]<-paste(sd(latempspring_tRH55C$Temp_Avg[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,9]<-paste(mean(latempspring_tRH55C$massflow_factor_1[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,10]<-paste(sd(latempspring_tRH55C$massflow_factor_1[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
  latempspring_regs2E[i,11]<-paste(mean(latempspring_tRH55C$vpd_1[latempspring_tRH55C$rtime>=latempspring_E$start_time[latempspring_E$segment==i]&latempspring_tRH55C$rtime<=latempspring_E$end_time[latempspring_E$segment==i]],na.rm=T))
}


#----

latempspring_regs2<-rbind(latempspring_regs2A,latempspring_regs2B,latempspring_regs2C,latempspring_regs2D,latempspring_regs2E)
#calculate water loss rate in mols per second
# divide by 18.01528 to convert from grams to moles
latempspring_regs2$slope<-as.numeric(latempspring_regs2$slope)
latempspring_regs2$D<-as.numeric(latempspring_regs2$D)
latempspring_regs2$mols_per_sec<-(-latempspring_regs2$slope/3600)/18.01528
latempspring_regs2$mf<-as.numeric(latempspring_regs2$mf)


##merge lm coefs /with stem dimensions
latempspring_dims<-read.csv("temperature x gbark traits_spring.csv")

df<-unique(latempspring[,c("segment","treatment")])
latempspring_dims<-merge(latempspring_dims,df,all.y = FALSE)
latempspring_dims$sampleID<-latempspring_dims$segment

##calculate surface area in m2 as cylinder
latempspring_dims$diam<-(latempspring_dims$basal_diam1+latempspring_dims$basal_diam2+latempspring_dims$distal_diam1+latempspring_dims$distal_diam2)/4
latempspring_dims$wood_diam<-(latempspring_dims$basal_xylem_diam1+latempspring_dims$basal_xylem_diam2+latempspring_dims$distal_xylem_diam1+latempspring_dims$distal_xylem_diam2)/4
#latempspring_dims$exp_length_try1<-(latempspring_dims$exp_length1+latempspring_dims$exp_length2)/2
#latempspring_dims$exp_length_try2<-ifelse(is.na(latempspring_dims$exp_length_try1),latempspring_dims$exp_length1,latempspring_dims$exp_length_try1)
#latempspring_dims$exp_length<-ifelse(is.na(latempspring_dims$exp_length_try2),latempspring_dims$exp_length2,latempspring_dims$exp_length_try2)
latempspring_dims$sa<-pi*((latempspring_dims$diam)/1000)*(latempspring_dims$exposed_length/1000)

##calculate bark thickness
latempspring_dims$bark_thickness<-(latempspring_dims$diam-latempspring_dims$wood_diam)/2

#calculate volume and initial water content, and stem density
latempspring_dims$vol<-latempspring_dims$total_length*(latempspring_dims$diam/2)^2
latempspring_dims$water_content<-(latempspring_dims$fresh_mass-latempspring_dims$dry_mass)/latempspring_dims$vol
latempspring_dims$density<-latempspring_dims$dry_mass/latempspring_dims$vol


#calculate lenticel density and size
latempspring_dims$lenticel_size<-rowMeans(latempspring_dims[,c(10:19)],na.rm = T)
latempspring_dims$LD<-latempspring_dims$n_lenticels/(latempspring_dims$diam*latempspring_dims$length_counted_lenticel)
##
#names(latempspring_regs2)
#names(latempspring_dims)
latempspring_gs<-merge(latempspring_regs2,latempspring_dims,all.x=T)


#calculate E and gs over entire drying time
latempspring_gs$E<-(latempspring_gs$mols_per_sec*1000)/latempspring_gs$sa #mmol/s/m2
latempspring_gs$gs<-latempspring_gs$E*latempspring_gs$mf/latempspring_gs$D #mmol/s/m2

#calculate bark thickness ratio
latempspring_gs$bark_thickness_ratio<-latempspring_gs$bark_thickness/latempspring_gs$diam

seg_obs_spring <- latempspring %>%
  group_by(segment) %>%
  summarise(
    start_rtime = min(rtime, na.rm = TRUE),
    end_rtime   = max(rtime, na.rm = TRUE),
    m_start_obs = mass_w_parafilm[which.min(rtime)],  # g, with parafilm
    m_end_obs   = mass_w_parafilm[which.max(rtime)],  # g, with parafilm
    .groups = "drop"
  )
library(dplyr)
endwater_tbl_spring <- latempspring %>%
  group_by(segment) %>%
  summarise(
    end_rtime   = max(rtime, na.rm = TRUE),
    m_start_obs = mass_w_parafilm[which.min(rtime)],
    m_end_obs   = mass_w_parafilm[which.max(rtime)],
    .groups = "drop"
  ) %>%
  left_join(latempspring_dims %>% dplyr::select(segment, fresh_mass, dry_mass), by = "segment") %>%
  mutate(
    parafilm_mass_est = m_start_obs - fresh_mass,
    end_mass_sample   = m_end_obs - parafilm_mass_est,
    water_mass_end_g  = end_mass_sample - dry_mass,
    water_mass_init_g = fresh_mass - dry_mass,
    RWC_end           = water_mass_end_g / water_mass_init_g,
    LossFrac_end      = 1 - RWC_end
  ) %>%
  dplyr::select(segment, RWC_end, LossFrac_end)   # keep ONLY what you want

# join ONLY those columns into your final dataset
latempspring_gs <- latempspring_gs %>%
  left_join(endwater_tbl_spring, by = "segment")




write.csv(latempspring_gs,"gbark_temp_Goff_1957_spring_w_mf.csv")
####look at R2  
par(mfrow=c(1,1))
plot(latempspring_gs$R2,latempspring_gs$gs)
text(latempspring_gs$R2,latempspring_gs$gs,latempspring_gs$sampleID,cex=.5)
dev.off()  





##plot at dry downs with regression lines and gbark values
par(mfrow=c(5,5))
par(mar=c(2,2,1,1))
for(i in unique(latempspring$segment)){
  plot(latempspring$dry_time[latempspring$segment==i],latempspring$mass_w_parafilm[latempspring$segment==i])
  points(latempspring$dry_time[latempspring$segment==i&latempspring$exclude==0],latempspring$mass_w_parafilm[latempspring$segment==i&latempspring$exclude==0],pch=16)
  abline(latempspring_gs$int[latempspring_gs$segment==i],latempspring_gs$slope[latempspring_gs$segment==i])
  text(adj=c(0,0),x=0.1,y=I(min(latempspring$mass_w_parafilm[latempspring$segment==i])+.1*(max(latempspring$mass_w_parafilm[latempspring$segment==i])-min(latempspring$mass_w_parafilm[latempspring$segment==i]))),bquote(italic(g)[bark]==.(round(latempspring_gs$gs[latempspring_gs$segment==i],2))))
  mtext(i)
}


#dev.off()  
#how to use par?

#plot(latempspring$dry_time[latempspring$segment=="001b"],latempspring$mass_w_parafilm[latempspring$segment=="001b"])
#points(latempspring$dry_time[latempspring$segment=="001b"&latempspring$exclude==0],latempspring$mass_w_parafilm[latempspring$segment=="001b"&latempspring$exclude==0],pch=16)
#abline(latempspring_gs$int[latempspring_gs$segment=="001b"],latempspring_gs$slope[latempspring_gs$segment=="001b"])
#text(adj=c(0,0),x=0.1,y=I(min(latempspring$mass_w_parafilm[latempspring$segment=="001b"])+.1*(max(latempspring$mass_w_parafilm[latempspring$segment=="001b"])-min(latempspring$mass_w_parafilm[latempspring$segment=="001b"]))),bquote(italic(g)[bark]==.(round(latempspring_gs$gs[latempspring_gs$segment==i],2))))
#mtext("001b")

library(ggplot2)
plot<-ggplot(latempspring,aes(x=dry_time,y=mass_w_parafilm))+
  geom_point()+
  geom_smooth(method = "lm",formula = y~x)+
  facet_wrap(~segment,scales = "free")+
  theme_bw()
plot
