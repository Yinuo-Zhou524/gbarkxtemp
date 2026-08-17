setwd("E://LSU//research//gbark_temperature")
library(MASS)
library(leaps)
library(dplyr)
latemp1<-read.csv("temperature x gbark dry down.csv")
#colnames(latemp1)<-c("tag","date","sp",	"treatment"	,"segment_location",	"segment"	,"mass_w_parafilm",	"hour",	"minute")

#latemp1$tag<-as.numeric(latemp1$tag)
#latemp1$sp[is.na(latemp1$sp)] = "NA"
#latemp1$ob<-seq(1,length(latemp1$tag),1)


#make date and time readable
latemp1$time<-paste(latemp1$hour,latemp1$minute,sep=":")
latemp1$date_and_time1<-paste(latemp1$measurement_date,latemp1$time,sep=" ")
latemp1$rtime<-as.POSIXct(latemp1$date_and_time1,format="%Y.%m.%d %H:%M")

###calculate hours of dry down 
###What is end included time?
latemp_dry1<-data.frame(matrix(nrow=0,ncol=3,data=NA))
names(latemp_dry1)<-c("segment","start_time","end_time")


for(i in unique(latemp1$segment[!is.na(latemp1$mass_w_parafilm)]))
{latemp_dry1[i,1]<-paste(i)
latemp_dry1[i,2]<-paste(min(latemp1$rtime[latemp1$segment==i]))
#latemp_dry1[i,3]<-paste(max(latemp1$rtime[latemp1$exclude==0 & latemp1$tag==i]))
latemp_dry1[i,3]<-paste(max(latemp1$rtime[latemp1$segment==i]))}

latemp_dry1$start_time<-as.POSIXct(latemp_dry1$start_time,format="%Y-%m-%d %H:%M")
#latemp_dry1$end_included_time<-as.POSIXct(latemp_dry1$end_included_time,format="%Y-%m-%d %H:%M")
latemp_dry1$end_time<-as.POSIXct(latemp_dry1$end_time,format="%Y-%m-%d %H:%M")
latemp_dry1$total_time<-as.numeric(latemp_dry1$end_time-latemp_dry1$start_time)


#
latemp<-merge(latemp1,latemp_dry1,by="segment")
latemp$dry_time<-(as.numeric(latemp$rtime-latemp$start_time))/3600


#get the R^2 for each segment
x<-latemp$mass_w_parafilm[latemp$segment=="023E"]
y<-as.numeric(latemp$dry_time[latemp$segment=="023E"])

for(i in 1:length(x))
{xi<-x[x<= c(sort(x)[length(x)+1-i])]
yi<-y[y>= c(sort(y,decreasing=T)[length(y)+1-i])]
lm<-lm(xi~yi)
print(summary(lm)$r.squared)
}

##Make it for each segment, a function?
r2_exclude=function(a){
  x<-latemp$mass_w_parafilm[latemp$segment==a]
  y<-as.numeric(latemp$dry_time[latemp$segment==a])
  
  for(i in 1:length(x))
  {xi<-x[x<= c(sort(x)[length(x)+1-i])]
  yi<-y[y>= c(sort(y,decreasing=T)[length(y)+1-i])]
  lm<-lm(xi~yi)
  print(summary(lm)$r.squared)
  }
}

#R2>0.99
  
###Input the excluded data
latemp2<-read.csv("temperature x gbark dry down_exclude.csv")
colnames(latemp2)<-c("tag","collection_date","round","date","sp",	"treatment"	,"treatment_ID",	"segment"	,"mass_w_parafilm",	"exclude","hour",	"minute")


#make date and time readable
latemp2$time<-paste(latemp2$hour,latemp2$minute,sep=":")
latemp2$date_and_time1<-paste(latemp2$date,latemp2$time,sep=" ")
latemp2$rtime<-as.POSIXct(latemp2$date_and_time1,format="%Y.%m.%d %H:%M")

###calculate hours of dry down 
###What is end included time?
latemp_dry2<-data.frame(matrix(nrow=0,ncol=4,data=NA))
names(latemp_dry2)<-c("segment","start_time","end_included_time","end_time")


for(i in unique(latemp2$segment[!is.na(latemp2$mass_w_parafilm)]))
{latemp_dry2[i,1]<-paste(i)
latemp_dry2[i,2]<-paste(min(latemp2$rtime[latemp2$segment==i]))
latemp_dry2[i,3]<-paste(max(latemp2$rtime[latemp2$exclude==0 & latemp2$segment==i]))
latemp_dry2[i,4]<-paste(max(latemp2$rtime[latemp2$segment==i]))}

latemp_dry2$start_time<-as.POSIXct(latemp_dry2$start_time,format="%Y-%m-%d %H:%M")
latemp_dry2$end_included_time<-as.POSIXct(latemp_dry2$end_included_time,format="%Y-%m-%d %H:%M")
latemp_dry2$end_time<-as.POSIXct(latemp_dry2$end_time,format="%Y-%m-%d %H:%M")
latemp_dry2$total_time<-as.numeric(latemp_dry2$end_time-latemp_dry2$start_time)


#
latemp<-merge(latemp2,latemp_dry2,by="segment")
latemp$dry_time<-(as.numeric(latemp$rtime-latemp$start_time))/3600

#exlude outlier points within stems that are probably due to wet stems or lack of water in the segment
latemp<-latemp[latemp$exclude!=1,]


 
#input air temp and RH
latemp_tRH10C<-read.csv("temp x gbark tRH 10C.csv")
latemp_tRH25C<-read.csv("temp x gbark tRH 25C.csv")
latemp_tRH35C<-read.csv("temp x gbark tRH 35C.csv")
latemp_tRH45C<-read.csv("temp x gbark tRH 45C.csv")
latemp_tRH55C<-read.csv("temp x gbark tRH 55C.csv")

#separate segment by temperature id 
latemp_A<-latemp[latemp$treatment_ID=="A",]
latemp_B<-latemp[latemp$treatment_ID=="B",]
latemp_C<-latemp[latemp$treatment_ID=="C",]
latemp_D<-latemp[latemp$treatment_ID=="D",]
latemp_E<-latemp[latemp$treatment_ID=="E",]

#Sys.setenv(TZ="CDT")
#Sys.getenv("TZ")
#options(tz="CDT")

###10C---A----
latemp_tRH10C$rtime<-as.POSIXct(latemp_tRH10C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference

#Log10 ew =  10.79574 (1-273.16/T)                                                          [3]
#- 5.02800 Log10(T/273.16)
#+ 1.50475 10-4 (1 - 10^(-8.2969*(T/273.16-1)))
#+ 0.42873 10-3 (10^(+4.76955*(1-273.16/T)) - 1)
#+ 0.78614
#with T in [K] and ew in [hPa] (Goff, 1957)

#latemp_tRH10C$svp_1<-0.6135*exp(17.502*latemp_tRH10C$Temp_Avg/(240.97+latemp_tRH10C$Temp_Avg)) #from LI610 manual
latemp_tRH10C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latemp_tRH10C$Temp_Avg+273.15)) - 5.02800*log10((latemp_tRH10C$Temp_Avg+273.15)/273.16)+
                           1.50475*10^(-4)*(1-10^(-8.2969*((latemp_tRH10C$Temp_Avg+273.15)/273.16-1)))+
                           0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_tRH10C$Temp_Avg+273.15))))+
                           0.78614)) # Golf 1957 unit:kPa
#latemp_tRH10C$svp_1<-0.1*(10^(-7.90298*(373.16/(latemp_tRH10C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latemp_tRH10C$Temp_Avg+273.15))-
#                           1.3816*10^(-7)*(10^(11.344*(1-(latemp_tRH10C$Temp_Avg+273.15)/373.16))-1)+
#                           8.1328*10^(-3)*(10^(-3.49149*(373.16/(latemp_tRH10C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latemp_tRH10C$vpd_1<-((100-latemp_tRH10C$RH_Avg)/100)*latemp_tRH10C$svp_1
latemp_tRH10C$D1<-(latemp_tRH10C$vpd_1/101.3)


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2A<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2A)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd")

for (i in unique(latemp_A$segment))
{
  lm1<-lm(latemp_A$mass_w_parafilm[latemp_A$segment==i]~as.numeric(latemp_A$dry_time[latemp_A$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2A[i,1]<-paste(i)
  latemp_regs2A[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2A[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2A[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latemp_tRH10C$rtime<-as.numeric(latemp_tRH10C$rtime)
for (i in unique(latemp_A$segment))
{
  latemp_regs2A[i,5]<-paste(mean(latemp_tRH10C$D1[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,6]<-paste(sd(latemp_tRH10C$D1[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,7]<-paste(mean(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,8]<-paste(sd(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
}


###25C---B----
latemp_tRH25C$rtime<-as.POSIXct(latemp_tRH25C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latemp_tRH25C$svp_1<-0.6135*exp(17.502*latemp_tRH25C$Temp_Avg/(240.97+latemp_tRH25C$Temp_Avg)) #from LI610 manual
latemp_tRH25C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latemp_tRH25C$Temp_Avg+273.15)) - 5.02800*log10((latemp_tRH25C$Temp_Avg+273.15)/273.16)+
                           1.50475*10^(-4)*(1-10^(-8.2969*((latemp_tRH25C$Temp_Avg+273.15)/273.16-1)))+
                           0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_tRH25C$Temp_Avg+273.15))))+
                           0.78614))#Goff 1957
#latemp_tRH25C$svp_1<-0.1*(10^(-7.90298*(373.16/(latemp_tRH25C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latemp_tRH25C$Temp_Avg+273.15))-
 #                          1.3816*10^(-7)*(10^(11.344*(1-(latemp_tRH25C$Temp_Avg+273.15)/373.16))-1)+
  #                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latemp_tRH25C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latemp_tRH25C$vpd_1<-((100-latemp_tRH25C$RH_Avg)/100)*latemp_tRH25C$svp_1
latemp_tRH25C$D1<-(latemp_tRH25C$vpd_1/101.3)


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2B<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2B)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd")

for (i in unique(latemp_B$segment))
{
  lm1<-lm(latemp_B$mass_w_parafilm[latemp_B$segment==i]~as.numeric(latemp_B$dry_time[latemp_B$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2B[i,1]<-paste(i)
  latemp_regs2B[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2B[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2B[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latemp_tRH25C$rtime<-as.numeric(latemp_tRH25C$rtime)
for (i in unique(latemp_B$segment))
{
  latemp_regs2B[i,5]<-paste(mean(latemp_tRH25C$D1[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,6]<-paste(sd(latemp_tRH25C$D1[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,7]<-paste(mean(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,8]<-paste(sd(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
}

###35C---C----
latemp_tRH35C$rtime<-as.POSIXct(latemp_tRH35C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latemp_tRH35C$svp_1<-0.6135*exp(17.502*latemp_tRH35C$Temp_Avg/(240.97+latemp_tRH35C$Temp_Avg)) #from LI610 manual
latemp_tRH35C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latemp_tRH35C$Temp_Avg+273.15)) - 5.02800*log10((latemp_tRH35C$Temp_Avg+273.15)/273.16)+
                           1.50475*10^(-4)*(1-10^(-8.2969*((latemp_tRH35C$Temp_Avg+273.15)/273.16-1)))+
                           0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_tRH35C$Temp_Avg+273.15))))+
                           0.78614))#Golf 1957
#latemp_tRH35C$svp_1<-0.1*(10^(-7.90298*(373.16/(latemp_tRH35C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latemp_tRH35C$Temp_Avg+273.15))-
 #                          1.3816*10^(-7)*(10^(11.344*(1-(latemp_tRH35C$Temp_Avg+273.15)/373.16))-1)+
  #                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latemp_tRH35C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latemp_tRH35C$vpd_1<-((100-latemp_tRH35C$RH_Avg)/100)*latemp_tRH35C$svp_1
latemp_tRH35C$D1<-(latemp_tRH35C$vpd_1/101.3)


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2C<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd")

for (i in unique(latemp_C$segment))
{
  lm1<-lm(latemp_C$mass_w_parafilm[latemp_C$segment==i]~as.numeric(latemp_C$dry_time[latemp_C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C[i,1]<-paste(i)
  latemp_regs2C[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latemp_tRH35C$rtime<-as.numeric(latemp_tRH35C$rtime)
for (i in unique(latemp_C$segment))
{
  latemp_regs2C[i,5]<-paste(mean(latemp_tRH35C$D1[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,6]<-paste(sd(latemp_tRH35C$D1[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,7]<-paste(mean(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,8]<-paste(sd(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
}

###45C---D----
latemp_tRH45C$rtime<-as.POSIXct(latemp_tRH45C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latemp_tRH45C$svp_1<-0.6135*exp(17.502*latemp_tRH45C$Temp_Avg/(240.97+latemp_tRH45C$Temp_Avg)) #from LI610 manual
latemp_tRH45C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latemp_tRH45C$Temp_Avg+273.15)) - 5.02800*log10((latemp_tRH45C$Temp_Avg+273.15)/273.16)+
                           1.50475*10^(-4)*(1-10^(-8.2969*((latemp_tRH45C$Temp_Avg+273.15)/273.16-1)))+
                           0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_tRH45C$Temp_Avg+273.15))))+
                           0.78614))
#latemp_tRH45C$svp_1<-0.1*(10^(-7.90298*(373.16/(latemp_tRH45C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latemp_tRH45C$Temp_Avg+273.15))-
 #                          1.3816*10^(-7)*(10^(11.344*(1-(latemp_tRH45C$Temp_Avg+273.15)/373.16))-1)+
  #                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latemp_tRH45C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latemp_tRH45C$vpd_1<-((100-latemp_tRH45C$RH_Avg)/100)*latemp_tRH45C$svp_1
latemp_tRH45C$D1<-(latemp_tRH45C$vpd_1/101.3)


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2D<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2D)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd")

for (i in unique(latemp_D$segment))
{
  lm1<-lm(latemp_D$mass_w_parafilm[latemp_D$segment==i]~as.numeric(latemp_D$dry_time[latemp_D$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2D[i,1]<-paste(i)
  latemp_regs2D[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2D[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2D[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latemp_tRH45C$rtime<-as.numeric(latemp_tRH45C$rtime)
for (i in unique(latemp_D$segment))
{
  latemp_regs2D[i,5]<-paste(mean(latemp_tRH45C$D1[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,6]<-paste(sd(latemp_tRH45C$D1[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,7]<-paste(mean(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,8]<-paste(sd(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
}


###55C---E----
latemp_tRH55C$rtime<-as.POSIXct(latemp_tRH55C$TIMESTAMP,format="%Y/%m/%d %H:%M")
#calculate D as mole fraction difference
#latemp_tRH55C$svp_1<-0.6135*exp(17.502*latemp_tRH55C$Temp_Avg/(240.97+latemp_tRH55C$Temp_Avg)) #from LI610 manual
latemp_tRH55C$svp_1<-0.1*(10^(10.79574*(1-273.16/(latemp_tRH55C$Temp_Avg+273.15)) - 5.02800*log10((latemp_tRH55C$Temp_Avg+273.15)/273.16)+
                           1.50475*10^(-4)*(1-10^(-8.2969*((latemp_tRH55C$Temp_Avg+273.15)/273.16-1)))+
                           0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_tRH55C$Temp_Avg+273.15))))+
                           0.78614))
#latemp_tRH55C$svp_1<-0.1*(10^(-7.90298*(373.16/(latemp_tRH55C$Temp_Avg+273.15)-1)+5.02808*log10(373.16/(latemp_tRH55C$Temp_Avg+273.15))-
 #                          1.3816*10^(-7)*(10^(11.344*(1-(latemp_tRH55C$Temp_Avg+273.15)/373.16))-1)+
  #                         8.1328*10^(-3)*(10^(-3.49149*(373.16/(latemp_tRH55C$Temp_Avg+273.15)-1))-1)+log10(1013.246))) #Goff Gratch 1946
latemp_tRH55C$vpd_1<-((100-latemp_tRH55C$RH_Avg)/100)*latemp_tRH55C$svp_1
latemp_tRH55C$D1<-(latemp_tRH55C$vpd_1/101.3)


#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2E<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2E)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd")

for (i in unique(latemp_E$segment))
{
  lm1<-lm(latemp_E$mass_w_parafilm[latemp_E$segment==i]~as.numeric(latemp_E$dry_time[latemp_E$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2E[i,1]<-paste(i)
  latemp_regs2E[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2E[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2E[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean D over the course of each stem dry down
latemp_tRH55C$rtime<-as.numeric(latemp_tRH55C$rtime)
for (i in unique(latemp_E$segment))
{
  latemp_regs2E[i,5]<-paste(mean(latemp_tRH55C$D1[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,6]<-paste(sd(latemp_tRH55C$D1[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,7]<-paste(mean(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,8]<-paste(sd(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
}


#----

latemp_regs2<-rbind(latemp_regs2A,latemp_regs2B,latemp_regs2C,latemp_regs2D,latemp_regs2E)
#calculate water loss rate in mols per second
# divide by 18.01528 to convert from grams to moles
latemp_regs2$slope<-as.numeric(latemp_regs2$slope)
latemp_regs2$D<-as.numeric(latemp_regs2$D)
latemp_regs2$mols_per_sec<-(-latemp_regs2$slope/3600)/18.01528


##merge lm coefs /with stem dimensions
latemp_dims<-read.csv("temperature x gbark traits.csv")

df<-unique(latemp[,c("segment","treatment")])
latemp_dims<-merge(latemp_dims,df,all.y = FALSE)
latemp_dims$sampleID<-latemp_dims$segment

##calculate surface area in m2 as cylinder
latemp_dims$diam<-(latemp_dims$basal_diam1+latemp_dims$basal_diam2+latemp_dims$distal_diam1+latemp_dims$distal_diam2)/4
latemp_dims$wood_diam<-(latemp_dims$basal_xylem_diam1+latemp_dims$basal_xylem_diam2+latemp_dims$distal_xylem_diam1+latemp_dims$distal_xylem_diam2)/4
#latemp_dims$exp_length_try1<-(latemp_dims$exp_length1+latemp_dims$exp_length2)/2
#latemp_dims$exp_length_try2<-ifelse(is.na(latemp_dims$exp_length_try1),latemp_dims$exp_length1,latemp_dims$exp_length_try1)
#latemp_dims$exp_length<-ifelse(is.na(latemp_dims$exp_length_try2),latemp_dims$exp_length2,latemp_dims$exp_length_try2)
latemp_dims$sa<-pi*((latemp_dims$diam)/1000)*(latemp_dims$exposed_length/1000)

##calculate bark thickness
latemp_dims$bark_thickness<-(latemp_dims$diam-latemp_dims$wood_diam)/2

#calculate volume and initial water content, and stem density
latemp_dims$vol<-latemp_dims$total_length*(latemp_dims$diam/2)^2
latemp_dims$water_content<-(latemp_dims$fresh_mass-latemp_dims$dry_mass)/latemp_dims$vol
latemp_dims$density<-latemp_dims$dry_mass/latemp_dims$vol


#calculate lenticel density and size
latemp_dims$lenticel_size<-rowMeans(latemp_dims[,c(10:19)],na.rm = T)
latemp_dims$LD<-latemp_dims$n_lenticels/(latemp_dims$diam*latemp_dims$length_counted_lenticel)
##
#names(latemp_regs2)
#names(latemp_dims)
latemp_gs<-merge(latemp_regs2,latemp_dims,all.x=T)


#calculate E and gs over entire drying time
latemp_gs$E<-(latemp_gs$mols_per_sec*1000)/latemp_gs$sa #mmol/s/m2
latemp_gs$gs<-latemp_gs$E/latemp_gs$D #mmol/s/m2

#calculate bark thickness ratio
latemp_gs$bark_thickness_ratio<-latemp_gs$bark_thickness/latemp_gs$diam

####look at R2  
par(mfrow=c(1,1))
plot(latemp_gs$R2,latemp_gs$gs,col=NA)
text(latemp_gs$R2,latemp_gs$gs,latemp_gs$sampleID,cex=.5)
dev.off()  





##plot at dry downs with regression lines and gbark values
par(mfrow=c(5,5))
par(mar=c(2,2,1,1))
for(i in unique(latemp$segment)){
  plot(latemp$dry_time[latemp$segment==i],latemp$mass_w_parafilm[latemp$segment==i])
  points(latemp$dry_time[latemp$segment==i&latemp$exclude==0],latemp$mass_w_parafilm[latemp$segment==i&latemp$exclude==0],pch=16)
  abline(latemp_gs$int[latemp_gs$segment==i],latemp_gs$slope[latemp_gs$segment==i])
  text(adj=c(0,0),x=0.1,y=I(min(latemp$mass_w_parafilm[latemp$segment==i])+.1*(max(latemp$mass_w_parafilm[latemp$segment==i])-min(latemp$mass_w_parafilm[latemp$segment==i]))),bquote(italic(g)[bark]==.(round(latemp_gs$gs[latemp_gs$segment==i],2))))
  mtext(i)
}
dev.off()  
#how to use par?

#plot(latemp$dry_time[latemp$segment=="001b"],latemp$mass_w_parafilm[latemp$segment=="001b"])
#points(latemp$dry_time[latemp$segment=="001b"&latemp$exclude==0],latemp$mass_w_parafilm[latemp$segment=="001b"&latemp$exclude==0],pch=16)
#abline(latemp_gs$int[latemp_gs$segment=="001b"],latemp_gs$slope[latemp_gs$segment=="001b"])
#text(adj=c(0,0),x=0.1,y=I(min(latemp$mass_w_parafilm[latemp$segment=="001b"])+.1*(max(latemp$mass_w_parafilm[latemp$segment=="001b"])-min(latemp$mass_w_parafilm[latemp$segment=="001b"]))),bquote(italic(g)[bark]==.(round(latemp_gs$gs[latemp_gs$segment==i],2))))
#mtext("001b")

library(ggplot2)
plot<-ggplot(latemp,aes(x=dry_time,y=mass_w_parafilm))+
  geom_point()+
  geom_smooth(method = "lm",formula = y~x)+
  facet_wrap(~segment,scales = "free")+
  theme_bw()
plot

#Linear model
lm1<-lm(gs~sp*treatment*segment_location*LD,latemp_gs)
summary(lm1)
anova(lm1)

lm1ls<-lm(gs~treatment*segment_location*LD,ls)
anova(lm1ls)

lm1na<-lm(gs~treatment*segment_location*LD,na)
anova(lm1na)


#lenticel size
lmls<-lm(lenticel_size~treatment*segment_location,ls)
anova(lmls)
lmna<-lm(lenticel_size~treatment*segment_location,na)
anova(lmna)

#lenticel density
lm2<-lm(LD~sp*treatment*segment_location,latemp_gs)
anova(lm2)
lm2ls<-lm(LD~treatment*segment_location,ls)
anova(lm2ls)
lm2na<-lm(LD~treatment*segment_location,na)
anova(lm2na)

#calculate relative gbark
latemp_relgs<-data.frame(matrix(nrow=0,ncol=4,data=NA))
names(latemp_relgs)<-c("PlantID","sp","treatment","relgs")


for(i in unique(latemp_gs$plantID))
{latemp_relgs[i,1]<-paste(i)
latemp_relgs[i,2]<-paste(unique(latemp_gs$sp[latemp_gs$plantID==i]))
latemp_relgs[i,3]<-paste(unique(latemp_gs$treatment[latemp_gs$plantID==i]))
latemp_relgs[i,4]<-paste(latemp_gs$gs[latemp_gs$plantID==i & latemp_gs$segment_location=="base"]/latemp_gs$gs[latemp_gs$plantID==i & latemp_gs$segment_location=="top"])
}

latemp_relgs$relgs<-as.numeric(latemp_relgs$relgs)
plot(latemp_relgs$PlantID,latemp_relgs$relgs)
plot(latemp_relgs$PlantID[latemp_relgs$sp=="NA"],latemp_relgs$relgs[latemp_relgs$sp=="NA"])

latemp_relgs$relgs<-as.numeric(latemp_relgs$relgs)
prelgs1<-ggplot(latemp_relgs,aes(sp,relgs))+
  geom_boxplot()+
  theme_classic()
prelgs1
prelgs2<-ggplot(latemp_relgs,aes(treatment,relgs,fill=sp))+
  geom_boxplot()+
  theme_classic()
prelgs2

latemp_gspos<-data.frame(matrix(nrow=0,ncol=5,data=NA))
names(latemp_gspos)<-c("PlantID","sp","treatment","topgs","basegs")


for(i in unique(latemp_gs$plantID))
{latemp_gspos[i,1]<-paste(i)
latemp_gspos[i,2]<-paste(unique(latemp_gs$sp[latemp_gs$plantID==i]))
latemp_gspos[i,3]<-paste(unique(latemp_gs$treatment[latemp_gs$plantID==i]))
latemp_gspos[i,4]<-paste(latemp_gs$gs[latemp_gs$plantID==i & latemp_gs$segment_location=="top"])
latemp_gspos[i,5]<-paste(latemp_gs$gs[latemp_gs$plantID==i & latemp_gs$segment_location=="base"])
}
latemp_gspos$topgs<-as.numeric(latemp_gspos$topgs)
latemp_gspos$basegs<-as.numeric(latemp_gspos$basegs)

t.test(latemp_gspos[latemp_gspos$sp=="LS",]$topgs,latemp_gspos[latemp_gspos$sp=="LS",]$basegs,paired = T)

t.test(latemp_gspos[latemp_gspos$sp=="NA",]$topgs,latemp_gspos[latemp_gspos$sp=="NA",]$basegs,paired = T)


t.test(latemp_gspos$topgs,latemp_gspos$basegs,paired = T)

var.test(relgs~sp,latemp_relgs)
t.test(relgs~sp,latemp_relgs)

p<-ggplot(latemp_gspos[latemp_gspos$sp=="LS",],aes(basegs,topgs,color=treatment))+
  geom_point()+
  geom_abline()+
  scale_color_grey()+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  
p
p<-ggplot(latemp_gspos[!latemp_gspos$PlantID=="145",],aes(basegs,topgs,color=sp))+
  geom_point()+
  #geom_abline()+
  scale_color_grey()+
  xlab("gbark_base")+
  ylab("gbark_top")+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p

p<-ggplot(latemp_gs[latemp_gs$segment_location=="base",],aes(LD,gs,color=sp))+
  geom_point()+
  #geom_abline()+
  theme_classic()

p

p<-ggplot(latemp_gs[latemp_gs$segment_location=="top",],aes(LD,gs,color=sp))+
  geom_point()+
  #geom_abline()+
  theme_classic()

p

p<-ggplot(latemp_gs[latemp_gs$segment_location=="top",],aes(lenticel_size,gs,color=sp))+
  geom_point()+
  #geom_abline()+
  theme_classic()

p
p<-ggplot(latemp_gs[latemp_gs$segment_location=="base",],aes(lenticel_size,gs,color=sp))+
  geom_point()+
  #geom_abline()+
  theme_classic()

p

p<-ggplot(latemp_gs[latemp_gs$sp=="NA",],aes(lenticel_size,gs))+
  geom_point()+
  geom_abline()+
  theme_classic()

p