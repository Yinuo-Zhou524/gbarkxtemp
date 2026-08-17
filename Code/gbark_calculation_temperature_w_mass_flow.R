setwd("E://LSU//research//gbark_temperature")
library(MASS)
library(leaps)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
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
x<-latemp$mass_w_parafilm[latemp$segment=="019A"]
y<-as.numeric(latemp$dry_time[latemp$segment=="019A"])
plot(y,x)
for(i in 1:length(x))
{xi<-x[x<= c(sort(x)[length(x)+1-i])]
yi<-y[y>= c(sort(y,decreasing=T)[length(y)+1-i])]
lm<-lm(xi~yi)
print(summary(lm)$r.squared)
}
p=ggplot(latemp[latemp$segment=="010E",], aes(x=dry_time, y=mass_w_parafilm))+
  geom_point()+facet_wrap(~segment,scales="free")+
  xlab("dry time (hr)")+
  ylab("mass (g)")+
  theme_Publication(base_size = 14)+
  theme(
    axis.title = element_text(face = "plain"),
    panel.border = element_rect(color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(face = "plain"),aspect.ratio = 1
  )



p

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
#mass flow correction
latemp_tRH10C$ws_1<-(latemp_tRH10C$svp_1/101.3)
latemp_tRH10C$wa_1<-(latemp_tRH10C$RH_Avg/100)*latemp_tRH10C$ws_1
latemp_tRH10C$wmean_1<- (latemp_tRH10C$ws_1+latemp_tRH10C$wa_1)/2
latemp_tRH10C$massflow_factor_1<- 1 - latemp_tRH10C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2A<-data.frame(matrix(nrow=0,ncol=23,data=NA))
names(latemp_regs2A)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd","RH","RH_sd",
                        "D_early","D_late","mf_early","mf_late",
                        "vpd_early","vpd_late","RH_early","RH_late","Tc_early","Tc_late")

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
  latemp_regs2A[i,9]<-paste(mean(latemp_tRH10C$massflow_factor_1[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,10]<-paste(sd(latemp_tRH10C$massflow_factor_1[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,11]<-paste(mean(latemp_tRH10C$vpd_1[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,12]<-paste(mean(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
  latemp_regs2A[i,13]<-paste(sd(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=latemp_A$start_time[latemp_A$segment==i]&latemp_tRH10C$rtime<=latemp_A$end_time[latemp_A$segment==i]],na.rm=T))
}
#seperate early and late
for (i in unique(latemp_A$segment))
{
  t0  <- latemp_A$start_time[latemp_A$segment==i][1]
  t1  <- latemp_A$end_time[latemp_A$segment==i][1]
  tcut <- t0 + 0.6*(t1 - t0)   # last 40% is "late"
  
  # EARLY
  latemp_regs2A[i,14] <- paste(mean(latemp_tRH10C$D1[latemp_tRH10C$rtime>=t0 & latemp_tRH10C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2A[i,16] <- paste(mean(latemp_tRH10C$massflow_factor_1[latemp_tRH10C$rtime>=t0 & latemp_tRH10C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2A[i,18] <- paste(mean(latemp_tRH10C$vpd_1[latemp_tRH10C$rtime>=t0 & latemp_tRH10C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2A[i,20] <- paste(mean(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=t0 & latemp_tRH10C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2A[i,22] <- paste(mean(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=t0 & latemp_tRH10C$rtime<=tcut], na.rm=TRUE))
  
  # LATE
  latemp_regs2A[i,15] <- paste(mean(latemp_tRH10C$D1[latemp_tRH10C$rtime>=tcut & latemp_tRH10C$rtime<=t1], na.rm=TRUE))
  latemp_regs2A[i,17] <- paste(mean(latemp_tRH10C$massflow_factor_1[latemp_tRH10C$rtime>=tcut & latemp_tRH10C$rtime<=t1], na.rm=TRUE))
  latemp_regs2A[i,19] <- paste(mean(latemp_tRH10C$vpd_1[latemp_tRH10C$rtime>=tcut & latemp_tRH10C$rtime<=t1], na.rm=TRUE))
  latemp_regs2A[i,21] <- paste(mean(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=tcut & latemp_tRH10C$rtime<=t1], na.rm=TRUE))
  latemp_regs2A[i,23] <- paste(mean(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=tcut & latemp_tRH10C$rtime<=t1], na.rm=TRUE))
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
#mass flow correction
latemp_tRH25C$ws_1<-(latemp_tRH25C$svp_1/101.3)
latemp_tRH25C$wa_1<-(latemp_tRH25C$RH_Avg/100)*latemp_tRH25C$ws_1
latemp_tRH25C$wmean_1<- (latemp_tRH25C$ws_1+latemp_tRH25C$wa_1)/2
latemp_tRH25C$massflow_factor_1<- 1 - latemp_tRH25C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2B<-data.frame(matrix(nrow=0,ncol=23,data=NA))
names(latemp_regs2B)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd","RH","RH_sd",
                        "D_early","D_late","mf_early","mf_late",
                        "vpd_early","vpd_late","RH_early","RH_late","Tc_early","Tc_late")

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
  latemp_regs2B[i,9]<-paste(mean(latemp_tRH25C$massflow_factor_1[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,10]<-paste(sd(latemp_tRH25C$massflow_factor_1[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,11]<-paste(mean(latemp_tRH25C$vpd_1[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,12]<-paste(mean(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
  latemp_regs2B[i,13]<-paste(sd(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=latemp_B$start_time[latemp_B$segment==i]&latemp_tRH25C$rtime<=latemp_B$end_time[latemp_B$segment==i]],na.rm=T))
}
#seperate early and late
for (i in unique(latemp_B$segment))
{
  t0  <- latemp_B$start_time[latemp_B$segment==i][1]
  t1  <- latemp_B$end_time[latemp_B$segment==i][1]
  tcut <- t0 + 0.6*(t1 - t0)   # last 40% is "late"
  
  # EARLY
  latemp_regs2B[i,14] <- paste(mean(latemp_tRH25C$D1[latemp_tRH25C$rtime>=t0 & latemp_tRH25C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2B[i,16] <- paste(mean(latemp_tRH25C$massflow_factor_1[latemp_tRH25C$rtime>=t0 & latemp_tRH25C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2B[i,18] <- paste(mean(latemp_tRH25C$vpd_1[latemp_tRH25C$rtime>=t0 & latemp_tRH25C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2B[i,20] <- paste(mean(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=t0 & latemp_tRH25C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2B[i,22] <- paste(mean(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=t0 & latemp_tRH25C$rtime<=tcut], na.rm=TRUE))
  
  # LATE
  latemp_regs2B[i,15] <- paste(mean(latemp_tRH25C$D1[latemp_tRH25C$rtime>=tcut & latemp_tRH25C$rtime<=t1], na.rm=TRUE))
  latemp_regs2B[i,17] <- paste(mean(latemp_tRH25C$massflow_factor_1[latemp_tRH25C$rtime>=tcut & latemp_tRH25C$rtime<=t1], na.rm=TRUE))
  latemp_regs2B[i,19] <- paste(mean(latemp_tRH25C$vpd_1[latemp_tRH25C$rtime>=tcut & latemp_tRH25C$rtime<=t1], na.rm=TRUE))
  latemp_regs2B[i,21] <- paste(mean(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=tcut & latemp_tRH25C$rtime<=t1], na.rm=TRUE))
  latemp_regs2B[i,23] <- paste(mean(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=tcut & latemp_tRH25C$rtime<=t1], na.rm=TRUE))
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
#mass flow correction
latemp_tRH35C$ws_1<-(latemp_tRH35C$svp_1/101.3)
latemp_tRH35C$wa_1<-(latemp_tRH35C$RH_Avg/100)*latemp_tRH35C$ws_1
latemp_tRH35C$wmean_1<- (latemp_tRH35C$ws_1+latemp_tRH35C$wa_1)/2
latemp_tRH35C$massflow_factor_1<- 1 - latemp_tRH35C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2C<-data.frame(matrix(nrow=0,ncol=23,data=NA))
names(latemp_regs2C)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd","RH","RH_sd",
                        "D_early","D_late","mf_early","mf_late",
                        "vpd_early","vpd_late","RH_early","RH_late","Tc_early","Tc_late")

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
  latemp_regs2C[i,9]<-paste(mean(latemp_tRH35C$massflow_factor_1[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,10]<-paste(sd(latemp_tRH35C$massflow_factor_1[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,11]<-paste(mean(latemp_tRH35C$vpd_1[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,12]<-paste(mean(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
  latemp_regs2C[i,13]<-paste(sd(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=latemp_C$start_time[latemp_C$segment==i]&latemp_tRH35C$rtime<=latemp_C$end_time[latemp_C$segment==i]],na.rm=T))
}

#separate early and late
for (i in unique(latemp_C$segment))
{
  t0  <- latemp_C$start_time[latemp_C$segment==i][1]
  t1  <- latemp_C$end_time[latemp_C$segment==i][1]
  tcut <- t0 + 0.6*(t1 - t0)   # last 40% is "late"
  
  # EARLY
  latemp_regs2C[i,14] <- paste(mean(latemp_tRH35C$D1[latemp_tRH35C$rtime>=t0 & latemp_tRH35C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2C[i,16] <- paste(mean(latemp_tRH35C$massflow_factor_1[latemp_tRH35C$rtime>=t0 & latemp_tRH35C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2C[i,18] <- paste(mean(latemp_tRH35C$vpd_1[latemp_tRH35C$rtime>=t0 & latemp_tRH35C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2C[i,20] <- paste(mean(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=t0 & latemp_tRH35C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2C[i,22] <- paste(mean(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=t0 & latemp_tRH35C$rtime<=tcut], na.rm=TRUE))
  
  # LATE
  latemp_regs2C[i,15] <- paste(mean(latemp_tRH35C$D1[latemp_tRH35C$rtime>=tcut & latemp_tRH35C$rtime<=t1], na.rm=TRUE))
  latemp_regs2C[i,17] <- paste(mean(latemp_tRH35C$massflow_factor_1[latemp_tRH35C$rtime>=tcut & latemp_tRH35C$rtime<=t1], na.rm=TRUE))
  latemp_regs2C[i,19] <- paste(mean(latemp_tRH35C$vpd_1[latemp_tRH35C$rtime>=tcut & latemp_tRH35C$rtime<=t1], na.rm=TRUE))
  latemp_regs2C[i,21] <- paste(mean(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=tcut & latemp_tRH35C$rtime<=t1], na.rm=TRUE))
  latemp_regs2C[i,23] <- paste(mean(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=tcut & latemp_tRH35C$rtime<=t1], na.rm=TRUE))
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
#mass flow correction
latemp_tRH45C$ws_1<-(latemp_tRH45C$svp_1/101.3)
latemp_tRH45C$wa_1<-(latemp_tRH45C$RH_Avg/100)*latemp_tRH45C$ws_1
latemp_tRH45C$wmean_1<- (latemp_tRH45C$ws_1+latemp_tRH45C$wa_1)/2
latemp_tRH45C$massflow_factor_1<- 1 - latemp_tRH45C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2D<-data.frame(matrix(nrow=0,ncol=23,data=NA))
names(latemp_regs2D)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd","RH","RH_sd",
                        "D_early","D_late","mf_early","mf_late",
                        "vpd_early","vpd_late","RH_early","RH_late","Tc_early","Tc_late")

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
  latemp_regs2D[i,9]<-paste(mean(latemp_tRH45C$massflow_factor_1[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,10]<-paste(sd(latemp_tRH45C$massflow_factor_1[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,11]<-paste(mean(latemp_tRH45C$vpd_1[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,12]<-paste(mean(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
  latemp_regs2D[i,13]<-paste(sd(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=latemp_D$start_time[latemp_D$segment==i]&latemp_tRH45C$rtime<=latemp_D$end_time[latemp_D$segment==i]],na.rm=T))
}

#separate early and late
for (i in unique(latemp_D$segment))
{
  t0  <- latemp_D$start_time[latemp_D$segment==i][1]
  t1  <- latemp_D$end_time[latemp_D$segment==i][1]
  tcut <- t0 + 0.6*(t1 - t0)   # last 40% is "late"
  
  # EARLY
  latemp_regs2D[i,14] <- paste(mean(latemp_tRH45C$D1[latemp_tRH45C$rtime>=t0 & latemp_tRH45C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2D[i,16] <- paste(mean(latemp_tRH45C$massflow_factor_1[latemp_tRH45C$rtime>=t0 & latemp_tRH45C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2D[i,18] <- paste(mean(latemp_tRH45C$vpd_1[latemp_tRH45C$rtime>=t0 & latemp_tRH45C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2D[i,20] <- paste(mean(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=t0 & latemp_tRH45C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2D[i,22] <- paste(mean(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=t0 & latemp_tRH45C$rtime<=tcut], na.rm=TRUE))
  
  # LATE
  latemp_regs2D[i,15] <- paste(mean(latemp_tRH45C$D1[latemp_tRH45C$rtime>=tcut & latemp_tRH45C$rtime<=t1], na.rm=TRUE))
  latemp_regs2D[i,17] <- paste(mean(latemp_tRH45C$massflow_factor_1[latemp_tRH45C$rtime>=tcut & latemp_tRH45C$rtime<=t1], na.rm=TRUE))
  latemp_regs2D[i,19] <- paste(mean(latemp_tRH45C$vpd_1[latemp_tRH45C$rtime>=tcut & latemp_tRH45C$rtime<=t1], na.rm=TRUE))
  latemp_regs2D[i,21] <- paste(mean(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=tcut & latemp_tRH45C$rtime<=t1], na.rm=TRUE))
  latemp_regs2D[i,23] <- paste(mean(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=tcut & latemp_tRH45C$rtime<=t1], na.rm=TRUE))
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
#mass flow correction
latemp_tRH55C$ws_1<-(latemp_tRH55C$svp_1/101.3)
latemp_tRH55C$wa_1<-(latemp_tRH55C$RH_Avg/100)*latemp_tRH55C$ws_1
latemp_tRH55C$wmean_1<- (latemp_tRH55C$ws_1+latemp_tRH55C$wa_1)/2
latemp_tRH55C$massflow_factor_1<- 1 - latemp_tRH55C$wmean_1

#for each stem, calculate water loss rate with an lm and put coeffients in a table
latemp_regs2E<-data.frame(matrix(nrow=0,ncol=23,data=NA))
names(latemp_regs2E)<-c("segment","int","slope","R2","D","D_sd","Tc","Tc_sd","mf","mf_sd","vpd","RH","RH_sd",
                        "D_early","D_late","mf_early","mf_late",
                        "vpd_early","vpd_late","RH_early","RH_late","Tc_early","Tc_late")

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
  latemp_regs2E[i,9]<-paste(mean(latemp_tRH55C$massflow_factor_1[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,10]<-paste(sd(latemp_tRH55C$massflow_factor_1[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,11]<-paste(mean(latemp_tRH55C$vpd_1[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,12]<-paste(mean(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  latemp_regs2E[i,13]<-paste(sd(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=latemp_E$start_time[latemp_E$segment==i]&latemp_tRH55C$rtime<=latemp_E$end_time[latemp_E$segment==i]],na.rm=T))
  
}
#seperate early and late
for (i in unique(latemp_E$segment))
{
  t0  <- latemp_E$start_time[latemp_E$segment==i][1]
  t1  <- latemp_E$end_time[latemp_E$segment==i][1]
  tcut <- t0 + 0.6*(t1 - t0)   # last 40% is "late"
  
  # EARLY
  latemp_regs2E[i,14] <- paste(mean(latemp_tRH55C$D1[latemp_tRH55C$rtime>=t0 & latemp_tRH55C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2E[i,16] <- paste(mean(latemp_tRH55C$massflow_factor_1[latemp_tRH55C$rtime>=t0 & latemp_tRH55C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2E[i,18] <- paste(mean(latemp_tRH55C$vpd_1[latemp_tRH55C$rtime>=t0 & latemp_tRH55C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2E[i,20] <- paste(mean(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=t0 & latemp_tRH55C$rtime<=tcut], na.rm=TRUE))
  latemp_regs2E[i,22] <- paste(mean(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=t0 & latemp_tRH55C$rtime<=tcut], na.rm=TRUE))
  
  # LATE
  latemp_regs2E[i,15] <- paste(mean(latemp_tRH55C$D1[latemp_tRH55C$rtime>=tcut & latemp_tRH55C$rtime<=t1], na.rm=TRUE))
  latemp_regs2E[i,17] <- paste(mean(latemp_tRH55C$massflow_factor_1[latemp_tRH55C$rtime>=tcut & latemp_tRH55C$rtime<=t1], na.rm=TRUE))
  latemp_regs2E[i,19] <- paste(mean(latemp_tRH55C$vpd_1[latemp_tRH55C$rtime>=tcut & latemp_tRH55C$rtime<=t1], na.rm=TRUE))
  latemp_regs2E[i,21] <- paste(mean(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=tcut & latemp_tRH55C$rtime<=t1], na.rm=TRUE))
  latemp_regs2E[i,23] <- paste(mean(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=tcut & latemp_tRH55C$rtime<=t1], na.rm=TRUE))
}

#----

latemp_regs2<-rbind(latemp_regs2A,latemp_regs2B,latemp_regs2C,latemp_regs2D,latemp_regs2E)
#calculate water loss rate in mols per second
# divide by 18.01528 to convert from grams to moles
latemp_regs2$slope<-as.numeric(latemp_regs2$slope)
latemp_regs2$D<-as.numeric(latemp_regs2$D)
latemp_regs2$mols_per_sec<-(-latemp_regs2$slope/3600)/18.01528
latemp_regs2$mf<-as.numeric(latemp_regs2$mf)
latemp_regs2$RH<-as.numeric(latemp_regs2$RH)
latemp_regs2$D_early  <- as.numeric(latemp_regs2$D_early)
latemp_regs2$D_late   <- as.numeric(latemp_regs2$D_late)
latemp_regs2$mf_early <- as.numeric(latemp_regs2$mf_early)
latemp_regs2$mf_late  <- as.numeric(latemp_regs2$mf_late)
latemp_regs2$vpd_early<- as.numeric(latemp_regs2$vpd_early)
latemp_regs2$vpd_late <- as.numeric(latemp_regs2$vpd_late)
latemp_regs2$RH_early <- as.numeric(latemp_regs2$RH_early)
latemp_regs2$RH_late  <- as.numeric(latemp_regs2$RH_late)
latemp_regs2$Tc_early <- as.numeric(latemp_regs2$Tc_early)
latemp_regs2$Tc_late  <- as.numeric(latemp_regs2$Tc_late)

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
latemp_gs$gs<-latemp_gs$E*latemp_gs$mf/latemp_gs$D #mmol/s/m2

#calculate bark thickness ratio
latemp_gs$bark_thickness_ratio<-latemp_gs$bark_thickness/latemp_gs$diam

seg_obs <- latemp %>%
  group_by(segment) %>%
  summarise(
    start_rtime = min(rtime, na.rm = TRUE),
    end_rtime   = max(rtime, na.rm = TRUE),
    m_start_obs = mass_w_parafilm[which.min(rtime)],  # g, with parafilm
    m_end_obs   = mass_w_parafilm[which.max(rtime)],  # g, with parafilm
    .groups = "drop"
  )
library(dplyr)
endwater_tbl <- latemp %>%
  group_by(segment) %>%
  summarise(
    end_rtime   = max(rtime, na.rm = TRUE),
    m_start_obs = mass_w_parafilm[which.min(rtime)],
    m_end_obs   = mass_w_parafilm[which.max(rtime)],
    .groups = "drop"
  ) %>%
  left_join(latemp_dims %>% dplyr::select(segment, fresh_mass, dry_mass), by = "segment") %>%
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
latemp_gs <- latemp_gs %>%
  left_join(endwater_tbl, by = "segment")

write.csv(latemp_gs,"gbark_temp_Goff_1957_w_mf.csv")

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









library(dplyr)

# 1) Build a per-segment table of start/end dry_time and predicted masses
#    Use exclude==0 window (your kept data) so "end" matches what contributed to slope.
seg_end <- latemp %>%
  group_by(segment) %>%
  summarise(
    t_start = min(dry_time, na.rm = TRUE),
    t_end_included = max(dry_time, na.rm = TRUE),  # after you filtered exclude!=1 this is "included end"
    m_start_obs = mass_w_parafilm[which.min(dry_time)],  # observed start mass (w/ parafilm)
    m_end_obs   = mass_w_parafilm[which.max(dry_time)],  # observed end mass (w/ parafilm)
    n_points = n(),
    .groups = "drop"
  )

# Helper: fit LM per segment and return predictions + diagnostics
fit_pred <- function(df) {
  df <- df[is.finite(df$dry_time) & is.finite(df$mass_w_parafilm), ]
  if (nrow(df) < 3) return(data.frame(m_start_pred = NA, m_end_pred = NA, r2_full = NA,
                                      slope_early = NA, slope_late = NA, r2_late = NA))
  
  lm1 <- lm(mass_w_parafilm ~ dry_time, data = df)
  
  t0 <- min(df$dry_time, na.rm = TRUE)
  t1 <- max(df$dry_time, na.rm = TRUE)
  
  # predicted masses at start/end of included period
  m_start_pred <- predict(lm1, newdata = data.frame(dry_time = t0))
  m_end_pred   <- predict(lm1, newdata = data.frame(dry_time = t1))
  
  r2_full <- summary(lm1)$adj.r.squared
  
  # --- shrinkage-ish diagnostics: early vs late slopes and late-fit R2 ---
  t_cut <- t0 + 0.6*(t1 - t0)   # last 40% as "late" window (adjust if needed)
  df_early <- df[df$dry_time <= t_cut, ]
  df_late  <- df[df$dry_time >= t_cut, ]
  
  slope_early <- if (nrow(df_early) >= 3) coef(lm(mass_w_parafilm ~ dry_time, data = df_early))[2] else NA
  slope_late  <- if (nrow(df_late)  >= 3) coef(lm(mass_w_parafilm ~ dry_time, data = df_late))[2]  else NA
  r2_late     <- if (nrow(df_late)  >= 3) summary(lm(mass_w_parafilm ~ dry_time, data = df_late))$adj.r.squared else NA
  
  data.frame(m_start_pred, m_end_pred, r2_full, slope_early, slope_late, r2_late)
}

seg_pred <- latemp %>%
  group_by(segment) %>%
  group_modify(~fit_pred(.x)) %>%
  ungroup()

seg_stats <- seg_end %>%
  left_join(seg_pred, by = "segment")

# 2) Join segment stats into latemp_gs (one row per segment in latemp_gs)
latemp_gs <- latemp_gs %>%
  left_join(seg_stats, by = "segment")

# 3) Estimate parafilm mass and compute end water content
#    Assumption: fresh_mass measured pre-wrap and close to start of drydown.
latemp_gs <- latemp_gs %>%
  mutate(
    parafilm_mass_est = m_start_pred - fresh_mass,                # g
    m_end_sample_est  = m_end_pred - parafilm_mass_est,           # g, sample only
    water_mass_end_g  = m_end_sample_est - dry_mass,              # g water remaining at end
    WC_end_g_per_vol  = water_mass_end_g / vol,                   # your same units as water_content
    RWC_end           = water_mass_end_g / (fresh_mass - dry_mass),# relative water content at end
    slope_ratio_late_early = slope_late / slope_early,
    slope_delta_lateminusearly = slope_late - slope_early,
    gbark_early = ((((-slope_early/3600)/18.01528)*1000/sa)* mf/D_early),
    gbark_late = ((((-slope_late/3600)/18.01528)*1000/sa)* mf/D_late),
    gbark_delta_late_early = gbark_late - gbark_early,
    gbark_ratio_late_early = gbark_late - gbark_early
  )

# Optional safety checks
# Flag weird values (helps catch parafilm/fresh_mass mismatch)
latemp_gs <- latemp_gs %>%
  mutate(
    flag_bad_WC = ifelse(is.na(WC_end_g_per_vol) | WC_end_g_per_vol < 0 | RWC_end < 0 | RWC_end > 1.2, 1, 0)
  )

# Now latemp_gs has WC_end_g_per_vol and RWC_end you can model against shrinkage diagnostics.
latemp_gs$Tc<-as.numeric(latemp_gs$Tc)
slopediffnum<-lmer(slope_ratio_late_early~sp*Tc +(1|TreeID), latemp_gs)
summary(slopediffnum)
anova(slopediffnum)
#summary(latemp_gs$Tc)

library(lme4)
library(lmerTest)
slopediff<-lmer(slope_ratio_late_early~sp*treatment +(1|TreeID), latemp_gs)
summary(slopediff)
anova(slopediff)
library(emmeans)
library(MuMIn)
library(multcomp)
emm_slope<-emmeans(slopediff,~treatment)
pairs(emm_slope, adjust = "tukey")

cld(emm_slope, Letters=letters, adjust="tukey")


slopediff_add <- lmer(slope_ratio_late_early ~ sp + treatment + (1|TreeID), latemp_gs)
emm_slope <- emmeans(slopediff_add, ~ treatment)
pairs(emm_slope, adjust="tukey")
cld(emm_slope, Letters=letters, adjust="tukey")

#emm_slope_num<-emmeans(slopediffnum,~Tc)

#cld(emm_slope_num, Letters=letters, adjust="tukey")

