setwd("F://LSU//research//gbark_temperature")
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
latemp_10C<-latemp[latemp$treatment_ID=="A",]
latemp_25C<-latemp[latemp$treatment_ID=="B",]
latemp_35C<-latemp[latemp$treatment_ID=="C",]
latemp_45C<-latemp[latemp$treatment_ID=="D",]
latemp_55C<-latemp[latemp$treatment_ID=="E",]


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


#calculate lenticel density and size, bark thickness ratio
latemp_dims$lenticel_size<-rowMeans(latemp_dims[,c(10:19)],na.rm = T)
latemp_dims$LD<-latemp_dims$n_lenticels/(latemp_dims$diam*latemp_dims$length_counted_lenticel)
latemp_dims$bark_thickness_ratio<-latemp_dims$bark_thickness/latemp_dims$diam

###10C---A----
latemp_tRH10C$rtime<-as.POSIXct(latemp_tRH10C$TIMESTAMP,format="%Y/%m/%d %H:%M")
latemp_regs2C10<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C10)<-c("segment","int","slope","R2","Tac","Tac_sd","RH","RH_sd")

for (i in unique(latemp_10C$segment))
{
  lm1<-lm(latemp_10C$mass_w_parafilm[latemp_10C$segment==i]~as.numeric(latemp_10C$dry_time[latemp_10C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C10[i,1]<-paste(i)
  latemp_regs2C10[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C10[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C10[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean T over the course of each stem dry down
latemp_tRH10C$rtime<-as.numeric(latemp_tRH10C$rtime)
for (i in unique(latemp_10C$segment))
{
  latemp_regs2C10[i,5]<-paste(mean(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=latemp_10C$start_time[latemp_10C$segment==i]&latemp_tRH10C$rtime<=latemp_10C$end_time[latemp_10C$segment==i]],na.rm=T))
  latemp_regs2C10[i,6]<-paste(sd(latemp_tRH10C$Temp_Avg[latemp_tRH10C$rtime>=latemp_10C$start_time[latemp_10C$segment==i]&latemp_tRH10C$rtime<=latemp_10C$end_time[latemp_10C$segment==i]],na.rm=T))
  latemp_regs2C10[i,7]<-paste(mean(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=latemp_10C$start_time[latemp_10C$segment==i]&latemp_tRH10C$rtime<=latemp_10C$end_time[latemp_10C$segment==i]],na.rm=T))
  latemp_regs2C10[i,8]<-paste(sd(latemp_tRH10C$RH_Avg[latemp_tRH10C$rtime>=latemp_10C$start_time[latemp_10C$segment==i]&latemp_tRH10C$rtime<=latemp_10C$end_time[latemp_10C$segment==i]],na.rm=T))
}
latemp_regs2C10$slope<-as.numeric(latemp_regs2C10$slope)
latemp_regs2C10$R2<-as.numeric(latemp_regs2C10$R2)
latemp_regs2C10$Tac<-as.numeric(latemp_regs2C10$Tac)
latemp_regs2C10$Tac_sd<-as.numeric(latemp_regs2C10$Tac_sd)
latemp_regs2C10$RH<-as.numeric(latemp_regs2C10$RH)
latemp_regs2C10$RH_sd<-as.numeric(latemp_regs2C10$RH_sd)


latemp_regs2C10$mols_per_sec<-(-latemp_regs2C10$slope/3600)/18.01528


latemp_gsC10<-merge(latemp_regs2C10,latemp_dims[latemp_dims$treatment=="10C",],all.x=T)
latemp_gsC10$E<-(latemp_gsC10$mols_per_sec*1000)/latemp_gsC10$sa #mmol/s/m2

Tcavg10C<-mean(as.numeric(latemp_regs2C10$Tac))
#L=44.6 KJ/mol
latemp_gsC10$deltaT<-44.6*latemp_gsC10$E/(29.3*0.135*sqrt(2/(latemp_gsC10$diam*10^(-3))))
latemp_gsC10$Tsc<-latemp_gsC10$Tac - latemp_gsC10$deltaT
latemp_gsC10$wi<-0.1*(10^(10.79574*(1-273.16/(latemp_gsC10$Tsc+273.15)) - 5.02800*log10((latemp_gsC10$Tsc+273.15)/273.16)+
                            1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC10$Tsc+273.15)/273.16-1)))+
                            0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC10$Tsc+273.15))))+
                            0.78614))/101.3
latemp_gsC10$w0<-latemp_gsC10$RH/100* 0.1*(10^(10.79574*(1-273.16/(latemp_gsC10$Tac+273.15)) - 5.02800*log10((latemp_gsC10$Tac+273.15)/273.16)+
                                                 1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC10$Tac+273.15)/273.16-1)))+
                                                 0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC10$Tac+273.15))))+
                                                 0.78614))/101.3
latemp_gsC10$deltaw<-latemp_gsC10$wi-latemp_gsC10$w0
latemp_gsC10$gs<-latemp_gsC10$E/latemp_gsC10$deltaw #mmol/s/m2

###25C---B----
latemp_tRH25C$rtime<-as.POSIXct(latemp_tRH25C$TIMESTAMP,format="%Y/%m/%d %H:%M")
latemp_regs2C25<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C25)<-c("segment","int","slope","R2","Tac","Tac_sd","RH","RH_sd")

for (i in unique(latemp_25C$segment))
{
  lm1<-lm(latemp_25C$mass_w_parafilm[latemp_25C$segment==i]~as.numeric(latemp_25C$dry_time[latemp_25C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C25[i,1]<-paste(i)
  latemp_regs2C25[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C25[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C25[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean T over the course of each stem dry down
latemp_tRH25C$rtime<-as.numeric(latemp_tRH25C$rtime)
for (i in unique(latemp_25C$segment))
{
  latemp_regs2C25[i,5]<-paste(mean(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=latemp_25C$start_time[latemp_25C$segment==i]&latemp_tRH25C$rtime<=latemp_25C$end_time[latemp_25C$segment==i]],na.rm=T))
  latemp_regs2C25[i,6]<-paste(sd(latemp_tRH25C$Temp_Avg[latemp_tRH25C$rtime>=latemp_25C$start_time[latemp_25C$segment==i]&latemp_tRH25C$rtime<=latemp_25C$end_time[latemp_25C$segment==i]],na.rm=T))
  latemp_regs2C25[i,7]<-paste(mean(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=latemp_25C$start_time[latemp_25C$segment==i]&latemp_tRH25C$rtime<=latemp_25C$end_time[latemp_25C$segment==i]],na.rm=T))
  latemp_regs2C25[i,8]<-paste(sd(latemp_tRH25C$RH_Avg[latemp_tRH25C$rtime>=latemp_25C$start_time[latemp_25C$segment==i]&latemp_tRH25C$rtime<=latemp_25C$end_time[latemp_25C$segment==i]],na.rm=T))
}
latemp_regs2C25$slope<-as.numeric(latemp_regs2C25$slope)
latemp_regs2C25$R2<-as.numeric(latemp_regs2C25$R2)
latemp_regs2C25$Tac<-as.numeric(latemp_regs2C25$Tac)
latemp_regs2C25$Tac_sd<-as.numeric(latemp_regs2C25$Tac_sd)
latemp_regs2C25$RH<-as.numeric(latemp_regs2C25$RH)
latemp_regs2C25$RH_sd<-as.numeric(latemp_regs2C25$RH_sd)


latemp_regs2C25$mols_per_sec<-(-latemp_regs2C25$slope/3600)/18.01528


latemp_gsC25<-merge(latemp_regs2C25,latemp_dims[latemp_dims$treatment=="25C",],all.x=T)
latemp_gsC25$E<-(latemp_gsC25$mols_per_sec*1000)/latemp_gsC25$sa #mmol/s/m2

Tcavg25C<-mean(as.numeric(latemp_regs2C25$Tac))
#L=44.03 KJ/mol
latemp_gsC25$deltaT<-44.03*latemp_gsC25$E/(29.3*0.135*sqrt(2/(latemp_gsC25$diam*25^(-3))))
latemp_gsC25$Tsc<-latemp_gsC25$Tac - latemp_gsC25$deltaT
latemp_gsC25$wi<-0.1*(10^(10.79574*(1-273.16/(latemp_gsC25$Tsc+273.15)) - 5.02800*log10((latemp_gsC25$Tsc+273.15)/273.16)+
                            1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC25$Tsc+273.15)/273.16-1)))+
                            0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC25$Tsc+273.15))))+
                            0.78614))/101.3
latemp_gsC25$w0<-latemp_gsC25$RH/100* 0.1*(10^(10.79574*(1-273.16/(latemp_gsC25$Tac+273.15)) - 5.02800*log10((latemp_gsC25$Tac+273.15)/273.16)+
                                                 1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC25$Tac+273.15)/273.16-1)))+
                                                 0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC25$Tac+273.15))))+
                                                 0.78614))/101.3
latemp_gsC25$deltaw<-latemp_gsC25$wi-latemp_gsC25$w0
latemp_gsC25$gs<-latemp_gsC25$E/latemp_gsC25$deltaw #mmol/s/m2

###35C---c----
latemp_tRH35C$rtime<-as.POSIXct(latemp_tRH35C$TIMESTAMP,format="%Y/%m/%d %H:%M")
latemp_regs2C35<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C35)<-c("segment","int","slope","R2","Tac","Tac_sd","RH","RH_sd")

for (i in unique(latemp_35C$segment))
{
  lm1<-lm(latemp_35C$mass_w_parafilm[latemp_35C$segment==i]~as.numeric(latemp_35C$dry_time[latemp_35C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C35[i,1]<-paste(i)
  latemp_regs2C35[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C35[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C35[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean T over the course of each stem dry down
latemp_tRH35C$rtime<-as.numeric(latemp_tRH35C$rtime)
for (i in unique(latemp_35C$segment))
{
  latemp_regs2C35[i,5]<-paste(mean(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=latemp_35C$start_time[latemp_35C$segment==i]&latemp_tRH35C$rtime<=latemp_35C$end_time[latemp_35C$segment==i]],na.rm=T))
  latemp_regs2C35[i,6]<-paste(sd(latemp_tRH35C$Temp_Avg[latemp_tRH35C$rtime>=latemp_35C$start_time[latemp_35C$segment==i]&latemp_tRH35C$rtime<=latemp_35C$end_time[latemp_35C$segment==i]],na.rm=T))
  latemp_regs2C35[i,7]<-paste(mean(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=latemp_35C$start_time[latemp_35C$segment==i]&latemp_tRH35C$rtime<=latemp_35C$end_time[latemp_35C$segment==i]],na.rm=T))
  latemp_regs2C35[i,8]<-paste(sd(latemp_tRH35C$RH_Avg[latemp_tRH35C$rtime>=latemp_35C$start_time[latemp_35C$segment==i]&latemp_tRH35C$rtime<=latemp_35C$end_time[latemp_35C$segment==i]],na.rm=T))
}
latemp_regs2C35$slope<-as.numeric(latemp_regs2C35$slope)
latemp_regs2C35$R2<-as.numeric(latemp_regs2C35$R2)
latemp_regs2C35$Tac<-as.numeric(latemp_regs2C35$Tac)
latemp_regs2C35$Tac_sd<-as.numeric(latemp_regs2C35$Tac_sd)
latemp_regs2C35$RH<-as.numeric(latemp_regs2C35$RH)
latemp_regs2C35$RH_sd<-as.numeric(latemp_regs2C35$RH_sd)


latemp_regs2C35$mols_per_sec<-(-latemp_regs2C35$slope/3600)/18.01528


latemp_gsC35<-merge(latemp_regs2C35,latemp_dims[latemp_dims$treatment=="35C",],all.x=T)
latemp_gsC35$E<-(latemp_gsC35$mols_per_sec*1000)/latemp_gsC35$sa #mmol/s/m2

Tcavg35C<-mean(as.numeric(latemp_regs2C35$Tac))
#L=43.64 KJ/mol
latemp_gsC35$deltaT<-43.64*latemp_gsC35$E/(29.3*0.135*sqrt(2/(latemp_gsC35$diam*25^(-3))))
latemp_gsC35$Tsc<-latemp_gsC35$Tac - latemp_gsC35$deltaT
latemp_gsC35$wi<-0.1*(10^(10.79574*(1-273.16/(latemp_gsC35$Tsc+273.15)) - 5.02800*log10((latemp_gsC35$Tsc+273.15)/273.16)+
                            1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC35$Tsc+273.15)/273.16-1)))+
                            0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC35$Tsc+273.15))))+
                            0.78614))/101.3
latemp_gsC35$w0<-latemp_gsC35$RH/100* 0.1*(10^(10.79574*(1-273.16/(latemp_gsC35$Tac+273.15)) - 5.02800*log10((latemp_gsC35$Tac+273.15)/273.16)+
                                                 1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC35$Tac+273.15)/273.16-1)))+
                                                 0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC35$Tac+273.15))))+
                                                 0.78614))/101.3
latemp_gsC35$deltaw<-latemp_gsC35$wi-latemp_gsC35$w0
latemp_gsC35$gs<-latemp_gsC35$E/latemp_gsC35$deltaw #mmol/s/m2

###45C---D----
latemp_tRH45C$rtime<-as.POSIXct(latemp_tRH45C$TIMESTAMP,format="%Y/%m/%d %H:%M")
latemp_regs2C45<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C45)<-c("segment","int","slope","R2","Tac","Tac_sd","RH","RH_sd")

for (i in unique(latemp_45C$segment))
{
  lm1<-lm(latemp_45C$mass_w_parafilm[latemp_45C$segment==i]~as.numeric(latemp_45C$dry_time[latemp_45C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C45[i,1]<-paste(i)
  latemp_regs2C45[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C45[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C45[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean T over the course of each stem dry down
latemp_tRH45C$rtime<-as.numeric(latemp_tRH45C$rtime)
for (i in unique(latemp_45C$segment))
{
  latemp_regs2C45[i,5]<-paste(mean(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=latemp_45C$start_time[latemp_45C$segment==i]&latemp_tRH45C$rtime<=latemp_45C$end_time[latemp_45C$segment==i]],na.rm=T))
  latemp_regs2C45[i,6]<-paste(sd(latemp_tRH45C$Temp_Avg[latemp_tRH45C$rtime>=latemp_45C$start_time[latemp_45C$segment==i]&latemp_tRH45C$rtime<=latemp_45C$end_time[latemp_45C$segment==i]],na.rm=T))
  latemp_regs2C45[i,7]<-paste(mean(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=latemp_45C$start_time[latemp_45C$segment==i]&latemp_tRH45C$rtime<=latemp_45C$end_time[latemp_45C$segment==i]],na.rm=T))
  latemp_regs2C45[i,8]<-paste(sd(latemp_tRH45C$RH_Avg[latemp_tRH45C$rtime>=latemp_45C$start_time[latemp_45C$segment==i]&latemp_tRH45C$rtime<=latemp_45C$end_time[latemp_45C$segment==i]],na.rm=T))
}
latemp_regs2C45$slope<-as.numeric(latemp_regs2C45$slope)
latemp_regs2C45$R2<-as.numeric(latemp_regs2C45$R2)
latemp_regs2C45$Tac<-as.numeric(latemp_regs2C45$Tac)
latemp_regs2C45$Tac_sd<-as.numeric(latemp_regs2C45$Tac_sd)
latemp_regs2C45$RH<-as.numeric(latemp_regs2C45$RH)
latemp_regs2C45$RH_sd<-as.numeric(latemp_regs2C45$RH_sd)


latemp_regs2C45$mols_per_sec<-(-latemp_regs2C45$slope/3600)/18.01528


latemp_gsC45<-merge(latemp_regs2C45,latemp_dims[latemp_dims$treatment=="45C",],all.x=T)
latemp_gsC45$E<-(latemp_gsC45$mols_per_sec*1000)/latemp_gsC45$sa #mmol/s/m2

Tcavg45C<-mean(as.numeric(latemp_regs2C45$Tac))
#L=43.31 KJ/mol
latemp_gsC45$deltaT<-43.31*latemp_gsC45$E/(29.3*0.135*sqrt(2/(latemp_gsC45$diam*25^(-3))))
latemp_gsC45$Tsc<-latemp_gsC45$Tac - latemp_gsC45$deltaT
latemp_gsC45$wi<-0.1*(10^(10.79574*(1-273.16/(latemp_gsC45$Tsc+273.15)) - 5.02800*log10((latemp_gsC45$Tsc+273.15)/273.16)+
                            1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC45$Tsc+273.15)/273.16-1)))+
                            0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC45$Tsc+273.15))))+
                            0.78614))/101.3
latemp_gsC45$w0<-latemp_gsC45$RH/100* 0.1*(10^(10.79574*(1-273.16/(latemp_gsC45$Tac+273.15)) - 5.02800*log10((latemp_gsC45$Tac+273.15)/273.16)+
                                                 1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC45$Tac+273.15)/273.16-1)))+
                                                 0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC45$Tac+273.15))))+
                                                 0.78614))/101.3
latemp_gsC45$deltaw<-latemp_gsC45$wi-latemp_gsC45$w0
latemp_gsC45$gs<-latemp_gsC45$E/latemp_gsC45$deltaw #mmol/s/m2

###55C---E----
latemp_tRH55C$rtime<-as.POSIXct(latemp_tRH55C$TIMESTAMP,format="%Y/%m/%d %H:%M")
latemp_regs2C55<-data.frame(matrix(nrow=0,ncol=8,data=NA))
names(latemp_regs2C55)<-c("segment","int","slope","R2","Tac","Tac_sd","RH","RH_sd")

for (i in unique(latemp_55C$segment))
{
  lm1<-lm(latemp_55C$mass_w_parafilm[latemp_55C$segment==i]~as.numeric(latemp_55C$dry_time[latemp_55C$segment==i]))
  
  #print(lm1$coefficients[2])
  latemp_regs2C55[i,1]<-paste(i)
  latemp_regs2C55[i,2]<-paste(lm1$coefficients[1])
  latemp_regs2C55[i,3]<-paste(lm1$coefficients[2])
  latemp_regs2C55[i,4]<-paste(round(summary(lm1)$adj.r.squared,5))
}

#put in mean T over the course of each stem dry down
latemp_tRH55C$rtime<-as.numeric(latemp_tRH55C$rtime)
for (i in unique(latemp_55C$segment))
{
  latemp_regs2C55[i,5]<-paste(mean(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=latemp_55C$start_time[latemp_55C$segment==i]&latemp_tRH55C$rtime<=latemp_55C$end_time[latemp_55C$segment==i]],na.rm=T))
  latemp_regs2C55[i,6]<-paste(sd(latemp_tRH55C$Temp_Avg[latemp_tRH55C$rtime>=latemp_55C$start_time[latemp_55C$segment==i]&latemp_tRH55C$rtime<=latemp_55C$end_time[latemp_55C$segment==i]],na.rm=T))
  latemp_regs2C55[i,7]<-paste(mean(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=latemp_55C$start_time[latemp_55C$segment==i]&latemp_tRH55C$rtime<=latemp_55C$end_time[latemp_55C$segment==i]],na.rm=T))
  latemp_regs2C55[i,8]<-paste(sd(latemp_tRH55C$RH_Avg[latemp_tRH55C$rtime>=latemp_55C$start_time[latemp_55C$segment==i]&latemp_tRH55C$rtime<=latemp_55C$end_time[latemp_55C$segment==i]],na.rm=T))
}
latemp_regs2C55$slope<-as.numeric(latemp_regs2C55$slope)
latemp_regs2C55$R2<-as.numeric(latemp_regs2C55$R2)
latemp_regs2C55$Tac<-as.numeric(latemp_regs2C55$Tac)
latemp_regs2C55$Tac_sd<-as.numeric(latemp_regs2C55$Tac_sd)
latemp_regs2C55$RH<-as.numeric(latemp_regs2C55$RH)
latemp_regs2C55$RH_sd<-as.numeric(latemp_regs2C55$RH_sd)


latemp_regs2C55$mols_per_sec<-(-latemp_regs2C55$slope/3600)/18.01528


latemp_gsC55<-merge(latemp_regs2C55,latemp_dims[latemp_dims$treatment=="55C",],all.x=T)
latemp_gsC55$E<-(latemp_gsC55$mols_per_sec*1000)/latemp_gsC55$sa #mmol/s/m2

Tcavg55C<-mean(as.numeric(latemp_regs2C55$Tac))
#L=42.73 KJ/mol
latemp_gsC55$deltaT<-42.73*latemp_gsC55$E/(29.3*0.135*sqrt(2/(latemp_gsC55$diam*10^(-3))))
latemp_gsC55$Tsc<-latemp_gsC55$Tac - latemp_gsC55$deltaT
latemp_gsC55$wi<-0.1*(10^(10.79574*(1-273.16/(latemp_gsC55$Tsc+273.15)) - 5.02800*log10((latemp_gsC55$Tsc+273.15)/273.16)+
                                1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC55$Tsc+273.15)/273.16-1)))+
                                0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC55$Tsc+273.15))))+
                                0.78614))/101.3
latemp_gsC55$w0<-latemp_gsC55$RH/100* 0.1*(10^(10.79574*(1-273.16/(latemp_gsC55$Tac+273.15)) - 5.02800*log10((latemp_gsC55$Tac+273.15)/273.16)+
                            1.50475*10^(-4)*(1-10^(-8.2969*((latemp_gsC55$Tac+273.15)/273.16-1)))+
                            0.42873*10^(-3)*(10^(4.76955*(1-273.16/(latemp_gsC55$Tac+273.15))))+
                            0.78614))/101.3
latemp_gsC55$deltaw<-latemp_gsC55$wi-latemp_gsC55$w0
latemp_gsC55$gs<-latemp_gsC55$E/latemp_gsC55$deltaw #mmol/s/m2

latemp_gs<-rbind(latemp_gsC10,latemp_gsC25,latemp_gsC35,latemp_gsC45,latemp_gsC55)
write.csv(latemp_gs,"gbark_temperature_evapheatadj.csv")

summary(latemp_gs[latemp_gs$treatment=="10C",]$deltaT)
std.error(latemp_gs[latemp_gs$treatment=="10C",]$deltaT)
summary(latemp_gs[latemp_gs$treatment=="25C",]$deltaT)
std.error(latemp_gs[latemp_gs$treatment=="25C",]$deltaT)
summary(latemp_gs[latemp_gs$treatment=="35C",]$deltaT)
std.error(latemp_gs[latemp_gs$treatment=="35C",]$deltaT)
summary(latemp_gs[latemp_gs$treatment=="45C",]$deltaT)
std.error(latemp_gs[latemp_gs$treatment=="45C",]$deltaT)
summary(latemp_gs[latemp_gs$treatment=="55C",]$deltaT)
std.error(latemp_gs[latemp_gs$treatment=="55C",]$deltaT)




#analysis
tempadj<-read.csv("gbark_temperature_evapheatadj.csv")
tempadj$sp<-as.factor(tempadj$sp)
tempadj$TreeID<-as.factor(tempadj$TreeID)
p1<-ggplot(tempadj,aes(x=treatment,y=gbark))+
  geom_point(aes(group=interaction(sp,treatment),color=sp))+
  facet_wrap(~TreeID)
p1
p1box<-ggplot(tempadj)+
  geom_boxplot(aes(x=treatment,y=gbark))+
  geom_point((aes(x=treatment,y=gbark,group=interaction(TreeID,treatment),color=TreeID)))+
  facet_wrap(~sp)+
  theme_classic()
p1box

lm1<-lmer(gbark~sp*treatment+(1|TreeID),tempadj)
summary(lm1)  
anova(lm1)
emmeans(lm1,list(pairwise~treatment|sp),adj="Tukey")
lm_means<-emmeans(object = lm1,specs =c("sp","treatment"))
lm_means_cld<-multcomp::cld(lm_means,adjust="Tukey",Letters = letters,alpha = 0.05)
lm_means_cld

CC<-tempadj[tempadj$sp=="CC",]
LS<-tempadj[tempadj$sp=="LS",]
MG<-tempadj[tempadj$sp=="MG",]
PE<-tempadj[tempadj$sp=="PE",]
TD<-tempadj[tempadj$sp=="TD",]

lmsp1<-lm(CC$gbark~CC$treatment)
anova(lmsp1)
LSD.test(aov(gbark~treatment,data=CC),"treatment")$groups

lmsp2<-lm(LS$gbark~LS$treatment)
anova(lmsp2)
LSD.test(aov(gbark~treatment,data=LS),"treatment")$groups

lmsp3<-lm(MG$gbark~MG$treatment)
anova(lmsp3)
lmsp4<-lm(PE$gbark~PE$treatment)
anova(lmsp4)
lmsp5<-lm(TD$gbark~TD$treatment)
anova(lmsp5)

#CC LS significant, MG marginal significant

#bark traits
lm3<-lmList(gbark~LD+lenticel_size+bark_thickness_ratio | sp,data = tempadj,pool = F)
lm3
summary(lm3)
#none of them are significant
treetrait1<-dplyr::select(treetrait,c("TreeID","DBH_cm","Height_m","Crown_Height_m"))
temptrait<-merge(tempadj,treetrait1,by="TreeID")
lm4<-lmList(gbark~DBH_cm+Height_m|sp,temptrait,pool = F)
summary(lm4)
#none significant

lm5<-lm(gbark~lichen*sp*treatment,tempadj)
summary(lm5)
anova(lm5)
lm5.1<-lmer(gbark~lichen+sp+treatment+(1|TreeID),tempadj)
anova(lm5.1)
lm5.2<-lmList(gbark~lichen|sp,tempadj,pool = F)
summary(lm5.2)