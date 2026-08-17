setwd("F://LSU//research//gbark_temperature")
library(agricolae)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(olsrr)
library(plotrix)
library(multcompView)
temp1spring<-read.csv("gbark_temp_Goff_1957_spring.csv")
temp1spring$TreeID<-as.factor(temp1spring$TreeID)
temp1spring$Round<-as.factor(temp1spring$Round)
temp1spring$treatment<-as.factor(temp1spring$treatment)
treetraitspring<-read.csv("tree traits spring.csv")
treetraitspring$sp<-as.factor(treetraitspring$sp)
treetraitspring1<-treetraitspring%>%dplyr::select(TreeID,DBH_cm,Height_m,sp)%>%na.omit()

summary(treetraitspring1[treetraitspring1$sp=="Liquidambar styraciflua",])
#std.error(na.omit(treetraitspring[treetraitspring$sp=="Liquidambar styraciflua",]$DBH_cm))
#std.error(na.omit(treetraitspring[treetraitspring$sp=="Liquidambar styraciflua",]$Height_m))
mean(treetraitspring1[treetraitspring1$sp=="Liquidambar styraciflua",]$DBH_cm)
mean(treetraitspring1[treetraitspring1$sp=="Liquidambar styraciflua",]$Height_m)
sd(treetraitspring1[treetraitspring1$sp=="Liquidambar styraciflua",]$DBH_cm)
sd(treetraitspring1[treetraitspring1$sp=="Liquidambar styraciflua",]$Height_m)

summary(treetraitspring1[treetraitspring1$sp=="Pinus elliottii",])
#std.error(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$DBH_cm)
#std.error(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$Height_m)
mean(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$DBH_cm)
mean(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$Height_m)
sd(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$DBH_cm)
sd(treetraitspring[treetraitspring$sp=="Pinus elliottii",]$Height_m)

plot(treetraitspring$sp,treetraitspring$DBH_cm)

temp1spring<-temp1spring[!temp1spring$treatment=="55C",]

p1<-ggplot(temp1spring,aes(x=treatment,y=gbark))+
  geom_point(aes(group=interaction(sp,treatment),color=sp))+
  facet_wrap(~TreeID)
p1
p1box<-ggplot(temp1spring)+
  geom_boxplot(aes(x=treatment,y=gbark))+
  geom_point((aes(x=treatment,y=gbark,group=interaction(TreeID,treatment),color=TreeID)))+
  facet_wrap(~sp)+
  theme_classic()
p1box
  
#in the order of median-----
temp_median<-temp1spring %>%group_by(sp,treatment)%>% summarize(median_gbark=median(gbark))%>%arrange(sp,desc(median_gbark))
temp2<-temp1spring
temp2$treatment<-factor(temp2$treatment,levels = temp_median$treatment)
p1boxmedord<-ggplot(temp2)+
  geom_boxplot(aes(x=treatment,y=gbark))+
  geom_point((aes(x=treatment,y=gbark,group=interaction(TreeID,treatment),color=TreeID)))+
  facet_wrap(~sp,scales = "free")
p1boxmedord

tempsp<-temp1spring%>%group_by(sp,treatment)%>%
  summarise(gbark_avg=mean(gbark))%>%
  arrange(sp,desc(gbark_avg))
p2<-ggplot(tempsp,aes(treatment,y=gbark_avg))+
  geom_point()+
  facet_wrap(~sp)
p2

#sp and trt effect-----
lm1<-lmer(gbark~sp*treatment+(1|TreeID),temp1spring)
summary(lm1)  
anova(lm1)
emmeans(lm1,list(pairwise~treatment|sp),adj="Tukey")
lm_means<-emmeans(object = lm1,specs =c("sp","treatment"))
lm_means_cld<-multcomp::cld(lm_means,adj="Tukey",Letters = letters,alpha = 0.05)
lm_means_cld

#with round info
lm2<-lmer(gbark~sp+treatment*Round+(1|TreeID),temp1spring)
anova(lm2)
lm2.1<-lmer(gbark~sp*treatment+(1|Round),temp1spring)
anova(lm2.1)
#?it's the same results as the previous one
#lm2.2<-lmList(gbark~treatment | sp,data = temp1spring)
#lm2.2
#summary(lm2.2)

LSspring<-temp1spring[temp1spring$sp=="LS",]
PEspring<-temp1spring[temp1spring$sp=="PE",]


lmsp2<-lm(LSspring$gbark~LSspring$treatment)
anova(lmsp2)
LSD.test(aov(gbark~treatment,data=LSspring),"treatment")$groups#inappropriate


lmsp4<-lm(PEspring$gbark~PEspring$treatment)
anova(lmsp4)
LSD.test(aov(gbark~treatment,data=PEspring),"treatment")$groups
#PE significant, LS marginal significant

lm1lsspring<-lmer(gbark~treatment+(1|TreeID),LSspring)
summary(lm1lsspring)
anova(lm1lsspring)
lm1lsspring_means<-emmeans(object = lm1lsspring,specs ="treatment")
lm1lsspring_means_cld<-multcomp::cld(lm1lsspring_means,Letters = letters,alpha = 0.05)
lm1lsspring_means_cld


lm1pespring<-lmer(gbark~treatment+(1|TreeID),PEspring)
summary(lm1pespring)
anova(lm1pespring)
lm1pespring_means<-emmeans(object = lm1pespring,specs ="treatment")
lm1pespring_means_cld<-multcomp::cld(lm1pespring_means,Letters = letters,alpha = 0.05)
lm1pespring_means_cld


##bark traits

lm3<-lm(gbark~LD+lenticel_size+bark_thickness_ratio,data = LSspring,pool = F)
lm3
summary(lm3)
#none of them are significant
treetraitspring1<-dplyr::select(treetraitspring,c("TreeID","DBH_cm","Height_m"))
temptraitspring<-merge(temp1spring,treetraitspring1,by="TreeID")
lm4<-lmList(gbark~DBH_cm+Height_m|sp,temptraitspring,pool = F)
summary(lm4)
#none significant

lm5<-lm(gbark~lichen*sp*treatment,temp1spring)
summary(lm5)
anova(lm5)
lm5.1<-lmer(gbark~lichen+sp+treatment+(1|TreeID),temp1spring)
anova(lm5.1)
lm5.2<-lmList(gbark~lichen|sp,temp1spring[temp1spring$lichen>0,],pool = F)
summary(lm5.2)

#anova(lm(CC$gbark~CC$lichen))
anova(lm(LSspring$gbark~LSspring$lichen))
#anova(lm(MG$gbark~MG$lichen))
anova(lm(PEspring$gbark~PEspring$lichen))
#anova(lm(TD$gbark~TD$lichen))


#none sp. related to lichen 

p3<-ggplot(temp1spring[temp1spring$lichen>0,],aes(x=lichen,y=gbark))+
  geom_point((aes(x=lichen,y=gbark,color=treatment)))+
  facet_wrap(~sp,scales = "free")
  p3
  
  ols_test_normality(lm(gbark~lichen,temp1spring))
  leveneTest(lm(gbark~lichen,temp1spring))
  
#Sensor Validation-----
s2523<-read.csv("F://LSU//research//gbark_temperature//Sensor Validation//2523.csv")
s2524<-read.csv("F://LSU//research//gbark_temperature//Sensor Validation//2524.csv")
s2528<-read.csv("F://LSU//research//gbark_temperature//Sensor Validation//2528.csv")
s4261<-read.csv("F://LSU//research//gbark_temperature//Sensor Validation//4261.csv")
s4307<-read.csv("F://LSU//research//gbark_temperature//Sensor Validation//4307.csv")

s2523$TIMESTAMP<-as.POSIXct(s2523$TIMESTAMP)
s2524$TIMESTAMP<-as.POSIXct(s2524$TIMESTAMP)
s2528$TIMESTAMP<-as.POSIXct(s2528$TIMESTAMP)
s4261$TIMESTAMP<-as.POSIXct(s4261$TIMESTAMP)
s4307$TIMESTAMP<-as.POSIXct(s4307$TIMESTAMP)

library(ggplot2)
psensortemp<-ggplot()+
  geom_line(data=s2523,mapping=aes(x=TIMESTAMP,y=Temp_Avg_C),color="red")+
  geom_line(data=s2524,mapping=aes(x=TIMESTAMP,y=Temp_Avg_C),color="orange")+
  geom_line(data=s2528,mapping=aes(x=TIMESTAMP,y=Temp_Avg_C),color="yellow")+
  geom_line(data=s4261,mapping=aes(x=TIMESTAMP,y=Temp_Avg_C),color="green")+
  geom_line(data=s4307,mapping=aes(x=TIMESTAMP,y=Temp_Avg_C),color="blue")+
  scale_x_datetime(date_breaks = "8 hour", date_labels =  "%Y/%m/%d %T")+
  scale_y_continuous(breaks =seq(0,60,2))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1))
  
  
psensortemp
  
psensorrh<-ggplot()+
  geom_line(data=s2523,mapping=aes(x=TIMESTAMP,y=RH_Avg_.),color="red")+
  geom_line(data=s2524,mapping=aes(x=TIMESTAMP,y=RH_Avg_.),color="orange")+
  geom_line(data=s2528,mapping=aes(x=TIMESTAMP,y=RH_Avg_.),color="yellow")+
  geom_line(data=s4261,mapping=aes(x=TIMESTAMP,y=RH_Avg_.),color="green")+
  geom_line(data=s4307,mapping=aes(x=TIMESTAMP,y=RH_Avg_.),color="blue")+
  scale_x_datetime(date_breaks = "8 hour", date_labels =  "%Y/%m/%d %T")+
  scale_y_continuous(breaks =seq(0,100,2))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1))


psensorrh
  



#

#Compare season effect-------
temp1merge<-read.csv("gbark_temp_Goff_1957_merge.csv")
temp1merge$TreeID<-as.factor(temp1merge$TreeID)
temp1merge$Round<-as.factor(temp1merge$Round)
temp1merge$treatment<-as.factor(temp1merge$treatment)
lm6<-lmer(gbark~Season+sp+Tc+(1|TreeID),temp1merge)
summary(lm6)

#lm7.1<-lmer(lenticel_size~sp+Tc+Season+(1|TreeID),temp1merge)
#summary(lm7.1)
#anova(lm7.1)
lmls7.1<-lmer(lenticel_size~Season+(1|TreeID),temp1merge[temp1merge$sp=="LS",])
summary(lmls7.1)
anova(lmls7.1)
##season significant (spring smaller lenticel_size)

#lmpe7.1<-lmer(lenticel_size~Tc+Season+(1|TreeID),temp1merge[temp1merge$sp=="PE",])
#anova(lmpe7.1)

lmls7.2<-lmer(LD~sp+Season+(1|TreeID),temp1merge[temp1merge$sp=="LS",])
summary(lmls7.2)
anova(lmls7.2)
##non sig

lm7.3<-lmer(bark_thickness_ratio~sp+Season+(1|TreeID),temp1merge)
summary(lm7.3)
anova(lm7.3)
lmls7.3<-lmer(bark_thickness_ratio~Season+(1|TreeID),temp1merge[temp1merge$sp=="LS",])
anova(lmls7.3)
lmpe7.3<-lmer(bark_thickness_ratio~Season+(1|TreeID),temp1merge[temp1merge$sp=="PE",])
anova(lmpe7.3)
#non sig

