setwd("E://LSU//research//gbark_temperature")
library(agricolae)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(olsrr)
library(plotrix)
temp1<-read.csv("gbark_temp_Goff_1957.csv")
temp1$TreeID<-as.factor(temp1$TreeID)
temp1$Round<-as.factor(temp1$Round)
temp1$treatment<-as.factor(temp1$treatment)
treetrait<-read.csv("tree traits.csv")
treetrait$sp<-as.factor(treetrait$sp)
treetrait1<-treetrait%>%dplyr::select(TreeID,DBH_cm,Height_m,Crown_Height_m)

summary(temp1[temp1$sp=="CC",])
std.error(temp1[temp1$sp=="CC",]$apex_dist)
sd(temp1[temp1$sp=="CC",]$apex_dist)
summary(temp1[temp1$sp=="LS",])
std.error(temp1[temp1$sp=="LS",]$apex_dist)
sd(temp1[temp1$sp=="LS",]$apex_dist)
summary(temp1[temp1$sp=="MG",])
std.error(temp1[temp1$sp=="MG",]$apex_dist)
sd(temp1[temp1$sp=="MG",]$apex_dist)
summary(temp1[temp1$sp=="TD",])
std.error(temp1[temp1$sp=="TD",]$apex_dist)
sd(temp1[temp1$sp=="TD",]$apex_dist)
summary(temp1[temp1$sp=="PE",])
std.error(temp1[temp1$sp=="PE",]$apex_dist)
sd(temp1[temp1$sp=="PE",]$apex_dist)



gbark_traits<-read.csv("temperature x gbark traits.csv")
gbark_traits$basal_diam = (gbark_traits$basal_diam1+gbark_traits$basal_diam2)/2
gbark_traits$distal_diam = (gbark_traits$distal_diam1+gbark_traits$distal_diam2)/2
summary(gbark_traits$basal_diam)
sd(gbark_traits$basal_diam)
summary(gbark_traits$distal_diam)
sd(gbark_traits$distal_diam)

summary(treetrait[treetrait$sp=="Carpinus caroliniana",])
std.error(treetrait[treetrait$sp=="Carpinus caroliniana",]$DBH_cm)
std.error(treetrait[treetrait$sp=="Carpinus caroliniana",]$Height_m)
sd(treetrait[treetrait$sp=="Carpinus caroliniana",]$DBH_cm)
sd(treetrait[treetrait$sp=="Carpinus caroliniana",]$Height_m)

summary(treetrait[treetrait$sp=="Liquidambar styraciflua",])
std.error(treetrait[treetrait$sp=="Liquidambar styraciflua",]$DBH_cm)
std.error(treetrait[treetrait$sp=="Liquidambar styraciflua",]$Height_m)
sd(treetrait[treetrait$sp=="Liquidambar styraciflua",]$DBH_cm)
sd(treetrait[treetrait$sp=="Liquidambar styraciflua",]$Height_m)

summary(treetrait[treetrait$sp=="Magnolia grandiflora",])
std.error(treetrait[treetrait$sp=="Magnolia grandiflora",]$DBH_cm)
std.error(treetrait[treetrait$sp=="Magnolia grandiflora",]$Height_m)
sd(treetrait[treetrait$sp=="Magnolia grandiflora",]$DBH_cm)
sd(treetrait[treetrait$sp=="Magnolia grandiflora",]$Height_m)

summary(treetrait[treetrait$sp=="Taxodium distichum",])
std.error(treetrait[treetrait$sp=="Taxodium distichum",]$DBH_cm)
std.error(treetrait[treetrait$sp=="Taxodium distichum",]$Height_m)
sd(treetrait[treetrait$sp=="Taxodium distichum",]$DBH_cm)
sd(treetrait[treetrait$sp=="Taxodium distichum",]$Height_m)

summary(treetrait[treetrait$sp=="Pinus elliottii",])
std.error(treetrait[treetrait$sp=="Pinus elliottii",]$DBH_cm)
std.error(treetrait[treetrait$sp=="Pinus elliottii",]$Height_m)
sd(treetrait[treetrait$sp=="Pinus elliottii",]$DBH_cm)
sd(treetrait[treetrait$sp=="Pinus elliottii",]$Height_m)


plot(treetrait$sp,treetrait$DBH_cm)

p1<-ggplot(temp1,aes(x=treatment,y=gbark))+
  geom_point(aes(group=interaction(sp,treatment),color=sp))+
  facet_wrap(~TreeID)
p1
p1box<-ggplot(temp1)+
  geom_boxplot(aes(x=treatment,y=gbark))+
  geom_point((aes(x=treatment,y=gbark,group=interaction(TreeID,treatment),color=TreeID)))+
  facet_wrap(~sp)+
  theme_classic()
p1box
  
#in the order of median-----
temp_median<-temp1 %>%group_by(sp,treatment)%>% summarize(median_gbark=median(gbark))%>%arrange(sp,desc(median_gbark))
temp2<-temp1
temp2$treatment<-factor(temp2$treatment,levels = temp_median$treatment)
p1boxmedord<-ggplot(temp2)+
  geom_boxplot(aes(x=treatment,y=gbark))+
  geom_point((aes(x=treatment,y=gbark,group=interaction(TreeID,treatment),color=TreeID)))+
  facet_wrap(~sp,scales = "free")
p1boxmedord

tempsp<-temp1%>%group_by(sp,treatment)%>%
  summarise(gbark_avg=mean(gbark))%>%
  arrange(sp,desc(gbark_avg))
p2<-ggplot(tempsp,aes(treatment,y=gbark_avg))+
  geom_point()+
  facet_wrap(~sp)
p2

#sp and trt effect-----
lm1<-lmer(gbark~sp*treatment+(1|TreeID),temp1)
summary(lm1)  
anova(lm1)
emmeans(lm1,list(pairwise~treatment|sp),adj="Tukey")
lm_means<-emmeans(object = lm1,specs =c("sp","treatment"))
lm_means_cld<-multcomp::cld(lm_means,adjust="Tukey",Letters = letters,alpha = 0.05)
lm_means_cld

#with round info
lm2<-lmer(gbark~sp+treatment*Round+(1|TreeID),temp1)
anova(lm2)
lm2.1<-lmer(gbark~sp*treatment+(1|Round),temp1)
anova(lm2.1)
#?it's the same results as the previous one
#lm2.2<-lmList(gbark~treatment | sp,data = temp1)
#lm2.2
#summary(lm2.2)

#make temperature effect continuous
lm1con<-lmer(gbark~sp*Tc + (1|TreeID),temp1)
summary(lm1con)
anova(lm1con)


CC<-temp1[temp1$sp=="CC",]
LS<-temp1[temp1$sp=="LS",]
MG<-temp1[temp1$sp=="MG",]
PE<-temp1[temp1$sp=="PE",]
TD<-temp1[temp1$sp=="TD",]

summary(CC[CC$treatment=="10C",])#3.059
std.error(CC[CC$treatment=="10C",]$gbark)#0.2383052
summary(CC[CC$treatment=="25C",])#4.407
std.error(CC[CC$treatment=="25C",]$gbark)#0.7349746
summary(CC[CC$treatment=="35C",])#2.778
std.error(CC[CC$treatment=="35C",]$gbark)#0.2086763
summary(CC[CC$treatment=="45C",])#2.979
std.error(CC[CC$treatment=="45C",]$gbark)#0.4896261
summary(CC[CC$treatment=="55C",])#2.306
std.error(CC[CC$treatment=="55C",]$gbark)#0.2467613

summary(LS[LS$treatment=="10C",])#6.543
std.error(LS[LS$treatment=="10C",]$gbark)#0.6571911
summary(LS[LS$treatment=="25C",])#8.149
std.error(LS[LS$treatment=="25C",]$gbark)#0.7757994
summary(LS[LS$treatment=="35C",])#4.835
std.error(LS[LS$treatment=="35C",]$gbark)#0.7572976
summary(LS[LS$treatment=="45C",])#3.891
std.error(LS[LS$treatment=="45C",]$gbark)#0.6209759
summary(LS[LS$treatment=="55C",])#5.584
std.error(LS[LS$treatment=="55C",]$gbark)#0.6910319

summary(MG[MG$treatment=="10C",])#5.599
std.error(MG[MG$treatment=="10C",]$gbark)#1.028109
summary(MG[MG$treatment=="25C",])#6.129
std.error(MG[MG$treatment=="25C",]$gbark)#0.9464785
summary(MG[MG$treatment=="35C",])#3.662
std.error(MG[MG$treatment=="35C",]$gbark)#0.3290681
summary(MG[MG$treatment=="45C",])#3.815
std.error(MG[MG$treatment=="45C",]$gbark)#0.8909456
summary(MG[MG$treatment=="55C",])#3.394
std.error(MG[MG$treatment=="55C",]$gbark)#0.3061158

summary(PE[PE$treatment=="10C",])#4.841
std.error(PE[PE$treatment=="10C",]$gbark)#0.5665478
summary(PE[PE$treatment=="25C",])#5.426
std.error(PE[PE$treatment=="25C",]$gbark)#0.8727509
summary(PE[PE$treatment=="35C",])#3.19
std.error(PE[PE$treatment=="35C",]$gbark)#0.2512251
summary(PE[PE$treatment=="45C",])#5.308
std.error(PE[PE$treatment=="45C",]$gbark)#0.7687821
summary(PE[PE$treatment=="55C",])#4.723
std.error(PE[PE$treatment=="55C",]$gbark)#0.4989006

summary(TD[TD$treatment=="10C",])#7.997
std.error(TD[TD$treatment=="10C",]$gbark)#2.1296
summary(TD[TD$treatment=="25C",])#8.13
std.error(TD[TD$treatment=="25C",]$gbark)#1.299471
summary(TD[TD$treatment=="35C",])#4.596
std.error(TD[TD$treatment=="35C",]$gbark)#0.7570776
summary(TD[TD$treatment=="45C",])#4.702
std.error(TD[TD$treatment=="45C",]$gbark)#0.6285674
summary(TD[TD$treatment=="55C",])#5.477
std.error(TD[TD$treatment=="55C",]$gbark)#1.041379











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
lm3<-lmList(gbark~LD+lenticel_size+bark_thickness_ratio | sp,data = temp1,pool = F)
lm3
summary(lm3)

lm3.2<-lmList(log(gbark)~LD+lenticel_size+bark_thickness_ratio | sp,data = temp1,pool = F)
lm3.2
summary(lm3.2)
#none of them are significant
#####? bark thickness ratio for LS is significant
summary(lm(gbark~bark_thickness_ratio,data = LS))
plot(LS$gbark~LS$bark_thickness_ratio)
abline(LS$gbark~LS$bark_thickness_ratio)

lm3.1<-lmList(gbark~LD+lenticel_size+bark_thickness_ratio | treatment,data = CC,pool = F)
lm3.1<-lmList(gbark~LD | treatment,data = CC,pool = F)

summary(lm3.1)

lm3.1.1<-lmer(gbark~bark_thickness_ratio*treatment+(1|TreeID),data = CC)
summary(lm3.1.1)

treetrait1<-dplyr::select(treetrait,c("TreeID","DBH_cm","Height_m","Crown_Height_m"))
temptrait<-merge(temp1,treetrait1,by="TreeID")
lm4<-lmList(gbark~DBH_cm+Height_m|sp,temptrait,pool = F)
summary(lm4)
#none significant
###?DBH and height for CC PE significant



lm5<-lm(gbark~lichen*sp*treatment,temp1)
summary(lm5)
anova(lm5)
lm5.1<-lmer(gbark~lichen+sp+treatment+(1|TreeID),temp1)
anova(lm5.1)
lm5.2<-lmList(gbark~lichen|sp,temp1,pool = F)
summary(lm5.2)

anova(lm(CC$gbark~CC$lichen))
anova(lm(LS$gbark~LS$lichen))
anova(lm(MG$gbark~MG$lichen))
anova(lm(PE$gbark~PE$lichen))
anova(lm(TD$gbark~TD$lichen))


#Only one sp. (MG) related to lichen P=0.029

p3<-ggplot(temp1,aes(x=lichen,y=gbark))+
  geom_point((aes(x=lichen,y=gbark,color=treatment)))+
  facet_wrap(~sp,scales = "free")
  p3
  temp1$lichen<-as.numeric(temp1$lichen)
  
  ols_test_normality(lm(gbark~lichen,temp1))
  leveneTest(lm(gbark~lichen,temp1))
  
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


#apex_length
apl1<-lm(CC$gbark~CC$apex_dist)
anova(apl1)

apl2<-lm(LS$gbark~LS$apex_dist)
anova(apl2)

apl3<-lm(MG$gbark~MG$apex_dist)
anova(apl3)

apl4<-lm(PE$gbark~PE$apex_dist)
anova(apl4)

apl5<-lm(TD$gbark~TD$apex_dist)
anova(apl5)
plot(apl5)

#polynomial
library(MuMIn)

#sp1 CC
plmCC2<-lm(CC$gbark~poly(CC$treatment,2))
summary(plmCC2)
plmCC3<-lm(CC$gbark~poly(CC$treatment,3))
summary(plmCC3)
plmCC4<-lm(CC$gbark~poly(CC$treatment,4))
summary(plmCC4)
print(anova(lmsp1,plmCC2,plmCC3,plmCC4))

AICc(lmsp1,plmCC2,plmCC3,plmCC4)
#lm and ply 4 best - lmsp1, plmCC4

#k-fold cross-validation with k = 10 folds
CC.shuffled<-CC[sample(nrow(CC)),]
K<-10
degree<-4
folds<-cut(seq(1,nrow(CC.shuffled)),breaks=K,labels=F)
mseCC<-matrix(data = NA,nrow = K,ncol = degree)

for(i in 1:K){
  testIndexes<-which(folds==i,arr.ind = T)
  testData<-CC.shuffled[testIndexes,]
  trainData<-CC.shuffled[-testIndexes,]
  
  for (j in 1:degree){
    fit.train=lm(gbark ~poly(treatment,j),data = trainData)
    fit.test=predict(fit.train,newdata = testData)
    mseCC[i,j] = mean((fit.test-testData$gbark)^2)
  }
}
colMeans(mseCC)
#plmCC4 best
CC$treatmentnum<-paste0(substr(CC$treatment,1,2))
pplmCC<-ggplot(CC,aes(x=as.numeric(treatmentnum),y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_bw()
pplmCC
 


#sp2 LS
lmsp2
plmLS2<-lm(LS$gbark~poly(LS$treatment,2))
summary(plmLS2)
plmLS3<-lm(LS$gbark~poly(LS$treatment,3))
summary(plmLS3)
plmLS4<-lm(LS$gbark~poly(LS$treatment,4))
summary(plmLS4)
print(anova(lmsp2,plmLS2,plmLS3,plmLS4))

AICc(lmsp2,plmLS2,plmLS3,plmLS4)
#poly 3 best - plmLS3

#k-fold cross-validation with k = 10 folds
LS.shuffled<-LS[sample(nrow(LS)),]
K<-10
degree<-4
folds<-cut(seq(1,nrow(LS.shuffled)),breaks=K,labels=F)
mseLS<-matrix(data = NA,nrow = K,ncol = degree)

for(i in 1:K){
  testIndexes<-which(folds==i,arr.ind = T)
  testData<-LS.shuffled[testIndexes,]
  trainData<-LS.shuffled[-testIndexes,]
  
  for (j in 1:degree){
    fit.train=lm(gbark ~poly(treatment,j),data = trainData)
    fit.test=predict(fit.train,newdata = testData)
    mseLS[i,j] = mean((fit.test-testData$gbark)^2)
  }
}
colMeans(mseLS)
#plmLS3 best
LS$treatmentnum<-paste0(substr(LS$treatment,1,2))
pplmLS<-ggplot(LS,aes(x=as.numeric(treatmentnum),y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,3),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_bw()
pplmLS



#sp3 MG
plmMG2<-lm(MG$gbark~poly(MG$treatment,2))
summary(plmMG2)
plmMG3<-lm(MG$gbark~poly(MG$treatment,3))
summary(plmMG3)
plmMG4<-lm(MG$gbark~poly(MG$treatment,4))
summary(plmMG4)
print(anova(lmsp3,plmMG2,plmMG3,plmMG4))

AICc(lmsp3,plmMG2,plmMG3,plmMG4)
#poly 2 best - plmMG2
#k-fold cross-validation with k = 10 folds
MG.shuffled<-MG[sample(nrow(MG)),]
K<-10
degree<-4
folds<-cut(seq(1,nrow(MG.shuffled)),breaks=K,labels=F)
mseMG<-matrix(data = NA,nrow = K,ncol = degree)

for(i in 1:K){
  testIndexes<-which(folds==i,arr.ind = T)
  testData<-MG.shuffled[testIndexes,]
  trainData<-MG.shuffled[-testIndexes,]
  
  for (j in 1:degree){
    fit.train=lm(gbark ~poly(treatment,j),data = trainData)
    fit.test=predict(fit.train,newdata = testData)
    mseMG[i,j] = mean((fit.test-testData$gbark)^2)
  }
}
colMeans(mseMG)
#plmMG2 best
MG$treatmentnum<-paste0(substr(MG$treatment,1,2))
pplmMG<-ggplot(MG,aes(x=as.numeric(treatmentnum),y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,2),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_bw()
pplmMG



#sp4 PE
plmPE2<-lm(PE$gbark~poly(PE$treatment,2))
summary(plmPE2)
plmPE3<-lm(PE$gbark~poly(PE$treatment,3))
summary(plmPE3)
plmPE4<-lm(PE$gbark~poly(PE$treatment,4))
summary(plmPE4)
print(anova(lmsp4,plmPE2,plmPE3,plmPE4))

AICc(lmsp4,plmPE2,plmPE3,plmPE4)
#lm best

#k-fold cross-validation with k = 10 folds
PE.shuffled<-PE[sample(nrow(PE)),]
K<-10
degree<-4
folds<-cut(seq(1,nrow(PE.shuffled)),breaks=K,labels=F)
msePE<-matrix(data = NA,nrow = K,ncol = degree)

for(i in 1:K){
  testIndexes<-which(folds==i,arr.ind = T)
  testData<-PE.shuffled[testIndexes,]
  trainData<-PE.shuffled[-testIndexes,]
  
  for (j in 1:degree){
    fit.train=lm(gbark ~poly(treatment,j),data = trainData)
    fit.test=predict(fit.train,newdata = testData)
    msePE[i,j] = mean((fit.test-testData$gbark)^2)
  }
}
colMeans(msePE)
#lmsp4 best
PE$treatmentnum<-paste0(substr(PE$treatment,1,2))
pplmPE<-ggplot(PE,aes(x=as.numeric(treatmentnum),y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_bw()
pplmPE



#sp5 TD
plmTD2<-lm(TD$gbark~poly(TD$treatment,2))
summary(plmTD2)
plmTD3<-lm(TD$gbark~poly(TD$treatment,3))
summary(plmTD3)
plmTD4<-lm(TD$gbark~poly(TD$treatment,4))
summary(plmTD4)
print(anova(lmsp5,plmTD2,plmTD3,plmTD4))

AICc(lmsp5,plmTD2,plmTD3,plmTD4)
#poly 2 best -plmTD2

#k-fold cross-validation with k = 10 folds
TD.shuffled<-TD[sample(nrow(TD)),]
K<-10
degree<-4
folds<-cut(seq(1,nrow(TD.shuffled)),breaks=K,labels=F)
mseTD<-matrix(data = NA,nrow = K,ncol = degree)

for(i in 1:K){
  testIndexes<-which(folds==i,arr.ind = T)
  testData<-TD.shuffled[testIndexes,]
  trainData<-TD.shuffled[-testIndexes,]
  
  for (j in 1:degree){
    fit.train=lm(gbark ~poly(treatment,j),data = trainData)
    fit.test=predict(fit.train,newdata = testData)
    mseTD[i,j] = mean((fit.test-testData$gbark)^2)
  }
}
colMeans(mseTD)
#lmsp4 best
TD$treatmentnum<-paste0(substr(TD$treatment,1,2))
pplmTD<-ggplot(TD,aes(x=as.numeric(treatmentnum),y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,2),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='TemTDrature Treatments')+
  theme_bw()
pplmTD

pplmCC#poly 1
pplmLS#poly 3
pplmMG#poly 2
pplmPE#poly 1
pplmTD#poly 2

library(dplyr)
CC_summary<-CC %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))
LS_summary<-LS %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))

LSspring_summary<-LSspring %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))

MG_summary<-MG %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))
PE_summary<-PE %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))
PEspring_summary<-PEspring %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))

TD_summary<-TD %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))

####with Tc (numeric)
library(gridExtra)
source('F://LSU//research//gbark_temperature//Rtheme//ggplot_theme_Publication-2.R')
#sp1 CC
lmCC1num<-lm(CC$gbark~poly(CC$Tc,1))
summary(lmCC1num)
plmCC2num<-lm(CC$gbark~poly(CC$Tc,2))
summary(plmCC2num)
plmCC3num<-lm(CC$gbark~poly(CC$Tc,3))
summary(plmCC3num)
plmCC4num<-lm(CC$gbark~poly(CC$Tc,4))
summary(plmCC4num)
print(anova(lmCC1num,plmCC2num,plmCC3num,plmCC4num))

AICc(lmCC1num,plmCC2num,plmCC3num,plmCC4num)
#ply 4 best - plmCC4num

pplmCCnum<-ggplot(CC,aes(x=Tc,y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,4),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_minimal()
pplmCCnum

pplmCCsum<-ggplot(CC,aes(x=Tc,y=gbark))+
  #geom_point()+
  stat_smooth(linetype = "dashed",method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_point(data=CC_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=CC_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                    ymin=mean_gbark - se_gbark,
                    ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=CC_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                    xmax = mean_Tc + se_Tc),
                height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Capinus~caroliniana))),tag = "A")+
  theme_Publication(base_size = 14)+
  theme(panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
pplmCCsum  
  
#sp2 LS
lmLS1num<-lm(LS$gbark~poly(LS$Tc,1))
summary(lmLS1num)
plmLS2num<-lm(LS$gbark~poly(LS$Tc,2))
summary(plmLS2num)
plmLS3num<-lm(LS$gbark~poly(LS$Tc,3))
summary(plmLS3num)
plmLS4num<-lm(LS$gbark~poly(LS$Tc,4))
summary(plmLS4num)
print(anova(lmLS1num,plmLS2num,plmLS3num,plmLS4num))

AICc(lmLS1num,plmLS2num,plmLS3num,plmLS4num)
#ply 3 best - plmLS3num

pplmLSnum<-ggplot(LS,aes(x=Tc,y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,3),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_minimal()
pplmLSnum

pplmLSsum<-ggplot(LS,aes(x=Tc,y=gbark))+
  #geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,3),linewidth=0.5,color="black",fill="darkgreen")+
  stat_smooth(data=LSspring, aes(x=Tc,y=gbark),method = 'lm',formula = y~poly(x,2),linewidth=0.5,color="black",fill="chocolate4")+
  geom_point(data=LS_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=LS_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_point(data=LSspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "chocolate4")+
  geom_errorbar(data=LSspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=LS_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  geom_errorbarh(data=LSspring_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Liquidamber~styraciflua))),tag = "B")+
  theme_Publication(base_size = 14)+
  theme(panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pplmLSsum  

#spring
lmLSspring1num<-lm(LSspring$gbark~poly(LSspring$Tc,1))
summary(lmLSspring1num)
plmLSspring2num<-lm(LSspring$gbark~poly(LSspring$Tc,2))
summary(plmLSspring2num)
plmLSspring3num<-lm(LSspring$gbark~poly(LSspring$Tc,3))
summary(plmLSspring3num)
plmLSspring4num<-lm(LSspring$gbark~poly(LSspring$Tc,4))
summary(plmLSspring4num)
print(anova(lmLSspring1num,plmLSspring2num,plmLSspring3num,plmLSspring4num))

AICc(lmLSspring1num,plmLSspring2num,plmLSspring3num,plmLSspring4num)
#ply 2 best


#sp3 MG
lmMG1num<-lm(MG$gbark~poly(MG$Tc,1))
summary(lmMG1num)
plmMG2num<-lm(MG$gbark~poly(MG$Tc,2))
summary(plmMG2num)
plmMG3num<-lm(MG$gbark~poly(MG$Tc,3))
summary(plmMG3num)
plmMG4num<-lm(MG$gbark~poly(MG$Tc,4))
summary(plmMG4num)
print(anova(lmMG1num,plmMG2num,plmMG3num,plmMG4num))

AICc(lmMG1num,plmMG2num,plmMG3num,plmMG4num)
#ply 1 best - lmMG1num

pplmMGnum<-ggplot(MG,aes(x=Tc,y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_minimal()
pplmMGnum

pplmMGsum<-ggplot(MG,aes(x=Tc,y=gbark))+
  #geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_point(data=MG_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=MG_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=MG_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Magnolia~grandifolia))),tag = "C")+
  theme_Publication(base_size = 14)+
  theme(panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pplmMGsum  

#sp4 PE
lmPE1num<-lm(PE$gbark~poly(PE$Tc,1))
summary(lmPE1num)
plmPE2num<-lm(PE$gbark~poly(PE$Tc,2))
summary(plmPE2num)
plmPE3num<-lm(PE$gbark~poly(PE$Tc,3))
summary(plmPE3num)
plmPE4num<-lm(PE$gbark~poly(PE$Tc,4))
summary(plmPE4num)
print(anova(lmPE1num,plmPE2num,plmPE3num,plmPE4num))

AICc(lmPE1num,plmPE2num,plmPE3num,plmPE4num)
#ply 1 best - lmPE4num (not significant)

pplmPEnum<-ggplot(PE,aes(x=Tc,y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_minimal()
pplmPEnum

pplmPEsum<-ggplot(PE,aes(x=Tc,y=gbark))+
  #geom_point()+
  stat_smooth(linetype = "dashed",method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  stat_smooth(data = PEspring,aes(x=Tc,y=gbark),method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="chocolate4")+
  geom_point(data=PE_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_point(data=PEspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "chocolate4")+
  geom_errorbar(data=PE_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=PE_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  geom_errorbar(data=PEspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=PEspring_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Pinus~elliottii))),tag = "D")+
  theme_Publication(base_size = 14)+
  theme(panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pplmPEsum  

#spring
lmPEspring1num<-lm(PEspring$gbark~poly(PEspring$Tc,1))
summary(lmPEspring1num)
plmPEspring2num<-lm(PEspring$gbark~poly(PEspring$Tc,2))
summary(plmPEspring2num)
plmPEspring3num<-lm(PEspring$gbark~poly(PEspring$Tc,3))
summary(plmPEspring3num)
plmPEspring4num<-lm(PEspring$gbark~poly(PEspring$Tc,4))
summary(plmPEspring4num)
print(anova(lmPEspring1num,plmPEspring2num,plmPEspring3num,plmPEspring4num))

AICc(lmPEspring1num,plmPEspring2num,plmPEspring3num,plmPEspring4num)
#ply 1 best

#sp5 TD
lmTD1num<-lm(TD$gbark~poly(TD$Tc,1))
summary(lmTD1num)
plmTD2num<-lm(TD$gbark~poly(TD$Tc,2))
summary(plmTD2num)
plmTD3num<-lm(TD$gbark~poly(TD$Tc,3))
summary(plmTD3num)
plmTD4num<-lm(TD$gbark~poly(TD$Tc,4))
summary(plmTD4num)
print(anova(lmTD1num,plmTD2num,plmTD3num,plmTD4num))

AICc(lmTD1num,plmTD2num,plmTD3num,plmTD4num)
#ply 1 best - lmTD1num (not significant)

pplmTDnum<-ggplot(TD,aes(x=Tc,y=gbark))+
  geom_point()+
  stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments')+
  theme_minimal()
pplmTDnum

pplmTDsum<-ggplot(TD,aes(x=Tc,y=gbark))+
  #geom_point()+
  stat_smooth(linetype = "dashed",method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_point(data=TD_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=TD_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=TD_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),x='Temperature Treatments (°C)',title = expression(paste(italic(Taxodium~distchum))),tag = "E")+
  theme_Publication(base_size = 14)+
  theme(panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pplmTDsum  

pplmCCsum#ply 4
pplmLSsum#ply 3
pplmMGsum#ply 1
pplmPEsum#ply 1 (not significant)
pplmTDsum#ply 1 (not significant)

library(patchwork)
wrap_plots(list(pplmCCsum,pplmLSsum,pplmMGsum,pplmPEsum,pplmTDsum),ncol=1)

