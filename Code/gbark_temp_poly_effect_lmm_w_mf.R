source("E:/LSU/research/gbark_temperature/gbark_temp_data_analyses_w_mf.R") 


####with Tc (numeric)
#####
library(gridExtra)
source('E://LSU//research//gbark_temperature//Rtheme//ggplot_theme_Publication-2.R')
#sp1 CC

lmCC1numlmmlog<-lmer(log(gbark)~poly(CC$Tc,1)+(1|TreeID),data = CC,REML = F)
summary(lmCC1numlmmlog)

plmCC2numlmmlog<-lmer(log(gbark)~poly(CC$Tc,2)+(1|TreeID),data = CC,REML = F)
summary(plmCC2numlog)

plmCC3numlmmlog<-lmer(log(gbark)~poly(CC$Tc,3)+(1|TreeID),data = CC,REML = F)
summary(plmCC3numlog)


print(anova(lmCC1numlmmlog,plmCC2numlmmlog,plmCC3numlmmlog))
AICc(lmCC1numlmmlog,plmCC2numlmmlog,plmCC3numlmmlog)



#select linear
newCC <- data.frame(Tc = seq(5, 60, length.out = 300))
newCC$Tc<-as.numeric(newCC$Tc)
#newCC <- newCC[order(newCC$Tc), ]
predCC <- predict(lmCC1numlog, newdata = newCC,se.fit = T)

newCC$fit  <- exp(predCC$fit)
newCC$lwr  <- exp(predCC$fit - 1.96 * predCC$se.fit)
newCC$upr  <- exp(predCC$fit + 1.96 * predCC$se.fit)

pplmCCsum<-ggplot(CC,aes(x=Tc,y=gbark))+
  #geom_point()+
  #stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_ribbon(data = newCC, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="darkgreen") +
  geom_line(data = newCC, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=CC_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=CC_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=CC_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Carpinus~caroliniana))),tag = "A")+
  theme_Publication(base_size = 14)+
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
pplmCCsum  

#sp2 LS
lmLS1numlmmlog<-lmer(log(gbark)~poly(LS$Tc,1)+(1|TreeID),data = LS,REML = F)
summary(lmLS1numlmmlog)

plmLS2numlmmlog<-lmer(log(gbark)~poly(LS$Tc,2)+(1|TreeID),data = LS,REML = F)
summary(plmLS2numlog)

plmLS3numlmmlog<-lmer(log(gbark)~poly(LS$Tc,3)+(1|TreeID),data = LS,REML = F)
summary(plmLS3numlog)


print(anova(lmLS1numlmmlog,plmLS2numlmmlog,plmLS3numlmmlog))
AICc(lmLS1numlmmlog,plmLS2numlmmlog,plmLS3numlmmlog)


#ply 3 best - plmLS3num

print(anova(lmLS1numlog,plmLS2numlog,plmLS3numlog))

AICc(lmLS1numlog,plmLS2numlog,plmLS3numlog)
#ply 3 best - plmLS3numlog

plmLS3numlog<-lm(log(gbark)~poly(Tc,3,raw = T),data = LS)
LS_Tc_min <- min(LS$Tc, na.rm = TRUE)
LS_Tc_max <- max(LS$Tc, na.rm = TRUE)
newLS <- data.frame(Tc = seq(LS_Tc_min, LS_Tc_max, length.out = 300))
#newLS <- data.frame(Tc = seq(5, 60, length.out = 300))
newLS$Tc<-as.numeric(newLS$Tc)
predLS <- predict(plmLS3numlog, newdata = newLS,se.fit = T)
newLS$fit  <- exp(predLS$fit)
newLS$lwr  <- exp(predLS$fit - 1.96 * predLS$se.fit)
newLS$upr  <- exp(predLS$fit + 1.96 * predLS$se.fit)



#spring
lmLSspring1numlmmlog<-lmer(log(gbark)~poly(LSspring$Tc,1)+(1|TreeID),data = LSspring,REML = F)
summary(lmLSspring1numlmmlog)

plmLSspring2numlmmlog<-lmer(log(gbark)~poly(LSspring$Tc,2)+(1|TreeID),data = LSspring,REML = F)
summary(plmLSspring2numlog)
anova(plmLSspring2numlog)

plmLSspring3numlmmlog<-lmer(log(gbark)~poly(LSspring$Tc,3)+(1|TreeID),data = LSspring,REML = F)
summary(plmLSspring3numlog)


print(anova(lmLSspring1numlmmlog,plmLSspring2numlmmlog,plmLSspring3numlmmlog))
AICc(lmLSspring1numlmmlog,plmLSspring2numlmmlog,plmLSspring3numlmmlog)

#ply 2 best

print(anova(lmLSspring1numlog,plmLSspring2numlog,plmLSspring3numlog))

AICc(lmLSspring1numlog,plmLSspring2numlog,plmLSspring3numlog)
#ply 2 best

plmLSspring2numlog<-lm(log(gbark)~poly(Tc,2,raw=T),data = LSspring)
LS_Tc_min2 <- min(LSspring$Tc, na.rm = TRUE)
LS_Tc_max2 <- max(LSspring$Tc, na.rm = TRUE)
newLSspring <- data.frame(Tc = seq(LS_Tc_min2, LS_Tc_max2, length.out = 300))
#newLSspring <- data.frame(Tc = seq(5, LS_Tc_max2, length.out = 300))
newLSspring$Tc<-as.numeric(newLSspring$Tc)
#newLS <- newLS[order(newLS$Tc), ]
predLSspring <- predict(plmLSspring2numlog, newdata = newLSspring,se.fit = T)

newLSspring$fit  <- exp(predLSspring$fit)
newLSspring$lwr  <- exp(predLSspring$fit - 1.96 * predLSspring$se.fit)
newLSspring$upr  <- exp(predLSspring$fit + 1.96 * predLSspring$se.fit)

newLS <- newLS %>%
  mutate(Tc = as.numeric(Tc)) %>%
  arrange(Tc) %>%
  filter(is.finite(fit), is.finite(lwr), is.finite(upr))

newLSspring <- newLSspring %>%
  mutate(Tc = as.numeric(Tc)) %>%
  arrange(Tc) %>%
  filter(is.finite(fit), is.finite(lwr), is.finite(upr))

pplmLSsum<-ggplot(LS,aes(x=Tc,y=gbark))+
  #geom_point()+
  #stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_ribbon(data = newLS, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="darkgreen") +
  geom_line(data = newLS, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=LS_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=LS_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=LS_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  geom_ribbon(data = newLSspring, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="chocolate4") +
  geom_line(data = newLSspring, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=LSspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "chocolate4")+
  geom_errorbar(data=LSspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                               ymin=mean_gbark - se_gbark,
                                                               ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=LSspring_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                                xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Liquidamber~styraciflua))),tag = "B")+
  theme_Publication(base_size = 14)+
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  scale_y_continuous(limits = c(0,12),breaks = seq(0,12,by = 4))+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))

pplmLSsum  





#sp3 MG
lmMG1numlmmlog<-lmer(log(gbark)~poly(MG$Tc,1)+(1|TreeID),data = MG,REML = F)
summary(lmMG1numlmmlog)

plmMG2numlmmlog<-lmer(log(gbark)~poly(MG$Tc,2)+(1|TreeID),data = MG,REML = F)
summary(plmMG2numlog)

plmMG3numlmmlog<-lmer(log(gbark)~poly(MG$Tc,3)+(1|TreeID),data = MG,REML = F)
summary(plmMG3numlog)


print(anova(lmMG1numlmmlog,plmMG2numlmmlog,plmMG3numlmmlog))
AICc(lmMG1numlmmlog,plmMG2numlmmlog,plmMG3numlmmlog)
#ply 1 best - lmMG1numlog


lmMG1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data=MG)
newMG <- data.frame(Tc = seq(5, 60, length.out = 300))
newMG$Tc<-as.numeric(newMG$Tc)
#newMG <- newMG[order(newMG$Tc), ]
predMG <- predict(lmMG1numlog, newdata = newMG,se.fit = T)

newMG$fit  <- exp(predMG$fit)
newMG$lwr  <- exp(predMG$fit - 1.96 * predMG$se.fit)
newMG$upr  <- exp(predMG$fit + 1.96 * predMG$se.fit)



pplmMGsum<-ggplot(MG,aes(x=Tc,y=gbark))+
  #geom_point()+
  #stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_ribbon(data = newMG, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="darkgreen") +
  geom_line(data = newMG, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=MG_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=MG_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=MG_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Magnolia~grandiflora))),tag = "C")+
  theme_Publication(base_size = 14)+
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
pplmMGsum  


#sp4 PE
lmPE1numlmmlog<-lmer(log(gbark)~poly(PE$Tc,1)+(1|TreeID),data = PE,REML = F)
summary(lmPE1numlmmlog)

plmPE2numlmmlog<-lmer(log(gbark)~poly(PE$Tc,2)+(1|TreeID),data = PE,REML = F)
summary(plmPE2numlog)

plmPE3numlmmlog<-lmer(log(gbark)~poly(PE$Tc,3)+(1|TreeID),data = PE,REML = F)
summary(plmPE3numlog)


print(anova(lmPE1numlmmlog,plmPE2numlmmlog,plmPE3numlmmlog))
AICc(lmPE1numlmmlog,plmPE2numlmmlog,plmPE3numlmmlog)


#ply 1 best  (not significant)

lmPE1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data=PE)
newPE <- data.frame(Tc = seq(5, 60, length.out = 300))
newPE$Tc<-as.numeric(newPE$Tc)
#newPE <- newPE[order(newPE$Tc), ]
predPE <- predict(lmPE1numlog, newdata = newPE,se.fit = T)

newPE$fit  <- exp(predPE$fit)
newPE$lwr  <- exp(predPE$fit - 1.96 * predPE$se.fit)
newPE$upr  <- exp(predPE$fit + 1.96 * predPE$se.fit)


#spring
lmPEspring1numlmmlog<-lmer(log(gbark)~poly(PEspring$Tc,1)+(1|TreeID),data = PEspring,REML = F)
summary(lmPEspring1numlmmlog)

plmPEspring2numlmmlog<-lmer(log(gbark)~poly(PEspring$Tc,2)+(1|TreeID),data = PEspring,REML = F)
summary(plmPEspring2numlog)

plmPEspring3numlmmlog<-lmer(log(gbark)~poly(PEspring$Tc,3)+(1|TreeID),data = PEspring,REML = F)
summary(plmPEspring3numlog)


print(anova(lmPEspring1numlmmlog,plmPEspring2numlmmlog,plmPEspring3numlmmlog))
AICc(lmPEspring1numlmmlog,plmPEspring2numlmmlog,plmPEspring3numlmmlog)

#ply 1 best

lmPEspring1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data = PEspring)
PE_Tc_min2 <- min(PEspring$Tc, na.rm = TRUE)
PE_Tc_max2 <- max(PEspring$Tc, na.rm = TRUE)
newPEspring <- data.frame(Tc = seq(PE_Tc_min2, PE_Tc_max2, length.out = 300))
#newPEspring <- data.frame(Tc = seq(5, PE_Tc_max2, length.out = 300))
newPEspring$Tc<-as.numeric(newPEspring$Tc)
#newPE <- newPE[order(newPE$Tc), ]
predPEspring <- predict(lmPEspring1numlog, newdata = newPEspring,se.fit = T)

newPEspring$fit  <- exp(predPEspring$fit)
newPEspring$lwr  <- exp(predPEspring$fit - 1.96 * predPEspring$se.fit)
newPEspring$upr  <- exp(predPEspring$fit + 1.96 * predPEspring$se.fit)


pplmPEsum<-ggplot(PE,aes(x=Tc,y=gbark))+
  #geom_point()+
  #stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_ribbon(data = newPE, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="darkgreen") +
  geom_line(data = newPE, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5,linetype = "dashed") +
  geom_point(data=PE_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=PE_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=PE_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  geom_ribbon(data = newPEspring, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="chocolate4") +
  geom_line(data = newPEspring, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=PEspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "chocolate4")+
  geom_errorbar(data=PEspring_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                               ymin=mean_gbark - se_gbark,
                                                               ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=PEspring_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                                xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Pinus~elliottii))),tag = "D")+
  theme_Publication(base_size = 14)+
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  #scale_y_continuous(limits = c(0,12),breaks = seq(0,12,by = 4))+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))

pplmPEsum  



#sp5 TD
lmTD1numlmmlog<-lmer(log(gbark)~poly(TD$Tc,1)+(1|TreeID),data = TD,REML = F)
summary(lmTD1numlmmlog)

plmTD2numlmmlog<-lmer(log(gbark)~poly(TD$Tc,2)+(1|TreeID),data = TD,REML = F)
summary(plmTD2numlog)

plmTD3numlmmlog<-lmer(log(gbark)~poly(TD$Tc,3)+(1|TreeID),data = TD,REML = F)
summary(plmTD3numlog)


print(anova(lmTD1numlmmlog,plmTD2numlmmlog,plmTD3numlmmlog))
AICc(lmTD1numlmmlog,plmTD2numlmmlog,plmTD3numlmmlog)

#ply 1 best - lmTD1num 

lmTD1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data=TD)
newTD <- data.frame(Tc = seq(5, 60, length.out = 300))
newTD$Tc<-as.numeric(newTD$Tc)
#newTD <- newTD[order(newTD$Tc), ]
predTD <- predict(lmTD1numlog, newdata = newTD,se.fit = T)

newTD$fit  <- exp(predTD$fit)
newTD$lwr  <- exp(predTD$fit - 1.96 * predTD$se.fit)
newTD$upr  <- exp(predTD$fit + 1.96 * predTD$se.fit)

pplmTDsum<-ggplot(TD,aes(x=Tc,y=gbark))+
  #geom_point()+
  #stat_smooth(method = 'lm',formula = y~poly(x,1),linewidth=0.5,color="black",fill="darkgreen")+
  geom_ribbon(data = newTD, aes(x = Tc, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.4,fill="darkgreen") +
  geom_line(data = newTD, aes(x = Tc, y = fit),
            inherit.aes = FALSE, linewidth = 0.5) +
  geom_point(data=TD_summary, inherit.aes = FALSE,aes(x=mean_Tc,y=mean_gbark),size = 3, color = "darkgreen")+
  geom_errorbar(data=TD_summary, inherit.aes = FALSE,aes(x=mean_Tc,
                                                         ymin=mean_gbark - se_gbark,
                                                         ymax=mean_gbark + se_gbark),
                width=0.5)+
  geom_errorbarh(data=TD_summary, inherit.aes = FALSE,aes(y=mean_gbark,xmin = mean_Tc - se_Tc,
                                                          xmax = mean_Tc + se_Tc),
                 height = 0.1)+
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y='',x='Air Temperature (°C)',title = expression(paste(italic(Taxodium~distichum))),tag = "E")+
  theme_Publication(base_size = 14)+
  theme(axis.title.x = element_text(face = "plain"))+
  theme(plot.tag.position = c(0.05,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(-0.5, 0.1, 0.1, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
pplmTDsum  



pplmCCsum#ply 4
pplmLSsum#ply 3
pplmMGsum#ply 1
pplmPEsum#ply 1 (not significant)
pplmTDsum#ply 1 (significant)

library(cowplot)

temp_poly_newmarg <- plot_grid(
  pplmCCsum, pplmLSsum, pplmMGsum, pplmPEsum, pplmTDsum,
  ncol = 1,
  align = "v",
  axis = "lr"
)

#library(patchwork)
#temp_poly_newmarg<-wrap_plots(list(pplmCCsum,pplmLSsum,pplmMGsum,pplmPEsum,pplmTDsum),ncol=1)


temp_poly_newmarg
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp poly log new margin.pdf",temp_poly_newmarg, width = 4, height = 14,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp poly log new margin.png",temp_poly_newmarg, width = 4, height = 14,dpi = 300)

