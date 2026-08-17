source("E:/LSU/research/gbark_temperature/gbark_calculation_temperature_w_mass_flow.R") 


library(dplyr)
library(broom)
latemp_dims2 <- latemp_dims %>%
  mutate(
    d_basal_mm  = rowMeans(across(c(basal_diam1, basal_diam2)), na.rm=TRUE),
    d_distal_mm = rowMeans(across(c(distal_diam1, distal_diam2)), na.rm=TRUE),
    r_basal_m   = (d_basal_mm / 1000) / 2,
    r_distal_m  = (d_distal_mm / 1000) / 2,
    L_m         = (total_length / 1000),   # or total_length, depending on what your mass represents
    # Frustum volume (m^3): V = (pi*L/3) * (r1^2 + r1*r2 + r2^2)
    V_m3        = (pi * L_m / 3) * (r_basal_m^2 + r_basal_m*r_distal_m + r_distal_m^2)
  )
endwater_tbl_full <- latemp %>%
  dplyr::group_by(segment) %>%
  dplyr::summarise(
    end_rtime   = max(rtime, na.rm = TRUE),
    m_start_obs = mass_w_parafilm[which.min(rtime)],
    m_end_obs   = mass_w_parafilm[which.max(rtime)],
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    latemp_dims %>% dplyr::select(segment, fresh_mass, dry_mass),
    by = "segment"
  ) %>%
  dplyr::left_join(
    latemp_gs %>% dplyr::select(segment,Tc,RH,treatment),
    by = "segment"
  ) %>%
  dplyr::mutate(
    parafilm_mass_est = m_start_obs - fresh_mass,
    end_mass_sample   = m_end_obs - parafilm_mass_est,
    water_mass_end_g  = end_mass_sample - dry_mass,
    water_mass_init_g = fresh_mass - dry_mass,
    RWC_end           = water_mass_end_g / water_mass_init_g,
    LossFrac_end      = 1 - RWC_end
  )
endmass_tbl <- endwater_tbl_full %>%   # <- use the version before dropped columns
  dplyr::select(segment, end_mass_sample)

sat_check <- latemp_dims %>%
  mutate(
    segment = as.character(segment),
    
    d_basal_mm  = rowMeans(across(c(basal_diam1, basal_diam2)), na.rm=TRUE),
    d_distal_mm = rowMeans(across(c(distal_diam1, distal_diam2)), na.rm=TRUE),
    
    r_basal_cm  = (d_basal_mm / 10) / 2,    # mm -> cm, then /2
    r_distal_cm = (d_distal_mm / 10) / 2,
    L_cm        = total_length / 10,        # mm -> cm
    
    # frustum volume in cm^3
    V_cm3 = (pi * L_cm / 3) * (r_basal_cm^2 + r_basal_cm*r_distal_cm + r_distal_cm^2),
    
    rho_g_cm3 = dry_mass / V_cm3            # dry density in g/cm^3
  ) 

R  <- 8.314462618      # J mol^-1 K^-1
Vw <- 18.01528e-6      # m^3 mol^-1 (partial molal volume of water ~ molar volume)
#T_K <- T_C + 273.15
#library(weathermetrics)
sat_check2 <- sat_check %>%
  left_join(endwater_tbl_full %>% dplyr::select(segment, end_mass_sample,Tc,RH), by="segment") %>%
  mutate(
    # Wolfe/Simpson saturation mass (g); 1.54 g/cm3 is cell-wall density
    sat_mass_g = V_cm3 * (rho_g_cm3 + (1 - rho_g_cm3/1.54)),
    
    RWC_end_sat = (end_mass_sample - dry_mass) / (sat_mass_g - dry_mass),
    deficit_end = 1 - RWC_end_sat,
    CWR_end = deficit_end *(sat_mass_g-dry_mass)*(rho_g_cm3/dry_mass*1000),#unit: kg m-3
    psi_end = -1.3*CWR_end/(490.4-CWR_end), #(MPa) Jupa et al. 2016
    
    RH_instem = 100*exp(psi_end/(((R*(Tc+273.15))/Vw)/1e6)) #Stewart and Broadbridge 1999
  )
 
endwater_tbl_full$RWC_end_sat<-sat_check2$RWC_end_sat
latemp_gs$RWC_end_sat<-sat_check2$RWC_end_sat

summary(sat_check2$V_cm3)
summary(sat_check2$rho_g_cm3)
summary(sat_check2$sat_mass_g)
summary(sat_check2$RWC_end_sat)
summary(sat_check2$deficit_end)
summary(sat_check2$psi_end)
quantile(sat_check2$psi_end,c(.01,.99))

summary(sat_check2$RH_instem)
quantile(sat_check2$RH_instem,c(.01,.99))

sat_check3 <- sat_check2 %>%
  mutate(
    # saturation vapor pressure at Tc (kPa) 
    svp_kPa = 0.1*(10^(10.79574*(1-273.16/(Tc+273.15)) -
                         5.02800*log10((Tc+273.15)/273.16)+
                         1.50475e-4*(1-10^(-8.2969*((Tc+273.15)/273.16-1)))+
                         0.42873e-3*(10^(4.76955*(1-273.16/(Tc+273.15))))+
                         0.78614)),# kPa #used in gbark calculation
    
    ws = svp_kPa/101.3,
    wa = (RH/100)*ws,   # actual air RH column here
    w_instem = (RH_instem/100)*ws,
    
    D_100 = ws - wa,
    D_unsat  = w_instem - wa,
    
    mf_unsat = 1-(w_instem+wa)/2
  )

gbark_unsat<-latemp_gs %>%
  dplyr::select(segment,TreeID,sp,treatment,Tc,E,D,mf,gs)%>%
  distinct()

gbark_unsat2<-gbark_unsat %>%
   left_join(sat_check3%>% dplyr::select(segment,RH_air=RH,RH_instem,D_100,D_unsat,mf_unsat),by="segment" )%>%
  mutate(
    gbark_unsat = E*mf_unsat/D_unsat)

gbark_unsat3<-gbark_unsat2%>%
  mutate(
    
gbark_pct_change=(gbark_unsat-gs)/gs
)

summary(gbark_unsat3$gbark_pct_change)
quantile(gbark_unsat3$gbark_pct_change,c(.01,.99))
quantile(gbark_unsat3$gbark_pct_change,c(.01,.95))

hist(gbark_unsat3$gbark_pct_change * 100,
     main="",xlab = expression(italic(g)[bark]~Percentage~Change~"(%)"))

hist(gbark_unsat3$RH_instem,
     main="",xlab = expression(RH[stem]~"(%)"))

#summary(sat_check3$gbark_pct_change)
#quantile(sat_check3$gbark_pct_change,c(.01,.99))

gbark_unsat3 %>%
  filter(gbark_pct_change > 0.10) %>%
  mutate(delta_RH = RH_instem - RH_air) %>%
  summarise(min_delta_RH = min(delta_RH, na.rm=TRUE),
            median_delta_RH = median(delta_RH, na.rm=TRUE))

temp1 %>%
  group_by(treatment) %>%
  summarise(mean_g=mean(gbark, na.rm=TRUE),
            se_g=sd(gbark, na.rm=TRUE)/sqrt(sum(!is.na(gbark))),
            .groups="drop")
#ggplot(sat_check3, aes(x = treatment, y = gbark_pct_change)) +
#  geom_boxplot() +
#  geom_point()



latemp_slope<-latemp_gs %>%
  dplyr::select(segment,sp,Tc,treatment,TreeID,slope_early,slope_late,slope_ratio_late_early,gbark_delta_late_early, gbark_ratio_late_early,slope_delta_lateminusearly,diam,RWC_end,RWC_end_sat,WC_end_g_per_vol) %>%
  left_join(gbark_unsat3%>%dplyr::select(segment,gbark_pct_change),by="segment" )
  

slopedeltadiff<-lmer(slope_delta_lateminusearly~sp*treatment +(1|TreeID), latemp_slope)
summary(slopedeltadiff)
anova(slopedeltadiff)

emm_slopedelta<-emmeans(slopedeltadiff,~treatment)
pairs(emm_slopedelta, adjust = "tukey")

cld(emm_slopedelta, Letters=letters, adjust="tukey")



gbarkratiodiff  <- lmer(gbark_ratio_late_early ~ sp * treatment + (1|TreeID), latemp_slope)
summary(gbarkratiodiff )
anova(gbarkratiodiff )

emm_gbarkratio  <- emmeans(gbarkratiodiff , ~ treatment)
pairs(emm_gbarkratio , adjust="tukey")
cld(emm_gbarkratio , Letters=letters, adjust="tukey")

gbarkdeltadiff <- lmer(gbark_delta_late_early ~ sp * treatment + (1|TreeID), latemp_slope)
summary(gbarkdeltadiff)
anova(gbarkdeltadiff)

emm_gbarkdelta <- emmeans(gbarkdeltadiff, ~ treatment)
pairs(emm_gbarkdelta, adjust="tukey")
cld(emm_gbarkdelta, Letters=letters, adjust="tukey")

plot(as.factor(latemp_slope$treatment),as.numeric(latemp_slope$gbark_delta_late_early))
points(as.factor(latemp_slope$treatment),as.numeric(latemp_slope$gbark_delta_late_early))

gbarkdeltadiff_diam <- lmer(gbark_delta_late_early ~treatment+diam  + (1|TreeID), latemp_slope)
anova(gbarkdeltadiff_diam)

#test if gbark_delta differ from 0 in each treatment
gbarkdeltadiffnosp <- lmer(gbark_delta_late_early ~ treatment + (1|TreeID), latemp_slope)
anova(gbarkdeltadiffnosp)

emm_gbarkdeltatest<-emmeans(gbarkdeltadiffnosp,~treatment)
test(emm_gbarkdeltatest,null=0)



##############################
p_gbarkdelta_data<-latemp_gs%>%
  dplyr::select(sp,segment,treatment,gbark_early,gbark_late)%>%
  pivot_longer(
    cols = c(gbark_early, gbark_late),
    names_to = "time",
    values_to = "gbark"
  ) %>%
  mutate(
    time = factor(time, levels = c("gbark_early", "gbark_late"),
                  labels = c("Early", "Late")),
    Treatment = factor(treatment, levels = c("10C", "25C", "35C", "45C", "55C"))
  )
trt_cols <- c(
  "10C" = "#2C7BB6",  # blue
  "25C" = "#00A6CA",  # blue-green
  "35C" = "#F9D057",  # yellow
  "45C" = "#F46D43",  # orange
  "55C" = "#D73027"   # red
)
sp_labs <- c(
  CC = "italic('C. caroliniana')",
  LS = "italic('L. styraciflua')",
  MG = "italic('M. grandiflora')",
  PE = "italic('P. elliottii')",
  TD = "italic('T. distichum')"
)
trt_labs <- c("10C"="10", "25C"="25", "35C"="35", "45C"="45", "55C"="55")

gbark_early_late_ttest<-latemp_gs %>%
  group_by(sp) %>%
  summarise(
    mean_early = mean(gbark_early, na.rm = TRUE),
    mean_late = mean(gbark_late, na.rm = TRUE),
    mean_delta = mean(gbark_late - gbark_early, na.rm = TRUE),
    p_value = t.test(gbark_late, gbark_early, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("italic(P) == ",
                     format.pval(p_value, digits = 1, eps = 0.001)),
    sp_strip = sp_labs[as.character(sp)]
  )

t.test(latemp_gs[latemp_gs$sp=="CC",]$gbark_early,latemp_gs[latemp_gs$sp=="CC",]$gbark_late,paired = T)
t.test(latemp_gs[latemp_gs$sp=="LS",]$gbark_early,latemp_gs[latemp_gs$sp=="LS",]$gbark_late,paired = T)
t.test(latemp_gs[latemp_gs$sp=="MG",]$gbark_early,latemp_gs[latemp_gs$sp=="MG",]$gbark_late,paired = T)
t.test(latemp_gs[latemp_gs$sp=="PE",]$gbark_early,latemp_gs[latemp_gs$sp=="PE",]$gbark_late,paired = T)
t.test(latemp_gs[latemp_gs$sp=="TD",]$gbark_early,latemp_gs[latemp_gs$sp=="TD",]$gbark_late,paired = T)


p_gbarkdelta_data$sp_strip <- sp_labs[as.character(p_gbarkdelta_data$sp)]

p_gbarkdelta <- ggplot(p_gbarkdelta_data,
                  aes(x = time, y = gbark,
                      group = segment, color = treatment)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(
    data = gbark_early_late_ttest,
    aes(x = 1.5, y = Inf, label = p_label),
    inherit.aes = FALSE,
    parse = TRUE,
    vjust = 1.4,
    size = 4
  ) +
  labs(x = NULL,y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  facet_wrap(~sp_strip,ncol=1,labeller = label_parsed)+theme_Publication(base_size = 14) +
  theme(
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "plain")
  )

p_gbarkdelta

ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/gbark late early line.pdf",p_gbarkdelta, width = 4, height = 12,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/gbark late early line.png",p_gbarkdelta, width = 4, height = 12,dpi = 300)

cld_gbark_delta<-c("ab","a","ab","b","a")
cld_gbark_delta_df<-as.data.frame(cld_gbark_delta)
cld_gbark_delta_df$trt_strip<-c("10", "25", "35", "45", "55")

latemp_slope$trt_strip<- trt_labs[as.character(latemp_slope$treatment)]
sp_long<- c(
  CC = "C. caroliniana",
  LS = "L. styraciflua",
  MG = "M. grandiflora",
  PE = "P. elliottii",
  TD = "T. distichum"
)
latemp_slope$sp_long<- sp_long[as.character(latemp_slope$sp)]


p_gbark_delta_box<-ggplot(data=latemp_slope,aes(x=trt_strip,y=gbark_delta_late_early))+
  geom_boxplot(outlier.shape = 1)+
  geom_point(aes(group = interaction(TreeID, treatment),color=sp_long), size = 2,shape=1)+
  labs(x = "Temperature Treatment (°C)", y = expression(paste(Delta,italic(g)[bark]~(mmol~m^-2~s^-1))),color="Species") +

  geom_text(data = cld_gbark_delta_df,
            aes(x = trt_strip, y = 2, label = cld_gbark_delta),
            inherit.aes = FALSE, size = 6)+theme_Publication(base_size = 14)+
  theme(
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "plain"),
        legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "plain"))

p_gbark_delta_box
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/gbark late early box.pdf",p_gbark_delta_box, width = 8, height = 6,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/gbark late early box.png",p_gbark_delta_box, width = 8, height = 6,dpi = 300)


cld_slope_delta<-c("a","a","a","a","b")
cld_slope_delta_df<-as.data.frame(cld_slope_delta)
cld_slope_delta_df$trt_strip<-c("10", "25", "35", "45", "55")

p_slope_delta_box<-ggplot(data=latemp_slope,aes(x=trt_strip,y=slope_delta_lateminusearly))+
  geom_boxplot(outlier.shape = 1)+
  geom_point(aes(group = interaction(TreeID, treatment),color=sp_long), size = 2,shape=1)+
  labs(x = "Temperature Treatment (°C)", y = expression(paste(Delta,Slope~(g~h^-1))),color="Species") +
  
  geom_text(data = cld_slope_delta_df,
            aes(x = trt_strip, y = 0.018, label = cld_slope_delta),
            inherit.aes = FALSE, size = 6)+theme_Publication(base_size = 14)+
  theme(
    axis.title = element_text(face = "plain"),
    panel.border = element_rect(color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "plain"),
    legend.text = element_text(face = "italic"),
    strip.text = element_text(face = "plain"))

p_slope_delta_box
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/slope delta early box.pdf",p_slope_delta_box, width = 8, height = 6,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/slope delta late early box.png",p_slope_delta_box, width = 8, height = 6,dpi = 300)




lm_endRWC_gbarkdelta<-lm(latemp_slope$gbark_delta_late_early~latemp_slope$RWC_end)
anova(lm_endRWC_gbarkdelta)

lm_endRWCsat_gbarkdelta<-lm(gbark_delta_late_early ~ RWC_end_sat,data=latemp_slope)
anova(lm_endRWCsat_gbarkdelta)

cor.test(latemp_slope$gbark_delta_late_early,latemp_slope$RWC_end_sat)

##############################
##############################
##############################
##############################
##############################


##############################
#####analyses without the >0.1
reduce_id<- gbark_unsat3 %>%
  filter(gbark_pct_change>0.1)%>%
  pull(segment) %>%
  unique()

reduce_id

temp1reduce <- temp1 %>%
  filter(!segment %in% reduce_id )
temp1reduce$loggbark<-log(temp1reduce$gbark)
lmm1logreduce<-lmer(loggbark ~ sp * treatment + (1|TreeID), data = temp1reduce)
plot(lmm1logreduce)
lmm1logreduce
summary(lmm1logreduce)
anova(lmm1logreduce)
#sp x trt 0.04 changed to 0.050

CCreduce<-temp1[temp1reduce$sp=="CC",]
LSreduce<-temp1[temp1reduce$sp=="LS",]
MGreduce<-temp1[temp1reduce$sp=="MG",]
PEreduce<-temp1[temp1reduce$sp=="PE",]
TDreduce<-temp1[temp1reduce$sp=="TD",]

library(dplyr)
CCreduce_summary<-CCreduce %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))
LSreduce_summary<-LSreduce %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))



MGreduce_summary<-MGreduce %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))
PEreduce_summary<-PEreduce %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))


TDreduce_summary<-TDreduce %>% 
  group_by(treatment) %>%
  summarise(mean_gbark = mean(gbark),
            se_gbark=sd(gbark)/sqrt(n()),
            mean_Tc = mean(Tc),
            se_Tc = sd(Tc)/sqrt(n()))


#sp1 CC

lmCC1numlmmlogreduce <- lmer(log(gbark) ~ poly(Tc, 1) + (1|TreeID),
                             data = CCreduce, REML = FALSE)
summary(lmCC1numlmmlogreduce)

plmCC2numlmmlogreduce<-lmer(log(gbark)~poly(Tc,2)+(1|TreeID),data = CCreduce,REML = F)
summary(plmCC2numlmmlogreduce)

plmCC3numlmmlogreduce<-lmer(log(gbark)~poly(Tc,3)+(1|TreeID),data = CCreduce,REML = F)
summary(plmCC3numlmmlogreduce)


print(anova(lmCC1numlmmlogreduce,plmCC2numlmmlogreduce,plmCC3numlmmlogreduce))
AICc(lmCC1numlmmlogreduce,plmCC2numlmmlogreduce,plmCC3numlmmlogreduce)

#delta AICc  = 2.05529 

###############
m1 <- lmer(log(gbark) ~ Tc + (1|TreeID), data=CCreduce, REML=FALSE)
m3 <- lmer(log(gbark) ~ poly(Tc,3,raw=TRUE) + (1|TreeID), data=CCreduce, REML=FALSE)

# Compare predicted curves visually
new <- data.frame(Tc = seq(min(CCreduce$Tc), max(CCreduce$Tc), length.out=200))
new$pred1 <- predict(m1, newdata=new, re.form=NA)
new$pred3 <- predict(m3, newdata=new, re.form=NA)
max(abs(new$pred3 - new$pred1))  # max difference on log scale



#sp2 LS

lmLS1numlmmlogreduce <- lmer(log(gbark) ~ poly(Tc, 1) + (1|TreeID),
                             data = LSreduce, REML = FALSE)
summary(lmLS1numlmmlogreduce)

plmLS2numlmmlogreduce<-lmer(log(gbark)~poly(Tc,2)+(1|TreeID),data = LSreduce,REML = F)
summary(plmLS2numlmmlogreduce)

plmLS3numlmmlogreduce<-lmer(log(gbark)~poly(Tc,3)+(1|TreeID),data = LSreduce,REML = F)
summary(plmLS3numlmmlogreduce)


print(anova(lmLS1numlmmlogreduce,plmLS2numlmmlogreduce,plmLS3numlmmlogreduce))
AICc(lmLS1numlmmlogreduce,plmLS2numlmmlogreduce,plmLS3numlmmlogreduce)


#sp3 MG

lmMG1numlmmlogreduce <- lmer(log(gbark) ~ poly(Tc, 1) + (1|TreeID),
                             data = MGreduce, REML = FALSE)
summary(lmMG1numlmmlogreduce)

plmMG2numlmmlogreduce<-lmer(log(gbark)~poly(Tc,2)+(1|TreeID),data = MGreduce,REML = F)
summary(plmMG2numlmmlogreduce)

plmMG3numlmmlogreduce<-lmer(log(gbark)~poly(Tc,3)+(1|TreeID),data = MGreduce,REML = F)
summary(plmMG3numlmmlogreduce)


print(anova(lmMG1numlmmlogreduce,plmMG2numlmmlogreduce,plmMG3numlmmlogreduce))
AICc(lmMG1numlmmlogreduce,plmMG2numlmmlogreduce,plmMG3numlmmlogreduce)

#sp4 PE

lmPE1numlmmlogreduce <- lmer(log(gbark) ~ poly(Tc, 1) + (1|TreeID),
                             data = PEreduce, REML = FALSE)
summary(lmPE1numlmmlogreduce)

plmPE2numlmmlogreduce<-lmer(log(gbark)~poly(Tc,2)+(1|TreeID),data = PEreduce,REML = F)
summary(plmPE2numlmmlogreduce)

plmPE3numlmmlogreduce<-lmer(log(gbark)~poly(Tc,3)+(1|TreeID),data = PEreduce,REML = F)
summary(plmPE3numlmmlogreduce)


print(anova(lmPE1numlmmlogreduce,plmPE2numlmmlogreduce,plmPE3numlmmlogreduce))
AICc(lmPE1numlmmlogreduce,plmPE2numlmmlogreduce,plmPE3numlmmlogreduce)


#sp5 TD

lmTD1numlmmlogreduce <- lmer(log(gbark) ~ poly(Tc, 1) + (1|TreeID),
                             data = TDreduce, REML = FALSE)
summary(lmTD1numlmmlogreduce)

plmTD2numlmmlogreduce<-lmer(log(gbark)~poly(Tc,2)+(1|TreeID),data = TDreduce,REML = F)
summary(plmTD2numlmmlogreduce)

plmTD3numlmmlogreduce<-lmer(log(gbark)~poly(Tc,3)+(1|TreeID),data = TDreduce,REML = F)
summary(plmTD3numlmmlogreduce)


print(anova(lmTD1numlmmlogreduce,plmTD2numlmmlogreduce,plmTD3numlmmlogreduce))
AICc(lmTD1numlmmlogreduce,plmTD2numlmmlogreduce,plmTD3numlmmlogreduce)


###################################################################################################################################################################################
#gbark calculation

latemp_slope_reduce <- latemp_slope %>%
  filter(!segment %in% reduce_id )

ggplot(latemp_slope_reduce,aes(as.factor(treatment),slope_ratio_late_early))+
  geom_boxplot()
 
ggplot(latemp_slope_reduce,aes(as.factor(treatment),gbark_delta_late_early))+
  geom_boxplot()


latemp_slope_reduce$Tc<-as.numeric(latemp_slope_reduce$Tc)
slopediffrednum<-lmer(slope_ratio_late_early~sp*Tc +(1|TreeID), latemp_slope_reduce)
summary(slopediffrednum)
anova(slopediffrednum)
#summary(latemp_gs$Tc)

library(lme4)
library(lmerTest)
sloperatiodiffred<-lmer(slope_ratio_late_early~sp*treatment +(1|TreeID), latemp_slope_reduce)
summary(sloperatiodiffred)
anova(sloperatiodiffred)
library(emmeans)
library(MuMIn)
library(multcomp)
emm_sloperatio_red<-emmeans(sloperatiodiffred,~treatment)
pairs(emm_sloperatio_red, adjust = "tukey")

cld(emm_sloperatio_red, Letters=letters, adjust="tukey")



slopedeltadiffred<-lmer(slope_delta_lateminusearly~sp*treatment +(1|TreeID), latemp_slope_reduce)
summary(slopedeltadiffred)
anova(slopedeltadiffred)
library(emmeans)
library(MuMIn)
library(multcomp)
emm_slopedelta_red<-emmeans(slopedeltadiffred,~treatment)
pairs(emm_slopedelta_red, adjust = "tukey")

cld(emm_slopedelta_red, Letters=letters, adjust="tukey")


gbarkratiodiff_red <- lmer(gbark_ratio_late_early ~ sp * treatment + (1|TreeID), latemp_slope_reduce)
summary(gbarkratiodiff_red)
anova(gbarkratiodiff_red)

emm_gbarkratio_red <- emmeans(gbarkratiodiff_red, ~ treatment)
pairs(emm_gbarkratio_red, adjust="tukey")
cld(emm_gbarkratio_red, Letters=letters, adjust="tukey")

gbarkdeltadiff_red <- lmer(gbark_delta_late_early ~ sp * treatment + (1|TreeID), latemp_slope_reduce)
summary(gbarkdeltadiff_red)
anova(gbarkdeltadiff_red)

emm_gbarkdelta_red <- emmeans(gbarkdeltadiff_red, ~ treatment)
pairs(emm_gbarkdelta_red, adjust="tukey")
cld(emm_gbarkdelta_red, Letters=letters, adjust="tukey")

gbarkdeltadiff_diam_red <- lmer(gbark_delta_late_early ~ treatment+diam + (1|TreeID), latemp_slope_reduce)
anova(gbarkdeltadiff_diam_red)

plot(as.factor(latemp_slope_reduce$treatment),latemp_slope_reduce$gbark_delta_late_early)

#


#check diameter
slopediff_diam<-lmer(slope_ratio_late_early~sp*treatment*diam +(1|TreeID), latemp_gs)
anova(slopediff_diam)
slopediff_diam<-lmer(log(slope_ratio_late_early)~sp*treatment*diam +(1|TreeID), latemp_gs)
slopediff_diam<-lmer(slope_ratio_late_early~diam +(1|TreeID), latemp_gs)

diam_check <- sat_check %>%
  left_join(endwater_tbl_full %>% dplyr::select(segment, end_mass_sample,Tc,RH,RWC_end,RWC_end), by="segment") 

diam_WC<-lmer(RWC_end~diam*sp*treatment+ (1|TreeID),diam_check)
anova(diam_WC)
 #diameter related to rwc_end
summary(lm(RWC_end~diam,data=diam_check))

diamdiff<-lmer(diam~treatment+(1|TreeID),latemp_slope)
anova(diamdiff)


plot(latemp_slope$gbark_delta_late_early~latemp_slope$RWC_end)
