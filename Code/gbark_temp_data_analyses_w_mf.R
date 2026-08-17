#with mf
#test with log transformation



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
library(MuMIn)
library(gridExtra)
library(performance)
library(DHARMa)
library(multcomp)

source('E://LSU//research//gbark_temperature//Rtheme//ggplot_theme_Publication-2.R')

temp1<-read.csv("data/gbark corrected w mf/gbark_temp_Goff_1957_corr_w_mf.csv")
temp1$TreeID<-as.factor(temp1$TreeID)
temp1$Round<-as.factor(temp1$Round)
temp1$treatment<-as.factor(temp1$treatment)
temp1spring<-read.csv("data/gbark corrected w mf/gbark_temp_Goff_1957_spring_corr_w_mf.csv")
temp1spring$TreeID<-as.factor(temp1spring$TreeID)
temp1spring$Round<-as.factor(temp1spring$Round)
temp1spring$treatment<-as.factor(temp1spring$treatment)
treetraitspring<-read.csv("tree traits spring.csv")
treetraitspring$sp<-as.factor(treetraitspring$sp)
treetraitspring1<-treetraitspring%>%dplyr::select(TreeID,DBH_cm,Height_m,sp)%>%na.omit()
temp1spring<-temp1spring[!temp1spring$treatment=="55C",]

lmm1<-lmer(gbark ~ sp * treatment + (1|TreeID), data = temp1)
performance::check_model(lmm1)
check_normality(lmm1)
#Warning: Non-normality of residuals detected (p < .001).
check_homogeneity(lmm1)
#Warning: Variances differ between groups (Bartlett Test, p = 0.000).
res <- resid(lmm1)
fit <- fitted(lmm1)

plot(fit, res); abline(h=0, lty=2)
qqnorm(res); qqline(res)

boxplot(res ~ temp1$sp, main="Residuals by species")
boxplot(res ~ temp1$treatment, main="Residuals by treatment")

temp1$loggbark<-log(temp1$gbark)
lmm1log<-lmer(loggbark ~ sp * treatment + (1|TreeID), data = temp1)
plot(lmm1log)
lmm1log
summary(lmm1log)
anova(lmm1log)

emm_sp_by_trt <- emmeans(lmm1log, ~ sp | treatment)         # on log scale
cld_sp_by_trt <- cld(emm_sp_by_trt, Letters=letters, adjust="none")
cld_sp_by_trt

library(influence.ME)
infl <- influence(lmm1log, group = "TreeID")
plot(infl, which = "cook")

# Back-transform EMMs if you want means on original scale
summary(emm_sp_by_trt, type="response")  # geometric means on original scale

joint_tests(emm_sp_by_trt)
test(emm_sp_by_trt, joint = TRUE)   # gives one omnibus test per treatment
pairs(emm_sp_by_trt, adjust = "none")   # LSD-style comparisons

eff <- contrast(emm_sp_by_trt, method = "eff")

gate_correct <- test(eff, joint = TRUE)
gate_correct

library(dplyr)
library(stringr)

# treatments that passed the correct gate
sig_trt <- gate_correct %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>%
  pull(treatment)

cld_df <- as.data.frame(cld_sp_by_trt) %>%
  mutate(.group = str_replace_all(.group, " ", "")) %>%
  filter(treatment %in% sig_trt)

letters_df <- cld_df %>%
  mutate(.group = str_replace_all(.group, " ", ""),
         y = 16.5)
ns_df <- as.data.frame(gate_correct) %>%
  filter(p.value >= 0.05) %>%
  mutate(sp = unique(temp1$sp)[1],
         y = 16.5,
         label = "ns")

sp_labs <- c(
  CC = "italic('C. caroliniana')",
  LS = "italic('L. styraciflua')",
  MG = "italic('M. grandiflora')",
  PE = "italic('P. elliottii')",
  TD = "italic('T. distichum')"
)


p1boxtemp <- ggplot(temp1) +
  geom_boxplot(aes(x = sp, y = gbark)) +
  geom_point(aes(x = sp, y = gbark, group = interaction(TreeID, sp)),
             size = 2) +
  facet_wrap(~treatment) +
  labs(y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))), x = "Species") +
  scale_x_discrete(labels = function(x) parse(text = sp_labs[x])) +
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",axis.text.x = element_text(size = 8,angle = 60, hjust = 1),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(data = letters_df,
            aes(x = sp, y = y, label = .group),
            inherit.aes = FALSE,
            size = 4)#+
 # geom_text(data = ns_df,
            #aes(x = sp, y = y, label = label),
            #inherit.aes = FALSE,
            #size = 4) 
p1boxtemp

ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box loggbark.pdf",p1boxtemp, width = 8, height = 7,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box loggbark.png",p1boxtemp, width = 8, height = 7,dpi = 300)


#gbark faceted by sp

#trt_labs <- c("10C"="10 °C", "25C"="25 °C", "35C"="35 °C", "45C"="45 °C", "55C"="55 °C")

trt_labs <- c("10C"="10", "25C"="25", "35C"="35", "45C"="45", "55C"="55")

#sp_labs_1 <- c(
  #CC = 'C. caroliniana',
  #LS = 'L. styraciflua',
  #MG = 'M. grandiflora',
  #PE = 'P. elliottii',
  #TD = 'T. distichum'
#)

emm_trt_by_sp <- emmeans(lmm1log, ~ treatment | sp)
eff2 <- contrast(emm_trt_by_sp, method = "eff")
gate_trt_by_sp <- test(eff2, joint = TRUE) 
gate_trt_by_sp

cld_trt_by_sp <- cld(emm_trt_by_sp, Letters = letters, adjust = "none") %>%
  as.data.frame() %>%
  mutate(.group = str_replace_all(.group, " ", ""))

sig_sp <- as.data.frame(gate_trt_by_sp) %>%
  filter(p.value < 0.05) %>%
  pull(sp)
letters_trt_df <- cld_trt_by_sp %>%
  filter(sp %in% sig_sp) %>%
  mutate(y = 16.5)

panel_label<-data.frame(label = c("A","B","C","D","E"),sp_strip = c(
  "italic('C. caroliniana')",
   "italic('L. styraciflua')",
  "italic('M. grandiflora')",
  "italic('P. elliottii')",
  "italic('T. distichum')"
))

temp1$sp_strip <- sp_labs[as.character(temp1$sp)]
letters_trt_df <- letters_trt_df %>%
  mutate(sp_strip = sp_labs[as.character(sp)]) %>%   
  distinct(sp_strip, treatment, .group, .keep_all = TRUE)



CC<-temp1[temp1$sp=="CC",]
LS<-temp1[temp1$sp=="LS",]
MG<-temp1[temp1$sp=="MG",]
PE<-temp1[temp1$sp=="PE",]
TD<-temp1[temp1$sp=="TD",]

LSspring<-temp1spring[temp1spring$sp=="LS",]
PEspring<-temp1spring[temp1spring$sp=="PE",]

p1box_trt_by_sp<-ggplot(temp1, aes(x = treatment, y = gbark)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  facet_wrap(~ sp_strip,labeller = label_parsed,ncol=1) +
  scale_x_discrete(labels = trt_labs) +
  labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  geom_text(data = letters_trt_df,
            aes(x = treatment, y = y, label = .group),
            inherit.aes = FALSE, size = 4)+theme_Publication(base_size = 14) +
  theme(strip.text.x = element_text(face = "italic"),
    legend.position = "none",
    axis.text.x = element_text(size = 8, #angle = 45, hjust = 1
                               ),
    axis.title = element_text(face = "plain"),
    panel.border = element_rect(color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(face = "plain")
  )

p1box_trt_CC<-ggplot(CC, aes(x = treatment, y = gbark)) +
  geom_boxplot() +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  scale_x_discrete(labels = trt_labs) +
  #labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Carpinus~caroliniana))),tag = "A")+
  geom_text(data = letters_trt_df[letters_trt_df$sp=="CC",],
            aes(x = treatment, y = 8.5, label = .group),
            inherit.aes = FALSE, size = 7)+theme_Publication(base_size = 14) +
  scale_y_continuous(limits = c(1,9),breaks = c(2,4,6,8))+
  theme(        legend.position = "none",axis.text.y = element_text(size = 16)   )  +
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.08,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
p1box_trt_CC

p1box_trt_LS<-ggplot(LS, aes(x = treatment, y = gbark)) +
  geom_boxplot() +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  scale_x_discrete(labels = trt_labs) +
  #labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Liquidamber~styraciflua))),tag = "B")+
  geom_text(data = letters_trt_df[letters_trt_df$sp=="LS",],
            aes(x = treatment, y = 12, label = .group),
            inherit.aes = FALSE, size = 7)+theme_Publication(base_size = 14) +
  #scale_y_continuous(limits = c(0,10),breaks = c(2,4,6,8))+
  scale_y_continuous(limits = c(1,12.5),breaks = c(2,4,6,8,10))+
  theme(        legend.position = "none",axis.text.y = element_text(size = 16) )  +
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.1,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
p1box_trt_LS

p1box_trt_MG<-ggplot(MG, aes(x = treatment, y = gbark)) +
  geom_boxplot() +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  scale_x_discrete(labels = trt_labs) +
  #labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  labs(y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
       x='',title = expression(paste(italic(Magnolia~grandiflora))),tag = "C")+
  geom_text(data = letters_trt_df[letters_trt_df$sp=="MG",],
            aes(x = treatment, y = 12, label = .group),
            inherit.aes = FALSE, size = 7)+theme_Publication(base_size = 14) +
  #scale_y_continuous(limits = c(0,10),breaks = c(2,4,6,8))+
  scale_y_continuous(limits = c(1,12.5),breaks = c(2,4,6,8,10))+
  theme(        legend.position = "none",axis.text.y = element_text(size = 16)   )  +
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.12,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
p1box_trt_MG

p1box_trt_PE<-ggplot(PE, aes(x = treatment, y = gbark)) +
  geom_boxplot() +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  scale_x_discrete(labels = trt_labs) +
  #labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),#x='Temperature Treatments (°C)',
    y='',x='',title = expression(paste(italic(Pinus~elliottii))),tag = "D")+
  geom_text(data = letters_trt_df[letters_trt_df$sp=="PE",],
            aes(x = treatment, y = 10, label = .group),
            inherit.aes = FALSE, size = 7)+theme_Publication(base_size = 14) +
  #scale_y_continuous(limits = c(0,10),breaks = c(2,4,6,8))+
  scale_y_continuous(limits = c(1,10.5),breaks = c(2,4,6,8))+
  theme(        legend.position = "none" ,axis.text.y = element_text(size = 16) )  +
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.08,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin=margin(-0.5, 0.1, -0.5, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
p1box_trt_PE


p1box_trt_TD<-ggplot(TD, aes(x = treatment, y = gbark)) +
  geom_boxplot() +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  scale_x_discrete(labels = trt_labs) +
  #labs(x = "Air Temperature (°C)", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  labs(#y=expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
  y='',x='Air Temperature (°C)',title = expression(paste(italic(Taxodium~distichum))),tag = "E")+
  geom_text(data = letters_trt_df[letters_trt_df$sp=="TD",],
            aes(x = treatment, y = 18.5, label = .group),
            inherit.aes = FALSE, size = 7)+theme_Publication(base_size = 14) +
  #scale_y_continuous(limits = c(0,10),breaks = c(2,4,6,8))+
 scale_y_continuous(limits = c(1,19.5),breaks = c(4,8,12,16))+
  theme(        legend.position = "none",axis.title.x = element_text(face = "plain"),axis.text = element_text(size = 16)  )  +
  theme(#axis.text.x = element_blank(),
    plot.tag.position = c(0.1,0.90),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin=margin(-0.5, 0.1, 0.1, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))
p1box_trt_TD
library(cowplot)

p1box_trt_by_sp <- plot_grid(
 p1box_trt_CC,p1box_trt_LS,p1box_trt_MG,p1box_trt_PE,p1box_trt_TD,
 ncol = 1,
  align = "v",
  axis = "lr"
)
p1box_trt_by_sp
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box loggbark facet by sp.pdf",p1box_trt_by_sp, width = 4, height = 12,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box loggbark facet by sp.png",p1box_trt_by_sp, width = 4, height = 12,dpi = 300)





#spring
lmm1spring<-lmer(gbark~sp*treatment+(1|TreeID),temp1spring)
plot(lmm1spring)
lmm1springlog<-lmer(log(gbark)~sp*treatment+(1|TreeID),temp1spring)
plot(lmm1springlog)
summary(lmm1springlog)
anova(lmm1springlog)
qqnorm(resid(lmm1springlog))
emm_trt_spring <- emmeans(lmm1springlog, ~ treatment)
joint_tests(emm_trt_spring)
test(emm_trt_spring, joint = TRUE)   # gives one omnibus test per treatment
pairs(emm_trt_spring, adjust = "none")   # LSD-style comparisons
cld_trt_spring <- cld(emm_trt_spring, Letters=letters, adjust="none")

library(dplyr)
library(stringr)

# treatments that passed the correct gate

sp_labs_spring <- c(
  LS = "italic('L. styraciflua')",
  PE = "italic('P. elliottii')"
)
trt_labs_spring <- c(
  "10C" = "10",
  "25C" = "25",
  "35C" = "35",
  "45C" = "45"
)

trt_letters_df_spring <- as.data.frame(cld_trt_spring) %>%
  mutate(.group = str_replace_all(.group, " ", ""),
         y = 16.5)

p1boxtempsp_spring <- ggplot(temp1spring) +
  geom_boxplot(aes(x = sp, y = gbark)) +
  geom_point(aes(x = sp, y = gbark, group = interaction(TreeID, sp)),
             size = 2) +
  #facet_wrap(~treatment) +
  labs(y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))), x = "Species",tag = "A") +
  scale_x_discrete(labels = function(x) parse(text = sp_labs_spring[x])) +
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",axis.text.x = element_text(size = 14),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

p1boxtempsp_spring

p1boxtemptrt_spring <- ggplot(temp1spring) +
  geom_boxplot(aes(x =treatment, y = gbark)) +
  geom_point(aes(x = treatment, y = gbark, group = interaction(TreeID, sp)),
             size = 2) +
geom_text(data = trt_letters_df_spring,
          aes(x = treatment, y = y, label = .group),
          inherit.aes = FALSE,
          size = 6) +

  labs(y = " ", x = "Air temperature (\u00B0C)",tag = "B") +
  scale_x_discrete(labels = trt_labs_spring) +
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",axis.text.x = element_text(size = 14),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p1boxtemptrt_spring
library(patchwork)

p_s1 <- (p1boxtempsp_spring + p1boxtemptrt_spring) + plot_annotation(tag_levels = "A")
p_s1
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box log spring.pdf",p_s1, width = 8, height = 4,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp box log spring.png",p_s1, width = 8, height = 4,dpi = 300)


emm_trt_by_sp_spring <- emmeans(lmm1springlog, ~ treatment | sp)
eff2 <- contrast(emm_trt_by_sp, method = "eff")
gate_trt_by_sp <- test(eff2, joint = TRUE) 
gate_trt_by_sp

cld_trt_by_sp <- cld(emm_trt_by_sp, Letters = letters, adjust = "none") %>%
  as.data.frame() %>%
  mutate(.group = str_replace_all(.group, " ", ""))

sig_sp <- as.data.frame(gate_trt_by_sp) %>%
  filter(p.value < 0.05) %>%
  pull(sp)
letters_trt_df <- cld_trt_by_sp %>%
  filter(sp %in% sig_sp) %>%
  mutate(y = 16.5)


temp1$sp_strip <- sp_labs[as.character(temp1$sp)]
letters_trt_df <- letters_trt_df %>%
  mutate(sp_strip = sp_labs[as.character(sp)]) %>%   
  distinct(sp_strip, treatment, .group, .keep_all = TRUE)
p1box_trt_by_sp<-ggplot(temp1, aes(x = treatment, y = gbark)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(TreeID, treatment)), size = 2) +  # no jitter
  facet_wrap(~ sp_strip,labeller = label_parsed) +
  scale_x_discrete(labels = trt_labs) +
  labs(x = "Air Temperature", y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))) +
  geom_text(data = letters_trt_df,
            aes(x = treatment, y = y, label = .group),
            inherit.aes = FALSE, size = 4)+theme_Publication(base_size = 14) +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "plain")
  )
p1box_trt_by_sp






###polynomial analyses
CC<-temp1[temp1$sp=="CC",]
LS<-temp1[temp1$sp=="LS",]
MG<-temp1[temp1$sp=="MG",]
PE<-temp1[temp1$sp=="PE",]
TD<-temp1[temp1$sp=="TD",]

LSspring<-temp1spring[temp1spring$sp=="LS",]
PEspring<-temp1spring[temp1spring$sp=="PE",]


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
source('E://LSU//research//gbark_temperature//Rtheme//ggplot_theme_Publication-2.R')
#sp1 CC
lmCC1num<-lm(CC$gbark~poly(CC$Tc,1))
summary(lmCC1num)
lmCC1numlog<-lm(log(CC$gbark)~poly(CC$Tc,1))
summary(lmCC1numlog)

plmCC2num<-lm(CC$gbark~poly(CC$Tc,2))
summary(plmCC2num)
plmCC2numlog<-lm(log(CC$gbark)~poly(CC$Tc,2))
summary(plmCC2numlog)

plmCC3num<-lm(CC$gbark~poly(CC$Tc,3))
summary(plmCC3num)
plmCC3numlog<-lm(log(CC$gbark)~poly(CC$Tc,3))
summary(plmCC3numlog)

plmCC4num<-lm(CC$gbark~poly(CC$Tc,4))
summary(plmCC4num)
plmCC4numlog<-lm(log(CC$gbark)~poly(CC$Tc,4))
summary(plmCC4numlog)

print(anova(lmCC1num,plmCC2num,plmCC3num,plmCC4num))
AICc(lmCC1num,plmCC2num,plmCC3num,plmCC4num)

print(anova(lmCC1numlog,plmCC2numlog,plmCC3numlog,plmCC4numlog))
AICc(lmCC1numlog,plmCC2numlog,plmCC3numlog,plmCC4numlog)
print(anova(lmCC1numlog,plmCC2numlog,plmCC3numlog))
AICc(lmCC1numlog,plmCC2numlog,plmCC3numlog)


#ply 4 best - plmCC4num
#select linear
lmCC1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data=CC)
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
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(1, 0.5, 0, 0.2, "cm"))
pplmCCsum  

#sp2 LS
lmLS1num<-lm(LS$gbark~poly(LS$Tc,1))
summary(lmLS1num)
lmLS1numlog<-lm(log(LS$gbark)~poly(LS$Tc,1))
summary(lmLS1numlog)

plmLS2num<-lm(LS$gbark~poly(LS$Tc,2))
summary(plmLS2num)
plmLS2numlog<-lm(log(LS$gbark)~poly(LS$Tc,2))
summary(plmLS2numlog)

plmLS3num<-lm(LS$gbark~poly(LS$Tc,3))
summary(plmLS3num)
plmLS3numlog<-lm(log(LS$gbark)~poly(LS$Tc,3))
summary(plmLS3numlog)

plmLS4num<-lm(LS$gbark~poly(LS$Tc,4))
summary(plmLS4num)
plmLS4numlog<-lm(log(LS$gbark)~poly(LS$Tc,4))
summary(plmLS4numlog)

print(anova(lmLS1num,plmLS2num,plmLS3num,plmLS4num))

AICc(lmLS1num,plmLS2num,plmLS3num,plmLS4num)
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
lmLSspring1num<-lm(LSspring$gbark~poly(LSspring$Tc,1))
summary(lmLSspring1num)
lmLSspring1numlog<-lm(log(LSspring$gbark)~poly(LSspring$Tc,1))
summary(lmLSspring1numlog)

plmLSspring2num<-lm(LSspring$gbark~poly(LSspring$Tc,2))
summary(plmLSspring2num)
plmLSspring2numlog<-lm(log(LSspring$gbark)~poly(LSspring$Tc,2))
summary(plmLSspring2numlog)

plmLSspring3num<-lm(LSspring$gbark~poly(LSspring$Tc,3))
summary(plmLSspring3num)
plmLSspring3numlog<-lm(log(LSspring$gbark)~poly(LSspring$Tc,3))
summary(plmLSspring3numlog)

plmLSspring4num<-lm(LSspring$gbark~poly(LSspring$Tc,4))
summary(plmLSspring4num)
plmLSspring4numlog<-lm(log(LSspring$gbark)~poly(LSspring$Tc,4))
summary(plmLSspring4numlog)

print(anova(lmLSspring1num,plmLSspring2num,plmLSspring3num,plmLSspring4num))

AICc(lmLSspring1num,plmLSspring2num,plmLSspring3num,plmLSspring4num)
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
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  scale_y_continuous(limits = c(0,12),breaks = seq(0,12,by = 4))+
  theme(plot.margin=margin(1, 0.5, 0, 0.2, "cm"))

pplmLSsum  





#sp3 MG
lmMG1num<-lm(MG$gbark~poly(MG$Tc,1))
summary(lmMG1num)
lmMG1numlog<-lm(log(MG$gbark)~poly(MG$Tc,1))
summary(lmMG1numlog)

plmMG2num<-lm(MG$gbark~poly(MG$Tc,2))
summary(plmMG2num)
plmMG2numlog<-lm(log(MG$gbark)~poly(MG$Tc,2))
summary(plmMG2numlog)

plmMG3num<-lm(MG$gbark~poly(MG$Tc,3))
summary(plmMG3num)
plmMG3numlog<-lm(log(MG$gbark)~poly(MG$Tc,3))
summary(plmMG3numlog)

plmMG4num<-lm(MG$gbark~poly(MG$Tc,4))
summary(plmMG4num)
plmMG4numlog<-lm(log(MG$gbark)~poly(MG$Tc,4))
summary(plmMG4numlog)

print(anova(lmMG1num,plmMG2num,plmMG3num,plmMG4num))

AICc(lmMG1num,plmMG2num,plmMG3num,plmMG4num)
#ply 1 best - lmMG1num

print(anova(lmMG1numlog,plmMG2numlog,plmMG3numlog))

AICc(lmMG1numlog,plmMG2numlog,plmMG3numlog)
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
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(1, 0.5, 0, 0.2, "cm"))
pplmMGsum  


#sp4 PE
lmPE1num<-lm(PE$gbark~poly(PE$Tc,1))
summary(lmPE1num)
lmPE1numlog<-lm(log(PE$gbark)~poly(PE$Tc,1))
summary(lmPE1numlog)
#anova(lmPE1numlog)

plmPE2num<-lm(PE$gbark~poly(PE$Tc,2))
summary(plmPE2num)
plmPE2numlog<-lm(log(PE$gbark)~poly(PE$Tc,2))
summary(plmPE2numlog)

plmPE3num<-lm(PE$gbark~poly(PE$Tc,3))
summary(plmPE3num)
plmPE3numlog<-lm(log(PE$gbark)~poly(PE$Tc,3))
summary(plmPE3numlog)

plmPE4num<-lm(PE$gbark~poly(PE$Tc,4))
summary(plmPE4num)
plmPE4numlog<-lm(log(PE$gbark)~poly(PE$Tc,4))
summary(plmPE4numlog)
print(anova(lmPE1num,plmPE2num,plmPE3num,plmPE4num))

AICc(lmPE1num,plmPE2num,plmPE3num,plmPE4num)
#ply 1 best - lmPE4num (not significant)
print(anova(lmPE1numlog,plmPE2numlog,plmPE3numlog))

AICc(lmPE1numlog,plmPE2numlog,plmPE3numlog)
#ply 1 best - lmPE4num (not significant)

lmPE1numlog<-lm(log(gbark)~poly(Tc,1,raw=T),data=PE)
newPE <- data.frame(Tc = seq(5, 60, length.out = 300))
newPE$Tc<-as.numeric(newPE$Tc)
#newPE <- newPE[order(newPE$Tc), ]
predPE <- predict(lmPE1numlog, newdata = newPE,se.fit = T)

newPE$fit  <- exp(predPE$fit)
newPE$lwr  <- exp(predPE$fit - 1.96 * predPE$se.fit)
newPE$upr  <- exp(predPE$fit + 1.96 * predPE$se.fit)


#spring
lmPEspring1num<-lm(PEspring$gbark~poly(PEspring$Tc,1))
summary(lmPEspring1num)
lmPEspring1numlog<-lm(log(PEspring$gbark)~poly(PEspring$Tc,1))
summary(lmPEspring1numlog)

plmPEspring2num<-lm(PEspring$gbark~poly(PEspring$Tc,2))
summary(plmPEspring2num)
plmPEspring2numlog<-lm(log(PEspring$gbark)~poly(PEspring$Tc,2))
summary(plmPEspring2numlog)

plmPEspring3num<-lm(PEspring$gbark~poly(PEspring$Tc,3))
summary(plmPEspring3num)
plmPEspring3numlog<-lm(log(PEspring$gbark)~poly(PEspring$Tc,3))
summary(plmPEspring3numlog)

plmPEspring4num<-lm(PEspring$gbark~poly(PEspring$Tc,4))
summary(plmPEspring4num)
plmPEspring4numlog<-lm(log(PEspring$gbark)~poly(PEspring$Tc,4))
summary(plmPEspring4numlog)

print(anova(lmPEspring1num,plmPEspring2num,plmPEspring3num,plmPEspring4num))

AICc(lmPEspring1num,plmPEspring2num,plmPEspring3num,plmPEspring4num)
#ply 1 best
print(anova(lmPEspring1numlog,plmPEspring2numlog,plmPEspring3numlog))

AICc(lmPEspring1numlog,plmPEspring2numlog,plmPEspring3numlog)
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
  theme(axis.text.x = element_blank(),plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  #scale_y_continuous(limits = c(0,12),breaks = seq(0,12,by = 4))+
  theme(plot.margin=margin(1, 0.5, 0, 0.2, "cm"))

pplmPEsum  



#sp5 TD
lmTD1num<-lm(TD$gbark~poly(TD$Tc,1))
summary(lmTD1num)
lmTD1numlog<-lm(log(TD$gbark)~poly(TD$Tc,1))
summary(lmTD1numlog)

plmTD2num<-lm(TD$gbark~poly(TD$Tc,2))
summary(plmTD2num)
plmTD2numlog<-lm(log(TD$gbark)~poly(TD$Tc,2))
summary(plmTD2numlog)

plmTD3num<-lm(TD$gbark~poly(TD$Tc,3))
summary(plmTD3num)
plmTD3numlog<-lm(log(TD$gbark)~poly(TD$Tc,3))
summary(plmTD3numlog)

plmTD4num<-lm(TD$gbark~poly(TD$Tc,4))
summary(plmTD4num)
plmTD4numlog<-lm(log(TD$gbark)~poly(TD$Tc,4))
summary(plmTD4numlog)

print(anova(lmTD1num,plmTD2num,plmTD3num,plmTD4num))

AICc(lmTD1num,plmTD2num,plmTD3num,plmTD4num)
#ply 1 best - lmTD1num 

print(anova(lmTD1numlog,plmTD2numlog,plmTD3numlog))

AICc(lmTD1numlog,plmTD2numlog,plmTD3numlog)
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
  theme(plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5,60),breaks = seq(0,60,by = 10))+
  theme(plot.margin=margin(1, 0.5, 0.2, 0.2, "cm"))
pplmTDsum  



pplmCCsum#ply 4
pplmLSsum#ply 3
pplmMGsum#ply 1
pplmPEsum#ply 1 (not significant)
pplmTDsum#ply 1 (significant)

#CC: ply 4 (but up to 3rd, so 1, sig)
#LS: summer (sig) ply 3 winter (sig) ply 2
#TD: summer (sig) ply1
#MG: ply 1(sig)
#PE: summer (non-sig) winter (sig) both ply 1


library(patchwork)
temp_poly<-wrap_plots(list(pplmCCsum,pplmLSsum,pplmMGsum,pplmPEsum,pplmTDsum),ncol=1)

temp_poly
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp poly log.pdf",temp_poly, width = 4, height = 15,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/temp poly log.png",temp_poly, width = 4, height = 15,dpi = 300)

#NOTES:
#for the 95% CI in this situation, draws a 95% confidence interval for the mean response of the fitted regression line:
#E(gbark∣T)
#Although some treatment means fall partially outside the 95% confidence band of the fitted regression, 
#this band represents uncertainty in the estimated mean response rather than a test of individual treatment differences. 
#The ANOVA evaluates the global temperature effect across all observations and therefore remains significant 

#bark trait for each trt
summary(lmer(gbark~treatment*bark_thickness_ratio+(1|TreeID),data = CC))
anova(lmer(gbark~treatment*bark_thickness_ratio+(1|TreeID),data = CC))
summary(lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = CC))
anova(lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = CC))

summary(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
anova(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
summary(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
anova(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
lmmLSlog<-lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS)
emtrends(lmmLSlog, ~ treatment, var = "bark_thickness_ratio")

summary(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = MG))
anova(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = MG))
summary(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = MG))
anova(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = MG))

summary(lmer(gbark~treatment*bark_thickness_ratio+(1|TreeID),data = PE))
anova(lmer(gbark~treatment*bark_thickness_ratio+(1|TreeID),data = PE))
summary(lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = PE))
anova(lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = PE))

summary(lmer(gbark~treatment+LD+lenticel_size+bark_thickness_ratio+(1|TreeID),data = TD))
anova(lmer(gbark~treatment+LD+lenticel_size+bark_thickness_ratio+(1|TreeID),data = TD))
summary(lmer(log(gbark)~treatment+LD+lenticel_size+bark_thickness_ratio+(1|TreeID),data = TD))
anova(lmer(log(gbark)~treatment+LD+lenticel_size+bark_thickness_ratio+(1|TreeID),data = TD))

lmmgbarkLD<-lmer(log(gbark)~sp*LD+(1|TreeID),data = temp1)
anova(lmmgbarkLD)

lmmgbarkbtr<-lmer(log(gbark)~sp*bark_thickness_ratio+(1|TreeID),data = temp1)
anova(lmmgbarkbtr)

lmmgbarklentsize<-lmer(log(gbark)~sp*lenticel_size+(1|TreeID),data = temp1)
anova(lmmgbarklentsize)


library(emmeans)

m_LS <- lmer(log(gbark) ~ treatment*(LD + lenticel_size + bark_thickness_ratio) + (1|TreeID),
             data = LS)

# treatment-specific slope for bark_thickness_ratio
tr_btr <- emtrends(m_LS, ~ treatment, var = "bark_thickness_ratio")
tr_btr
test(tr_btr)   


summary(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
lmmLSspring<-lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring)
emtrends(lmmLSspring, ~ treatment, var = "bark_thickness_ratio")
summary(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
anova(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
lmmLSspringlog<-lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring)
emtrends(lmmLSspringlog, ~ treatment, var = "bark_thickness_ratio")





summary(lmer(gbark~treatment*(bark_thickness_ratio)+(1|TreeID),data = PEspring))
lmmPEspring<-lmer(gbark~treatment*bark_thickness_ratio+(1|TreeID),data=PEspring)
emtrends(lmmPEspring, ~treatment,var = "bark_thickness_ratio")
summary(lmer(log(gbark)~treatment*(bark_thickness_ratio)+(1|TreeID),data = PEspring))
anova(lmer(log(gbark)~treatment*(bark_thickness_ratio)+(1|TreeID),data = PEspring))
lmmPEspringlog<-lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data=PEspring)
emtrends(lmmPEspringlog, ~treatment,var = "bark_thickness_ratio")

#gbark vs. bark trait for each trait separately for LS

anova(lmer(log(gbark) ~ treatment * bark_thickness_ratio + (1 | TreeID), data = LS))
anova(lmer(log(gbark) ~ bark_thickness_ratio + (1 | TreeID), data = LS))
anova(lmer(log(gbark) ~ treatment + bark_thickness_ratio + (1 | TreeID), data = LS))

anova(lmer(log(gbark) ~ treatment * LD + (1 | TreeID), data = LS))

anova(lmer(log(gbark) ~ treatment * lenticel_size + (1 | TreeID), data = LS))

anova(lmer(log(gbark) ~ treatment * bark_thickness_ratio + (1 | TreeID), data = LSspring))

anova(lmer(log(gbark) ~ treatment * LD + (1 | TreeID), data = LSspring))

anova(lmer(log(gbark) ~ treatment * lenticel_size + (1 | TreeID), data = LSspring))
   
lmmLSbtrspringlog<-lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = LSspring)
anova(lmmLSbtrspringlog)
emtrends(lmmLSbtrspringlog, ~ treatment,var = "bark_thickness_ratio")



#Compare season effect-------
temp1merge<-read.csv("gbark_temp_Goff_1957_merge_w_mf.csv")
temp1merge$TreeID<-as.factor(temp1merge$TreeID)
temp1merge$Round<-as.factor(temp1merge$Round)
temp1merge$treatment<-as.factor(temp1merge$treatment)
lm6<-lmer(gbark~Season+sp+Tc+(1|TreeID),temp1merge)
summary(lm6)
anova(lm6)
lm6log<-lmer(log(gbark)~Season+sp+Tc+(1|TreeID),temp1merge)
summary(lm6log)
anova(lm6log)

lm6.1<-lmer(gbark~Season*sp*Tc+(1|TreeID),temp1merge)
summary(lm6.1)
anova(lm6.1)
lm6.1log<-lmer(log(gbark)~Season*sp*Tc+(1|TreeID),temp1merge)
summary(lm6.1log)
anova(lm6.1log)

lm6.1logpoly<-lmer(log(gbark)~Season*sp*poly(Tc,2)+(1|TreeID),temp1merge)
summary(lm6.1logpoly)
anova(lm6.1logpoly)

lm6.1logcat<-lmer(log(gbark)~Season*sp*treatment+(1|TreeID),temp1merge)
summary(lm6.1logcat)
anova(lm6.1logcat)
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
#sp effect

anova(lm(gbark ~ Season*sp*Tc, data=temp1merge))

#fig S3

#shape_map <- c("Lenticel density" = 16, "Lenticel size (mm)" = 17, "Bark thickness ratio" = 15)
trt_shape <- c("10C" = 16, "25C" = 17, "35C" = 15, "45C" = 18, "55C" = 3)
trt_cols <- c(
  "10C" = "#2C7BB6",  # blue
  "25C" = "#00A6CA",  # blue-green
  "35C" = "#F9D057",  # yellow
  "45C" = "#F46D43",  # orange
  "55C" = "#D73027"   # red
)

#CC
ylimprimCC<-range(CC$lenticel_size,na.rm = T)
ylimsecCC<-range(CC$bark_thickness_ratio,na.rm = T)
bCC<-diff(ylimprimCC)/diff(ylimsecCC)
aCC<-ylimprimCC[1]-bCC*ylimsecCC[1]

 
barktraitCCopt<-ggplot(CC)+
  geom_point(aes(y=LD,x=gbark),size=2,shape = 16,color = "darkgreen")+
  scale_y_continuous(
    "Lenticel Size (mm) or Bark Thickness Ratio", 
    sec.axis = sec_axis(~ (. - aCC)/bCC, name = "Lenticel Density")
  )+
  geom_point(aes(y=lenticel_size,x=gbark),size=2,shape = 17,color = "darkgreen")+
  geom_point(aes(y=bark_thickness_ratio,x=gbark),size=2,shape = 18,color = "darkgreen")+
  labs(x = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1)))
) +
  theme_Publication(base_size = 14)+
  theme(legend.position = c(0.5,0.5), strip.text = element_text(face = "plain"),,plot.tag.position = c(0.05,1),  panel.border = element_rect(color="black",fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

barktraitCCopt

library(dplyr)
library(tidyr)
CC_long <- CC %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitCC <- ggplot(CC_long, aes(x = value, y = gbark,
                                   color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
   labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",color = "Treatment",
    #tag = "A"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
 scale_color_manual(values = trt_cols, name = "Treatment")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
    strip.background = element_rect(fill = NA,colour = NA),
    strip.text = element_text(face = "plain"),
    strip.text.x.bottom = element_text(
      margin = margin(t = -2, b = 0)
    ),#plot.tag.position = c(0.005,0.98),
    panel.border = element_rect(color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")


barktraitCC



#LS
LS_long <- LS %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )

library(ggh4x)

annotation1<-annotate("text", x = 0.04, y = 9, label = "italic(P)==0.03", parse=TRUE, size = 6)

barktraitLS <- ggplot(LS_long, aes(x = value, y = gbark,
                                   color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",
    color = "Treatment",
   # tag = "B"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),#plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")#+
  #at_panel(annotation1, PANEL == 3)




barktraitLS




#MG
MG_long <- MG %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitMG <- ggplot(MG_long, aes(x = value, y = gbark,
                                   color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    color = "Treatment",
   # tag = "C"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),#plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")


barktraitMG




#PE
PE_long <- PE %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitPE <- ggplot(PE_long, aes(x = value, y = gbark,
                                   color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",color = "Treatment",
  #  tag = "D"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),#plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")


barktraitPE




#TD
TD_long <- TD %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitTD <- ggplot(TD_long, aes(x = value, y = gbark,
                                   color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",color = "Treatment",
   # tag = "E"
  ) +
  guides(shape = "none") + 
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  theme_Publication(base_size = 14) +
  # Source - https://stackoverflow.com/a/78915470
  # Posted by Federica Gazzelloni
  # Retrieved 2026-03-25, License - CC BY-SA 4.0
  
  guides(color = guide_legend("Treatment"),
         shape = guide_legend("Treatment")) +

  theme(#legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),# plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")

barktraitTD

aligned_plots <- cowplot::align_plots(
  barktraitCC, barktraitLS, barktraitMG, barktraitPE, barktraitTD,
  align = "v",
  axis = "lr"
)

p_s3 <- plot_grid(  plotlist = aligned_plots,
                  ncol = 1,
                  labels = c("A", "B", "C", "D", "E"),
                  label_size = 14,
                  label_fontface = "plain",
                  label_x = 0.005,
                  label_y = 0.995,
                  hjust = 0,
                  vjust = 1)

#p_s3 <- wrap_plots(list(barktraitCC , barktraitLS,barktraitMG,barktraitPE,barktraitTD),ncol = 1)
p_s3
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/bark trait stack fig s3.pdf",p_s3, width = 10, height = 16,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/bark trait stack fig s3.png",p_s3, width = 10, height = 16,dpi = 300)




#LSspring
LSspring_long <- LSspring %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )

annotation2<-annotate("text", x = 0.09, y = 10, label = "italic(P)==0.005", parse=TRUE, size = 6)

barktraitLSspring <- ggplot(LSspring_long, aes(x = value, y = gbark,
                                               color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",
    color = "Treatment",
    #tag = "A"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),#plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.2, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")+
  at_panel(annotation2, PANEL == 3)



barktraitLSspring

#PEspring
PEspring_long <- PEspring %>%
  dplyr::select(gbark, treatment, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    treatment = factor(treatment, levels = c("10C","25C","35C","45C","55C")),
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitPEspring <- ggplot(PEspring_long, aes(x = value, y = gbark,
                                               color = treatment, shape = treatment)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",
    color = "Treatment",
    #tag = "B"
  ) +
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  guides(color = guide_legend("Treatment"),
         shape = guide_legend("Treatment")) +
  theme_Publication(base_size = 14) +
  theme(#legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        strip.text.x.bottom = element_text(
          margin = margin(t = -2, b = 0)
        ),
        #plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  theme(plot.margin=margin(0.1, 0.1, 0, 0.2, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"),strip.placement = "outside")


barktraitPEspring

p_s4 <- wrap_plots(list( barktraitLSspring,barktraitPEspring),ncol = 1)
p_s4 <-p_s4+plot_annotation(
  theme = theme(
    plot.margin = margin(5.5, 5.5, 5.5, 40)
  )
) &
  labs(y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))))

library(cowplot)

p_s4 <-  barktraitLSspring/barktraitPEspring

aligned_plots_spring <- cowplot::align_plots(
  barktraitLSspring, barktraitPEspring, 
  align = "v",
  axis = "lr"
)

p_s4 <- plot_grid(  plotlist = aligned_plots_spring,
                    ncol = 1,
                    labels = c("A", "B"),
                    label_size = 14,
                    label_fontface = "plain",
                    label_x = 0.02,
                    label_y = 0.998,
                    hjust = 0,
                    vjust = 1)

p_s4<-ggdraw(p_s4) +
  draw_label(
    expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    x = 0.01, y = 0.5, angle = 90, vjust = 0.5
  )
p_s4
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/bark trait stack spring fig s4.pdf",p_s4, width = 10, height = 8,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/bark trait stack spring fig s4.png",p_s4, width = 10, height = 8,dpi = 300)


#check end relative water content for shrinkage theory

lmmRWC<-lmer(RWC_end~sp*treatment+(1|TreeID),temp1)
plot(lmmRWC)
summary(lmmRWC)
anova(lmmRWC)

lmmRWCsat<-lmer(RWC_end_sat~sp*treatment+(1|TreeID),temp1)
plot(lmmRWCsat)
summary(lmmRWCsat)
anova(lmmRWCsat)

emm_RWC <- emmeans(lmmRWC, ~ treatment | sp)
pairs(emm_RWC, adjust = "tukey")

cld(emm_RWC, Letters=letters, adjust="none")

emm_RWCsat <- emmeans(lmmRWCsat, ~ treatment | sp)
pairs(emm_RWCsat, adjust = "tukey")

cld(emm_RWCsat, Letters=letters, adjust="none")


m0 <- lmer(log(gbark) ~ sp * treatment + (1|TreeID), data = temp1)
m1 <- lmer(log(gbark) ~ sp * treatment + RWC_end + (1|TreeID), data = temp1)
m1sat <- lmer(log(gbark) ~ sp * treatment + RWC_end_sat + (1|TreeID), data = temp1)

anova(m0, m1)
summary(m1)

anova(m0,m1sat)

anova(m1sat)

ggplot(temp1, aes(RWC_end, log(gbark))) +
  geom_point() +
  facet_wrap(~treatment) +
  geom_smooth(method="lm", se=FALSE)

ggplot(temp1, aes(RWC_end, log(gbark))) +
  geom_point(aes(color=treatment)) +
  facet_wrap(~sp) +
  geom_smooth(method="lm", se=FALSE)

RWC_Tc<-ggplot(temp1, aes(Tc, RWC_end)) +
  geom_point(aes(color=gbark)) +
  scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~sp) +
  geom_smooth(method="lm", se=FALSE)

RWC_Tc


p_S5<-ggplot(temp1, aes(Tc, RWC_end)) +
  geom_point() +
  #scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~sp_strip,labeller = label_parsed) +
  labs(y = "Relative Water Content", x = "Air Temperature(°C)") +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "plain")
  )
p_S5


ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWCend Tc figS5.pdf",p_S5, width = 8, height = 8,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWCend Tc figS5.png",p_S5, width = 8, height = 8,dpi = 300)


emm_RWC_tbl <- as.data.frame(emm_RWC) %>%
  rename(
    mean = emmean,
    lwr  = lower.CL,
    upr  = upper.CL
  )

cld_RWC_tbl <- cld(emm_RWC, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  mutate(.group = gsub("\\s+", "", .group)) %>%   # remove spaces
  dplyr::select(sp, treatment, .group)

# Merge letters onto emmeans table
emm_RWC_tbl <- emm_RWC_tbl %>%
  left_join(cld_RWC_tbl, by = c("sp","treatment"))


# Per species, put letters a bit above the max CI
letter_RWC_pos <- emm_RWC_tbl %>%
  group_by(sp) %>%
  summarise(y_top = max(upr, na.rm = TRUE),.groups="drop")


emm_RWC_tbl <- emm_RWC_tbl %>%
  left_join(letter_RWC_pos, by = "sp") %>%
  mutate(label_y = 1.04)   # adjust this offset if needed (RWC is ~0–1)
temp1$sp_strip <- sp_labs[as.character(temp1$sp)]
emm_RWC_tbl$sp_strip<- sp_labs[as.character(emm_RWC_tbl$sp)]
p_s5_trt <- ggplot(temp1, aes(x = treatment, y = RWC_end)) +
  geom_point(color="grey") +
  geom_errorbar(data = emm_RWC_tbl,
                aes(x=treatment,ymin = mean-SE, ymax = mean+SE),
                width = 0.18, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = emm_RWC_tbl,aes(x = treatment,y = mean),
             size = 2, inherit.aes = FALSE) +
  geom_text(data = emm_RWC_tbl,
            aes(x = treatment,y = label_y, label = .group),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  scale_x_discrete(labels = trt_labs) +facet_wrap(~ sp_strip, labeller = label_parsed) +
  labs(
    x = "Temperature treatment (°C)",
    y = "End Relative Water Content") +
  theme_Publication(base_size = 14) +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "plain"))
        

p_s5_trt
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWC_end vs temp trt facet by sp.pdf",p_s5_trt, width = 8, height = 7,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWC_end vs temp trt facet by sp.png",p_s5_trt, width = 8, height = 7,dpi = 300)

emm_RWCsat_tbl <- as.data.frame(emm_RWCsat) %>%
  rename(
    mean = emmean,
    lwr  = lower.CL,
    upr  = upper.CL
  )

cld_RWCsat_tbl <- cld(emm_RWCsat, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  mutate(.group = gsub("\\s+", "", .group)) %>%   # remove spaces
  dplyr::select(sp, treatment, .group)

# Merge letters onto emmeans table
emm_RWCsat_tbl <- emm_RWCsat_tbl %>%
  left_join(cld_RWCsat_tbl, by = c("sp","treatment"))


# Per species, put letters a bit above the max CI
letter_RWCsat_pos <- emm_RWCsat_tbl %>%
  group_by(sp) %>%
  summarise(y_top = max(upr, na.rm = TRUE),.groups="drop")


emm_RWCsat_tbl <- emm_RWCsat_tbl %>%
  left_join(letter_RWCsat_pos, by = "sp") %>%
  mutate(label_y = 1.04)   # adjust this offset if needed (RWC is ~0–1)
temp1$sp_strip <- sp_labs[as.character(temp1$sp)]
emm_RWCsat_tbl$sp_strip<- sp_labs[as.character(emm_RWCsat_tbl$sp)]
p_s5_trt_sat <- ggplot(temp1, aes(x = treatment, y = RWC_end_sat)) +
  geom_point(color="grey") +
  geom_errorbar(data = emm_RWCsat_tbl,
                aes(x=treatment,ymin = mean-SE, ymax = mean+SE),
                width = 0.18, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = emm_RWCsat_tbl,aes(x = treatment,y = mean),
             size = 2, inherit.aes = FALSE) +
  geom_text(data = emm_RWCsat_tbl,
            aes(x = treatment,y = label_y, label = .group),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  scale_x_discrete(labels = trt_labs) +facet_wrap(~ sp_strip, labeller = label_parsed) +
  labs(
    x = "Air Temperature (°C)",
    y = "End Relative Water Content") +
  theme_Publication(base_size = 14) +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title = element_text(face = "plain"),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "plain"))


p_s5_trt_sat
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWC_end_sat vs temp trt facet by sp.pdf",p_s5_trt_sat, width = 8, height = 7,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect poly/RWC_end_sat vs temp trt facet by sp.png",p_s5_trt_sat, width = 8, height = 7,dpi = 300)


  








temp1RWC <- temp1 %>%
  group_by(sp, treatment) %>%
  mutate(RWC_wST = RWC_end - mean(RWC_end, na.rm=TRUE)) %>%
  ungroup()

m1st <- lmer(log(gbark) ~ sp * treatment + RWC_wST + (1|TreeID), data=temp1RWC)
summary(m1st)
anova(m1st)



m1gbark_RWC_end_sat <- lmer(log(gbark) ~ sp * treatment + RWC_wST + (1|TreeID), data=temp1RWC)
summary(m1st)
anova(m1st)



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

ggplot(latemp, aes(x=dry_time, y=mass_w_parafilm))+
  geom_point()+facet_wrap(~segment,scales="free")


latemp2<-read.csv("temperature x gbark dry down_exclude.csv")
colnames(latemp2)<-c("tag","collection_date","round","date","sp",	"treatment"	,"treatment_ID",	"segment"	,"mass_w_parafilm",	"exclude","hour",	"minute")


#make date and time readable
latemp2$time<-paste(latemp2$hour,latemp2$minute,sep=":")
latemp2$date_and_time1<-paste(latemp2$date,latemp2$time,sep=" ")
latemp2$rtime<-as.POSIXct(latemp2$date_and_time1,format="%Y.%m.%d %H:%M")

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

latemp001015<-latemp[(1:910),]
latemp016025<-latemp[(911:1540),]

pSdrydown1<-ggplot(latemp001015, aes(x=dry_time, y=mass_w_parafilm))+
  geom_point(aes(color=as.factor(exclude)))+facet_wrap(~segment,scales="free",ncol=5)+
  scale_color_manual(values = c("black","red3"))+
  labs(x="Dry Time (Hour)",y="Mass (g)", color = "Exclusion")
pSdrydown2<-ggplot(latemp016025, aes(x=dry_time, y=mass_w_parafilm))+
  geom_point(aes(color=as.factor(exclude)))+facet_wrap(~segment,scales="free",ncol=5)+
  scale_color_manual(values = c("black","red3"))+
  labs(x="Dry Time (Hour)",y="Mass (g)", color = "Exclusion")

pSdrydown1
pSdrydown2


ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/temp drydown 1.pdf",pSdrydown1, width = 14, height = 18,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/temp drydown 1.png",pSdrydown1, width = 14, height = 18,dpi = 300)

ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/temp drydown 2.pdf",pSdrydown2, width = 14, height = 12,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp effect drydown/temp drydown 2.png",pSdrydown2, width = 14, height = 12,dpi = 300)


######bark traits vs sp

lmmbtsp<-lmer(bark_thickness~sp+(1|TreeID),temp1)
summary(lmmbtsp)
anova(lmmbtsp)
emm_bt_sp <- emmeans(lmmbtsp, ~ sp )         # on log scale
cld_bt_sp <- cld(emm_bt_sp, Letters=letters, adjust="none")

emm_bt_sp
cld_bt_sp

lmmbtrsp<-lmer(bark_thickness_ratio~sp+(1|TreeID),temp1)
summary(lmmbtrsp)
anova(lmmbtrsp)
emm_btr_sp <- emmeans(lmmbtrsp, ~ sp )         # on log scale
cld_btr_sp <- cld(emm_btr_sp, Letters=letters, adjust="none")

emm_btr_sp
cld_btr_sp

lmmbtrspSpring<-lmer(bark_thickness_ratio~sp+(1|TreeID),temp1spring)
summary(lmmbtrspSpring)
anova(lmmbtrspSpring)

emm_btr_spSpring <- emmeans(lmmbtrspSpring, ~ sp )         # on log scale
emm_btr_spSpring

lmmdiamsp<-lmer(diam~sp+(1|TreeID),temp1)
summary(lmmdiamsp)
anova(lmmdiamsp)

lmmlentsizesp<-lmer(lenticel_size~sp+(1|TreeID),temp1)
summary(lmmlentsizesp)
anova(lmmlentsizesp)
emm_lentsize_sp <- emmeans(lmmlentsizesp, ~ sp )         # on log scale
#cld_lentsize_sp <- cld(emm_lentsize_sp, Letters=letters, adjust="none")

emm_lentsize_sp
#cld_lentsize_sp

lmmlentdensitysp<-lmer(LD~sp+(1|TreeID),temp1)
summary(lmmlentdensitysp)
anova(lmmlentdensitysp)

emm_lentdensity_sp <- emmeans(lmmlentdensitysp, ~ sp )         # on log scale
cld_lentdensity_sp <- cld(emm_lentdensity_sp, Letters=letters, adjust="none")

emm_lentdensity_sp
cld_lentdensity_sp




### covariance between temp and vpd 
cov(temp1$Tc,temp1$vpd)
cor(temp1$Tc,temp1$vpd)
par(pty = "s")

plot(temp1$Tc, temp1$vpd)
abline(lm(vpd ~ Tc, data = temp1), lty = 1)






##=========================================temp1##==============================================
##########discarded
CC_long <- CC %>%
  dplyr::select(gbark, LD, lenticel_size, bark_thickness_ratio) %>%
  pivot_longer(
    cols = c(LD, lenticel_size, bark_thickness_ratio),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(gbark), is.finite(value)) %>%
  mutate(
    trait = factor(trait,
                   levels = c("LD","lenticel_size","bark_thickness_ratio"),
                   labels = c("Lenticel density", "Lenticel size (mm)", "Bark thickness ratio"))
  )


barktraitCC <- ggplot(CC_long, aes(x = value, y = gbark, shape = trait)) +
  geom_point(size = 2, color = "darkgreen") +
  facet_wrap(~ trait, scales = "free_x", nrow = 1) +  # horizontal; free trait x scales
  scale_shape_manual(values = shape_map) +
  labs(
    x = NULL,
    y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    tag = "A"
  ) +
  theme_Publication(base_size = 14) +
  theme(
    
    axis.title = element_text(face = "plain"),
    plot.tag.position = c(0.02, 1.02),
    panel.border = element_rect(color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

barktraitCC