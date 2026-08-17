library(tidyverse)

# summarize from raw data
env_sum <- temp1 %>%
  group_by(treatment) %>%
  summarise(
    Tc_mean  = mean(Tc, na.rm = TRUE),
    Tc_se    = sd(Tc, na.rm = TRUE) / sqrt(sum(!is.na(Tc))),
    VPD_mean = mean(vpd, na.rm = TRUE),
    VPD_se   = sd(vpd, na.rm = TRUE) / sqrt(sum(!is.na(vpd))),
    .groups = "drop"
  )

# long format for raw points
env_raw_long <- temp1 %>%
  dplyr::select(treatment, Tc, vpd) %>%
  pivot_longer(
    cols = c(Tc, vpd),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = recode(
      variable,
      Tc = "Air temperature (°C)",
      vpd = "VPD (kPa)"
    )
  )

# long format for means ± SE
env_sum_long <- env_sum %>%
  pivot_longer(
    cols = c(Tc_mean, VPD_mean, Tc_se, VPD_se),
    names_to = c("variable", ".value"),
    names_pattern = "(Tc|VPD)_(mean|se)"
  ) %>%
  mutate(
    variable = recode(
      variable,
      Tc = "Air temperature (°C)",
      VPD = "VPD (kPa)"
    )
  )


p_env <- ggplot() +
  geom_point(
    data = env_raw_long,
    aes(x = factor(treatment), y = value),
    alpha = 0.35,
    size = 1.8,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_errorbar(
    data = env_sum_long,
    aes(x = factor(treatment), ymin = mean - se, ymax = mean + se),
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_point(
    data = env_sum_long,
    aes(x = factor(treatment), y = mean),
    size = 2.8
  ) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  labs(
    x = "Target temperature treatment (°C)",
    y = NULL
  ) +
  scale_x_discrete(labels = c("10","25","35","45","55"))+

  theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.title.x  = element_text(face = "plain")
    
  )

p_env

env_check<-temp1 %>%
  dplyr::select(treatment,Tc,Round,vpd)
par(pty="s")
plot(temp1$Tc,temp1$vpd)
model<-lm(log(temp1$vpd)~temp1$Tc)
curve(exp(coef(model)[1] + coef(model)[2] * x), add = TRUE, col = "red", lwd = 2)
dev.off()


summary(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
anova(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
summary(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
anova(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS))
lmmLSlog<-lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LS)
emtrends(lmmLSlog, ~ treatment, var = "bark_thickness_ratio")

anova(lmer(log(gbark)~bark_thickness_ratio+(1|TreeID),data = LS))

m_LS_add <- lmer(log(gbark) ~ treatment + bark_thickness_ratio + (1 | TreeID),
  data = LS
   )
anova(m_LS_add, type = 3)

summary(m_LS_add)

isSingular(m_LS_add)




LS_barkthicknessratio<-LS_long[LS_long$trait=="Bark thickness ratio",]
barkthicknessratioLS <- ggplot(LS_barkthicknessratio, aes(x = value, y = gbark
                                   )) +
  geom_point(aes(color = treatment, shape = treatment),size = 2) +
  #facet_wrap(~ trait, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",
    color = "Treatment",
    #tag = "B"
  ) +geom_smooth(aes(x = value, y = gbark),method = "lm", color="black",size=1,se=T,alpha = 0.2)+
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) #+
  #theme(plot.margin=margin(0, 0.1, 0, 0.1, "cm"),plot.title = element_text(vjust = -4), panel.spacing = unit(0, "cm"))+
  #at_panel(annotation1, PANEL == 3)




barkthicknessratioLS


summary(lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
lmmLSspring<-lmer(gbark~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring)
emtrends(lmmLSspring, ~ treatment, var = "bark_thickness_ratio")
summary(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
anova(lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring))
lmmLSspringlog<-lmer(log(gbark)~treatment*(LD+lenticel_size+bark_thickness_ratio)+(1|TreeID),data = LSspring)
emtrends(lmmLSspringlog, ~ treatment, var = "bark_thickness_ratio")

anova(lmer(log(gbark)~bark_thickness_ratio+(1|TreeID),data = LSspring))
summary(lmer(log(gbark)~treatment*bark_thickness_ratio+(1|TreeID),data = LSspring))

LSspring_barkthicknessratio<-LSspring_long[LSspring_long$trait=="Bark thickness ratio",]

barkthicknessratioLSspring <- ggplot(LSspring_barkthicknessratio, aes(x = value, y = gbark
,color = treatment, shape = treatment)) +
  geom_point(size = 2.5) +
  #facet_wrap(~ trait, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    #y = expression(paste(italic(g)[bark]~(mmol~m^-2~s^-1))),
    y="",
    color = "Treatment",
    #tag = "B"
  ) +geom_smooth(aes(x = value, y = gbark,color=treatment),method = "lm",size=1,se=F)+
  scale_shape_manual(values = trt_shape,name = "Treatment") +
  scale_color_manual(values = trt_cols, name = "Treatment")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme_Publication(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA,colour = NA),
        strip.text = element_text(face = "plain"),
        plot.tag.position = c(0.05,0.98),
        panel.border = element_rect(color="black", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) #+
barkthicknessratioLSspring
