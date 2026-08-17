setwd("E://LSU//research//gbark_temperature/")

# Packages
library(dplyr)
library(lubridate)
library(ggplot2)

# Read data
dat <- read.csv("data/Baton Rouge weather data 2000-2024.csv", stringsAsFactors = FALSE)
dat$Temp_mean<-as.numeric(dat$Temp_mean)
# Fahrenheit to Celsius
f_to_c <- function(x) (x - 32) * 5/9

dat <- dat %>%
  mutate(
    date  = mdy(date),
    year  = year(date),
    month = month(date),
    Temp_mean_C = f_to_c(Temp_mean)
  )

# ---- Reference climatology: 2000–2024 monthly mean ± SD (from daily values) ----
clim_2000_2024 <- dat %>%
  filter(year >= 2000, year <= 2024) %>%
  group_by(month) %>%
  summarise(
    mean_C = mean(Temp_mean_C, na.rm = TRUE),
    sd_C   = sd(Temp_mean_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(series = "2000–2024 mean")

# ---- 2023 monthly mean ----
y2023 <- dat %>%
  filter(year == 2023) %>%
  group_by(month) %>%
  summarise(mean_C = mean(Temp_mean_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(series = "2023")

# ---- 2024 monthly mean ----
y2024 <- dat %>%
  filter(year == 2024) %>%
  group_by(month) %>%
  summarise(mean_C = mean(Temp_mean_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(series = "2024")

library(gridExtra)
source('E://LSU//research//gbark_temperature//Rtheme//ggplot_theme_Publication-2.R')

latemp1<-read.csv("temperature x gbark dry down.csv")
collectiondate2023<-latemp1%>%
  dplyr::select("segment","sp","collection_date")%>%
  distinct(segment, collection_date, sp) %>%   # remove duplicates
  mutate(collection_date=ymd(collection_date),
         month=month(collection_date),
         day=day(collection_date),
         x_pos = month + (day - 1) /31)# convert to continuous month scale
latempspring1<-read.csv("temperature x gbark dry down 202403.csv")
collectiondate2024<-latempspring1%>%
  dplyr::select("segment","sp","collection_date")%>%
  distinct(segment, collection_date, sp) %>%   # remove duplicates
  mutate(collection_date=ymd(collection_date),
         month=month(collection_date),
         day=day(collection_date),
         x_pos = month + (day - 1) /31)# convert to continuous month scale


# ---- Plot ----
temp_comparison <- ggplot() +
  geom_ribbon(
    data = clim_2000_2024,
    aes(x = month, ymin = mean_C - sd_C, ymax = mean_C + sd_C),
    fill = "grey70", alpha = 1
  ) +
  geom_line(
    data = clim_2000_2024,
    aes(x = month, y = mean_C, color = series),
    linewidth = 1
  ) +
  geom_line(
    data = y2023,
    aes(x = month, y = mean_C, color = series),
    linewidth = 1
  ) +
  geom_line(
    data = y2024,
    aes(x = month, y = mean_C, color = series),
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Series",
    values = c(
      "2000–2024 mean" = "black",
      "2023" = "red3",
      "2024" = "cornflowerblue"
    )
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    x = "Month",
    y = "Average temperature (°C)"
  ) +
  geom_point(data = collectiondate2023,aes(x=x_pos,y=min(clim_2000_2024$mean_C)-5),shape=124,size=6,color="red3")+
  geom_point(data = collectiondate2024,aes(x=x_pos,y=min(clim_2000_2024$mean_C)-5),shape=124,size=6,color="cornflowerblue")+
  theme_Publication(base_size = 14) +
  theme(
    axis.title = element_text(face = "plain"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.55, 0.15),
    legend.title = element_blank(),
    legend.background = element_blank()
  )

temp_comparison

ggsave("E:/LSU/research/gbark_temperature/figure/temp comparison.pdf",
       temp_comparison, width = 5, height = 4, dpi = 300)

ggsave("E:/LSU/research/gbark_temperature/figure/temp comparison.png",
       temp_comparison, width = 5, height = 4, dpi = 300)

summer_ambient <- dat %>%
  filter(year == 2023, month == 06) %>%
  summarise(
    mean_C = mean(Temp_mean_C, na.rm = TRUE),
    sd_C   = sd(Temp_mean_C, na.rm = TRUE),
  )

winter_ambient <- dat %>%
  filter(year == 2024, month == 03) %>%
  summarise(
    mean_C = mean(Temp_mean_C, na.rm = TRUE),
    sd_C   = sd(Temp_mean_C, na.rm = TRUE),
  )
