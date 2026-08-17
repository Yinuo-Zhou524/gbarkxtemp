setwd("E://LSU//research//gbark_temperature/")


# Packages
library(dplyr)
library(lubridate)
library(ggplot2)

# Read data
dat <- read.csv("data/Baton Rouge weather data.csv", stringsAsFactors = FALSE)


# Fahrenheit to Celsius
f_to_c <- function(x) (x - 32) * 5/9

dat <- dat %>%
  mutate(
    date  = mdy(date),
    year  = year(date),
    month = month(date),
    Temp_mean_C = f_to_c(Temp_mean)
  )

# Function to convert p-value to stars
p_to_star <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

# Month-wise t-tests: 2023 vs 2000–2022
sig_df <- dat %>%
  filter(year >= 2000) %>%
  group_by(month) %>%
  summarise(
    p_value = tryCatch(
      t.test(
        Temp_mean_C[year == 2023],
        Temp_mean_C[year >= 2000 & year <= 2022]
      ,alternative = "greater")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    star = sapply(p_value, p_to_star)
  )



# Merge stars with 2023 data
clim_23yr <- dat %>%
  filter(year >= 2000, year <= 2022) %>%
  group_by(month) %>%
  summarise(
    mean_C = mean(Temp_mean_C, na.rm = TRUE),
    se_C   = sd(Temp_mean_C, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(series = "2000–2022 mean")

y2023 <- dat %>%
  filter(year == 2023) %>%
  group_by(month) %>%
  summarise(mean_C = mean(Temp_mean_C, na.rm = TRUE), .groups = "drop")

# Merge stars with 2023 data
y2023 <- left_join(y2023, sig_df, by = "month")

y2023 <- y2023 %>%
  mutate(series = "2023")

temp_comparison<-ggplot() +
  geom_ribbon(
    data = clim_23yr,
    aes(x = month, ymin = mean_C - se_C, ymax = mean_C + se_C),
    fill = "grey70", alpha = 1
  ) +
  geom_line(data = clim_23yr,
            aes(x = month, y = mean_C, color = series),
            linewidth = 1) +
  #geom_point(data = clim_23yr,
   #          aes(x = month, y = mean_C, color = series),
    #         size = 2) +
  
  geom_line(data = y2023,
            aes(x = month, y = mean_C, color = series),
            linewidth = 1) +
  #geom_point(data = y2023,
   #          aes(x = month, y = mean_C, color = series),
    #         size = 2) +
  
  # significance stars (do NOT map color here)
  geom_text(data = y2023,
            aes(x = month,
                y = mean_C + 0.6,
                label = star),
            color = "red",
            size = 5) +
  
  scale_color_manual(
    name = "Series",
    values = c("2000–2022 mean" = "black",
               "2023" = "red")
  ) +
  
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    x = "Month",
    y = "Average temperature (°C)",
  ) +
  labs(
    x = "Month",
    y = "Average temperature (°C)",
  ) +
  theme_Publication(base_size = 14)+
  theme(axis.title.x = element_text(face = "plain"), 
        panel.border = element_rect(color="black",fill=NA),
        panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(
    legend.position = c(0.5, 0.15),
    legend.title = element_blank(),
    legend.background = element_blank()
  )

temp_comparison

ggsave("E:/LSU/research/gbark_temperature/figure/temp comparison.pdf",temp_comparison, width = 5, height = 4,dpi = 300)
ggsave("E:/LSU/research/gbark_temperature/figure/temp comparison.png",temp_comparison, width = 5, height = 4,dpi = 300)


