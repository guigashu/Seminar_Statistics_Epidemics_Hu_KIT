library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(Metrics)

# 1) Data manipulation
data_set <- read_csv("Downloads/National-level data_SINGAPORE_20010101_20221224.csv", 
                     col_types = cols(adm_0_name = col_skip(), 
                                      adm_1_name = col_skip(), adm_2_name = col_skip(), 
                                      full_name = col_skip(), ISO_A0 = col_skip(), 
                                      FAO_GAUL_code = col_skip(), RNE_iso_code = col_skip(), 
                                      IBGE_code = col_skip(), calendar_start_date = col_character(), 
                                      calendar_end_date = col_character(), 
                                      case_definition_standardised = col_skip(), 
                                      S_res = col_skip(), 
                                      UUID = col_skip()))

# Filter data set for only weekly values
data_set <- data_set[data_set$T_res == "Week", ]

data_set <- data_set %>%
  mutate(
    calendar_start_date = ymd(calendar_start_date),
    calendar_end_date   = ymd(calendar_end_date)
  )

typeof(data_set$calendar_start_date)
typeof(data_set$calendar_end_date)

# Ensure your date is in Date format 
data_set$calendar_end_date <- as.Date(data_set$calendar_end_date, format = "%Y%m%d")

# Split data set between 2002-2012 and 2006-2017
year_count <- data_set %>% filter(Year %in% 2001:2022) %>% group_by(Year) %>% summarise(count = n())

data_set <- data_set[data_set$Year %in% 2006:2017, ]

# Plot for 2006-2017
p1 <- ggplot(data_set, aes(x = calendar_end_date, y = dengue_total)) +
  geom_line(size = 1, color = "darkblue") +
  labs(
    title = "Dengue Total in Singapore (2006-2017)",
    subtitle = "Weekly Dengue Case Counts",
    x = "Date",
    y = "Dengue Total Cases"
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
  theme_classic() +
  theme(
    plot.title = element_text(lineheight = 16, face = "bold"),
    axis.title = element_text(lineheight  = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p1)

X_t_data <- data_set$dengue_total

# Mean
mean(X_t_data)

# Variance
var(X_t_data)

# ACF
acf(X_t_data)











