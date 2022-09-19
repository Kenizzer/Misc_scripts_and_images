# KS precipitation Map

# Code by: Joel F. Swift 6/23/2022

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(ggpubr); packageVersion('ggpubr')
library(urbnmapr)

# Theme set 
theme_set(theme_pubr())

### 30 year average (1990-2021) of percipitation for KS counties
county_precip_KS <- read.csv('KS_precipitation1895to2021.csv', sep = ',', header = TRUE)
which(colnames(county_precip_KS)=="X1990") # index 99
which(colnames(county_precip_KS)=="X2021") # index 130
# Average the rainfall annually for 1990-2021
county_precip_KS$AVG_30_years <- rowMeans(county_precip_KS[,99:130])
# Make a simpler data frame for plotting
county_fip_30yrAVG <- data.frame('county_fips' = county_precip_KS$FIP_code ,'AVG_30_years'= county_precip_KS$AVG_30_years)
# Fips # to character to play nice with plotting function
county_fip_30yrAVG$county_fips <- as.character(county_fip_30yrAVG$county_fips)
# Sanity check
county_fip_30yrAVG[1:5,]
summary(county_fip_30yrAVG)


# Color palette
# Taken from the colors in https://climate.k-state.edu/basics/
colfunc<-colorRampPalette(c("#fffe7a",
                            "#bcf856",
                            "#80ed38",
                            "#32e007",
                            "#33c74c",
                            "#3baa7d",
                            "#2393aa",
                            "#27639d",
                            "#1b388d",
                            "#130f78"))


map_to_save <- county_fip_30yrAVG %>% 
  left_join(counties, by = "county_fips") %>% 
  filter(state_name =="Kansas") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = AVG_30_years)) +
  geom_polygon(color = "#000000", size = .25) +
  scale_fill_gradientn( guide = guide_colorbar(title.position = "top", frame.colour = "black", ticks.colour = "black"), colours = colfunc(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  theme(legend.title = element_text(),
        legend.key.width = unit(.5, "in")) +
  labs(fill = "Normal Annual Preciptation (in.) \n 1990-2021")

ggsave("KS_1990-2021_percip_map.svg",map_to_save, height = 8, width = 12)