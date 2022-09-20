# TX precipitation Map

# Code by: Joel F. Swift 09/19/2022

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(ggpubr); packageVersion('ggpubr')
library(urbnmapr)

# Theme set 
theme_set(theme_pubr())

### 30 year average (1990-2021) of percipitation for KS counties
county_precip_TX <- read.csv('1971_2000_texas_counties_precip_normal.csv', sep = ',', header = TRUE)
# Make a simpler data frame for plotting
county_precip_TX <- data.frame('county_fips' = county_precip_TX$FIP ,'AVG_30_years'= county_precip_TX$Normal)
# Fips # to character to play nice with plotting function
county_precip_TX$county_fips <- as.character(county_precip_TX$county_fips)
# Sanity check
county_precip_TX[1:5,]
summary(county_precip_TX)


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


map_to_save <- county_precip_TX %>% 
  left_join(counties, by = "county_fips") %>% 
  filter(state_name =="Texas") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = AVG_30_years)) +
  geom_polygon(color = "#000000", size = .25) +
  scale_fill_gradientn( guide = guide_colorbar(title.position = "top", frame.colour = "black", ticks.colour = "black"), colours = colfunc(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  theme(legend.title = element_text(),
        legend.key.width = unit(.5, "in")) +
  labs(fill = "Normal Annual Preciptation (in.) \n 1971-2000")

ggsave("TX_1971-2000_percip_map.svg",map_to_save, height = 8, width = 12)