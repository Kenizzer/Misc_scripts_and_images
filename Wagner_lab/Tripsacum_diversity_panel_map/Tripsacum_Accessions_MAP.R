#Tripsacum Accession Map
#10/13/2022
#Code by: Joel Swift

# Depens
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggpubr)

# Theme and data import
# Select plot theme globally
theme_set(theme_pubr())
# Load base map data for USA states
states <- map_data("state")
# Load tripsacum dactyloides data with GPS points (Decimal degrees)
Trip_gps_points <- read.csv("Tripsacum_Accessions_w_GPS.tsv", sep = "\t", header = TRUE)

# full USA map
ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = 'white', color = "black") + 
  coord_fixed(1.3) +
  guides(fill="none") + 
  geom_point(data = Trip_gps_points, aes(x = Long, y = Lat), color = "black", size = 4) +
  geom_point(data = Trip_gps_points, aes(x = Long, y = Lat), color = "green", size = 3)

# Select states where we have collections (or border collection states)
select_states <- subset(states, region %in% c("texas", "oklahoma", "arkansas", "missouri", "illinois",
                                              "indiana", "kansas", "mississippi", "georgia", "alabama",
                                              "south carolina", "north carolina", "tennessee", "kentucky",
                                              "louisiana", "florida"))

Plot_to_save <- ggplot(data = select_states) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = 'white', color = "black") + 
  coord_fixed(1.3) +
  guides(fill="none") + 
  geom_point(data = Trip_gps_points, aes(x = Long, y = Lat), pch = 21, color = "black", fill = "firebrick1", size = 4, alpha = 0.85, stroke=1.1) +
  ggsn::scalebar(select_states, dist = 200, dist_unit = 'km', st.size=3, height=0.05, transform  = TRUE, model = 'WGS84') +
  ggsn::north(select_states, symbol = 10) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle(expression(italic("Tripsacum dactyloides")~'accessions in common garden'))

ggsave("select_states_tripsacum_map.svg", Plot_to_save, height = 6, width = 6)
ggsave("select_states_tripsacum_map.png", Plot_to_save, height = 6, width = 6)