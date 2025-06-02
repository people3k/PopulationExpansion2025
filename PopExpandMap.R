
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(plyr)
library(tidyverse)

##Load the 33 radiocarbon data sets with spatial locations

cyc<- read.csv("data/expansion.csv")

##Subset and create data frame for innovation 1
inovate1<-subset(cyc, Innovation=="1")
inovate1<-cbind(Longitude=inovate1$Longitude, Latitude=inovate1$Latitude)
inovate1<-data.frame(inovate1)

# Convert to sf object and transform to Equal Earth
inovate1points <- st_as_sf(inovate1, coords = c("Longitude", "Latitude"), crs = 4326)
inovate1_equal_area <- st_transform(inovate1points, crs = 8857)

##Innovation 2
inovate2<-subset(cyc, Innovation=="2")
inovate2<-cbind(Longitude=inovate2$Longitude, Latitude=inovate2$Latitude)
inovate2<-data.frame(inovate2)

# Convert to sf object and transform to Equal Earth
inovate2points <- st_as_sf(inovate2, coords = c("Longitude", "Latitude"), crs = 4326)
inovate2_equal_area <- st_transform(inovate2points, crs = 8857)

##Innovation 3
inovate3<-subset(cyc, Innovation=="3")
inovate3<-cbind(Longitude=inovate3$Longitude, Latitude=inovate3$Latitude)
inovate3<-data.frame(inovate3)

# Convert to sf object and transform to Equal Earth
inovate3points <- st_as_sf(inovate3, coords = c("Longitude", "Latitude"), crs = 4326)
inovate3_equal_area <- st_transform(inovate3points, crs = 8857)


##World map
world <- ne_countries(scale = "medium", returnclass = "sf")
# Transform to Equal Earth projection (EPSG:8857)
world_equal_area <- st_transform(world, crs = 8857)

# Ensure it's using WGS84 (EPSG:4326)
#st_crs(world)  # Should return EPSG:4326, which is WGS 84

# Plot using ggplot
maptest<-ggplot(data = world_equal_area) +
  geom_sf(fill = "NA", color = "gray40") +
  #coord_sf(crs = st_crs(4326)) +  # Explicitly use WGS84
  labs(title = "Archaeological Case Studies") +
  theme_bw() +
  geom_sf(data = inovate1_equal_area, color = "red", size = 2) +
  geom_sf(data = inovate2_equal_area, color = "green", size = 2) +
  geom_sf(data = inovate3_equal_area, color = "blue", size = 2) 
maptest

##Create a more publication ready map
library(ggspatial)  # For scale bar and north arrow
library(RColorBrewer)

# Combine your three datasets into one with a new column indicating group
inovate1_equal_area$group <- "Internal Innovation"
inovate2_equal_area$group <- "Mixed Innovation"
inovate3_equal_area$group <- "External Innovation"
all_cases <- rbind(inovate1_equal_area, inovate2_equal_area, inovate3_equal_area)

# Map with improvements
map <- ggplot(data = world_equal_area) +
  geom_sf(fill = "gray95", color = "gray70", size = 0.3) +
  geom_sf(data = all_cases, aes(color = group), size = 2) +
  scale_color_brewer(palette = "Set1", name = "Innovation Rank") +
  labs(
    title = "Archaeological Case Studies"
    #subtitle = "Visualized on an equal-area projection",
    #caption = "Data sources: Your research data"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) 
  #annotation_scale(location = "bl", width_hint = 0.2) 
  #annotation_north_arrow(location = "bl", which_north = "true", 
                       #  style = north_arrow_fancy_orienteering)
map