## script for making rasters for radius species
## Tropics to Tundra project
## Erica Fischer


# 0 set up workspace --------------------------------------------------------

## load packages
library(fasterize)
library(tidyverse)
library(sf)
library(terra)
library(data.table)
library(raster)


## load in data
final <- read_csv("~/Desktop/T2T_data/full_T2T_database_10_22_24.csv")
final <- subset(final, removed_outlier == 0)
final <- subset(final, in_carribean == 0)

## species list
species_list <- read_csv("~/Desktop/T2T_data/polygon runs/radius_spp.csv")
binomials <- species_list$species
final2 <- subset(final, species %in% binomials)

## load in shapefile for relevant study area
globe <- readRDS("~/Desktop/PhD data/shapefiles/World_Continents.rds")
s_america <- globe[globe$CONTINENT=="South America",]
n_america <- globe[globe$CONTINENT=="North America",]
americas <- raster::union(n_america, s_america)
americas <- raster::crop(americas, extent(c(-200, -28, -65, 85))) # this is bigger than needed on purpose

## make base raster
r <- raster(ncol = 36000, nrow = 18000) #1km x 1km
r <- aggregate(r, fact = 5) #5km x 5km
r[r > 0] = 0
base <- crop(r, americas)

## convert data to simplefeatures
bee_lat_lon <- final2 %>% 
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    agr = "constant",
    crs = 4326,  ##WGS84   
    #crs = 4269,
    stringsAsFactors = FALSE,
    remove = FALSE) 



# 1 make radius rasters --------------------------------------------------------

for (i in 1:length(binomials)){ 
  
  sp <- binomials[i] 
  
  one_bee_data <- bee_lat_lon %>% 
    filter(species == sp)
  
  ## buffer 30km around points
  sp_circles <- one_bee_data %>%
    vect() %>%
    terra::buffer(width = 30000) %>% # 30 kilometers
    st_as_sf()
  
  ## make into a raster
  sp_raster <- fasterize(sp_circles, r)
  
  ## crop to shapefile
  sp_raster <- crop(sp_raster, americas)
  
  sp_raster[sp_raster > 0] = 1
  sp_raster[is.na(sp_raster)] = 0
  
  ## save as geoTIFF
  writeRaster(sp_raster, format = "GTiff", overwrite = TRUE, file = paste(paste("~/Desktop/T2T_data/polygon runs/radius_spp",paste(binomial),sep="/"),"raster.tiff",sep="_"))
  
  
}



# 2 generate pngs if needed --------------------------------------------------------

## read in the ranges
names <- list.files("~/Desktop/T2T_data/polygon runs/radius_spp")
names <- as.data.frame(names)
names <- str_split_fixed(names$names, "_", 2)
names <- as.data.frame(names)
names$V2 <- NULL
colnames(names)[1] ="names"

## base raster in ~1km resolution
r <- raster(ncol = 36000, nrow = 18000) #1km x 1km
r <- aggregate(r, fact = 5) #5km x 5km
r[r > 0] = 0
base <- crop(r, americas)


for (i in c(1:length(names$names))) { 
  
  file <- paste(paste("~/Desktop/T2T_data/polygon runs/radius_spp",paste(names$names[i]), sep = "/"),"_raster.tiff", sep = "")
  range <- raster(file)
  
  range <- mask(range, americas)
  
  png(file=paste(paste("~/Desktop/T2T_data/radius_pngs", paste(names$names[i], ".png"), sep="/")), height = 1600, width = 2000)
  
  ## plot range
  plot(base)
  box = c(-177, -27, -55, 82.8)
  cols_q <- colorRampPalette(c("#ffffff", "#e3ae1b"))
  plot(range, legend = FALSE, col = cols_q(5))
  plot(americas, add = T)
  
  dev.off()
}


