## script for making polygon species rasters
## Tropics to Tundra project
## Erica Fischer


# 0 set up workspace ----------------------------------------------------

## load packages
library(dplyr)
library(sf)
library(ggplot2)
library(stars)
library(raster)
library(tidyr)
library(readr)
library(fasterize)
library(stringr)
library(data.table)

## continent data for checking polygon
globe <- readRDS("~/Documents/grad school stuff/PhD stuff/data tables figures/shapefiles/World_Continents.rds")
s_america <- globe[globe$CONTINENT=="South America",]
n_america <- globe[globe$CONTINENT=="North America",]
americas <- raster::union(n_america, s_america)
americas <- raster::crop(americas, extent(c(-200, -28, -65, 85))) # this is bigger than needed on purpose
countries <- st_read("~/Downloads/world-administrative-boundaries/world-administrative-boundaries.shp")

## read in bee data 
bee_data <- read_csv("~/Desktop/usa checklist/full_T2T_database.csv")
bee_data <- subset(bee_data, bee_data$removed_outlier == 0)
bee_data <- subset(bee_data, bee_data$in_caribbean == 0)

## species list
all_species <- read_csv("~/Desktop/T2T_species_list.csv")
all_species <- subset(all_species, type != "SDM")
all_species <- all_species$species
all_species <- unique(all_species)

done_species <- read_csv("~/Desktop/areas_at_different_alphas.csv")
done_species <- done_species$species

all_species <- setdiff(all_species, done_species)

## base raster in ~5km resolution
r <- raster(ncol = 36000, nrow = 18000) #1km x 1km
r <- aggregate(r, fact = 5) #5km x 5km
r[r > 0] = 0
base <- crop(r, americas)

names <- species_list
all_species <- all_species$missing



# 1 set up data ---------------------------------------------------------

## filter to just relevant species 
bee_data2 <- filter(bee_data, species %in% all_species)

bee_data3 <- subset(bee_data2, decimalLatitude >= -67)
bee_data3 <- subset(bee_data3, decimalLongitude < 0)

bee_data4 <- data.frame(species = bee_data3$species, 
                           decimalLatitude = bee_data3$decimalLatitude, 
                           decimalLongitude = bee_data3$decimalLongitude)

bee_data4 <- na.omit(bee_data4)

## thin and see what's still >3

bee_data4$decimalLatitude <- round(bee_data4$decimalLatitude, 5)
bee_data4$decimalLongitude <- round(bee_data4$decimalLongitude, 5)

bee_data5 <- data.table(species = bee_data4$species,
                     decimalLatitude = bee_data4$decimalLatitude, 
                     decimalLongitude = bee_data4$decimalLongitude)

bee_data5 <- distinct(bee_data4)


count_occurrences <- bee_data4 %>% group_by(species) %>% summarise(count = length(species))

more_than_3 <- subset(count_occurrences, count >= 3)

bee_data4 <- subset(bee_data4, species %in% more_than_3$species)

thinned <- data.table(Longitude = double(),
                      Latitude = double(), 
                      species = character())


for (i in (1:length(more_than_3))){
  
  sp <- subset(bee_data5, species == more_than_3$species[i])
  
  sp$species <- NULL
  sp <- sp[,c(2, 1)]
  
  thin2 = spThin::thin.algorithm(sp, thin.par = 5, reps = 100)
  
  no <- sample(1:100, 1)
  
  thin2 <- thin2[[no]]
  
  thin2$species <- more_than_3$species[i]
  
  thinned <- rbind(thinned, thin2)
  
  thin <- NULL
  thin2 <- NULL
  
}


## edit the species list

thinned_count <- thinned %>% group_by(species) %>% summarise(count = length(species))
more_than_2 <- rbind(more_than_3, other_more_than_3)
more_than_2 <- distinct(more_than_2)

species_list <- as.data.frame(all_species)

other_methods <- subset(species_list, !(species_list %in% more_than_2$species))

## subset main data to other more than 3
more_than_data <- subset(bee_data4, species %in% other_more_than_3$species)

## change column names in thinned
thinned2 <- thinned[,c(3, 2, 1)]
colnames(thinned2)<- c("species", "decimalLatitude", "decimalLongitude")

missing <- subset(bee_data5, species %in% more_than_2$species)

more_than_data <- rbind(missing, thinned2)
more_than_data <- na.omit(more_than_data)

more_than_spp <- unique(more_than_data$species)


# 2 run polygons --------------------------------------------------------

for (i in 1:length(more_than_spp)){
  
  sp = filter(more_than_data, species == paste(more_than_spp[i]))
  sp = dplyr::select(sp, decimalLongitude, decimalLatitude)
  sp = na.omit(sp)
  
  ch <- chull(sp)
  coords <- sp[c(ch, ch[1]), ]  # closed polygon
  
  pols <- SpatialPolygons(list(Polygons(list(Polygon(coords)), 1)), proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  pols <- raster::buffer(pols, width = 0.3)
  
  sp_raster <- rasterize(pols, r)
  sp_raster[sp_raster > 0] = 1 
  
  sp_raster <- crop(sp_raster, americas)
  
  sp_raster[sp_raster > 0] = 1
  sp_raster[is.na(sp_raster)] = 0
  
  sp_raster <- mask(sp_raster, americas)
  
  ## save convex hull as geoTIFF file
  writeRaster(sp_raster, format = "GTiff", overwrite = TRUE, file = paste(paste("~/Desktop/T2T_data/polygon_runs",paste(more_than_spp[i]),sep="/"),"raster.tiff",sep="_"))
    
}



# 3 make png maps for polygon species -------------------------------------------

## read in the ranges
names <- list.files("~/Desktop/T2T_data/polygon_runs")
names <- as.data.frame(names)
names <- str_split_fixed(names$names, "_", 2)
names <- as.data.frame(names)
names$V2 <- NULL
colnames(names)[1] ="names"


## base raster in ~1km resolution
r <- raster(ncol = 36000, nrow = 18000) 
r[r > 0] = 0
base <- crop(r, americas)


## loop to make figures with occurrence dots
for (i in c(1:length(names$names))) { 
  
  file <- paste(paste("~/Desktop/T2T_data/polygon_runs",paste(names$names[i]), sep = "/"),"_raster.tiff", sep = "")
  range <- raster(file)
      
  range <- mask(range, americas)
      
  species_data <- subset(bee_data3, species == names$names[i])
  species_data <- species_data[,7:8]
   
      
  # write new png --change write-to folder as needed
  png(file=paste(paste("~/Desktop/T2T_data/png ranges", paste(names$names[i], ".png"), sep="/")), height = 1600, width = 2000)
      
  # plot range
  plot(base)
  box = c(-177, -27, -55, 82.8)
  cols_q <- colorRampPalette(c("#ffffff", "gold1"))
  plot(range, legend = FALSE, col = cols_q(5))
  plot(countries$geometry, add = T)
      
  points(species_data$decimalLongitude, species_data$decimalLatitude, col = "mediumblue", pch = 20, cex = 0.5)
  
  title(main = paste(names$names[i]), cex.main = 4)
      
  dev.off()
      
}



