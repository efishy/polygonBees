## script for making alpha hull rasters 
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
library(alphahull)
library(sparch)
library(spThin)
library(dismo)

## read in dataset
bee_data <- read_csv("~/Desktop/T2T_data/full_T2T_database.csv")
bee_data <- subset(bee_data, bee_data$removed_outlier == 0)
bee_data <- subset(bee_data, bee_data$in_carribean == 0)

## read in full species list and simplify
species <- read_csv("~/Downloads/Full Taxonomy T2T species 6 23 25.csv")
species <- species %>% group_by(species) %>% summarise(count = length(species))

## table for doing different species at different alpha values
types <- read_csv("~/Desktop/updated_area_centroid.csv")
alphas <- subset(types, type == "alpha hull")

alpha_level <- read_csv("~/Desktop/species list inc.csv")
level14 <- subset(alpha_level, alpha_level$`Alpha level` == 14)
level14 <- level14$species

species <- subset(species, species %in% level14)
species_final <- species$species

## read in/make relevant rasters and shapefiles
globe <- readRDS("~/Desktop/PhD data/shapefiles/World_Continents.rds")
s_america <- globe[globe$CONTINENT=="South America",]
n_america <- globe[globe$CONTINENT=="North America",]
americas <- raster::union(n_america, s_america)
americas <- raster::crop(americas, extent(c(-200, -28, -65, 85))) # this is bigger than needed on purpose

countries <- st_read("~/Downloads/world-administrative-boundaries/world-administrative-boundaries.shp")

r <- raster(ncol = 36000, nrow = 18000) #1km x 1km
r <- aggregate(r, fact = 5) #5km x 5km
r[r > 0] = 0
r2 <- r
base <- crop(r, americas)
base2 <- base

st_as_raster <- function(rstars){
  rext <- st_bbox(rstars)
  raster(t(rstars[[1]]), xmn = rext[1], xmx = rext[3],
         ymn = rext[2], ymx=rext[4],
         crs = st_crs(rstars)$proj4string)
}



# 1 data processing -------------------------------------------------------

## count occurrences per species, remove duplicates
final2 <- bee_data

final2$decimalLatitude <- round(final2$decimalLatitude, 5)
final2$decimalLongitude <- round(final2$decimalLongitude, 5)

final3 <- data.table(species = final2$species,
                     decimalLatitude = final2$decimalLatitude, 
                     decimalLongitude = final2$decimalLongitude)

final3 <- distinct(final3)
final3 <- subset(final3, species %in% listy)

count_occurrences <- final3 %>% group_by(species) %>% summarise(count = length(species))

more_than_5 <- subset(count_occurrences, count > 5)

final4 <- subset(final3, species %in% more_than_5$species)

thinned <- data.table(Longitude = double(),
                      Latitude = double(), 
                      species = character())

## spatially thin records for species with enough records to do so
for (i in (1:length(more_than_5$species))){
  
    sp <- subset(final4, species == more_than_5$species[i])
    
    sp$species <- NULL
    sp <- sp[,c(2, 1)]
    
    thin2 = spThin::thin.algorithm(sp, thin.par = 5, reps = 100)
    
    no <- sample(1:100, 1)
    
    thin2 <- thin2[[no]]

    thin2$species <- more_than_5$species[i]

    thinned <- rbind(thinned, thin2)

    thin <- NULL
    thin2 <- NULL
}


## make sure still enough records after thinning to generate hull
thinned_count <- thinned %>% group_by(species) %>% summarise(count = length(species))
more_than_3 <- subset(thinned_count, count >= 3)

other_more_than_3 <- subset(count_occurrences, count >= 3 & count <= 5)
more_than_2 <- rbind(more_than_3, other_more_than_3)

## subset main data to other more than 3
more_than_data <- subset(final3, species %in% other_more_than_3$species)


## change column names in thinned
thinned2 <- thinned[,c(3, 2, 1)]

colnames(thinned2)<- c("species", "decimalLatitude", "decimalLongitude")

more_than_data <- rbind(thinned2, more_than_data)
more_than_data <- distinct(more_than_data)

#more_than_2 <- unique(final3$species)



# 2 generate hulls (and optional pngs) ------------------------------------

problem_species <- data.frame(species = character())

for(i in c(1:length(more_than_2))){
  
  sp <- more_than_2$species[i]
  species_data <- subset(more_than_data, species == sp)
  
  sp_hull <- ahull(species_data$decimalLongitude, species_data$decimalLatitude, alpha = 9)
  
  sp_poly <- ah2sp(sp_hull)
  
  if(is.null(sp_poly)){
    
    problem_species <- rbind(problem_species, sp)
    
  }else{
    
    sp_poly2 <- st_as_sf(sp_poly)
    
    sp_poly2 <- sp_poly2 %>% st_buffer(dist = 0.3)
    
    sp_raster <- fasterize(sp_poly2, r)
    
    sp_raster <- crop(sp_raster, americas)
    
    sp_raster <- mask(sp_raster, americas)
    
    sp_raster[sp_raster > 0] = 1
    sp_raster[is.na(sp_raster)] = 0
    
    ## save as geoTIFF
    writeRaster(sp_raster, format = "GTiff", overwrite = TRUE, file = paste(paste("~/Desktop/usa checklist/new ranges/alpha hull/alpha 9",paste(more_than_2$species[i]),sep="/"),"raster.tiff",sep="_"))
    
    ## make png
    # sp_raster <- mask(sp_raster, americas)
    # 
    # png(file=paste(paste("~/Desktop/alpha pngs", paste(more_than_2[i], ".png"), sep="/")), height = 1600, width = 2000)
    # 
    # # plot range
    # plot(base)
    # box = c(-177, -27, -55, 82.8)
    # cols_q <- colorRampPalette(c("#ffffff", "#e3ae1b"))
    # plot(sp_raster, legend = FALSE, col = cols_q(5))
    # plot(americas, add = T)
    # 
    # # add occurrence points to map
    # points(species_data$decimalLongitude, species_data$decimalLatitude, col = "blue3", pch = 20, cex = 0.5)
    # 
    # title(main = paste(more_than_2[i]), cex.main = 4)
    # 
    # 
    # dev.off()
    
    base <- base2
    
  }
  
}

colnames(problem_species) <- "species"

write_csv(problem_species, "~/Desktop/alpha hull/alpha 14/impossible_hull_spp.csv")




# 3 hulls for spp with few points -----------------------------------------

## 5 species with 4 points
four_points <- subset(count_occurrences, count == 4)
rows <- sample.int(202, 5)

four_names <- four_points[rows, ]
four_data <- subset(final4, species %in% four_points$species)


for(i in c(1:5)){
  
  sp <- four_names$species[i]
  species_data <- subset(four_data, species == sp)
  
  sp_hull <- ahull(species_data$decimalLongitude, species_data$decimalLatitude, alpha = 15)
  
  sp_poly <- ahull2poly(sp_hull)
  
  sp_poly2 <- st_as_sf(sp_poly)
  
  sp_raster <- fasterize(sp_poly2, r)
  
  sp_raster <- crop(sp_raster, americas)
  
  sp_raster[sp_raster > 0] = 1
  sp_raster[is.na(sp_raster)] = 0
  
  ## save as geoTIFF
  writeRaster(sp_raster, format = "GTiff", overwrite = TRUE, file = paste(paste("~/Desktop/T2T_data/alpha hull test/tiffs",paste(four_points$species[i]),sep="/"),"raster.tiff",sep="_"))
  
  
  # ## make png
  # sp_raster <- mask(sp_raster, americas)
  # 
  # png(file=paste(paste("~/Desktop/T2T_data/alpha hull test/pngs", paste(nine_names$species[i], ".png"), sep="/")), height = 1600, width = 2000)
  # 
  # ## plot range
  # plot(base)
  # box = c(-177, -27, -55, 82.8)
  # cols_q <- colorRampPalette(c("#ffffff", "#e3ae1b"))
  # plot(sp_raster, legend = FALSE, col = cols_q(5))
  # plot(americas, add = T)
  # 
  # #sp <- sp[!is.na(inout),]
  # 
  # # add occurrence points to map
  # title(main = paste(four_points$species[i]), cex.main = 4)
  # 
  # 
  # dev.off()
  
}


