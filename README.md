R scripts for making range files for species without adequate data for generating species distribution models 



alpha_hulls.R generates alpha convex hulls for species with enough occurrence points to make valid shape

buffer_species.R makes 30-km radius circles around occurrence points for species with less than 3 valid occurrence points or for which polygons make no biological sense

polygons_bee.R generates minimum convex hulls for species with 3 or more occurrence points but for which alpha hull method did not produce a result


