rm(list=ls())
setwd("~/Documents/14 U of T/Mahler lab/Michelle/Analysis/bioclim wrangling")

# takes occurrence data (each specimen as lat and long coordinates)
# using occurrence data from Velasco 2020
# calculates bioclim values at each of those lat and long coordinates
# (note that you need bioclim variables downloaded)
# then takes mean of bioclim variables for each sp
# code originally written by Dan Greenberg in 2017
# and then modified/updated by Gavia 2022 (now using terra instead of raster)
# https://cran.r-project.org/web/packages/terra/terra.pdf 

library(terra)
library(sp)
#library('spocc')
require(maptools)  
library(ggplot2)
library(dplyr)
library(sf)

##### read in occurrence data ####
# there are 389 species here
# pulled from Velasco 2020 sup mat
lat.long <- read.csv("ALL_RECORDS_COMBINED_31052017.csv")[, -c(1,2)]

# compress gaigei into one record and rename others?
lat.long$species <- ifelse(lat.long$species=="gaigei_1" | lat.long$species=="gaigei_2" | lat.long$species=="gaigei_3", "gaigei", lat.long$species)

lat.long$species <- gsub("/", ".", lat.long$species)
sp <- unique(lat.long$species)


##### read in Bioclim data ####
# bioclim data downloaded from https://www.worldclim.org/data/worldclim21.html, Jan 2022
# updated to use terra instead
rasts <- rast(list.files("./wc2.1_30s_bio", full.names=T)) # in terra
name1 <- as.character(names(rasts)) #name of bioclim rasters

##### read in shapefiles ####
ranges <- sf::st_read ("./Velasco et al 2020_sup mat/RangesAnolis 2/ANOLE_RANGE_MAPS_POE_etal_TAXONOMY_27032017.shp")
ranges$BINOMIAL <- gsub("/", ".", ranges$BINOMIAL)
# notes that there are fewer shape files than there are sp occurences (378 vs 389)
uni <- unique(ranges$BINOMIAL)
setdiff(sp,uni)


##### Bioclim values per specimen #####
# using Velasco 2020 occurrence records 
for(i in sp){
  print(i)
  # pulls out lat and long of each sp
  i_lat.long <- subset(lat.long, species==i)[, 2:3] 
  rownames(i_lat.long) <- 1:nrow(i_lat.long)
  i_ranges <- subset(ranges, BINOMIAL==i)
  
  # convert i_lat.long into simple feature data
  i_lat.longSP <- st_as_sf(i_lat.long, coords = c("Longitude", "Latitude"), 
                   crs = 4326, agr = "constant")
  st_crs(i_ranges) <- CRS("+init=epsg:4326")  # add same coordinate system
  
  # error handling for if no range map for sp
  default <- "missing range map"  
  try({default<-st_intersects(i_lat.longSP, i_ranges, sparse=F)[,1]}, silent=TRUE)
  intersections <- default
  
  # removes points outside of range (note that seems to be overly conservative, removing points conceivably within range)
  keep <- which(intersections == TRUE)
  i_lat.long_trim <- i_lat.long[keep,]
  
  values<-c()
  for (j in 1:nlyr(rasts)){  
    rast.value<-lapply(terra::extract(rasts[[j]], i_lat.long_trim), na.exclude)[[2]]
    length=length(unlist(rast.value))
    values<-append(values,unlist(rast.value))
  }
  
  # associates bioclim names with each species specific file
  assign(gsub(" ",".",i), as.data.frame(matrix(values,nrow=length,ncol=19)))
  y=get(gsub(" ",".",i))
  colnames(y)=name1
  assign(gsub(" ",".",i),y)
  write.csv(y,paste("./Species_bioclimfiles_ranges/",gsub(" ",".",i),".csv",sep=""))

}
  
# read in sp files
setwd("Species_bioclimfiles_ranges")
files <- list.files()

# rename variables so look nicer
get_name <- function(x){
  name_temp <- x[3:4]
  paste(name_temp[1], name_temp[2])
}

bio_names <- do.call(rbind, lapply(strsplit(name1, "_"), get_name))

##### calculate mean for each bioclim variable over each occurrence record #####
# read in indiv sp files, extract sp name, rename bioclim columns, calculate mean for each bioclim variable
mean_bioclim = function(path) {
  temp.dat = read.csv(path, header = T)[, -1]
  names(temp.dat)[1:19] <- bio_names
  name <- gsub('.csv','', path)
  sample.length <- nrow(temp.dat)
  bio.means <- data.frame(t(data.frame(colMeans(temp.dat))))
  bio.means$species <- paste(name)
  bio.means$sample.length <- sample.length
  return(bio.means)
}

# bind all together
bioclim.data_list = lapply(files, mean_bioclim)
bioclim.data_df = data.frame(do.call(rbind, bioclim.data_list))
rownames(bioclim.data_df) <- 1:nrow(bioclim.data_df)

# there are 3 species for which there are no bioclim records
setwd("~/Documents/14 U of T/Mahler lab/Michelle/Analysis/bioclim wrangling")
#write.csv(bioclim.data_df, "anolis_bioclim_ranges.csv")

# differences between bioclim avgs?
means1 <- read.csv(file="anolis_bioclim_ranges.csv")
means2 <- read.csv(file="anolis_bioclim.csv")

means.join <- left_join(means1, means2, by="species")

dev.off()
plot(means.join$bio.2.x, means.join$bio.2.y, 
     xlab="occurences trimmed to range map",
     ylab="all cccurences")
abline(a=0, b=1, col="darkred")

 ##### plot occurrences ####
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html 
world <- map_data("world")
# only keep parts of map where actually have occurrences (+ a little bit of buffer)
anolis_world <- subset(world, long<=max(lat.long$Longitude)+10 & long>=min(lat.long$Longitude)-10 &
                         lat<=max(lat.long$Latitude)+10 & lat>=min(lat.long$Latitude)-10 )

# a reminder of which bioclim variables are which: https://www.worldclim.org/data/bioclim.html 
ggplot() + 
  geom_polygon(data = anolis_world, aes(x=long, y = lat, group = group), 
               fill="grey80", color="grey20") + 
  theme_bw() +
  geom_point(data=na.omit(left_join(lat.long, bioclim.data_df, by="species")),            
             aes(Longitude, Latitude, colour=bio.1), size=2, alpha=.5) +
  labs(y= "Latitude", x = "Longitude", size=3)


# any occurrences really far from others?
outliers <- lat.long %>% 
  group_by(species) %>% 
  filter(n()>1) %>% 
  mutate(mean.long=mean(Longitude), sd.long=4*sd(Longitude),
         mean.lat=mean(Latitude), sd.lat=4*sd(Latitude))

outliers$out.long <- ifelse(outliers$Longitude>=outliers$mean.long+outliers$sd.long, 1, 0)
outliers$out.lat <- ifelse(outliers$Latitude>=outliers$mean.lat+outliers$sd.lat, 1, 0)

outliers.trim <- subset(outliers, out.long!=0 | out.lat!=0)
#write.csv(outliers.trim, "potential outlier occurrences.csv")

ggplot() + 
  geom_polygon(data = anolis_world, aes(x=long, y = lat, group = group), 
               fill="grey80", color="grey20") + 
  theme_bw() +
  geom_point(data=subset(outliers, out.lat==1),            
             aes(Longitude, Latitude, colour=species), size=2, alpha=.5) +
  labs(y= "Latitude", x = "Longitude", size=3)


