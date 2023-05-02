# created Jan18 2017 by Dan Greenberg
# modified by Gavia 27Jan2022

# code to extract maximum and minimum latitudes for each species
# requires lat and long coordinates for specimens per sp

setwd("~/Documents/14 U of T/Mahler lab/Michelle/Analysis/bioclim wrangling")

#read in data set
#these are lat long coordinates of specimens from Velasco 2020 data
lat.long <- read.csv("ALL_RECORDS_COMBINED_31052017.csv")

#get max and min lat values for each species and find largest absolute 
max_y=data.frame(tapply(lat.long$Latitude, lat.long$species, max))  
min_y=data.frame(tapply(lat.long$Latitude, lat.long$species, min))
lat=cbind(max_y, min_y) 
names(lat)[1] <- "max_y"
names(lat)[2] <- "min_y"
lat$abs_lat=ifelse(abs(lat$min_y) > abs(lat$max_y), abs(lat$min_y), abs(lat$max_y))

head(lat)

write.csv(lat, file="Anolis latitudinal limits.csv")
