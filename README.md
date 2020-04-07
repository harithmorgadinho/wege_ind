# WEGE index
WEGE is an R package that allows the user to calculate the WEGE index for a particular area. Additionally it also calculates rasters of KBA criteria (A1a, A1b, A1e, and B1) Weighted endemism, the EDGE score, Evolutionary Distintiveness and Extintion risk.

# Usage
A simple example:
```r
require(sp)
require(sf)

species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)){
  temp  <-  Polygon(cbind(runif(3,1,50),runif(3,1,50)))
  range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))}
input <- st_as_sf(SpatialPolygons(range_list))
categories <- c('LC','NT','VU','EN','CR')
input$binomial <- species
input$category <- sample(size = nrow(input),x = categories,replace = TRUE)

target_area <- Polygon(cbind(runif(3,1,50),runif(3,1,50)))
target_area <- Polygons(list(target_area), ID = 'Target area')
target_area <- st_as_sf(SpatialPolygons(list(target_area)))
spat_ras(target_area,input,species = 'binomial',res=0.2)
```
