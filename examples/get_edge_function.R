library(WEGE) 
library(sp)
library(sf)

species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)) {
  temp0 <- cbind(runif(3,1,50),runif(3,1,50))
  temp <- Polygon(rbind(temp0,temp0[1,]))
  range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))
}
input <- st_as_sf(SpatialPolygons(range_list))
categories <- c('LC','NT','VU','EN','CR')
input$binomial <- species
input$category <- sample(size = nrow(input),x = categories,replace = TRUE)
input$ED <- runif(nrow(input),1,30)
temp0 <- cbind(runif(3,1,50),runif(3,1,50))
target_area <- Polygon(rbind(temp0,temp0[1,]))
target_area <- Polygons(list(target_area), ID = 'Target area')
target_area <- st_as_sf(SpatialPolygons(list(target_area)))
get_edge(target_area = target_area,input = input,species = 'binomial',category = 'category')
