require(WEGE)
require(sp)
require(sf)

# Dom: this example can fail randomly.
species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)){
  temp  <-  suppressWarnings(Polygon(cbind(runif(3,1,50),runif(3,1,50))))
  range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))
}
input <- st_as_sf(SpatialPolygons(range_list))
categories <- c('LC','NT','VU','EN','CR')
input$binomial <- species
input$category <- sample(size = nrow(input),x = categories,replace = TRUE)

target_area <- suppressWarnings(Polygon(cbind(runif(3,1,50),runif(3,1,50))))
target_area <- Polygons(list(target_area), ID = 'Target area')
target_area <- st_as_sf(SpatialPolygons(list(target_area)))
suppressWarnings(get_kba_criteria(target_area,input))