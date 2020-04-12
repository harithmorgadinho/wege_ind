require(WEGE) # Dom: I think it's better to use library, use require within functions
require(sp)
require(sf)

species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)) {
  # Dom: do expect warnings? If so, maybe add a comment to that effect?
  temp  <-  suppressWarnings(Polygon(cbind(runif(3,1,50),runif(3,1,50))))
  range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))
}
input <- st_as_sf(SpatialPolygons(range_list))
categories <- c('LC','NT','VU','EN','CR')
input$binomial <- species
input$category <- sample(size = nrow(input),x = categories,replace = TRUE)
input$ED <- runif(nrow(input),1,30)
target_area <- suppressWarnings(Polygon(cbind(runif(3,1,50),runif(3,1,50))))
target_area <- Polygons(list(target_area), ID = 'Target area')
target_area <- st_as_sf(SpatialPolygons(list(target_area)))
get_edge(target_area,input,species = 'binomial',category = 'category')