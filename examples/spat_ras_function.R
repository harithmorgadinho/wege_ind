
library(WEGE)
library(sp)
library(sf)
library(raster)

 species <- sample(letters, 10)
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
 temp0 <- cbind(runif(3,1,50),runif(3,1,50))
 target_area <- Polygon(rbind(temp0,temp0[1,]))
 target_area <- Polygons(list(target_area), ID = 'Target area')
 target_area <- st_as_sf(SpatialPolygons(list(target_area)))
 spat_ras(target_area, input,species = 'binomial', res=0.3)

