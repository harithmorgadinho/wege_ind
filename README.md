
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- devtools::rmarkdown::render("README.Rmd") -->
<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->
WEGE index
==========

WEGE is an R package that allows the user to calculate the WEGE index for a particular area. Additionally it also calculates rasters of KBA criteria (A1a, A1b, A1e, and B1) Weighted endemism, the EDGE score, Evolutionary Distintiveness and Extintion risk.

Install
=======

The package can currently only be installed through GitHub:

``` r
# install.packages("remotes")
remotes::install_github("harithmorgadinho/wege_ind")
```

Usage
=====

A spatial raster example:

``` r
library(WEGE)
library(sp)
library(sf)

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

A get\_edge example:

``` r
library(WEGE) 
library(sp)
#> Warning: package 'sp' was built under R version 3.5.2
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.1.3, PROJ 4.9.3

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
#> [[1]]
#> [1] 17.94026
```

Citation:
=========

Farooq, H.
