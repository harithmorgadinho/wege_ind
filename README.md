
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- devtools::rmarkdown::render("README.Rmd") -->

<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->

# WEGE index

WEGE is an R package that allows the user to calculate the WEGE index
for a particular area. Additionally it also calculates rasters of KBA
criteria (A1a, A1b, A1e, and B1) Weighted endemism, the EDGE score,
Evolutionary Distintiveness and Extintion risk.

# Install

The package can currently only be installed through GitHub:

``` r
# install.packages("remotes")
remotes::install_github("harithmorgadinho/wege_ind")
```

# Usage

A simple example:

``` r
library(WEGE)
library(sp)
library(sf)

# Dom: Hmmm doesn't work for me this example ...not running as part of the .Rmd
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

# Example data

``` r
library(WEGE)
# example data can be accessed...
data("amph_afr_df")
head(amph_afr_df)
#>                BINOMIAL decimallongitude decimallatitude category       ED
#> 1   Acanthixalus sonjae         -7.33333         5.83333       NT 19.64515
#> 2   Acanthixalus sonjae         -7.33333         5.83333       NT 19.64515
#> 3 Acanthixalus spinosus         11.96117         2.98450       LC 19.67318
#> 4 Acanthixalus spinosus          9.68000         4.83000       LC 19.67318
#> 5 Acanthixalus spinosus         11.96670         2.98330       LC 19.67318
#> 6 Acanthixalus spinosus         11.96670         2.98330       LC 19.67318

data("kruger_park")
(kruger_park)
#> Simple feature collection with 1 feature and 28 fields
#> geometry type:  POLYGON
#> dimension:      XY
#> bbox:           xmin: 30.88982 ymin: -25.5303 xmax: 32.032 ymax: -22.32804
#> epsg (SRID):    4326
#> proj4string:    +proj=longlat +datum=WGS84 +no_defs
#>   WDPAID WDPA_PID PA_DEF                 NAME            ORIG_NAME
#> 1    873      873      1 Kruger National Park Kruger National Park
#>           DESIG     DESIG_ENG DESIG_TYPE     IUCN_CAT       INT_CRIT MARINE
#> 1 National Park National Park   National Not Assigned Not Applicable      0
#>   REP_M_AREA GIS_M_AREA REP_AREA GIS_AREA        NO_TAKE NO_TK_AREA     STATUS
#> 1          0          0 19169.15 19171.61 Not Applicable          0 Designated
#>   STATUS_YR                               GOV_TYPE OWN_TYPE
#> 1      1926 Federal or national ministry or agency    State
#>                      MANG_AUTH    MANG_PLAN          VERIF METADATAID SUB_LOC
#> 1 South African National Parks Not Reported State Verified       1840   ZA-MP
#>   PARENT_ISO ISO3                       geometry
#> 1        ZAF  ZAF POLYGON ((31.38601 -24.4843...
```
