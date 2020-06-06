
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- devtools::rmarkdown::render("README.Rmd") -->
<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->
WEGE index
==========

WEGE is an R package that allows the user to calculate the WEGE index for a particular area. Additionally it also calculates rasters of KBA criteria (A1a, A1b, A1e, and B1) Weighted endemism, the EDGE score, Evolutionary Distinctiveness and Extinction risk.

Install
=======

The package can currently only be installed through GitHub:

``` r
# install.packages("remotes")
remotes::install_github("harithmorgadinho/wege_ind")
```

Usage
=====

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
#> [1] 19.20848
```

A get\_wege example:

``` r
library(sp)
library(sf)
library(WEGE)

species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)){
  temp  <-  Polygon(cbind(runif(4,1,50),runif(4,1,50)))
  range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))}
input <- st_as_sf(SpatialPolygons(range_list))
categories <- c('LC','NT','VU','EN','CR')
input$binomial <- species
input$category <- sample(size = nrow(input),x = categories,replace = TRUE)

target_area <- Polygon(cbind(runif(4,1,50),runif(4,1,50)))
target_area <- Polygons(list(target_area), ID = 'Target area')
target_area <- st_as_sf(SpatialPolygons(list(target_area)))
get_wege(target_area,input,species = 'binomial',category = 'category')
#> [[1]]
#> [1] 872.4039
```

A get\_kba-criteria example:

``` r
library(WEGE)
library(sp)
library(sf)

species <- letters[1:26]
range_list <- list()
for (i in seq_along(species)){
  temp0 <- cbind(runif(3,1,50),runif(3,1,50))
  temp  <-  Polygon(rbind(temp0,temp0[1,]))
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
get_kba_criteria(target_area,input)
#>   species         area category     area_kba    perc_kba A1a A1b A1e  B1
#> 1       d 3.044359e-04       CR 2.628431e-05  8.63377278 yes  no  no  no
#> 2       k 6.708239e-05       EN 5.435383e-06  8.10254778 yes  no  no  no
#> 3       o 3.379152e-06       CR 5.523783e-07 16.34665554 yes  no  no yes
#> 4       s 5.944492e-04       CR 3.284138e-07  0.05524673 yes  no  no  no
#> 5       x 2.204984e-04       VU 3.454541e-05 15.66696904  no yes  no yes
#> 6       z 2.528314e-05       NT 7.475781e-06 29.56824960  no  no  no yes
```

A raster example example:

``` r
library(WEGE)
library(sp)
library(sf)
library(raster)
#> Warning: package 'raster' was built under R version 3.5.2

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
 
 input$ed <- runif(runif(10,1,50))
 temp0 <- cbind(runif(3,1,50),runif(3,1,50))
 target_area <- Polygon(rbind(temp0,temp0[1,]))
 target_area <- Polygons(list(target_area), ID = 'Target area')
 target_area <- st_as_sf(SpatialPolygons(list(target_area)))
 spat_ras(target_area,input,species = 'binomial',ed='ed', res = 1)
#> 1 6722 6723 6724 6725 6726 6727 6728 6729 67210 67211 67212 67213 67214 67215 67216 67217 67218 67219 67220 67221 67222 67223 67224 67225 67226 67227 67228 67229 67230 67231 67232 67233 67234 67235 67236 67237 67238 67239 67240 67241 67242 67243 67244 67245 67246 67247 67248 67249 67250 67251 67252 67253 67254 67255 67256 67257 67258 67259 67260 67261 67262 67263 67264 67265 67266 67267 67268 67269 67270 67271 67272 67273 67274 67275 67276 67277 67278 67279 67280 67281 67282 67283 67284 67285 67286 67287 67288 67289 67290 67291 67292 67293 67294 67295 67296 67297 67298 67299 672100 672101 672102 672103 672104 672105 672106 672107 672108 672109 672110 672111 672112 672113 672114 672115 672116 672117 672118 672119 672120 672121 672122 672123 672124 672125 672126 672127 672128 672129 672130 672131 672132 672133 672134 672135 672136 672137 672138 672139 672140 672141 672142 672143 672144 672145 672146 672147 672148 672149 672150 672151 672152 672153 672154 672155 672156 672157 672158 672159 672160 672161 672162 672163 672164 672165 672166 672167 672168 672169 672170 672171 672172 672173 672174 672175 672176 672177 672178 672179 672180 672181 672182 672183 672184 672185 672186 672187 672188 672189 672190 672191 672192 672193 672194 672195 672196 672197 672198 672199 672200 672201 672202 672203 672204 672205 672206 672207 672208 672209 672210 672211 672212 672213 672214 672215 672216 672217 672218 672219 672220 672221 672222 672223 672224 672225 672226 672227 672228 672229 672230 672231 672232 672233 672234 672235 672236 672237 672238 672239 672240 672241 672242 672243 672244 672245 672246 672247 672248 672249 672250 672251 672252 672253 672254 672255 672256 672257 672258 672259 672260 672261 672262 672263 672264 672265 672266 672267 672268 672269 672270 672271 672272 672273 672274 672275 672276 672277 672278 672279 672280 672281 672282 672283 672284 672285 672286 672287 672288 672289 672290 672291 672292 672293 672294 672295 672296 672297 672298 672299 672300 672301 672302 672303 672304 672305 672306 672307 672308 672309 672310 672311 672312 672313 672314 672315 672316 672317 672318 672319 672320 672321 672322 672323 672324 672325 672326 672327 672328 672329 672330 672331 672332 672333 672334 672335 672336 672337 672338 672339 672340 672341 672342 672343 672344 672345 672346 672347 672348 672349 672350 672351 672352 672353 672354 672355 672356 672357 672358 672359 672360 672361 672362 672363 672364 672365 672366 672367 672368 672369 672370 672371 672372 672373 672374 672375 672376 672377 672378 672379 672380 672381 672382 672383 672384 672385 672386 672387 672388 672389 672390 672391 672392 672393 672394 672395 672396 672397 672398 672399 672400 672401 672402 672403 672404 672405 672406 672407 672408 672409 672410 672411 672412 672413 672414 672415 672416 672417 672418 672419 672420 672421 672422 672423 672424 672425 672426 672427 672428 672429 672430 672431 672432 672433 672434 672435 672436 672437 672438 672439 672440 672441 672442 672443 672444 672445 672446 672447 672448 672449 672450 672451 672452 672453 672454 672455 672456 672457 672458 672459 672460 672461 672462 672463 672464 672465 672466 672467 672468 672469 672470 672471 672472 672473 672474 672475 672476 672477 672478 672479 672480 672481 672482 672483 672484 672485 672486 672487 672488 672489 672490 672491 672492 672493 672494 672495 672496 672497 672498 672499 672500 672501 672502 672503 672504 672505 672506 672507 672508 672509 672510 672511 672512 672513 672514 672515 672516 672517 672518 672519 672520 672521 672522 672523 672524 672525 672526 672527 672528 672529 672530 672531 672532 672533 672534 672535 672536 672537 672538 672539 672540 672541 672542 672543 672544 672545 672546 672547 672548 672549 672550 672551 672552 672553 672554 672555 672556 672557 672558 672559 672560 672561 672562 672563 672564 672565 672566 672567 672568 672569 672570 672571 672572 672573 672574 672575 672576 672577 672578 672579 672580 672581 672582 672583 672584 672585 672586 672587 672588 672589 672590 672591 672592 672593 672594 672595 672596 672597 672598 672599 672600 672601 672602 672603 672604 672605 672606 672607 672608 672609 672610 672611 672612 672613 672614 672615 672616 672617 672618 672619 672620 672621 672622 672623 672624 672625 672626 672627 672628 672629 672630 672631 672632 672633 672634 672635 672636 672637 672638 672639 672640 672641 672642 672643 672644 672645 672646 672647 672648 672649 672650 672651 672652 672653 672654 672655 672656 672657 672658 672659 672660 672661 672662 672663 672664 672665 672666 672667 672668 672669 672670 672671 672672 672
```

![](README-raster-example-1.png)

    #> class      : RasterStack 
    #> dimensions : 14, 48, 672, 10  (nrow, ncol, ncell, nlayers)
    #> resolution : 1, 1  (x, y)
    #> extent     : 1.424993, 49.42499, 8.060887, 22.06089  (xmin, xmax, ymin, ymax)
    #> crs        : NA 
    #> names      :       A1a,       A1b,       A1e,        B1,        GE,        ED,      EDGE,      WEGE,        WE,      KBAs 
    #> min values :  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000, -4.037394,  0.000000,  0.000000,  0.000000 
    #> max values :  1.000000,  1.000000,  1.000000,  1.000000,  0.441800,  1.733206,  0.000000,  0.004418,  0.000300,  1.000000

Citation:
=========

Farooq, H., Azevedo, J., Belluardo F., Nanvonamuquitxo, C., Bennett, D., Moat, J. Soares, A., Faurby, S., Antonelli, A. WEGE: A new metric for ranking locations for biodiversity conservation. in prep.
