#' @name get_kba_criteria
#' @title Function to get the KBA criteria
#' @description A function to get the WEGE index value for a provided polygon.
#' @param target_area  Either a sp or sf object to which to calculate the WEGE index.
#' @param input Species ranges, either from a shapefile or from a georeferenced species list with a column for species, two columns for coordinates and one for 
#' the IUCN category.
#' @param x name of the longitude column.
#' @param y name of the latitude column.
#' @param  species name of the species column.
#' @param  category name of IUCN the category column. Terminology must be as follows: DD for Data Deficient, LC for Least Concern, NT for Near Threatened, 
#' EN, for Endangered, CR for Critically Endangered, EW for Extinct in the wild 
#' and EX for Extinct.
#' @param  res grid-cell size to use to calculate the range of the species in 
#' case a georeferenced species list was provided.
#' @examples
#' require(sp)
#' require(sf)
#' 
#' species <- letters[1:26]
#' range_list <- list()
#' for (i in seq_along(species)){
#'   temp  <-  Polygon(cbind(runif(3,1,50),runif(3,1,50)))
#'   range_list[[i]] <- Polygons(list(temp), ID = c(species[i]))}
#' input <- st_as_sf(SpatialPolygons(range_list))
#' categories <- c('LC','NT','VU','EN','CR')
#' input$binomial <- species
#' input$category <- sample(size = nrow(input),x = categories,replace = TRUE)

#' target_area <- Polygon(cbind(runif(3,1,50),runif(3,1,50)))
#' target_area <- Polygons(list(target_area), ID = 'Target area')
#' target_area <- st_as_sf(SpatialPolygons(list(target_area)))
#' get_kba_criteria(target_area,input)
#'
#'
#'@importFrom sf st_as_sf
#'@importFrom sf st_geometry
#'@importFrom sf st_crs
#'@importFrom sf st_intersection
#'@importFrom raster extent
#'@importFrom sp CRS
#'@importFrom raster raster
#'@importFrom raster ncell
#'@importFrom graphics plot
#'@importFrom utils stack
#' @export

get_kba_criteria <- function(target_area,input,x,y,species='binomial', category 
                             = 'category',res = 1) {

  
  if (any(class(input) %in% "sf")) {
    input_cl <- 'sf_ob'}else {input_cl <- 'df_ob'}

  if (any(class(target_area) %in% 'SpatialPolygonsDataFrame')) {
    target_area <- sf::st_as_sf(target_area)
  }
  if (any(class(input) %in% 'SpatialPolygonsDataFrame')) {
    input <- sf::st_as_sf(input)
    if (!sf::st_crs(target_area) == sf::st_crs(input)) {
      stop("Inputs have a different projection")
    }
  }
  if (any(class(input) %in% 'SpatialPolygonsDataFrame')) {
    input <- sf::st_as_sf(input)
    if (!sf::st_crs(target_area) == sf::st_crs(input)) {
      stop("Inputs have a different projection")
    }
  }
  if (any(class(input) %in% "data.frame")) {
    crs_ta <- sf::st_crs(target_area)
    input <-  sf::st_as_sf(x = input,coords = c(x,y),crs = crs_ta)
  }

  sps_grid <- sf::st_intersects(input,target_area)
  intersected_object_t <- t(sps_grid)
  sp_numbers <- unlist(intersected_object_t[1:nrow(target_area)])
  sp <- unique(input[[species]][sp_numbers])
  if (identical(sp, character(0))) {
    stop("No species found in selected area")}

  if (input_cl == 'df_ob') {

    rgrid  <-  raster::raster(raster::extent(input), resolution = res,crs = 
                                sp::CRS(crs_ta$proj4string))
    rgrid[] <- 1:raster::ncell(rgrid)
    rgrid <- sf::st_as_sf(raster::rasterToPolygons(rgrid))

    iucn_to_grid_range <- function(iucn_shp,grid_to_use) {
      r_grid_sf <- grid_to_use
      sf_to_intersect <- iucn_shp
      sf::st_crs(r_grid_sf) <- sf::st_crs(sf_to_intersect)
      sps_grid <- sf::st_intersects(sf_to_intersect,r_grid_sf)
      intersected_object <- sps_grid



      area <- unlist(lapply(intersected_object,length))

      sp_range_df <- cbind.data.frame(species = iucn_shp[[species]],area = area)

      return(sp_range_df)
    }
    input <- input[input[[species]] %in% sp,]
    input_combined <- stats::aggregate(input,
                                by = list(input[[species]]),
                                FUN = mean)
    input_combined <- input_combined[,c('Group.1','geometry')]
    colnames(input_combined)[1] <- species

    tmp <- iucn_to_grid_range(iucn_shp = input_combined,grid_to_use = rgrid)
    tmp <- merge(tmp,input[,c(species,'category')],by.x = 'species',by.y = 
                   species)
    tmp <- tmp[,-4]
    tmp <- unique(tmp)

    min_rangeA1a <- (200*res)+1
    min_rangeA1b <- (100*res)+1
    min_rangeB1 <- (1000*res)+1

    KBA_A1a <- tmp[tmp$category == "CR"|tmp$category == "EN" & tmp$area < 
                     min_rangeA1a,]
    if(nrow(KBA_A1a) == 0){
      KBA_A1a_2 <- data.frame()
    }else{KBA_A1a_2 <- cbind.data.frame(species = KBA_A1a$species,A1a = 'yes')}


    KBA_A1b <- tmp[tmp$category == "VU" & tmp$area < min_rangeA1b,]
    if(nrow(KBA_A1b) == 0){
      KBA_A1b_2 <- data.frame()
      
    }else{KBA_A1b_2 <- cbind.data.frame(species = KBA_A1b$species,A1b = 'yes')}

    KBA_A1e <- tmp[tmp$perc_kba == 100,]
    if(nrow(KBA_A1e) == 0){
      KBA_A1e_2 <- data.frame()
    }else{
      KBA_A1e_2 <- cbind.data.frame(species = KBA_A1e$species,A1e = 'yes')}

    KBA_B1 <- tmp[tmp$perc_kba < min_rangeB1,]
    if(nrow(KBA_B1) == 0){
      KBA_B1_2 <- data.frame()
    }else{
      KBA_B1_2 <- cbind.data.frame(species = KBA_B1$species,B1 = 'yes')}

  }else {

    area_input <- function(input,sp) {
      temp <- sf::st_area(input[input[[species]] %in% sp,])
      attributes(temp) <- NULL
      temp <- sum(temp)
      return(temp <- temp/1000000)
    }
    all_area <- lapply(sp,area_input,input = input)
    tmp <- cbind.data.frame(species = sp,area = unlist(all_area))
    input_subset <- st_intersection(input,target_area)
    area_kba <- lapply(sp,area_input,input = input_subset)
    tmp_2 <- cbind.data.frame(species = sp,area_kba = unlist(area_kba))
    tmp <- merge(tmp,input[,c(species,'category')],by.x = 'species',by.y = 
                   species)
    tmp <- tmp[,-4]
    tmp <- unique(tmp)
    tmp <- merge(tmp,tmp_2,by.x = 'species',by.y = 'species')
    tmp$perc_kba <- 100*(tmp$area_kba/tmp$area)
    min_rangeA1a <- 1
    min_rangeA1b <- 0.5
    min_rangeB1 <- 10


    KBA_A1a <- tmp[tmp$category == "CR"|tmp$category == "EN" & tmp$perc_kba > min_rangeA1a,]
    if(nrow(KBA_A1a) == 0){
      KBA_A1a_2=data.frame()
    }else{KBA_A1a_2 <- cbind.data.frame(species = KBA_A1a$species,A1a = 'yes')}


    KBA_A1b=tmp[tmp$category =="VU" & tmp$perc_kba > min_rangeA1b,]
    if(nrow(KBA_A1b) == 0){
      KBA_A1b_2=data.frame()
    }else{KBA_A1b_2 <- cbind.data.frame(species = KBA_A1b$species,A1b = 'yes')}

    KBA_A1e=tmp[tmp$perc_kba == 100,]
    if(nrow(KBA_A1e) == 0){
      KBA_A1e_2 <- data.frame()
    }else{
      KBA_A1e_2 <- cbind.data.frame(species = KBA_A1e$species,A1e = 'yes')}

    KBA_B1 <- tmp[tmp$perc_kba>min_rangeB1,]
    if(nrow(KBA_B1) == 0){
      KBA_B1_2 <- data.frame()
    }else{
      KBA_B1_2 <- cbind.data.frame(species = KBA_B1$species,B1 = 'yes')}

  }

  kba_df_tmp=unique(rbind(KBA_A1a,KBA_A1b,KBA_A1e,KBA_B1))
  if (nrow(kba_df_tmp) == 0){
    return(cat('No species found to trigger KBA status\n'))
  }else{
    if(nrow(KBA_A1a_2) == 0){
      kba_df_tmp$A1a <- 'no'}else{
        kba_df_tmp <- merge(kba_df_tmp,KBA_A1a_2,by.x = 'species',by.y 
                            = 'species',all.x=TRUE)}
    if(nrow(KBA_A1b_2) == 0){
      kba_df_tmp$A1b <- 'no'}else{
        kba_df_tmp <- merge(kba_df_tmp,KBA_A1b_2,by.x = 'species',by.y 
                            = 'species',all.x=TRUE)}
    if(nrow(KBA_A1e_2) == 0){
      kba_df_tmp$A1e <- 'no'}else{
        kba_df_tmp <- merge(kba_df_tmp,KBA_A1e_2,by.x = 'species',by.y 
                            = 'species',all.x=TRUE)}
    if(nrow(KBA_B1_2) == 0){
      kba_df_tmp$B1 <- 'no'}else{
        kba_df_tmp <- merge(kba_df_tmp,KBA_B1_2,by.x = 'species',by.y 
                            = 'species',all.x=TRUE)}
    
    kba_df_tmp$A1a <- as.character(kba_df_tmp$A1a)
    kba_df_tmp$A1b <- as.character(kba_df_tmp$A1b)
    kba_df_tmp$A1e <- as.character(kba_df_tmp$A1e)
    kba_df_tmp$B1 <- as.character(kba_df_tmp$B1)
    kba_df_tmp[is.na(kba_df_tmp)]='no'

    return(kba_df_tmp)
  }
}
