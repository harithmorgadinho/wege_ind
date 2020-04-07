#' @title spat_ras
#' @description A function to get the WEGE index value for a provided polygon.
#' @param target_area  Either a sp or sf object to which to calculate the WEGE index.
#' @param input Species ranges, either from a shapefile or from a georeferenced species list with a column for species, two columns for coordinates and one for 
#' the IUCN category.
#' @param x name of the longitude column.
#' @param y name of the latitude column.
#' @param  species name of the species column.
#' @param  ed name of the evolutionary distinctiveness column.
#' @param  category name of IUCN the category column. Terminology must be as 
#' follows: DD for Data Deficient, LC for Least Concern, NT for Near Threatened, 
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
#' spat_ras(target_area,input,species = 'binomial',res=0.2)
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
#'
#' @export

spat_ras <- function(target_area, input, x, y, species='binomial',
                     category = 'category', ed = FALSE, res = 1) {


  if (is.null(input[[species]])) {
    stop(paste0("No column found with the name - ",paste(species)))
  }
  if (is.null(input[[category]])) {
    stop(paste0("No column found with the name - ",paste(category)))
  }
  if (!missing(ed)) {
    if (is.null(input[[ed]])) {
      stop(paste0("No column found with the name - ",paste(ed)))
    }}

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
    stop("No species found in selected area")
  }

  if (input_cl == 'df_ob') {

    rgrid  <-  raster(raster::extent(input), resolution = res,crs = sp::CRS(crs_ta$proj4string))
    rgrid[] <- 1:raster::ncell(rgrid)
    rgrid <- sf::st_as_sf(raster::rasterToPolygons(rgrid))

    iucn_to_grid_range <- function(iucn_shp,grid_to_use) {
      r_grid_sf <- grid_to_use
      sf_to_intersect <- iucn_shp
      sf::st_crs(r_grid_sf) <- sf::st_crs(sf_to_intersect)
      sps_grid <- sf::st_intersects(sf_to_intersect,r_grid_sf)
      intersected_object <- sps_grid



      area <- unlist(lapply(intersected_object,length))

      sp_range_df <- cbind.data.frame(species = iucn_shp[[species]],area=area/(res*res*10000))

      return(sp_range_df)
    }
    input <- input[input[[species]] %in% sp,]
    input_combined <- stats::aggregate(input,
                                by = list(input$BINOMIAL),
                                FUN = mean)
    input_combined <- input_combined[,c('Group.1','geometry')]
    colnames(input_combined)[1] <- species

    tmp <- iucn_to_grid_range(iucn_shp = input_combined,grid_to_use = rgrid)


  }else {

    area_input <- function(input,sp) {
      temp <- sf::st_area(input[input[[species]] %in% sp,])
      attributes(temp) <- NULL
      temp <- sum(temp)
      return(temp <- temp/1000000)
    }
    all_area <- lapply(sp,area_input,input = input)
    tmp <- cbind.data.frame(species = sp,area = unlist(all_area))
  }

  tmp <- merge(tmp,input[,c(species,'category')],by.x = 'species',by.y = 
                 species)
  tmp <- tmp[,-4]
  er_df <- cbind.data.frame(status = c('DD','LC','NT','VU','EN','CR','EW','EX'),
                            ER = c(0.0513,0.0009,0.0071,0.0513,0.4276,0.9688,1,1
                                   ))

  tmp <- unique(tmp)
  tmp <- merge(tmp, er_df, by.x = category, by.y = 'status', sort = TRUE)

  if (!missing(ed)) {
    tmp <- merge(tmp, input[,c(species,ed)], by.x = 'species', by.y = species, 
                 sort = TRUE)
    tmp <- tmp[,-6]
  }
  
  rgrid  <-  raster(raster::extent(target_area), resolution = res,crs = sp::CRS(crs_ta$proj4string))
  rgrid[] <- 1:raster::ncell(rgrid)
  rgrid <- sf::st_as_sf(raster::rasterToPolygons(rgrid))

  iucn_to_grid <- function(iucn_shp,grid_to_use) {

    sf::st_crs(grid_to_use) <- sf::st_crs(iucn_shp)
    sps_grid <- sf::st_intersects(iucn_shp,grid_to_use)
    intersected_object <- sps_grid
    intersected_object_t <- t(intersected_object)

    list_final <- list()
    for (i in seq_along(intersected_object_t)) {
      cat(i,length(intersected_object_t),'\n')
      list_final[i] <- intersected_object_t[i]
    }

    names(list_final) <- 1:nrow(grid_to_use)

    li_2 <- lapply(seq_along(list_final), function(i) {
      list_final[[i]] <- as.character(iucn_shp[[species]][list_final[[i]]])
    })

    li_3 <- li_2
    names(li_3) <- 1:nrow(grid_to_use)

    return(li_3)
  }
  sp_grid <- iucn_to_grid(iucn_shp = input,grid_to_use = rgrid)

  kba_size <- res*res*10000
  tmp[tmp$area < kba_size,]$area =  kba_size
  tmp$perc_kba <- 100*(kba_size/tmp$area)
  min_rangeA1a <- 0.5
  min_rangeA1b <- 1
  min_rangeB1 <- 10


  df_final = data.frame()

  for (i in seq_along(sp_grid)) {
    #cat(i,length(sp_grid),'\n')
    temp_df <- tmp[tmp$species %in% sp_grid[[i]],]

    if (nrow(temp_df) == 0) {
      we_temp <- 0
      wege_temp <- 0
      ed_temp <- 0
      edge_temp <- 0
      ge_temp <- 0
      KBA_A1a_temp <- 0
      KBA_A1b_temp <- 0
      KBA_A1e_temp <- 0
      KBA_B1_temp <- 0

    }else{

      we_temp <- lapply(1, function(x) sum(1/temp_df$area))

      wege_temp <- lapply(1, function(x) sum(sqrt(1/temp_df$area)*temp_df$ER))



      if (!missing(ed)) {
        edge_temp <- lapply(1, function(x) sum(log((temp_df$ED)*(temp_df$ER + 1
                                                                 ))))
        ed_temp <- lapply(1, function(x) sum(temp_df$ED))}else{
          ed_temp <- 0
          edge_temp <- 0
        }


      ge_temp <- lapply(1, function(x) sum(temp_df$ER))

      KBA_A1a_temp <- lapply(1, function(x) if (any(temp_df[temp_df$category == "CR" | temp_df$category == "EN",]$perc_kba > min_rangeA1a)) {1} else {0})

      KBA_A1b_temp <- lapply(1, function(x) if (any(temp_df[temp_df$category == "VU",]$perc_kba > min_rangeA1b)) {1} else {0})

      KBA_A1e_temp <- lapply(1, function(x) if (any(temp_df$perc_kba == 100)) {1
        } else {0})

      KBA_B1_temp <- lapply(1, function(x) if (any(temp_df$perc_kba > 
                                                   min_rangeB1)) {1} else {0})



    }
    df_temp <- cbind.data.frame(i,we=unlist(we_temp),wege=unlist(wege_temp),GE = unlist(ge_temp),ED = unlist(ed_temp),EDGE = unlist(edge_temp),kba_A1a = unlist(KBA_A1a_temp),kba_A1b = unlist(KBA_A1b_temp),kba_A1e = unlist(KBA_A1e_temp),
                                kba_B1 = unlist(KBA_B1_temp))
    df_final <- rbind(df_final,df_temp)
  }

  r  <-  raster(extent(target_area), resolution = res,crs = 
                  CRS(crs_ta$proj4string))

  r_GE <- r
  r_GE[] <- df_final$GE

  r_ED <- r
  r_ED[] <- df_final$ED

  r_EDGE <- r
  r_EDGE[] <- df_final$EDGE


  r_wege <- r
  r_wege[] <- df_final$wege


  r_we <- r
  r_we[] <- df_final$we

  r_A1a <- r
  r_A1a[] <- df_final$kba_A1a

  r_A1b <- r
  r_A1b[] <- df_final$kba_A1b

  r_A1e <- r
  r_A1e[] <- df_final$kba_A1e

  r_B1 <- r
  r_B1[] <- df_final$kba_B1

  r_Kbas <- sum(r_A1a,r_A1b,r_A1e,r_B1)
  r_Kbas[r_Kbas > 0] <- 1


  raster_stack <- stack(r_A1a,r_A1b,r_A1e,r_B1,r_GE,r_ED,r_EDGE,r_wege,r_we,
                        r_Kbas)
  names(raster_stack) <- c('A1a','A1b','A1e','B1','GE','ED','EDGE','WEGE','WE','KBAs')
  plot(r_Kbas)
  return(raster_stack)
}
