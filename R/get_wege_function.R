#' @title get_wege
#' @description A function to get the WEGE index value for a provided polygon.
#' @param target_area  Either a sp or sf object to which to calculate the WEGE index.
#' @param input Species ranges, either from a shapefile or from a georeferenced species list with a column for species, two columns for coordinates and one for the IUCN category.
#' @param x name of the longitude column.
#' @param y name of the latitude column.
#' @param  species name of the species column.
#' @param  category name of IUCN the category column. Terminology must be as follows: DD for Data Deficient, LC for Least Concern, NT for Near Threatened, EN, for Endangered, CR for Critically Endangered, EW for Extinct in the wild and EX for Extinct.
#' @param  res grid-cell size to use to calculate the range of the species in case a georeferenced species list was provided.
#' @examples
#' get_wege(target_area,input,species = 'binomial',category = 'category')
#' get_wege(target_area,input,x = 'decimallongitude',y = 'decimallatitude',species = 'BINOMIAL',category = 'category')
#'
#'
#'
#' @export
#' @importFrom sf st_as_sf
#' @importFrom sf st_crs
#' @importFrom sf st_intersects
#' @importFrom sf st_area
#' @importFrom sf st_transform
#' @importFrom sf st_geometry
#' @importFrom stats aggregate
#' @importFrom raster rasterToPolygons

get_wege <- function(target_area,input,x,y,species='binomial',category = 'category',res = 1) {
  require(sf)
  require(raster)
if(is.null(input[[species]])){
  stop(paste0("No column found with the name - ",paste(species)))
}
  if(is.null(input[[category]])){
    stop(paste0("No column found with the name - ",paste(category)))
  }

  if (any(class(input) %in% "sf")) {
    input_cl <- 'sf_ob'}else {input_cl <- 'df_ob'}

  if (any(class(target_area) %in% 'SpatialPolygonsDataFrame')) {
    target_area <- st_as_sf(target_area)
  }
  if (any(class(input) %in% 'SpatialPolygonsDataFrame')) {
    input <- st_as_sf(input)
    if (!st_crs(target_area) == st_crs(input)) {
      stop("Inputs have a different projection")
    }
  }
  if (any(class(input) %in% 'SpatialPolygonsDataFrame')) {
    input <- st_as_sf(input)
  if (!st_crs(target_area) == st_crs(input)) {
    stop("Inputs have a different projection")
  }
  }
  if (any(class(input) %in% "data.frame")) {
    crs_ta <- st_crs(target_area)
    input <-  st_as_sf(x = input,coords = c(x,y),crs = crs_ta)
  }

    sps_grid <- st_intersects(input,target_area)
    intersected_object_t <- t(sps_grid)
    sp_numbers <- intersected_object_t[1]
    sp <- unique(input[[species]][sp_numbers[[1]]])
    if (identical(sp, character(0))) {
      stop("No species found in selected area")
    }

    if (input_cl == 'df_ob') {

      rgrid  <-  raster(extent(input), resolution = res,crs = CRS(crs_ta$proj4string))
      rgrid[] <- 1:ncell(rgrid)
      rgrid <- st_as_sf(rasterToPolygons(rgrid))

      iucn_to_grid_range <- function(iucn_shp,grid_to_use) {
        r_grid_sf <- grid_to_use
        sf_to_intersect <- iucn_shp
        st_crs(r_grid_sf) <- st_crs(sf_to_intersect)
        sps_grid <- st_intersects(sf_to_intersect,r_grid_sf)
        intersected_object <- sps_grid



        area <- unlist(lapply(intersected_object,length))

        sp_range_df <- cbind.data.frame(species = iucn_shp[[species]],area/(res*res*10000))

        return(sp_range_df)
      }
      input <- input[input[[species]] %in% sp,]
      input_combined <- aggregate(input,
                     by = list(input$BINOMIAL),
                     FUN = mean)
      input_combined <- input_combined[,c('Group.1','geometry')]
      colnames(input_combined)[1] <- species

      tmp <- iucn_to_grid_range(iucn_shp = input_combined,grid_to_use = rgrid)


    }else {

    area_input <- function(input,sp) {
      temp <- st_area(input[input[[species]] %in% sp,])
      attributes(temp) <- NULL
      temp <- sum(temp)
      return(temp <- temp/1000000)
    }
    all_area <- lapply(sp,area_input,input = input)
    tmp <- cbind.data.frame(species = input$binomial[sp_numbers[[1]]],area = unlist(all_area))
    }

    tmp <- merge(tmp,input[,c(species,'category')],by.x = 'species',by.y = species)
    tmp <- tmp[,-4]
    er_df <- cbind.data.frame(status = c('DD','LC','NT','VU','EN','CR','EW','EX'),ER = c(0.0513,0.0009,0.0071,0.0513,0.4276,0.9688,1,1))

    tmp <- unique(tmp)
    tmp <- merge(tmp, er_df, by.x = category, by.y = 'status', sort = TRUE)

    return(lapply(1, function(x) sum(sqrt(1/tmp$area)*tmp$ER)))
}

