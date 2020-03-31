#' @title get_edge
#' @description A function to get the WEGE index value for a provided polygon.
#' @param target_area  Either a sp or sf object to which to calculate the WEGE index.
#' @param input Species ranges, either from a shapefile or from a georeferenced species list with a column for species, two columns for coordinates and one for the IUCN category.
#' @param x name of the longitude column.
#' @param y name of the latitude column.
#' @param  species name of the species column.
#' @param  ed name of the evolutionary distinctiveness column.
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
get_edge <- function(target_area,input,x,y,species='binomial', ed = 'ED', category = 'category', res = 1) {
  require(sf)
  require(raster)

  if(is.null(input[[species]])){
    stop(paste0("No column found with the name - ",paste(species)))
  }
  if(is.null(input[[category]])){
    stop(paste0("No column found with the name - ",paste(category)))
  }
  if(is.null(input[[ed]])){
    stop(paste0("No column found with the name - ",paste(ed)))
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
  tmp <- input[input[[species]] %in% sp,]
  st_geometry(tmp) <- NULL
  tmp <- tmp[,c(species,category,ed)]
  tmp <- unique(tmp)
  er_df <- cbind.data.frame(status = c('DD','LC','NT','VU','EN','CR','EW','EX'),
                            ER = c(0.0513,0.0009,0.0071,0.0513,0.4276,0.9688,1,1))

  tmp <- unique(tmp)
  tmp <- merge(tmp, er_df, by.x = category, by.y = 'status', sort = TRUE)

  return(lapply(1, function(x) sum(log10((tmp[[ed]])*(tmp$ER + 1)))))
}

