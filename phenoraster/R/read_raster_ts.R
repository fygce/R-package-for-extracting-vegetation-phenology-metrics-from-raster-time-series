#' Read raster time series from a folder
#'
#' @param indir Character, folder containing monthly/regular TIFFs
#' @param pattern Regex pattern to select files (default: "\\.tif$")
#' @return A list containing SpatRaster and time vector
#' @export
read_raster_ts <- function(indir, pattern="\\.tif$") {
  library(terra)

  files <- list.files(indir, pattern = pattern, full.names = TRUE)
  stopifnot(length(files) > 0)

  # extract YYYYMM
  datestr <- gsub("\\D", "", basename(files))
  dates <- as.Date(paste0(substr(datestr,1,4), "-", substr(datestr,5,6), "-15"))

  ord <- order(dates)
  files <- files[ord]
  dates <- dates[ord]

  r <- rast(files)
  crs(r) <- "EPSG:4326"

  list(raster=r, dates=dates)
}
