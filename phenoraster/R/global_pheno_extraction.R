#' Extract phenology for all pixels and years
#'
#' @param r SpatRaster, raster stack (multi-layer GPP)
#' @param doy Numeric, DOY vector
#' @param years Integer, year vector
#' @param outdir Character, folder to save yearly single-band TIFFs
#' @param threshold Numeric, relative amplitude threshold fraction
#' @param years_sel Numeric, subset of years (optional)
#' @export
global_pheno_extraction <- function(r, doy, years, outdir, threshold=0.5, years_sel=NULL) {
  library(terra)

  ncell <- ncell(r)
  vals <- values(r)

  if (is.null(years_sel)) years_sel <- sort(unique(years))

  SOS_mat <- EOS_mat <- POS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(years_sel))

  # ------------------------
  # Step 1: Loop over pixels
  # ------------------------
  for(i in 1:ncell) {
    y_all <- vals[i, ]

    fits <- fit_pixel_DL_yearly(y_all, doy, years, years_sel=years_sel)
    thr <- calc_relative_threshold_pixel(fits, threshold=threshold)

    if(!is.na(thr)) {
      pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr=thr))
      SOS_mat[i, ] <- pheno_pixel[,"SOS"]
      EOS_mat[i, ] <- pheno_pixel[,"EOS"]
      POS_mat[i, ] <- pheno_pixel[,"POS"]
    }

    if(i %% 5000 == 0) cat("Processed pixel:", i, "/", ncell, "\n")
  }

  # ------------------------
  # Step 2: Create output folder
  # ------------------------
  if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

  # ------------------------
  # Step 3: Save each year as single-band TIFF
  # ------------------------
  for(j in seq_along(years_sel)) {
    yr <- years_sel[j]

    # SOS
    single_r <- rast(nrows=nrow(r), ncols=ncol(r), nlyrs=1,
                     ext=ext(r), crs=crs(r))
    single_r <- setValues(single_r, SOS_mat[, j])
    writeRaster(single_r, filename=file.path(outdir, paste0("SOS_", yr, ".tif")),
                overwrite=TRUE)

    # EOS
    single_r <- rast(nrows=nrow(r), ncols=ncol(r), nlyrs=1,
                     ext=ext(r), crs=crs(r))
    single_r <- setValues(single_r, EOS_mat[, j])
    writeRaster(single_r, filename=file.path(outdir, paste0("EOS_", yr, ".tif")),
                overwrite=TRUE)

    # POS
    single_r <- rast(nrows=nrow(r), ncols=ncol(r), nlyrs=1,
                     ext=ext(r), crs=crs(r))
    single_r <- setValues(single_r, POS_mat[, j])
    writeRaster(single_r, filename=file.path(outdir, paste0("POS_", yr, ".tif")),
                overwrite=TRUE)
  }

  cat("All phenology TIFFs saved to:", outdir, "\n")

  # Return matrices
  invisible(list(SOS=SOS_mat, EOS=EOS_mat, POS=POS_mat))
}
