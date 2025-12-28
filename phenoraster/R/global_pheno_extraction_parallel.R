#' Extract phenology for all pixels and years (parallel version, Windows-safe)
#'
#' @param r SpatRaster, raster stack (multi-layer GPP)
#' @param doy Numeric, DOY vector
#' @param years Integer, year vector
#' @param outdir Character, output directory
#' @param threshold Numeric, relative amplitude threshold fraction
#' @param years_sel Numeric, subset of years (optional)
#' @param n_workers Integer, number of parallel workers
#'
#' @export
global_pheno_extraction_parallel <- function(
    r, doy, years, outdir,
    threshold = 0.5,
    years_sel = NULL,
    n_workers = 4
) {
  library(terra)
  library(future.apply)

  # ------------------------
  # Check years_sel
  # ------------------------
  if (is.null(years_sel)) years_sel <- sort(unique(years))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # ------------------------
  # Raster info in main process
  # ------------------------
  nr <- nrow(r)
  nc <- ncol(r)
  ex <- ext(r)
  cr <- crs(r)
  vals <- values(r)
  ncell <- ncell(r)

  # ------------------------
  # Parallel plan
  # ------------------------
  plan(multisession, workers = n_workers)

  # ------------------------
  # Parallel by year
  # ------------------------
  future_lapply(years_sel, function(yr) {
    cat("Processing year:", yr, "\n")

    # Create raster template using info only
    template <- rast(nrows = nr, ncols = nc, ext = ex, crs = cr)

    SOS_vec <- EOS_vec <- POS_vec <- rep(NA_real_, ncell)

    for (i in seq_len(ncell)) {
      y_all <- vals[i, ]
      fits <- fit_pixel_DL_yearly(y_all, doy, years, years_sel = yr)
      thr <- calc_relative_threshold_pixel(fits, threshold)
      if (!is.na(thr)) {
        ph <- extract_pheno_DL_year(fits[[as.character(yr)]], thr)
        SOS_vec[i] <- ph["SOS"]
        EOS_vec[i] <- ph["EOS"]
        POS_vec[i] <- ph["POS"]
      }
    }

    # Write single-band rasters
    writeRaster(setValues(template, SOS_vec),
                file.path(outdir, paste0("SOS_", yr, ".tif")),
                overwrite = TRUE)
    writeRaster(setValues(template, EOS_vec),
                file.path(outdir, paste0("EOS_", yr, ".tif")),
                overwrite = TRUE)
    writeRaster(setValues(template, POS_vec),
                file.path(outdir, paste0("POS_", yr, ".tif")),
                overwrite = TRUE)

    cat("Finished year:", yr, "\n")
    invisible(NULL)
  }, future.seed = TRUE)

  cat("All phenology TIFFs saved to:", outdir, "\n")
}
