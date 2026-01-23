#' Global phenology extraction using TIMESAT method 1 (dynamic threshold)
#'
#' @param r SpatRaster, raster time series
#' @param doy Numeric vector, day of year
#' @param years Integer vector, year for each layer
#' @param outdir Character, output folder
#' @param threshold Numeric, relative amplitude fraction
#' @param years_sel Optional integer vector, subset of years
#'
#' @export
global_pheno_extraction_TS1 <- function(
    r, doy, years, outdir,
    threshold = 0.5,
    years_sel = NULL
) {

  library(terra)

  ncell <- ncell(r)
  vals <- values(r)

  if (is.null(years_sel)) years_sel <- sort(unique(years))

  SOS_mat <- EOS_mat <- POS_mat <-
    matrix(NA_real_, nrow=ncell, ncol=length(years_sel))

  for (i in 1:ncell) {
    y_all <- vals[i, ]

    fits <- fit_pixel_DL_yearly(
      y_all, doy, years, years_sel = years_sel)

    ph <- extract_pheno_pixel_TS1(
      fits,
      threshold = threshold,
      years_sel = years_sel
    )

    SOS_mat[i, ] <- ph[, "SOS"]
    EOS_mat[i, ] <- ph[, "EOS"]
    POS_mat[i, ] <- ph[, "POS"]


    if (i %% 5000 == 0)
      cat("Processed pixel:", i, "/", ncell, "\n")
  }

  ## --- 保存 TIFF（完全复用你现有逻辑） ---
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  for (j in seq_along(years_sel)) {
    yr <- years_sel[j]

    for (v in c("SOS","EOS","POS")) {
      mat <- get(paste0(v, "_mat"))
      single_r <- rast(
        nrows=nrow(r), ncols=ncol(r),
        ext=ext(r), crs=crs(r))
      single_r <- setValues(single_r, mat[, j])

      writeRaster(
        single_r,
        filename=file.path(outdir, paste0(v, "_", yr, "_TS1.tif")),
        overwrite=TRUE
      )
    }
  }

  cat("TIMESAT method 1 phenology saved to:", outdir, "\n")

  invisible(list(SOS=SOS_mat, EOS=EOS_mat, POS=POS_mat))
}
