#' Fit double-logistic (DL) curve for a single pixel across years
#'
#' This function fits the Zhang double-logistic curve to a pixel time series
#' on a per-year basis. It can be used for both single-pixel debugging (with plotting)
#' and for global raster processing (without plotting).
#'
#' @param y_all Numeric vector. Time series of the pixel (e.g., GPP).
#' @param doy_all Numeric vector. Day-of-year corresponding to each observation in y_all.
#' @param year_vec Integer vector. Year corresponding to each observation in y_all.
#' @param pixel_index Optional integer. If provided, enables plotting of per-year fits for debugging.
#' @param years_sel Optional integer vector. Subset of years to fit. Defaults to all years in year_vec.
#' @param min_n Integer. Minimum number of valid points required to fit a year. Default is 6.
#'
#' @return A named list of fitted results for each year. Each element contains:
#' \itemize{
#'   \item \code{par} - fitted DL parameters,
#'   \item \code{t_obs} - original DOY values,
#'   \item \code{y_obs} - original values,
#'   \item \code{t_dense} - daily interpolated DOY for fitted curve,
#'   \item \code{fhat} - fitted values corresponding to \code{t_dense}.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(phenofit)
#' r <- rast("H:/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180/gpp_195001.tif")
#' vals <- values(r)
#' pixel_index <- 1222
#' v_all <- vals[pixel_index, ]
#' doy <- as.numeric(format(dates, "%j"))
#' years <- as.integer(format(dates, "%Y"))
#' fits <- fit_pixel_DL_yearly(v_all, doy, years, pixel_index = pixel_index)
#' }
#'
#' @export
fit_pixel_DL_yearly <- function(y_all, doy_all, year_vec,
                                pixel_index = NULL,
                                years_sel = NULL,
                                min_n = 6) {
  
  if (is.null(years_sel)) years_sel <- sort(unique(year_vec))
  fits <- list()
  
  # Optional debugging plot
  if (!is.null(pixel_index)) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(2,3), mar = c(4,4,3,1))
  }
  
  for (yr in years_sel) {
    
    idx <- which(year_vec == yr & !is.na(y_all))
    if (length(idx) < min_n) {
      if (!is.null(pixel_index)) {
        plot.new()
        title(paste("Pixel", pixel_index, yr, "\n(not enough data)"))
      }
      next
    }
    
    y <- y_all[idx]
    t <- doy_all[idx]
    
    fit <- try(phenofit::FitDL.Zhang(y, t, tout = t), silent = TRUE)
    
    if (!is.null(pixel_index)) {
      plot(t, y, pch = 19, col = "blue",
           main = paste("Pixel", pixel_index, "Year", yr),
           xlab = "DOY", ylab = "Value")
    }
    
    if (!inherits(fit, "try-error") && !is.null(fit$par)) {
      par_opt <- fit$par[1, ]
      t_dense <- seq(min(t), max(t), by = 1)
      fhat <- phenofit::doubleLog.Zhang(par_opt, t_dense)
      
      if (!is.null(pixel_index)) {
        lines(t_dense, fhat, col = "red", lwd = 2)
      }
      
      fits[[as.character(yr)]] <- list(
        par = par_opt,
        t_obs = t,
        y_obs = y,
        t_dense = t_dense,
        fhat = fhat
      )
    } else {
      if (!is.null(pixel_index)) {
        mtext("DL fit failed", side = 3, col = "red", cex = 0.8)
      }
      fits[[as.character(yr)]] <- NULL
    }
  }
  
  return(fits)
}
