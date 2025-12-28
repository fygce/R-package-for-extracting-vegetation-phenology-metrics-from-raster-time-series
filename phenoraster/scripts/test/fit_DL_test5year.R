debug_pixel_DL_years_FINAL <- function(
    y_all,
    doy_all,
    year_vec,
    pixel_index,
    years_sel = 1995:2000
) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  
  for (yr in years_sel) {
    
    idx <- which(year_vec == yr & !is.na(y_all))
    
    if (length(idx) < 6) {
      plot.new()
      title(paste("Pixel", pixel_index, yr, "\n(no data)"))
      next
    }
    
    y <- y_all[idx]
    t <- doy_all[idx]
    
    # —— 关键：这里直接重新拟合 ——
    fit <- try(FitDL.Zhang(y, t, tout = t), silent = TRUE)
    
    plot(t, y, pch = 19, col = "blue",
         main = paste("Pixel", pixel_index, "Year", yr),
         xlab = "DOY", ylab = "Value")
    
    if (!inherits(fit, "try-error") && !is.null(fit$par)) {
      par_opt <- fit$par[1, ]
      t_dense <- seq(min(t), max(t), by = 1)
      fhat    <- doubleLog.Zhang(par_opt, t_dense)
      lines(t_dense, fhat, col = "red", lwd = 2)
    } else {
      mtext("DL fit failed", side = 3, col = "red", cex = 0.8)
    }
  }
  
  invisible(TRUE)
}

debug_pixel_DL_years_FINAL(
  y_all       = v_orig2,
  doy_all     = doy,
  year_vec    = years,
  pixel_index = 1222,
  years_sel   = 1995:2000
)
