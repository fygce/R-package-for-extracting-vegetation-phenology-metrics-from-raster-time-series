#' Calculate relative amplitude threshold for a pixel (Timesat method 3)
#'
#' @param fits List of yearly DL fits
#' @param threshold Fraction of relative amplitude (default 0.5)
#' @return Numeric, relative threshold for all years of this pixel
#' @export
calc_relative_threshold_pixel <- function(fits, threshold=0.5) {
  all_vals <- unlist(lapply(fits, function(x) if(!is.null(x)) x$fhat))
  all_vals <- all_vals[is.finite(all_vals)]
  if(length(all_vals) < 20) return(NA_real_)

  v_sorted <- sort(all_vals)
  n <- length(v_sorted)
  k <- floor(n*0.1) # trim 10% top/bottom
  if(n - 2*k <= 0) return(NA_real_)

  v_trim <- v_sorted[(k+1):(n-k)]
  m <- max(1, floor(length(v_trim)*0.1))

  base <- mean(v_trim[1:m])
  maxv <- mean(v_trim[(length(v_trim)-m+1):length(v_trim)])

  if(!is.finite(base) || !is.finite(maxv) || maxv <= base) return(NA_real_)

  base + (maxv-base)*threshold
}
