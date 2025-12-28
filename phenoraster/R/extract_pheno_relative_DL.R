#' Extract phenology (SOS/EOS/POS) from DL curve using threshold
#'
#' @param fit_year List, single year DL fit
#' @param thr Numeric, threshold value
#' @return Named vector c(SOS, EOS, POS)
#' @export
extract_pheno_DL_year <- function(fit_year, thr) {
  t <- fit_year$t_dense
  f <- fit_year$fhat
  diff <- f - thr

  # SOS
  idx_sos <- which(diff[-length(diff)] < 0 & diff[-1] >= 0)
  sos <- if(length(idx_sos)) {
    i <- idx_sos[1]
    approx(diff[c(i,i+1)], t[c(i,i+1)], xout=0)$y
  } else NA_real_

  # EOS
  idx_eos <- which(diff[-length(diff)] >= 0 & diff[-1] < 0)
  eos <- if(length(idx_eos)) {
    i <- idx_eos[length(idx_eos)]
    approx(diff[c(i,i+1)], t[c(i,i+1)], xout=0)$y
  } else NA_real_

  # POS
  pos <- t[which.max(f)]

  c(SOS=sos, EOS=eos, POS=pos)
}
