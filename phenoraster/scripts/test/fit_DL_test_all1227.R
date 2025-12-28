# -------------------------------
# Step 0: 加载必要包
# -------------------------------
library(terra)
library(raster)
library(stats)

# TIMESAT / Zhang 双Logistic函数
# 确保你已经 source 过，或者在搜索路径中
# source("FitDL.Zhang.R")
# source("doubleLog.Zhang.R")
# -------------------------------
# Step 1: 读取 ACCESS-ESM1-5 GPP 栅格
# -------------------------------

indir <- "H:/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180"

files <- list.files(indir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(files) > 0)

# 从文件名中提取 YYYYMM（假定类似 gpp_199501.tif）
datestr <- gsub("\\D", "", basename(files))

dates <- as.Date(
  paste0(substr(datestr, 1, 4), "-",
         substr(datestr, 5, 6), "-15")
)

ord <- order(dates)
files <- files[ord]
dates <- dates[ord]

# 读成 SpatRaster
r <- rast(files)

crs(r) <- "EPSG:4326"  # 强制 CRS，避免后面报错

cat("Raster loaded:\n")
print(r)




# -------------------------------
# Step 2: 构造时间轴
# -------------------------------

time  <- dates
years <- as.integer(format(time, "%Y"))
doy   <- as.numeric(format(time, "%j"))

table(years)






# -------------------------------
# Step 3: 提取单像元时间序列
# -------------------------------

vals <- values(r)

pixel_index <- 1222   # 如果这个不行，后面可以换

v_all <- vals[pixel_index, ]

summary(v_all)

# 快速画原始时间序列
plot(time, v_all, type = "l",
     main = paste("Raw GPP - pixel", pixel_index),
     xlab = "Time", ylab = "GPP")





# -------------------------------
# -------------------------------
# Step 4: 单像元按年份拟合 DL (基于 phenofit)
# -------------------------------
fit_pixel_DL_yearly <- function(y_all, doy_all, year_vec, pixel_index=NULL, years_sel=NULL) {
  
  if (is.null(years_sel)) years_sel <- sort(unique(year_vec))
  
  fits <- list()
  
  # 可选绘图调试
  if (!is.null(pixel_index)) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(2,3), mar = c(4,4,3,1))
  }
  
  for (yr in years_sel) {
    
    # 提取当年有效观测值
    idx <- which(year_vec == yr & !is.na(y_all))
    
    if (length(idx) < 6) {  # 至少6个点才能拟合
      if (!is.null(pixel_index)) {
        plot.new()
        title(paste("Pixel", pixel_index, yr, "\n(no enough data)"))
      }
      next
    }
    
    y <- y_all[idx]
    t <- doy_all[idx]
    
    # DL 拟合
    fit <- try(phenofit::FitDL.Zhang(y, t, tout=t), silent = TRUE)
    
    # 绘图调试
    if (!is.null(pixel_index)) {
      plot(t, y, pch=19, col="blue",
           main=paste("Pixel", pixel_index, "Year", yr),
           xlab="DOY", ylab="Value")
    }
    
    if (!inherits(fit, "try-error") && !is.null(fit$par)) {
      par_opt <- fit$par[1, ]
      t_dense <- seq(min(t), max(t), by = 1)  # 每日加密
      fhat    <- phenofit::doubleLog.Zhang(par_opt, t_dense)
      
      if (!is.null(pixel_index)) {
        lines(t_dense, fhat, col="red", lwd=2)
      }
      
      fits[[as.character(yr)]] <- list(
        par     = par_opt,
        t_obs   = t,
        y_obs   = y,
        t_dense = t_dense,
        fhat    = fhat
      )
    } else {
      if (!is.null(pixel_index)) {
        mtext("DL fit failed", side=3, col="red", cex=0.8)
      }
      fits[[as.character(yr)]] <- NULL  # 失败年份保留空元素
    }
    
  }
  
  fits
}


library(phenofit)

pixel_index <- 1222
fits <- fit_pixel_DL_yearly(v_all, doy, years, pixel_index=pixel_index)

length(fits)
names(fits)





# -------------------------------
# Step 5: 画一个年份检查 DL 拟合
# -------------------------------

yr <- names(fits)[1]  # 先随便看一个
fit <- fits[[yr]]

plot(fit$t_obs, fit$y_obs, pch = 19, col = "blue",
     main = paste("DL fit check - pixel", pixel_index, yr),
     xlab = "DOY", ylab = "GPP")
lines(fit$t_dense, fit$fhat, col = "red", lwd = 2)

legend("topleft",
       legend = c("Original", "DL fitted"),
       col = c("blue", "red"),
       pch = c(19, NA), lty = c(NA, 1))




# -------------------------------
# Step 6: 计算相对振幅阈值（像元级）
# -------------------------------

calc_relative_threshold_pixel <- function(fits,
                                          threshold = 0.5,
                                          trim = 0.1) {
  
  all_vals <- unlist(lapply(fits, function(x) x$fhat))
  all_vals <- all_vals[is.finite(all_vals)]
  
  if (length(all_vals) < 20) return(NA_real_)
  
  v_sorted <- sort(all_vals)
  n <- length(v_sorted)
  
  k <- floor(n * trim)
  if (n - 2 * k <= 0) return(NA_real_)
  
  v_trim <- v_sorted[(k + 1):(n - k)]
  
  m <- max(1, floor(length(v_trim) * 0.1))
  
  base <- mean(v_trim[1:m])
  maxv <- mean(v_trim[(length(v_trim) - m + 1):length(v_trim)])
  
  if (!is.finite(base) || !is.finite(maxv) || maxv <= base)
    return(NA_real_)
  
  base + (maxv - base) * threshold
}

thr_pixel <- calc_relative_threshold_pixel(fits, threshold = 0.5)
thr_pixel






# -------------------------------
# Step 7: 从 DL 曲线中解 SOS / EOS / POS
# -------------------------------

extract_pheno_DL_year <- function(fit_year, thr) {
  
  t <- fit_year$t_dense
  f <- fit_year$fhat
  diff <- f - thr
  
  # SOS
  idx_sos <- which(diff[-length(diff)] < 0 & diff[-1] >= 0)
  sos <- if (length(idx_sos)) {
    i <- idx_sos[1]
    approx(diff[c(i, i + 1)], t[c(i, i + 1)], xout = 0)$y
  } else NA_real_
  
  # EOS
  idx_eos <- which(diff[-length(diff)] >= 0 & diff[-1] < 0)
  eos <- if (length(idx_eos)) {
    i <- idx_eos[length(idx_eos)]
    approx(diff[c(i, i + 1)], t[c(i, i + 1)], xout = 0)$y
  } else NA_real_
  
  pos <- t[which.max(f)]
  
  c(SOS = sos, EOS = eos, POS = pos)
}

pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr = thr_pixel))
pheno_pixel







# -------------------------------
# Step 8: 物候点可视化验证
# -------------------------------

yr <- rownames(pheno_pixel)[1]

fit <- fits[[yr]]
ph  <- pheno_pixel[yr, ]

plot(fit$t_dense, fit$fhat, type = "l", lwd = 2,
     main = paste("Phenology check - pixel", pixel_index, yr),
     xlab = "DOY", ylab = "GPP")

abline(h = thr_pixel, lty = 2, col = "grey40")
abline(v = ph["SOS"], col = "green", lwd = 2)
abline(v = ph["POS"], col = "red",   lwd = 2)
abline(v = ph["EOS"], col = "brown", lwd = 2)

legend("topleft",
       legend = c("DL", "Threshold", "SOS", "POS", "EOS"),
       col = c("black", "grey40", "green", "red", "brown"),
       lty = c(1, 2, 1, 1, 1))

















