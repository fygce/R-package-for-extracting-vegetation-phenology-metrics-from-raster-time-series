library(terra)

# -----------------------------
# 1. 模拟单像元时间序列或使用已有向量
ntime <- 12
time <- seq(as.Date("1995-01-01"), as.Date("1995-12-31"), length.out = ntime)
doy  <- as.numeric(format(time, "%j"))
years <- as.integer(format(time, "%Y"))

# 假设 v_orig2 是已有的原始像元序列
# y <- v_orig2[idx_1995]  # 如果已有真实数据
# 这里模拟数据
set.seed(42)
y <- c(0,0,0,0,0.26,0.80,0.59,0.09,0.003,0,0,0)

pixel_index <- 1222

# -----------------------------
# 2. 拟合 DL 曲线并绘制
fit <- FitDL.Zhang(y, doy, tout = doy)
fit_par <- fit$par[1, ]

# 生成密集横轴平滑曲线
t_dense <- seq(min(doy), max(doy), by = 1)
fhat_dense <- doubleLog.Zhang(fit_par, t_dense)

# 绘制拟合
plot(doy, y, type="p", col="blue", pch=19,
     main=paste("Pixel", pixel_index, "Year", unique(years)),
     xlab="DOY", ylab="Value")
lines(t_dense, fhat_dense, col="red", lwd=2)
legend("topleft", legend=c("Original", "DL Fitted"),
       col=c("blue","red"), pch=c(19,NA), lty=c(NA,1))

# -----------------------------
# 3. 生成单像元 SpatRaster 并包装成 phenoRaster
r_single <- rast(nrows=1, ncols=ntime, nlyrs=1, vals=fhat_dense[seq(1, ntime)],
                 crs=NA)
ts_dl_single <- list(
  rast = r_single,
  time = time
)
class(ts_dl_single) <- "phenoRaster"

# 4. 计算 TIMESAT 风格相对阈值
# 单像元版本，不使用 terra::app
calc_relative_threshold_single <- function(y, threshold=0.5, trim=0.1) {
  # y: 向量 (该像元完整时间序列)
  y <- y[!is.na(y)]
  if (length(y) < 6) return(NA_real_)

  n <- length(y)
  k <- floor(n * trim)
  if (n - 2*k <= 0) return(NA_real_)

  y_trim <- sort(y)[(k+1):(n-k)]

  m <- max(1, floor(length(y_trim) * 0.1))
  base <- mean(y_trim[1:m])
  maxv <- mean(y_trim[(length(y_trim)-m+1):length(y_trim)])

  if (!is.finite(base) || !is.finite(maxv) || maxv <= base) return(NA_real_)

  base + (maxv - base) * threshold
}

# 使用示例
thr_val <- calc_relative_threshold_single(fhat_dense, threshold=0.5, trim=0.1)
thr_val


# -----------------------------
# 5. 提取 SOS/EOS/POS
pheno_single <- extract_pheno_relative_DL(ts_dl_single, thr_single)

# 6. 可视化 SOS/EOS/POS
sos <- terra::values(pheno_single$SOS)[1]
eos <- terra::values(pheno_single$EOS)[1]
pos <- terra::values(pheno_single$POS)[1]

plot(doy, y, type="p", col="blue", pch=19,
     main="Single Pixel Phenology",
     xlab="DOY", ylab="Value")
lines(t_dense, fhat_dense, col="red", lwd=2)
abline(h=terra::values(thr_single), col="green", lty=2)
points(sos, approx(t_dense, fhat_dense, xout=sos)$y, col="orange", pch=17, cex=1.5)
points(eos, approx(t_dense, fhat_dense, xout=eos)$y, col="purple", pch=17, cex=1.5)
points(pos, max(fhat_dense), col="brown", pch=19, cex=1.5)
legend("topleft", legend=c("Original","DL Fitted","Threshold","SOS","EOS","POS"),
       col=c("blue","red","green","orange","purple","brown"),
       pch=c(19,NA,NA,17,17,19), lty=c(NA,1,2,NA,NA,NA))
