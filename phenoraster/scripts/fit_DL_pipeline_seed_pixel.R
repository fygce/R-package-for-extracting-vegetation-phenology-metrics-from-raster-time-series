# ==============================================================
# phenoraster_example.R
# Example workflow for extracting SOS/EOS/POS from CMIP6 GPP
# ==============================================================

# -------------------------------
# Step 0: Load packages
# -------------------------------
library(phenoraster)  # your package
library(terra)
library(raster)
library(phenofit)

# -------------------------------
# Step 1: Define input/output paths
# -------------------------------
# Input: example CMIP6 GPP data included in the package
indir <- system.file(
  "extdata/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180",
  package = "phenoraster"
)

stopifnot(dir.exists(indir))

# Output: temporary directory for demonstration
outdir <- file.path(tempdir(), "phenoraster_output_TS1")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Step 2: Read raster time series
# -------------------------------
ts <- read_raster_ts(indir)
r <- ts$raster
dates <- ts$dates

years <- as.integer(format(dates, "%Y"))
doy   <- as.numeric(format(dates, "%j"))

cat("Raster stack loaded:", ncell(r), "pixels,", length(dates), "layers\n")

# -------------------------------
# -------------------------------
# Step 3: Single pixel DL fit test (随机或手动)
# -------------------------------
set.seed(123)  # 保证可重复

vals <- values(r)

# 找到有至少一些有效值的像元
good_pixels <- which(rowSums(!is.na(vals) & vals != 0) > 0)

# 用户可以手动指定像元 index，否则随机
manual_pixel_index <- NULL  # 随机时为"NULL" 手动时例如 1222，如果想固定就改成数字
pixel_index <- if(!is.null(manual_pixel_index)) {
  manual_pixel_index
} else {
  sample(good_pixels, 1)
}
cat("Selected test pixel index:", pixel_index, "\n")

v_all <- vals[pixel_index, ]

# 拟合 DL
fits <- fit_pixel_DL_yearly(v_all, doy, years)

# 计算相对阈值
thr <- calc_relative_threshold_pixel(fits, threshold=0.5)

# 提取物候点
pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr=thr))

# 绘图
yr <- rownames(pheno_pixel)[1]
fit <- fits[[yr]]

plot(fit$t_dense, fit$fhat, type="l", lwd=2, col="red",
     main=paste("Pixel", pixel_index, "Year", yr),
     xlab="DOY", ylab="GPP")
points(fit$t_obs, fit$y_obs, pch=19, col="blue")
abline(h=thr, lty=2, col="grey40")
abline(v=pheno_pixel[yr, "SOS"], col="green", lwd=2)
abline(v=pheno_pixel[yr, "POS"], col="red", lwd=2)
abline(v=pheno_pixel[yr, "EOS"], col="brown", lwd=2)
legend("topleft", legend=c("DL fit","Original","Threshold","SOS","POS","EOS"),
       col=c("red","blue","grey40","green","red","brown"),
       lty=c(1,NA,2,1,1,1), pch=c(NA,19,NA,NA,NA,NA))


# -------------------------------
# Step 4: Global subset test (3 years)
# -------------------------------

test_years <- 2000:2001  # subset of years to test

# ---- Serial run ----
cat("Running serial extraction...\n")
res_serial <- global_pheno_extraction(
  r = r,
  doy = doy,
  years = years,
  outdir = outdir,
  threshold = 0.5,
  years_sel = test_years
)

cat("Serial phenology extraction finished. Results saved to:", outdir, "\n")


