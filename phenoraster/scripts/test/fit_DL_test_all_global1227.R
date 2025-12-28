# -------------------------------
# Step 0: 加载必要包
# -------------------------------
library(terra)
library(stats)
library(phenofit)

# -------------------------------
# Step 1: 读取 ACCESS-ESM1-5 GPP 栅格
# -------------------------------
indir <- "H:/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180"
files <- list.files(indir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(files) > 0)

datestr <- gsub("\\D", "", basename(files))
dates <- as.Date(paste0(substr(datestr,1,4), "-", substr(datestr,5,6), "-15"))
ord <- order(dates)
files <- files[ord]
dates <- dates[ord]

r <- rast(files)
crs(r) <- "EPSG:4326"
cat("Raster loaded:\n")
print(r)





# -------------------------------
# Step 2: 构造时间轴
# -------------------------------
time  <- dates
years <- as.integer(format(time, "%Y"))
doy   <- as.numeric(format(time, "%j"))





# -------------------------------
# Step 3: 准备输出栅格
# -------------------------------
ncell <- ncell(r)
ny    <- length(unique(years))
SOS_r <- rast(r); EOS_r <- rast(r); POS_r <- rast(r)

# -------------------------------
# Step 4: 定义单像元 DL 拟合函数（按年份）
# -------------------------------
fit_pixel_DL_yearly <- function(y_all, doy_all, year_vec, years_sel=NULL) {
  if (is.null(years_sel)) years_sel <- sort(unique(year_vec))
  fits <- list()
  for (yr in years_sel) {
    idx <- which(year_vec==yr & !is.na(y_all))
    if (length(idx)<6) next
    y <- y_all[idx]; t <- doy_all[idx]
    fit <- try(phenofit::FitDL.Zhang(y,t,tout=t), silent=TRUE)
    if (!inherits(fit,"try-error") && !is.null(fit$par)) {
      par_opt <- fit$par[1,]
      t_dense <- seq(min(t), max(t), by=1)
      fhat    <- phenofit::doubleLog.Zhang(par_opt, t_dense)
      fits[[as.character(yr)]] <- list(par=par_opt, t_dense=t_dense, fhat=fhat)
    } else {
      fits[[as.character(yr)]] <- NULL
    }
  }
  fits
}





# -------------------------------
# Step 5: 相对阈值计算函数（像元级）
# -------------------------------
calc_relative_threshold_pixel <- function(fits, threshold=0.5, trim=0.1) {
  all_vals <- unlist(lapply(fits, function(x) if(!is.null(x)) x$fhat))
  all_vals <- all_vals[is.finite(all_vals)]
  if(length(all_vals)<20) return(NA_real_)
  v_sorted <- sort(all_vals)
  n <- length(v_sorted)
  k <- floor(n*trim)
  if(n-2*k<=0) return(NA_real_)
  v_trim <- v_sorted[(k+1):(n-k)]
  m <- max(1,floor(length(v_trim)*0.1))
  base <- mean(v_trim[1:m])
  maxv <- mean(v_trim[(length(v_trim)-m+1):length(v_trim)])
  if(!is.finite(base) || !is.finite(maxv) || maxv<=base) return(NA_real_)
  base + (maxv-base)*threshold
}





# -------------------------------
# Step 6: 解 SOS/EOS/POS
# -------------------------------
extract_pheno_DL_year <- function(fit_year, thr) {
  if(is.null(fit_year)) return(c(SOS=NA, EOS=NA, POS=NA))
  t <- fit_year$t_dense
  f <- fit_year$fhat
  diff <- f - thr
  idx_sos <- which(diff[-length(diff)]<0 & diff[-1]>=0)
  sos <- if(length(idx_sos)) approx(diff[c(idx_sos[1], idx_sos[1]+1)], t[c(idx_sos[1], idx_sos[1]+1)], xout=0)$y else NA
  idx_eos <- which(diff[-length(diff)]>=0 & diff[-1]<0)
  eos <- if(length(idx_eos)) approx(diff[c(idx_eos[length(idx_eos)], idx_eos[length(idx_eos)]+1)], t[c(idx_eos[length(idx_eos)], idx_eos[length(idx_eos)]+1)], xout=0)$y else NA
  pos <- t[which.max(f)]
  c(SOS=sos, EOS=eos, POS=pos)
}




# -------------------------------
# Step 7: 循环部分像元（子集年份）
# -------------------------------

# 测试用：只处理 3 年（例如 2000-2002）
test_years <- 2000:2002
unique_years <- intersect(sort(unique(years)), test_years)

ncell <- ncell(r)
vals <- values(r)

SOS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))
EOS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))
POS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))

for(i in 1:ncell) {
  y_all <- vals[i,]
  fits <- fit_pixel_DL_yearly(y_all, doy, years, years_sel=unique_years)
  thr  <- calc_relative_threshold_pixel(fits, threshold=0.5)
  if(!is.na(thr)) {
    pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr=thr))
    SOS_mat[i,] <- pheno_pixel[,"SOS"]
    EOS_mat[i,] <- pheno_pixel[,"EOS"]
    POS_mat[i,] <- pheno_pixel[,"POS"]
  }
  if(i %% 5000 == 0) cat("Processed pixel:", i, "/", ncell, "\n")
}

cat("Test run finished for years:", paste(unique_years, collapse=", "), "\n")


# -------------------------------
# Step 7: 循环处理部分像元（测试用年份）
# -------------------------------

# 测试用：只处理 3 年（例如 2000-2002）
test_years <- 2000:2002
unique_years <- intersect(sort(unique(years)), test_years)

ncell <- ncell(r)
vals <- values(r)

SOS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))
EOS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))
POS_mat <- matrix(NA_real_, nrow=ncell, ncol=length(unique_years))

for(i in 1:ncell) {
  y_all <- vals[i,]
  fits <- fit_pixel_DL_yearly(y_all, doy, years, years_sel=unique_years)
  thr  <- calc_relative_threshold_pixel(fits, threshold=0.5)
  if(!is.na(thr)) {
    pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr=thr))
    SOS_mat[i,] <- pheno_pixel[,"SOS"]
    EOS_mat[i,] <- pheno_pixel[,"EOS"]
    POS_mat[i,] <- pheno_pixel[,"POS"]
  }
  if(i %% 5000 == 0) cat("Processed pixel:", i, "/", ncell, "\n")
}

cat("Test run finished for years:", paste(unique_years, collapse=", "), "\n")


# -------------------------------
# Step 8: 构建单层栅格列表
# -------------------------------
SOS_r <- vector("list", length(unique_years))
EOS_r <- vector("list", length(unique_years))
POS_r <- vector("list", length(unique_years))

for(j in seq_along(unique_years)) {
  SOS_r[[j]] <- setValues(r[[1]], SOS_mat[,j])  # 单层模板
  EOS_r[[j]] <- setValues(r[[1]], EOS_mat[,j])
  POS_r[[j]] <- setValues(r[[1]], POS_mat[,j])
}

names(SOS_r) <- paste0("SOS_", unique_years)
names(EOS_r) <- paste0("EOS_", unique_years)
names(POS_r) <- paste0("POS_", unique_years)

cat("All done. Outputs are SOS_r, EOS_r, POS_r\n")


# -------------------------------
# Step 9: 保存每年单独 TIFF
# -------------------------------
outdir <- "H:/cmip6_gpp/output_pheno_test"  # 测试输出文件夹
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

for(j in seq_along(unique_years)) {
  yr <- unique_years[j]
  
  # SOS
  SOS_file <- file.path(outdir, paste0("SOS_", yr, ".tif"))
  writeRaster(SOS_r[[j]], filename = SOS_file, overwrite = TRUE)
  
  # EOS
  EOS_file <- file.path(outdir, paste0("EOS_", yr, ".tif"))
  writeRaster(EOS_r[[j]], filename = EOS_file, overwrite = TRUE)
  
  # POS
  POS_file <- file.path(outdir, paste0("POS_", yr, ".tif"))
  writeRaster(POS_r[[j]], filename = POS_file, overwrite = TRUE)
  
  if(j %% 5 == 0) cat("Saved year:", yr, "\n")  # 每5年提示一次
}

cat("All yearly SOS/EOS/POS rasters saved to:", outdir, "\n")

