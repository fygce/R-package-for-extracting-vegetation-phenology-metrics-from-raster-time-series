library(terra)
library(devtools)

# 加载你的包
devtools::load_all("E:/GlobalFromWhiteToGreen/phenoraster")

# 1. 读取时间序列
ts <- read_raster_ts(
  indir = "H:/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180"
)

stopifnot(inherits(ts, "phenoRaster"))
stopifnot(nlyr(ts$rast) == length(ts$time))

# 2. 年度 double logistic 拟合
ts_dl <- fit_DL_phenoRaster(ts)

stopifnot(inherits(ts_dl, "phenoRaster"))
stopifnot(nlyr(ts_dl$rast) == nlyr(ts$rast))

# 3. 计算相对阈值（TIMESAT method 3）
thr <- calc_relative_threshold(
  ts_dl,
  threshold = 0.5,  # 50% 相对振幅
  trim      = 0.1   # TIMESAT robust mean
)

plot(thr)  # 可视化检查

# 4. 提取 SOS/EOS/POS（基于 DL 拟合曲线）
pheno <- extract_pheno_relative_DL(
  ts_dl,
  thr
)

str(pheno)  # 检查输出结构

# 5. 输出到 GeoTIFF
outdir <- "H:/cmip6_gpp/historical/ACCESS-ESM1-5/phenology"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(pheno$years)) {
  yr <- pheno$years[i]

  writeRaster(
    pheno$SOS[[i]],
    file.path(outdir, paste0("SOS_", yr, ".tif")),
    overwrite = TRUE
  )

  writeRaster(
    pheno$EOS[[i]],
    file.path(outdir, paste0("EOS_", yr, ".tif")),
    overwrite = TRUE
  )

  writeRaster(
    pheno$POS[[i]],
    file.path(outdir, paste0("POS_", yr, ".tif")),
    overwrite = TRUE
  )
}

message("Pipeline 完成")
