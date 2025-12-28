# ==============================================================
# phenoRaster_example.R
# Example workflow for extracting SOS/EOS/POS from CMIP6 GPP
# ==============================================================

# -------------------------------
# Step 0: Load packages
# -------------------------------
library(phenoRaster)  # your package
library(terra)
library(raster)
library(phenofit)

# -------------------------------
# Step 1: Define input/output paths
# -------------------------------
indir <- "H:/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180"
outdir <- "H:/cmip6_gpp/output_pheno_test"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

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
# Step 3: Single pixel DL fit test
# -------------------------------
pixel_index <- 1222
v_all <- values(r)[pixel_index,]

fits <- fit_pixel_DL_yearly(v_all, doy, years, pixel_index=pixel_index)

thr <- calc_relative_threshold_pixel(fits, threshold=0.5)
pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr=thr))

# Plot example for one year
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
test_years <- 2000:2001
res <- global_pheno_extraction(r, doy, years, outdir, threshold=0.5, years_sel=test_years)

cat("Test extraction finished. Results saved to:", outdir, "\n")
