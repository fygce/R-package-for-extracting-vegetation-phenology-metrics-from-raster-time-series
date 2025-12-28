# phenoRaster

`phenoRaster` is an R package for extracting vegetation phenology metrics (Start of Season, SOS; End of Season, EOS; Peak of Season, POS) from raster time series, either globally or regionally.

## Features

- Read monthly/temporal raster stacks as SpatRaster objects.
- Fit double-logistic (DL) curves for each pixel.
- Calculate a robust relative threshold for phenology extraction.
- Extract SOS, POS, EOS for each pixel and each year.
- Export results as single-band TIFF files per year.

## Installation

```r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("yourusername/phenoraster")

## Scope and assumptions

`phenoRaster` is designed for regularly sampled raster time series (e.g., monthly GPP, NDVI, LAI).
The example workflow uses CMIP6 GPP data, but the core functions are not limited to CMIP6 products.

Users should ensure that:
- Raster layers are ordered in time
- File names or metadata allow correct reconstruction of observation dates
- The time series within each year is sufficiently sampled for DL fitting




Example Workflow




# ==============================================================
# phenoRaster_example.R
# Example workflow for extracting SOS/EOS/POS from CMIP6 GPP
# ==============================================================

# -------------------------------
# Step 0: Load required packages
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
# Step 3: Single pixel DL fit test (random or manual)
# -------------------------------
set.seed(123)  # ensure reproducibility

vals <- values(r)

# Select pixels with at least some valid values
good_pixels <- which(rowSums(!is.na(vals) & vals != 0) > 0)

# User can specify pixel index manually, otherwise sample randomly
manual_pixel_index <- NULL  # set to a number (e.g., 1222) to fix manually
pixel_index <- if(!is.null(manual_pixel_index)) {
  manual_pixel_index
} else {
  sample(good_pixels, 1)
}
cat("Selected test pixel index:", pixel_index, "\n")

v_all <- vals[pixel_index, ]

# Fit DL curves
fits <- fit_pixel_DL_yearly(v_all, doy, years)

# Calculate relative threshold
thr <- calc_relative_threshold_pixel(fits, threshold=0.5)

# Extract phenology points
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
# Step 4: Global phenology extraction (subset of years)
# -------------------------------
test_years <- 2000:2001
res <- global_pheno_extraction(r, doy, years, outdir, threshold=0.5, years_sel=test_years)

cat("Test extraction finished. Results saved to:", outdir, "\n")
