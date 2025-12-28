# ==============================================================
# phenoraster_example.R
# Example workflow for extracting SOS/EOS/POS from raster
# time series (demonstrated using CMIP6 GPP data)
# ==============================================================

# -------------------------------
# Step 0: Load required packages
# -------------------------------
library(phenoraster)  # your package
library(terra)
library(raster)
library(phenofit)

# -------------------------------
# Step 1: Define input/output paths
# -------------------------------
indir <- system.file("extdata/cmip6_gpp/historical/ACCESS-ESM1-5/tif_lon180",
                     package = "phenoraster")

outdir <- file.path(tempdir(), "pheno_output")
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
# Step 3: Single-pixel DL fit test
#        (randomly selected or manually specified)
# -------------------------------
set.seed(123)  # ensure reproducibility

vals <- values(r)

# Identify pixels with at least some valid values
good_pixels <- which(rowSums(!is.na(vals) & vals != 0) > 0)

# User can manually specify a pixel index;
# otherwise, one pixel is randomly sampled
manual_pixel_index <- NULL  # set to a number (e.g., 1222) to fix a pixel
pixel_index <- if (!is.null(manual_pixel_index)) {
  manual_pixel_index
} else {
  sample(good_pixels, 1)
}
cat("Selected test pixel index:", pixel_index, "\n")

v_all <- vals[pixel_index, ]

# Fit double-logistic (DL) curves
fits <- fit_pixel_DL_yearly(v_all, doy, years)

# Calculate relative amplitude threshold
thr <- calc_relative_threshold_pixel(fits, threshold = 0.5)

# Extract phenological metrics
pheno_pixel <- t(sapply(fits, extract_pheno_DL_year, thr = thr))

# Plot example for one year
yr <- rownames(pheno_pixel)[1]
fit <- fits[[yr]]

plot(fit$t_dense, fit$fhat, type = "l", lwd = 2, col = "red",
     main = paste("Pixel", pixel_index, "Year", yr),
     xlab = "DOY", ylab = "GPP")
points(fit$t_obs, fit$y_obs, pch = 19, col = "blue")
abline(h = thr, lty = 2, col = "grey40")
abline(v = pheno_pixel[yr, "SOS"], col = "green", lwd = 2)
abline(v = pheno_pixel[yr, "POS"], col = "red", lwd = 2)
abline(v = pheno_pixel[yr, "EOS"], col = "brown", lwd = 2)
legend("topleft",
       legend = c("DL fit", "Original", "Threshold", "SOS", "POS", "EOS"),
       col = c("red", "blue", "grey40", "green", "red", "brown"),
       lty = c(1, NA, 2, 1, 1, 1),
       pch = c(NA, 19, NA, NA, NA, NA))

# -------------------------------
# Step 4: Global phenology extraction (subset of years)
# -------------------------------
test_years <- 2000:2001  # subset of years for testing

# ---- Serial run ----
cat("Running serial phenology extraction...\n")
res_serial <- global_pheno_extraction(
  r = r,
  doy = doy,
  years = years,
  outdir = outdir,
  threshold = 0.5,
  years_sel = test_years
)

cat("Phenology extraction finished. Results saved to:", outdir, "\n")
