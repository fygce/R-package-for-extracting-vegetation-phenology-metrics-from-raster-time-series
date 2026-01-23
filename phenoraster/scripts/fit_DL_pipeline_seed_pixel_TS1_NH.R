# ==============================================================
# phenoraster_example_TS1_NH.R
# TIMESAT Method 1 (Dynamic threshold)
# Northern Hemisphere only (lat >= 0)
# ==============================================================

# -------------------------------
# Step 0: Load packages
# -------------------------------
library(phenoraster)
library(terra)
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
r  <- ts$raster
dates <- ts$dates

years <- as.integer(format(dates, "%Y"))
doy   <- as.numeric(format(dates, "%j"))

cat("Original raster:",
    ncell(r), "pixels,",
    length(dates), "layers\n")

# -------------------------------
# Step 2.5: Crop + mask to Northern Hemisphere (lat >= 0)
# -------------------------------

# 1) crop extent to lat >= 0
e_nh <- ext(r)
ymin(e_nh) <- 0
r <- crop(r, e_nh)

# 2) mask (safety)
lat <- init(r, "y")
r <- mask(r, lat >= 0)

cat("Cropped & masked to NH:",
    ncell(r), "pixels\n")
cat("NH extent:", as.character(ext(r)), "\n")

# -------------------------------
# Step 3: Single pixel test (TS1, NH only)
# -------------------------------
set.seed(123)

vals <- values(r)

# select pixels with valid data (NH only)
good_pixels <- which(rowSums(!is.na(vals) & vals != 0) > 0)

pixel_index <- sample(good_pixels, 1)
cat("Selected test pixel index:", pixel_index, "\n")

v_all <- vals[pixel_index, ]

# ---- DL fitting (yearly) ----
fits <- fit_pixel_DL_yearly(
  y_all = v_all,
  doy_all = doy,
  year_vec = years,
  pixel_index = pixel_index
)

# ---- Extract phenology (TS1) ----
pheno_pixel_TS1 <- extract_pheno_pixel_TS1(
  fits,
  threshold = 0.5
)

print(pheno_pixel_TS1)

# ---- Plot one example year ----
yr <- rownames(pheno_pixel_TS1)[1]
fit_year <- fits[[yr]]

thr_year <- calc_dynamic_threshold_year(
  fit_year,
  threshold = 0.5
)

plot(fit_year$t_dense, fit_year$fhat,
     type = "l", lwd = 2, col = "red",
     main = paste("Pixel", pixel_index, "Year", yr, "(TS1, NH)"),
     xlab = "DOY", ylab = "GPP")

points(fit_year$t_obs, fit_year$y_obs,
       pch = 19, col = "blue")

abline(h = thr_year, lty = 2, col = "grey40")
abline(v = pheno_pixel_TS1[yr, "SOS"], col = "green", lwd = 2)
abline(v = pheno_pixel_TS1[yr, "POS"], col = "red",   lwd = 2)
abline(v = pheno_pixel_TS1[yr, "EOS"], col = "brown", lwd = 2)

legend("topleft",
       legend = c("DL fit", "Original", "Dynamic threshold",
                  "SOS", "POS", "EOS"),
       col = c("red", "blue", "grey40",
               "green", "red", "brown"),
       lty = c(1, NA, 2, 1, 1, 1),
       pch = c(NA, 19, NA, NA, NA, NA),
       bty = "n")

# -------------------------------
# Step 4: Global phenology extraction (TS1, NH only)
# -------------------------------
test_years <- 2000:2001

cat("Running global TS1 extraction (NH only)...\n")

res_TS1 <- global_pheno_extraction_TS1(
  r = r,          # already cropped + masked
  doy = doy,
  years = years,
  outdir = outdir,
  threshold = 0.5,
  years_sel = test_years
)

cat("TS1 phenology extraction finished.\n")
cat("Results saved to:", outdir, "\n")

