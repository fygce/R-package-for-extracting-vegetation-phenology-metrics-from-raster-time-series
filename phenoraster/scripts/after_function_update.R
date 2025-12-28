# ============================================================
# Reinstall phenoraster after source code updates
# ============================================================

# --------- 1. Set package root path ---------
pkg_path <- "E:/GlobalFromWhiteToGreen/phenoraster"

stopifnot(dir.exists(pkg_path))
setwd(pkg_path)

cat("Working directory set to:\n", pkg_path, "\n")

# --------- 2. Regenerate documentation & NAMESPACE ---------
cat("Running devtools::document() ...\n")
devtools::document()

# --------- 3. Remove old installed version (if exists) -----
if ("phenoRaster" %in% rownames(installed.packages())) {
  cat("Removing existing phenoraster package ...\n")
  remove.packages("phenoraster")
}

# --------- 4. Install updated package ----------------------
cat("Installing updated phenoraster ...\n")
devtools::install(pkg_path, dependencies = TRUE, upgrade = "never")

# --------- 5. Load and verify --------------------------------
library(phenoRaster)

cat("Installed phenoraster version:\n")
print(packageVersion("phenoraster"))

cat("Available exported functions:\n")
print(ls("package:phenoraster"))

cat("Reinstall finished successfully.\n")
