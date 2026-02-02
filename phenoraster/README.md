## phenoraster: R package for phenology extraction from raster time series

`phenoraster` is an R package for extracting vegetation phenology metrics
(Start of Season, Peak of Season, End of Season) from raster time series.

The package is based on **Double-Logistic (DL) curve fitting** and currently
supports **three phenology extraction strategies**, corresponding to widely
used methods in the phenology literature and the TIMESAT framework.

## Installation
install.packages("remotes")
remotes::install_github("fygce/phenoraster")

## Supported Phenology Extraction Methods
1. TS3 — Relative Threshold Method (Multi-year Climatological Threshold)
Corresponds to TIMESAT Method 3 (Relative Amplitude Threshold)

This method is inspired by the relative amplitude approach implemented in
TIMESAT 3.3.

For each pixel, a single climatological threshold is calculated from the
entire multi-year fitted DL time series.
Start of Season (SOS) and End of Season (EOS) are identified when the yearly
fitted DL curve crosses a fixed fraction of this climatological amplitude.

Key characteristics
Threshold is derived from the full multi-year growth curve

A single, fixed threshold is applied consistently across all years

Designed for long-term phenology and climatological analysis

Advantages
Robust to noise and inter-annual variability

Extreme values in individual years do not bias the threshold

Preserves spatial heterogeneity while maintaining temporal consistency

Phenological metrics are directly comparable across years

Well suited for long-term ecosystem and climate change studies

Important usage note
Important:
Because the threshold is defined from the full multi-year record,
all years must be processed together in a single run.

Running this method on segmented subsets of years will result in different
thresholds and introduce artificial discontinuities in SOS/EOS.

Example workflow script

fit_DL_pipeline_seed_pixel_TS3.R
Typical input data
CMIP6 monthly GPP

Other regularly sampled long-term raster time series

2. TS1 — Dynamic Threshold Method (Year-specific Threshold)
Corresponds to TIMESAT Method 1 (Dynamic Relative Threshold)

In this approach, the threshold is calculated independently for each year
based on the fitted DL curve of that year.

SOS and EOS are determined when the yearly DL curve reaches a given fraction
of the annual amplitude.

Key characteristics
Threshold varies from year to year

Fully adaptive to inter-annual variability

Each year is treated independently

Advantages
Sensitive to year-to-year phenological shifts

Well suited for climate variability and extreme event analysis

Allows segmented or incremental processing

Computationally flexible for large datasets

Example workflow script

fit_DL_pipeline_seed_pixel_TS1_NH.R
Typical input data
CMIP6 monthly GPP

Other regularly sampled raster time series

3. D2 — Second-Order Derivative Method
Second-derivative phenology extraction based on DL curves

This method extracts phenological metrics using the second-order derivative
of the fitted Double-Logistic curve.

SOS: first local maximum of the second derivative in the first half of the year

POS: timing of the maximum fitted value

EOS: symmetric definition based on curvature in the second half of the year

This approach is widely used for high-temporal-resolution vegetation indices,
where curvature-based detection is more appropriate than threshold-based methods.

Key characteristics
Does not rely on amplitude thresholds

Based on curve shape and inflection behavior

Particularly suitable for dense or quasi-dense time series

Example workflow script

fit_DL_pipeline_seed_pixel_D2_NH.R
Typical input data
GIMMS NDVI3g+

MODIS NDVI / EVI

Other biweekly or high-frequency vegetation indices

## Scripts and Workflows
The package provides both reusable functions and end-to-end workflow scripts.

Method	Workflow script
TS3 (Relative threshold, climatological)	fit_DL_pipeline_seed_pixel_TS3.R
TS1 (Dynamic threshold, yearly)	fit_DL_pipeline_seed_pixel_TS1_NH.R
D2 (Second derivative)	fit_DL_pipeline_seed_pixel_D2_NH.R

Each workflow script includes:

Single-pixel diagnostic testing

Visualization of fitted DL curves and phenological points

Global raster-scale extraction

Optional Northern Hemisphere masking

## Development Status and Extensibility
phenoraster is actively developed for extracting vegetation phenology metrics
from raster time series.

The current implementation focuses on Double-Logistic curve fitting, but
the package is designed to be extensible. Possible future extensions include:

Parallel and block-wise optimization

Additional curve fitting methods

Alternative phenology definitions

Support for additional vegetation indices

Contributions are welcome via GitHub Pull Requests.
For major methodological extensions, users are encouraged to contact the
package author for discussion.
