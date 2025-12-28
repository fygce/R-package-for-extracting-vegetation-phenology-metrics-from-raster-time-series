# R-package-for-extracting-vegetation-phenology-metrics-from-raster-time-series

`phenoraster` is inspired by the **relative amplitude method** from TIMESAT3.3. 
For each pixel, a **threshold is calculated from the full time series**, and the start of season (SOS) 
and end of season (EOS) are determined when the fitted **Double-Logistic (DL) curve** reaches a specified 
fraction of this threshold. While following the same logic, the implementation is specific to this R package.

## Installation

install.packages("remotes")
remotes::install_github("fygce/phenoraster")

## Note

This approach reflects the **climatological background** of each pixel, providing several advantages:

- The threshold is **robust to noise and inter-annual variability**, because it uses the full time series rather than a single year.  
- Extreme values in individual years (e.g., an unusually high summer peak) **do not distort the threshold**, avoiding bias in SOS/EOS extraction.  
- Each pixel has a **local threshold**, preserving spatial heterogeneity while maintaining consistent yearly comparisons.  
- Because the threshold is based on the long-term climatology, **phenological metrics across different years are directly comparable**, enabling meaningful inter-annual analysis.  
- The method is flexible and can be applied to any regularly sampled raster time series, not limited to CMIP6 GPP data.

## Development Status and Extensibility

`phenoraster` is primarily developed for extracting vegetation phenology metrics from **raster time series**. 
While the current implementation uses the **Double-Logistic (DL) fitting method** and serial processing, 
the package is **not limited to this approach** and is actively being improved. 
Future versions may include **parallel processing** and support for **additional fitting methods**, 
which can be contributed by both the author and other users.

## Contributions are welcome! 
Please submit changes via GitHub Pull Requests. 
For major changes or new fitting methods, you may contact the package author for discussion.
