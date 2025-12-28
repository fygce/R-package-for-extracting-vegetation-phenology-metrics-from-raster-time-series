# R-package-for-extracting-vegetation-phenology-metrics-from-raster-time-series

`phenoraster` is inspired by the **relative amplitude method** from Timesat. 
For each pixel, a **threshold is calculated from the full time series**, and the start of season (SOS) 
and end of season (EOS) are determined when the fitted **Double-Logistic (DL) curve** reaches a specified 
fraction of this threshold. While following the same logic, the implementation is specific to this R package.

This approach reflects the **climatological background** of each pixel, providing several advantages:

- The threshold is **robust to noise and inter-annual variability**, because it uses the full time series rather than a single year.  
- Extreme values in individual years (e.g., an unusually high summer peak) **do not distort the threshold**, avoiding bias in SOS/EOS extraction.  
- Each pixel has a **local threshold**, preserving spatial heterogeneity while maintaining consistent yearly comparisons.  
- The method is flexible and can be applied to any regularly sampled raster time series, not limited to CMIP6 GPP data.


