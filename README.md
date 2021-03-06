Brightness Temperature Collation
================================


Hourly Product Algorithm

1) Cloud Mask

Two pass cloud filtering. First pass uses an assortment of masks to remove cloudy pixels from
granules. Second pass takes average over large time window and small time window, if the small time
window average < large time window average that pixel is marked as cloudy

2) Smoothing

Cloud filtered granules are then smoothed using a large (101,101,19) window in order to fill in 
missing pixels as well as get a general idea about what is happening in region of pixels

3) Approximation

Project clear values onto smooth.  Basically generate a coefficient for each pixel and multiply the 
smooth time series at that pixel by the coefficient to adjust curve.

4) Collation

Current attempt involves averaging clear and approx values in a window (5,5,13)

Issues
--------

	-Cloud leakage obviously hurts all curves
		- lowers smooth time series
		- worse approximation (projecting cloudy pixels)
		- Cloud leakage creeps into collation averaging
	-Some granules are missing due to calibration of satellite
		-new granules must be interpolated somehow

Current Ideas
--------------

	-Apply some more masks for IR bands to try to remove more clouds
	-Adjust approx the ignore outlying pixels
	-smoother collation in order to prevent bad splotches