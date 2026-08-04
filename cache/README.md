# Generated cache

This directory is intentionally empty in version control. The first run creates the required `Prepared_settings*` subdirectories and regenerates node coordinates, neighbor tables, B-spline bases, and covariance data when a required file is absent.

The generated `.dat` files are formatted Fortran data and can be very large. Keep them locally to accelerate later runs; do not commit them. A cache is tied to the grid dimensions, neighbor counts, kernel bandwidths, and other parameters that created it. Delete the affected cache directory after changing those parameters.

