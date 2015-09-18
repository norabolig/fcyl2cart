## -gas_in *file name* : For example, gasDensity.dat

Binary data file containing surface density information from FARGO simulation

## -tk_in *file name* : For example, gasTemperature.dat

Binary data file containing 2D temperature information from FARGO simulation

## -units *file name* : For example, units.dat

Ascii file produced by FARGO containing units of the simulation.  Can still be
scaled by user with the file defined by -scale. 

The file lists the following values for 1 code unit using a single line in the
following order: *mass length time temperature

## -param *file name* : For example, param.par

Ascii file produced by FARGO containing parameters for the simulation.
Variables are defined in the file itself.  See also FARGO documentation. (Yes,
I know this is bad, but there is a lot to list). 

## -scale *file name* : For example, scale.par

Ascii file supplied by user.  Contains the mass of the star (in code units),
the scaling for length, the mean molecular weight (in proton masses), the
dust to gas ratio, and a flag where '1' means to convert from MKS to cgs 
(FARGO units.dat files are in MKS, but cgs is needed for Rad3DMC). Set the flag
to '0' if no conversion is desired -- this does not affect the scaling for
length.

Values should be written on one line:
* 1.0 100. 2.35 0.01 1

## -grid *file name* : FOr example, cart.par

Ascii file supplied by user.  Contains grid information for creating 2D images
and the full 3D extrapolation to a data cube for producing Rad3DMC files.  Even
if only 2D maps are desired, the code must include values for the full 3D case.
When in doubt, make it a box. 

Values are given in a single line, giving the number of cells in x y z (nx, ny,
nz), the cell width in each direction (dx, dy, dz), and the number of
subdivisions in the FARGO cylindrical cells in radius (nrad) and azimuth
(nsec).  This is used for integrating the cylindrical cell mass into each
Cartesian cell. The higher the number, the more accurate the integration, but
it becomes much more costly.  16 and 16 are good values to use unless the 2D
maps show discreteness effects.  For example:

* nx ny nz dx dy dz nrad nsec

## -gas_out *file name* : For example, gout.dat

Produces surface density gnuplot image table.  Three columns are produced: x,
y, and surface density.  x and y are in scaled code units (scaling used in
-scale file), while the surface density is in the user-defined units (e.g.,
MKS). 

## -tk_out *file name*  : For example, tout.dat


## -amr_grid, -dust_density, -dust_tk

These are the gridding structure, the dust density, and the dust
temperature (assumed to be gas temperature).  Files are produced by 3D
extrapolation of interpolated 2D grid. All units are in physical units and must
cgs for rad3dmc. Please see Rad3DMC documentation for details of these files
(sorry!).
