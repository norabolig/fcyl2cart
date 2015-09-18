# fcyl2crt

A tool for taking 2D FARGO output and converting to:

* Density Maps
* Temperature Maps

and/or

* Rad3DMC input files

The code can be scripted using command line input to process multiple files at once.
The density and temperature maps can easily be plotted and customized in **gnuplot**.

## Basic Usage

The code takes output from a 2D hydrodynamics simulation run by the code FARGO.  
The user must include the gas and temperature files for the code, as well as the parameters file (usually param.par) and a units file (units.dat). 
In addition, the user must supply the following files (these names need not be unique), with the format described in the **code book**:

* scale.par -- used to scale the units (based on the units.dat) to the desired length scale. Stellar mass is also set here.
* cart.par -- sets the basic Cartesian grid dimensions and structure.
* fname\_default.par -- a list of input and output file names (expected FARGO and user-supplied can be overrided). However, fname\_default.par itself cannot be changed.

The following options are recognized when reading fname\_default.par. Note, the "-" must be included as if it were a command line option.

*  -gas\_in  *name*      // File name for surface density input
*  -gas\_out *name*      // File name for gridded surface density output
*  -tk\_in   *name*      // File name for temperature input (2D)
*  -tk\_out  *name*      // File name for gridded temperature output
*  -units    *name*      // File name for units file from FARGO
*  -scale    *name*      // File name for scaling file used to rescale units
*  -param    *name*      // File name for parameters file from FARGO
*  -grid     *name*      // File name for grid setup for making rad3dmc files
*  -amr\_grid  *name*    // File for amr\_grid.inp in rad3dmc (don't change unless you are prepared to accept the consequences
*  -dust\_density *name* // File for dust\_density,inp in rad3dmc (don't change unless you are prepared to accept the consequences
*  -dust\_tk    *name*   // File for dust\_temparture.inp in rad3dmc (don't change unless you are prepared to accept the consequences

The command line takes precedence over the list in fname\-default.par.  This allows easy scripting for different files. 
The same options are used on the command line to change the names as used in the fname\_default.par. There are two additional 
command line options:

* -h                 // Display options.
* -product *selection*  // Define products of program. Options include 'xy', 'rad3dmc', or 'both' (default).

## Assumptions for 3D construction.

* Only a static grid is used at the moment. This may be re-evaluated to include static mesh refinement for input into Rad3DMC.  
* The vertical direction is included assumption an isothermal vertical temperature structure.  This may be relaxed to include
other structural relations at a later time. 
* The 2D interpolation is done to conserve mass and energy density (for fixed adiabatic indices).  No attempt at conserving mass
is made for the extrapolation to 3D.  This may be re-evaluated.  Differences will be small as long as the structure is resolved. 
* Units are given according to the unit conversion given by the user.

## Density and Temperature Maps

The output for both files have three columns: x, y, and surface density (temperature).  The x and y values are in scaled code units,
while the density or temperature is given in the desired user-defined units (e.g., MKS). 
The images are easily plotted in gnuplot.  Minimum code to produce an image is as follows (assuming you are in gnuplot with the file in your working directory):

>gnuplot> plot "gasoutNNN.dat" w image

## Vertical Slices

The Extrapolate3D class has a method for producing gnuplot-friendly image tables for given x or y slices in the resulting 3D extrapolation. An example of producing a slice is commented out in the source code fcyl2crt.cpp. Uncommenting will produce a slice, and an option for making this more accessible to the user will be done shortly. 

