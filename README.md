libxfoil
================================================================================

libxfoil is a library to analyze the aerodynamics of 2-dimensional airfoil
sections programmatically using the Xfoil aerodynamics engine. It consists of
parts of the Xfoil source code wrapped with methods that take airfoil section
coordinates, operating points, options, and (if desired) flap settings as inputs
and calculate lift, drag, pitching moment coefficients, and transition
locations at each operating point. Functions are also included to generate 4-
and 5-digit NACA airfoils. The library is parallelized so that multiple analyses
can be performed at once. It also provides Fortran90, C, and Python bindings as
well as example programs in each language.
