libxfoil
================================================================================

libxfoil is a library to analyze the aerodynamics of 2-dimensional airfoil
sections programmatically using the Xfoil aerodynamics engine. The main focus
is on supporting the functionality of Xfoil's OPER direct analysis mode, but
some geometry functions including spline interpolation, flap deflections,
trailing edge gap modifications, and NACA 4- and 5-digit airfoil generation
are also included. In addition to integrated forces and moments, pressure
distributions and boundary layer variables (including skin friction coefficient,
displacement thickness, log of amplification ratio, etc.) can also be accessed.

The API employs a pseudo object-oriented approach; global variables in the 
original Xfoil code have been moved to a struct, allowing multiple instances of
the Xfoil engine, representing different airfoils, to be saved in memory 
simultaneously and executed in parallel. Memory for each instance is managed
dynamically through the xfoil_init and xfoil_cleanup methods. libxfoil provides
C, Fortran, and Python bindings as well as example programs in each language. 
