/* This file is part of libxfoil.
 
   libxfoil is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   libxfoil is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.
 
   Copyright (C) 2018 Daniel Prosser */

%module xfoil_interface
%include "cstring.i"
%include "cpointer.i"
%include "carrays.i"
%{
#include "xfoil_interface.h"
%}

/* Type handling */
%cstring_bounded_mutable(char *des, 1024);
%pointer_functions(int, intp);
%pointer_functions(double, doublep);
%pointer_functions(bool, boolp);
%array_functions(double, doublea);
%array_functions(int, inta);

/* Structs */
typedef struct
{
  double ncrit;
  double xtript, xtripb;
  bool viscous_mode;
  bool silent_mode;
  int maxit;
  double vaccel;
  bool fix_unconverged;
  bool reinitialize;
} xfoil_options_type;

typedef struct
{
  int npan;
  double cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2;
} xfoil_geom_options_type;

/* Interface functions */
extern void xfoil_init(void);
extern void xfoil_defaults(xfoil_options_type *xfoil_options);
extern void xfoil_set_paneling(xfoil_geom_options_type *geom_opts);
extern void xfoil_set_airfoil(double *xin, double *zin, int *npointin);
extern void xfoil_smooth_paneling(void);
extern void xfoil_apply_flap_deflection(double *xflap, double *yflap,
                                        int *y_flap_spec, double *degrees,
                                        int *npointout);
extern void xfoil_modify_tegap(double *gap, double *blendloc, int *npointout);
extern void xfoil_get_airfoil(double *xout, double *zout, int *npoint);
extern void xfoil_geometry_info(double *maxt, double *xmaxt, double *maxc,
                                double *xmaxc);
extern void xfoil_specal(double *alpha_spec, double *alpha, double *lift,
                         double *drag, double *moment);
extern void xfoil_speccl(double *cl_spec, double *alpha, double *lift,
                         double *drag, double *moment);
extern void xfoil_cleanup(void);
extern void naca_4_digit(char *des, int *npointside, double *xout,
                         double *zout, int *nout);
extern void naca_5_digit(char *des, int *npointside, double *xout,
                         double *zout, int *nout, int *stat);
extern void xfoil_spline_coordinates(double *x, double *z, int *npt, double *s,
                                     double *xs, double *zs);
extern void xfoil_eval_spline(double *x, double *z, double *s, double *xs,
                              double *zs, int *npt, double *sc, double *xc,
                              double *zc);
extern void xfoil_lefind(double *x, double *z, double *s, double *xs,
                         double *zs, int *npt, double *sle, double *xle,
                         double *zle);
extern void run_xfoil(int *npointin, double *xin, double *zin,
                      xfoil_geom_options_type *geom_options, int *noppoint,
                      double *operating_points, int *op_modes,
                      double *reynolds_numbers, double *mach_numbers,
                      bool *use_flap, double *x_flap, double *y_flap,
                      int *y_flap_spec, double *flap_degrees,
                      xfoil_options_type *xfoil_options, double *lift,
                      double *drag, double *moment, double *viscrms,
                      double *alpha, double *xtrt, double *xtrb,
                      double *ncript_per_point=NULL);
