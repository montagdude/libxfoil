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
 
   Copyright (C) 2018 Daniel Prosser

   See xfoil_interface.f90 for descriptions of inputs and outputs. */

#pragma once

#include <stdbool.h>

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

extern void xfoil_init(void);
extern void xfoil_defaults(const xfoil_options_type *xfoil_options);
extern void xfoil_set_paneling(const xfoil_geom_options_type *geom_opts);
extern void xfoil_set_airfoil(const double xin[], const double zin[],
                              const int *npointin);
extern void xfoil_smooth_paneling(void);
extern void xfoil_apply_flap_deflection(const double *xflap,
                                        const double *yflap,
                                        const int *y_flap_spec,
                                        const double *degrees, int *npointout);
extern void xfoil_modify_tegap(const double *gap, const double *blendloc,
                               int *npointout);
extern void xfoil_get_airfoil(double xout[], double zout[], const int *npoint);
extern void xfoil_geometry_info(double *maxt, double *xmaxt, double *maxc,
                                double *xmaxc);
extern void xfoil_specal(const double *alpha_spec, double *alpha, double *lift,
                         double *drag, double *moment);
extern void xfoil_speccl(const double *cl_spec, double *alpha, double *lift,
                         double *drag, double *moment);
extern void xfoil_cleanup(void);

/* The following methods utilize xfoil functionality, performing some
   calculations and returning a result, without needing to initialize xfoil,
   set a buffer airfoil, or smooth paneling first. */
extern void naca_4_digit(const char des[4], const int *npointside,
                         double xout[], double zout[], int *nout);
extern void naca_5_digit(const char des[5], const int *npointside,
                         double xout[], double zout[], int *nout, int *stat);
extern void xfoil_spline_coordinates(const double x[], const double z[],
                                     const int *npt, double s[], double xs[],
                                     double zs[]);
extern void xfoil_eval_spline(const double x[], const double z[],
                              const double s[], const double xs[],
                              const double zs[], const int *npt,
                              const double *sc, double *xc, double *zc);
extern void xfoil_lefind(const double x[], const double z[], const double s[],
                         const double xs[], const double zs[], const int *npt,
                         double *sle, double *xle, double *zle);

/* Convenience method to run xfoil at a bunch of different operating points,
   optionally changing flap deflections and ncrit values */
extern void run_xfoil(const int *npointin, const double xin[],
                      const double zin[],
                      const xfoil_geom_options_type *geom_options,
                      const int *noppoint, const double operating_points[],
                      const int op_modes[], const double reynolds_numbers[],
                      const double mach_numbers[], const bool *use_flap,
                      const double *x_flap, const double *y_flap,
                      const int *y_flap_spec, const double flap_degrees[],
                      const xfoil_options_type *xfoil_options, double lift[],
                      double drag[], double moment[], double viscrms[],
                      double alpha[], double xtrt[], double xtrb[],
                      const double ncrit_per_point[]);
