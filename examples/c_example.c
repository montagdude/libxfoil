#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "xfoil_interface.h"

int main()
{
  xfoil_options_type opts;
  xfoil_geom_options_type geom_opts;
  int nhalf = 100;
  int npoint = nhalf*2-1;
  double x[npoint], z[npoint];
  int noppoint = 5;
  double oppoints[5] = {-0.3, 0.0, 0.3, 0.6, 0.9};
  int opmodes[5] = {1, 1, 1, 1, 1};
  double re[5] = {1.E+05, 1.E+05, 1.E+05, 1.E+05, 1.E+05};
  double mach[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  bool use_flap = false;
  bool reinitialize = false;
  bool fix_unconverged = false;
  bool converged;
  double lift[noppoint], drag[noppoint], moment[noppoint], viscrms[noppoint];
  double alpha[noppoint], xtrt[noppoint], xtrb[noppoint];
  int i, stat;
  xfoil_data_group xdg;

  opts.ncrit = 9.;
  opts.xtript = 1.;
  opts.xtripb = 1.;
  opts.viscous_mode = true;
  opts.silent_mode = true;
  opts.maxit = 100;
  opts.vaccel = 0.01;

  geom_opts.npan = 200;
  geom_opts.cvpar = 1.0;
  geom_opts.cterat = 0.15;
  geom_opts.ctrrat = 0.2;
  geom_opts.xsref1 = 1.;
  geom_opts.xsref2 = 1.;
  geom_opts.xpref1 = 1.;
  geom_opts.xpref2 = 1.;

  // Generate an airfoil
  naca_5_digit("25012", &nhalf, x, z, &npoint, &stat);
  if (stat != 0) { return stat; }

  // Set up inputs
  xfoil_defaults(&xdg, &opts);
  xfoil_set_paneling(&xdg, &geom_opts);
  xfoil_set_airfoil(&xdg, x, z, &npoint);
  xfoil_smooth_paneling(&xdg, &stat);
  if (stat != 0) { return stat; }

  // Run it at all operating points
  printf("Running Xfoil ...\n");
  for ( i = 0; i < noppoint; i++ )
  {
    xfoil_set_reynolds_number(&xdg, &re[i]);
    xfoil_set_mach_number(&xdg, &mach[i]);
    xfoil_speccl(&xdg, &oppoints[i], &alpha[i], &lift[i], &drag[i], &moment[i],
                 &converged, &stat);
    if (stat != 0) { return stat; }
    if (! converged)
      printf("Point %d did not converge.\n", i);
    else
    {
      printf("Point %d converged.\n", i+1);
      printf("Point %d: AoA = %.4f, Cl = %.4f, Cd = %.4f, Cm = %.4f\n", i+1,
             alpha[i], lift[i], drag[i], moment[i]);
    }
  }

  // Now, do the same thing but with the run_xfoil method
  printf("\nRunning Xfoil using the run_xfoil method ...\n");
  run_xfoil(&npoint, x, z, &geom_opts, &noppoint, oppoints, opmodes, re, mach,
            &use_flap, NULL, NULL, NULL, NULL, &opts, &reinitialize,
            &fix_unconverged, lift, drag, moment, viscrms, alpha, xtrt, xtrb,
            &stat);
  if (stat != 0) { return stat; }

  for ( i = 0; i < noppoint; i++ )
  {
    printf("Point %d: AoA = %.4f, Cl = %.4f, Cd = %.4f, Cm = %.4f\n", i+1,
           alpha[i], lift[i], drag[i], moment[i]);
  }

  return 0;
}
