#!/usr/bin/env python

import sys
import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type
from matplotlib import pyplot as plt

def plot_spline(x, z, xspline, zspline, name):

  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.set_title(name)
  ax.plot(x, z)
  ax.plot(xspline, zspline, 'o')
  ax.grid()
  ax.set_aspect('equal', 'datalim')
  ax.legend(['Airfoil', 'Some spline interp points'])

  plt.show() 
  plt.clf()
  plt.close('all')

def plot_oldnew(x, z, xnew, znew):

  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.plot(x, z)
  ax.plot(xnew, znew)
  ax.grid()
  ax.set_aspect('equal', 'datalim')
  ax.legend(['Original', 'Modified'])

  plt.show()

def plot_cps(x_noflap, cp_noflap, x_withflap, cp_withflap):

  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('cp')
  ax.plot(x_noflap, cp_noflap)
  ax.plot(x_withflap, cp_withflap)
  ax.set_ylim(1.1, -2.7)
  ax.grid()
  ax.legend(['No flap', 'With flap'])

  plt.show()

if __name__ == "__main__":

  # Generate an airfoil
  digits = '2312'
  npointside = 100
  x, z, npoint = xiw.naca_4_digit('2312', npointside)

  # Test spline fitting and eval
  s, xs, zs = xiw.xfoil_spline_coordinates(x, z, npoint)
  smax = s[len(s)-1]
  x0, z0 = xiw.xfoil_eval_spline(x, z, s, xs, zs, npoint, 0.0)
  x14, z14 = xiw.xfoil_eval_spline(x, z, s, xs, zs, npoint, 0.25*smax)
  sle, xle, zle = xiw.xfoil_lefind(x, z, s, xs, zs, npoint)
  x34, z34 = xiw.xfoil_eval_spline(x, z, s, xs, zs, npoint, 0.75*smax)
  x1, z1 = xiw.xfoil_eval_spline(x, z, s, xs, zs, npoint, smax)

  plot_spline(x, z, [x0, x14, xle, x34, x1], [z0, z14, zle, z34, z1],
              'NACA ' + digits)

  # Xfoil settings
  opts = xfoil_options_type()
  opts.ncrit = 9.
  opts.xtript = 1.
  opts.xtripb = 1.
  opts.viscous_mode = True
  opts.silent_mode = True
  opts.maxit = 200
  opts.vaccel = 0.01
  opts.fix_unconverged = True
  opts.reinitialize = False

  geom_opts = xfoil_geom_options_type()
  geom_opts.npan = 160
  geom_opts.cvpar = 1.
  geom_opts.cterat = 0.15
  geom_opts.ctrrat = 0.2
  geom_opts.xsref1 = 1.
  geom_opts.xsref2 = 1.
  geom_opts.xpref1 = 1.
  geom_opts.xpref2 = 1.

  # Modify the trailing edge gap
  xiw.xfoil_init()
  xiw.xfoil_defaults(opts)
  if (xiw.xfoil_set_airfoil(x, z, npoint) != 0):
    print("Error setting airfoil: xfoil_init must be called first.")
    sys.exit(1)
  xiw.xfoil_set_paneling(geom_opts)
  if (xiw.xfoil_smooth_paneling() != 0):
    print("Error smoothing paneling: xfoil_set_airfoil must be called first.")
    sys.exit(1)
  npointnew, stat = xiw.xfoil_modify_tegap(0., 0.9) 
  if (stat != 0):
    print("Error modifying TE gap: xfoil_set_airfoil must be called first.")
    sys.exit(1)
  x_noflap, z_noflap = xiw.xfoil_get_airfoil(npointnew)
  plot_oldnew(x, z, x_noflap, z_noflap)

  # Set operating point
  re = 1.E+05
  mach = 0.1
  xiw.xfoil_set_reynolds_number(re)
  xiw.xfoil_set_mach_number(mach)

  # Run xfoil for current airfoil
  print("Running Xfoil without flap...")
  cl_spec = 1.0
  alpha, cl, cd, cm, converged, stat = xiw.xfoil_speccl(cl_spec)
  if (stat == 0):
    if not converged:
      print("Warning: Xfoil calculations did not converge.")
    print("Alpha: {:.4f}, Cl: {:.4f}, Cd: {:.4f}, Cm: {:.4f}"\
          .format(alpha, cl, cd, cm))
  elif stat == 1:
    print("Error running xfoil: xfoil_init must be called first.")
    sys.exit(1)

  # Get surface cp
  cp_noflap = xiw.xfoil_get_cp(npointnew)
   
  # Apply a flap deflection
  x_flap = 0.7
  y_flap = 0.0
  y_flap_spec = 0
  xiw.xfoil_set_paneling(geom_opts)
  if (xiw.xfoil_smooth_paneling() != 0):
    print("Error smoothing paneling: xfoil_set_airfoil must be called first.")
    sys.exit(1)
  npointnew, stat = xiw.xfoil_apply_flap_deflection(x_flap, y_flap,
                                                    y_flap_spec, 10.)
  if (stat != 0):
    print("Error applying flap deflection: " +
          "xfoil_set_airfoil must be called first.")
    sys.exit(1)
  x_withflap, z_withflap = xiw.xfoil_get_airfoil(npointnew)
  plot_oldnew(x, z, x_withflap, z_withflap)

  # Run xfoil again with flap
  print("Running Xfoil with flap...")
  xiw.xfoil_reinitialize_bl()
  alpha, cl, cd, cm, converged, stat = xiw.xfoil_speccl(cl_spec)
  if (stat == 0):
    if not converged:
      print("Warning: Xfoil calculations did not converge.")
    print("Alpha: {:.4f}, Cl: {:.4f}, Cd: {:.4f}, Cm: {:.4f}"\
          .format(alpha, cl, cd, cm))
  elif stat == 1:
    print("Error running xfoil: xfoil_init must be called first.")
    sys.exit(1)
   
  # Get surface cp
  cp_withflap = xiw.xfoil_get_cp(npointnew)

  # Plot cps
  plot_cps(x_noflap, cp_noflap, x_withflap, cp_withflap)
   
  xiw.xfoil_cleanup()
