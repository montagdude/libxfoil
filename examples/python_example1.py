#!/usr/bin/env python

import sys
import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type, \
                            xfoil_data_group
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

def plot_surface_var(x_noflap, var_noflap, x_withflap, var_withflap, varname,
                     reverse_y=False):

  fig, ax = plt.subplots()
  ax.set_xlabel('x/c')
  ax.set_ylabel(varname)
  ax.plot(x_noflap, var_noflap)
  ax.plot(x_withflap, var_withflap)
  if reverse_y:
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[1], ylim[0])
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

  geom_opts = xfoil_geom_options_type()
  geom_opts.npan = 160
  geom_opts.cvpar = 1.
  geom_opts.cterat = 0.15
  geom_opts.ctrrat = 0.2
  geom_opts.xsref1 = 1.
  geom_opts.xsref2 = 1.
  geom_opts.xpref1 = 1.
  geom_opts.xpref2 = 1.

  # Xfoil data
  xdg = xfoil_data_group()

  # Modify the trailing edge gap
  xiw.xfoil_init(xdg)
  xiw.xfoil_defaults(xdg, opts)
  xiw.xfoil_set_buffer_airfoil(xdg, x, z, npoint)
  xiw.xfoil_set_paneling(xdg, geom_opts)
  if (xiw.xfoil_smooth_paneling(xdg) != 0):
    print("Error smoothing paneling: xfoil_set_buffer_airfoil must be called " +
          "first.")
    sys.exit(1)
  npointnew, stat = xiw.xfoil_modify_tegap(xdg, 0., 0.9)
  if (stat != 0):
    print("Error modifying TE gap: xfoil_set_buffer_airfoil must be called " +
          "first.")
    sys.exit(1)
  x_noflap, z_noflap, stat = xiw.xfoil_get_current_airfoil(xdg, npointnew)
  if (stat != 0):
    print("Error getting airfoil: xfoil_smooth_paneling must be called first.")
    sys.exit(1)
  plot_oldnew(x, z, x_noflap, z_noflap)

  # Set operating point
  re = 1.E+05
  mach = 0.1
  xiw.xfoil_set_reynolds_number(xdg, re)
  xiw.xfoil_set_mach_number(xdg, mach)

  # Run xfoil for current airfoil
  print("Running Xfoil without flap...")
  cl_spec = 1.0
  alpha, cl, cd, cm, converged, stat = xiw.xfoil_speccl(xdg, cl_spec)
  if (stat == 0):
    if not converged:
      print("Warning: Xfoil calculations did not converge.")
    print("Alpha: {:.4f}, Cl: {:.4f}, Cd: {:.4f}, Cm: {:.4f}"\
          .format(alpha, cl, cd, cm))
  elif stat == 1:
    print("Error running xfoil: xfoil_smooth_paneling must be called first.")
    sys.exit(1)

  # Get surface cp, cf, transition location, ampl. ratio
  cp_noflap = xiw.xfoil_get_cp(xdg, npointnew)
  cf_noflap = xiw.xfoil_get_cf(xdg, npointnew)
  xtranst_noflap, _, xtransb_noflap, _ = xiw.xfoil_get_transloc(xdg)
  ampl_noflap = xiw.xfoil_get_ampl(xdg, npointnew)

  # Apply a flap deflection
  x_flap = 0.7
  y_flap = 0.0
  y_flap_spec = 0
  if (xiw.xfoil_smooth_paneling(xdg) != 0):
    print("Error smoothing paneling: xfoil_set_buffer_airfoil must be called " +
          "first.")
    sys.exit(1)
  npointnew, stat = xiw.xfoil_apply_flap_deflection(xdg, x_flap, y_flap,
                                                    y_flap_spec, 10.)
  if (stat != 0):
    print("Error applying flap deflection: " +
          "xfoil_set_buffer_airfoil must be called first.")
    sys.exit(1)
  x_withflap, z_withflap, stat = xiw.xfoil_get_current_airfoil(xdg, npointnew)
  if (stat != 0):
    print("Error getting airfoil: xfoil_smooth_paneling must be called first.")
    sys.exit(1)
  plot_oldnew(x, z, x_withflap, z_withflap)

  # Run xfoil again with flap
  print("Running Xfoil with flap...")
  xiw.xfoil_reinitialize_bl(xdg)
  alpha, cl, cd, cm, converged, stat = xiw.xfoil_speccl(xdg, cl_spec)
  if (stat == 0):
    if not converged:
      print("Warning: Xfoil calculations did not converge.")
    print("Alpha: {:.4f}, Cl: {:.4f}, Cd: {:.4f}, Cm: {:.4f}"\
          .format(alpha, cl, cd, cm))
  elif stat == 1:
    print("Error running xfoil: xfoil_smooth_paneling must be called first.")
    sys.exit(1)

  # Get surface cp, cf, transition location, ampl. ratio
  cp_withflap = xiw.xfoil_get_cp(xdg, npointnew)
  cf_withflap = xiw.xfoil_get_cf(xdg, npointnew)
  xtranst_withflap, _, xtransb_withflap, _ = xiw.xfoil_get_transloc(xdg)
  ampl_withflap = xiw.xfoil_get_ampl(xdg, npointnew)
  xiw.xfoil_cleanup(xdg)

  print("Transition locations:")
  print("No flap:   xtranstop = {:.4f}, xtransbot = {:.4f}"\
        .format(xtranst_noflap, xtransb_noflap))
  print("With flap: xtranstop = {:.4f}, xtransbot = {:.4f}"\
        .format(xtranst_withflap, xtransb_withflap))

  # Plot cp, cf, and amplification ratio
  plot_surface_var(x_noflap, cp_noflap, x_withflap, cp_withflap,
                   "Pressure coefficient", reverse_y=True)
  plot_surface_var(x_noflap, cf_noflap, x_withflap, cf_withflap,
                   "Skin friction coefficient")
  plot_surface_var(x_noflap, ampl_noflap, x_withflap, ampl_withflap,
                   "Amplitude ratio")
