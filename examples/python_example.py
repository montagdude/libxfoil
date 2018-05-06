#!/usr/bin/env python

import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type
from matplotlib import pyplot as plt

def plot_spline(x, z, xspline, zspline, name):

  plt.clf()
  plt.close()

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

  plt.clf()
  plt.close()
  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.plot(x, z)
  ax.plot(xnew, znew)
  ax.grid()
  ax.set_aspect('equal', 'datalim')
  ax.legend(['Original', 'Modified'])

  plt.show()

def plot_polars(alpha, cl, cd, cm, xtrt, xtrb, title=None):

  plt.clf()
  plt.close()

  fig, axarr = plt.subplots(2,2)
  if title is not None:
    fig.suptitle(title)

  axarr[0,0].set_xlabel('Angle of attack')
  axarr[0,0].set_ylabel('Lift coefficient')
  axarr[0,0].plot(alpha, cl)
  axarr[0,0].grid()

  axarr[0,1].set_xlabel('Drag coefficient')
  axarr[0,1].set_ylabel('Lift coefficient')
  axarr[0,1].plot(cd, cl)
  axarr[0,1].grid()

  axarr[1,0].set_xlabel('Angle of attack')
  axarr[1,0].set_ylabel('Pitching moment coefficient')
  axarr[1,0].plot(alpha, cm)
  axarr[1,0].grid()

  axarr[1,1].set_xlabel('Angle of attack')
  axarr[1,1].set_ylabel('Transition locations')
  axarr[1,1].plot(alpha, xtrt)
  axarr[1,1].plot(alpha, xtrb)
  axarr[1,1].grid()
  axarr[1,1].legend(['Top', 'Bottom'], loc='center left')

  plt.show()

if __name__ == "__main__":

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

  noppoint = 10
  oppoints = [-0.5, -0.25, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
  opmodes = noppoint*[1]
  re = noppoint*[3.E+06]
  mach = noppoint*[0.4]
  use_flap = False
  flapang = noppoint*[0]

  opts = xfoil_options_type()
  opts.ncrit = 9.
  opts.xtript = 1.
  opts.xtripb = 1.
  opts.viscous_mode = True
  opts.silent_mode = True
  opts.maxit = 100
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

  # Test modifying TE gap
  xiw.xfoil_init()
  xiw.xfoil_defaults(opts)
  xiw.xfoil_set_airfoil(x, z, npoint)
  npointnew = xiw.xfoil_modify_tegap(0., 0.9) 
  xnew, znew = xiw.xfoil_get_airfoil(npointnew)
  plot_oldnew(x, z, xnew, znew)

  # Test applying a flap deflection
  x_flap = 0.7
  y_flap = 0.0
  y_flap_spec = 0
  xiw.xfoil_set_paneling(geom_opts)
  xiw.xfoil_smooth_paneling()
  npointnew = xiw.xfoil_apply_flap_deflection(x_flap, y_flap, y_flap_spec, 10.)
  xnew, znew = xiw.xfoil_get_airfoil(npointnew)
  plot_oldnew(x, z, xnew, znew)

  print("Calculating aerodynamics with Xfoil (no flap) ...")
  lift, drag, moment, viscrms, alpha, xtrt, xtrb = \
    xiw.run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,
                  mach, use_flap, 0., 0., 0, flapang, opts)
  plot_polars(alpha, lift, drag, moment, xtrt, xtrb, 'No flap deflection')

  print("Calculating aerodynamics with Xfoil (with flap) ...")
  flapang = [0., 0., 0., 0., 0., 0., 3., 6., 9., 12., 15.]
  use_flap = True
  lift, drag, moment, viscrms, alpha, xtrt, xtrb = \
    xiw.run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,
                  mach, use_flap, x_flap, y_flap, y_flap_spec, flapang, opts)
  plot_polars(alpha, lift, drag, moment, xtrt, xtrb, 'With flap deflection')

  xiw.xfoil_cleanup()
