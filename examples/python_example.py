#!/usr/bin/env python

import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type
from matplotlib import pyplot as plt

def plot_airfoil(x, z, name):

  plt.clf()
  plt.close()

  fig, ax = plt.subplots()
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  ax.set_title(name)
  ax.plot(x, z)
  ax.grid()
  ax.set_aspect('equal', 'datalim')

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
  x, z, npoint = xiw.naca_4_digit('2312', 100)
  plot_airfoil(x, z, 'NACA ' + digits)

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
  geom_opts.npan = 200
  geom_opts.cvpar = 1.
  geom_opts.cterat = 0.15
  geom_opts.ctrrat = 0.2
  geom_opts.xsref1 = 1.
  geom_opts.xsref2 = 1.
  geom_opts.xpref1 = 1.
  geom_opts.xpref2 = 1.

  xiw.xfoil_init()

  print("Calculating aerodynamics with Xfoil (no flap) ...")
  lift, drag, moment, viscrms, alpha, xtrt, xtrb = \
    xiw.run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,
                  mach, use_flap, 0., 0., 0, flapang, opts)
  plot_polars(alpha, lift, drag, moment, xtrt, xtrb, 'No flap deflection')

  print("Calculating aerodynamics with Xfoil (with flap) ...")
  flapang = [0., 0., 0., 0., 0., 0., 3., 6., 9., 12., 15.]
  use_flap = True
  x_flap = 0.7
  y_flap = 0.0
  y_flap_spec = 0
  lift, drag, moment, viscrms, alpha, xtrt, xtrb = \
    xiw.run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,
                  mach, use_flap, x_flap, y_flap, y_flap_spec, flapang, opts)
  plot_polars(alpha, lift, drag, moment, xtrt, xtrb, 'With flap deflection')

  xiw.xfoil_cleanup()
