#!/usr/bin/env python

import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type
from matplotlib import pyplot as plt
import sys

def read_airfoil(filename):

  try:
    f = open(filename)
  except IOError:
    print("Error: unable to open {:s}.".format(filename))
    sys.exit(1)

  x = []
  z = []
  npoint = 0
  for line in f:
    x.append(float(line.split()[0]))
    z.append(float(line.split()[1]))
    npoint += 1

  return x, z, npoint

if __name__ == "__main__":

  x, z, npoint = read_airfoil("clarky.dat")

  noppoint = 8
  oppoints = [-6., -3., 0., 3., 6., 9., 12., 15.]
  opmodes = noppoint*[0]
  re = noppoint*[3.E+05]
  mach = noppoint*[0.3]
  flapang = noppoint*[0.]
  use_flap = False

  opts = xfoil_options_type()
  opts.ncrit = 9.
  opts.xtript = 1.
  opts.xtripb = 1.
  opts.viscous_mode = True
  opts.silent_mode = False
  opts.maxit = 100
  opts.vaccel = 0.01
  opts.fix_unconverged = False
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
  lift, drag, moment, viscrms, alpha, xtrt, xtrb = \
    xiw.run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,
                  mach, use_flap, 0., 0., 0, flapang, opts)

  for i in range(noppoint):
    print("Point {:d}: AoA = {:<.4f}, Cl = {:<.4f}, Cd = {:<.4f}, \
Cm = {:6.4f}".format(i, alpha[i], lift[i], drag[i], moment[i]))

  xiw.xfoil_cleanup()
