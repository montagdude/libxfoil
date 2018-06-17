#!/usr/bin/env python

import sys
import xfoil_interface_wrap as xiw
from xfoil_interface import xfoil_options_type, xfoil_geom_options_type, \
        xfoil_data_group
from matplotlib import pyplot as plt
from math import sqrt

def read_airfoil(fname):

    try:
        f = open(fname)
    except IOError:
        print("Unable to open {}.".format(fname))
        sys.exit(1)

    x = []
    z = []
    for line in f:
        if line.startswith("MH"):
            continue

        x.append(float(line.split()[0]))
        z.append(float(line.split()[1]))

    f.close()

    return x, z, len(x)

def centered_difference(yp1, y0, ym1, xp1, xm1):

    a = (yp1-y0 - (ym1-y0)*xp1/xm1) / (xp1*xp1 - xp1*xm1);
    b = (ym1-y0 - a*xm1*xm1) / xm1;

    return b

def surface_derivative(var, x, z):

    n = len(x)
    s = [0.]*n
    for i in range(1,n):
        dx = x[i] - x[i-1]
        dz = z[i] - z[i-1]
        ds = sqrt(dx**2. + dz**2.)
        s[i] = s[i-1] + ds

    dvar = [0.]*n
    ds = s[1] - s[0]
    dvar[0] = (var[1] - var[0])/ds
    for i in range(1,n-1):
        splus = s[i+1] - s[i]
        sminus = s[i-1] - s[i]
        dvar[i] = centered_difference(var[i+1], var[i], var[i-1], splus, sminus)
    ds = s[n-1] - s[n-2]
    dvar[n-1] = (var[n-1]-var[n-2])/ds

    return dvar

def plot_geometry(x, z, xwake, zwake, title=None):

    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    if title is not None:
        ax.set_title(title)
    ax.plot(x, z)
    ax.plot(xwake, zwake, marker='.')
    ax.set_aspect('equal', 'datalim')
    ax.grid()

    plt.show()

def plot(x, var, varname):

    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel(varname)
    ax.plot(x, var)
    ax.grid()

    plt.show()

def plot_with_wake(x, xwake, var, varwake, varname):

    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel(varname)
    ax.plot(x, var)
    ax.plot(xwake, varwake)
    ax.grid()
    ax.legend(['Surface', 'Wake'])

    plt.show()

if __name__ == "__main__":

    x, z, npoint = read_airfoil("mh32_smoothed.dat")

    opts = xfoil_options_type()
    opts.ncrit = 9.
    opts.xtript = 1.
    opts.xtripb = 1.
    opts.viscous_mode = True
    opts.silent_mode = True
    opts.maxit = 100
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

    xiw.xfoil_init(xdg)
    xiw.xfoil_defaults(xdg, opts)
    xiw.xfoil_set_buffer_airfoil(xdg, x, z, npoint)
    xiw.xfoil_set_paneling(xdg, geom_opts)
    if (xiw.xfoil_smooth_paneling(xdg) != 0):
        print("Err 1")
        sys.exit(1)
    xnew, znew, stat = xiw.xfoil_get_current_airfoil(xdg, geom_opts.npan)
    if (stat != 0):
        print("Err 2")
        sys.exit(1)

    xiw.xfoil_set_reynolds_number(xdg, 127000.)
    xiw.xfoil_set_mach_number(xdg, 0.0)
    alpha, cl, cd, cm, converged, stat = xiw.xfoil_speccl(xdg, 0.72)
    if (stat != 0):
        print("Err 3")
        sys.exit(1)

    cp = xiw.xfoil_get_cp(xdg, geom_opts.npan)
    uedge = xiw.xfoil_get_uedge(xdg, geom_opts.npan)
    dstar = xiw.xfoil_get_deltastar(xdg, geom_opts.npan)
    mass = []
    for i in range(geom_opts.npan):
        mass.append(uedge[i]*dstar[i])

    nwake = xiw.xfoil_get_wakepoints(xdg)
    xw, zw = xiw.xfoil_get_wake_geometry(xdg, nwake)
    cpw = xiw.xfoil_get_wake_cp(xdg, nwake)
    uedgew = xiw.xfoil_get_wake_uedge(xdg, nwake)

    xiw.xfoil_cleanup(xdg)

    # Mass defect derivative

    dmass = surface_derivative(mass, x, z)
    for i in range(81):
        dmass[i] *= -1.;    # Sign is flipped on top surface

    plot_geometry(xnew, znew, xw, zw, "Airfoil + wake")
    plot_with_wake(xnew, xw, cp, cpw, "Pressure coefficient")
    plot_with_wake(xnew, xw, uedge, uedgew, "Edge velocity")

    plot(x, dstar, "Displacement thickness")
    plot(x, mass, "Mass defect")
    plot(x, dmass, "Mass defect derivative")
