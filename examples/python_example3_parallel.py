#!/usr/bin/env python

import sys
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import time
import multiprocessing as mp

import libxfoil_wrap as xiw
from libxfoil import xfoil_options_type, xfoil_geom_options_type, \
                     xfoil_data_group

# Global Xfoil settings
opts = xfoil_options_type()
opts.ncrit = 7.5
opts.xtript = 1.
opts.xtripb = 1.
opts.viscous_mode = True
opts.silent_mode = True
opts.maxit = 300
opts.vaccel = 0.01

geom_opts = xfoil_geom_options_type()
geom_opts.npan = 200
geom_opts.cvpar = 1.
geom_opts.cterat = 0.15
geom_opts.ctrrat = 0.2
geom_opts.xsref1 = 1.
geom_opts.xsref2 = 1.
geom_opts.xpref1 = 1.
geom_opts.xpref2 = 1.

def analyze_all_re(digits, re, mach):

    xdg = xfoil_data_group()

    # Set up airfoil
    npointside = 100
    x, z, npoint = xiw.naca_4_digit(digits, npointside)
    xiw.xfoil_init(xdg)
    xiw.xfoil_defaults(xdg, opts)
    xiw.xfoil_set_paneling(xdg, geom_opts)
    xiw.xfoil_set_buffer_airfoil(xdg, x, z, npoint)
    if xiw.xfoil_smooth_paneling(xdg) != 0:
        print("Error smoothing paneling: xfoil_set_buffer_airfoil must be " +
              "called first.")
        sys.exit(1)

    xiw.xfoil_set_mach_number(xdg, mach)
    noper = 61
    alpha = np.linspace(-10., 10., noper)
    results = []
    for reynolds_number in re:
        xiw.xfoil_set_reynolds_number(xdg, reynolds_number)
        xiw.xfoil_reinitialize_bl(xdg)
        res = {'re': reynolds_number, 'mach': mach,
               'alpha': np.linspace(-10., 10., noper), 'cl': np.zeros((noper)),
               'cd': np.zeros((noper)), 'cm': np.zeros((noper))}
        for i in range(noper):
            alpha_out, res['cl'][i], res['cd'][i], res['cm'][i], converged, \
                stat = xiw.xfoil_specal(xdg, res['alpha'][i])
            if stat == 0:
                if not converged:
                    print("Warning: convergence failed for alpha {:.1f}." \
                          .format(res['alpha'][i]))
            else:
                print("Error: xfoil_smooth_paneling must be called first.")
                sys.exit(1)
        results.append(res)

    return results

def analyze_one_re(digits, re, mach):

    xdg = xfoil_data_group()

    # Set up airfoil
    npointside = 100
    x, z, npoint = xiw.naca_4_digit(digits, npointside)
    xiw.xfoil_init(xdg)
    xiw.xfoil_defaults(xdg, opts)
    xiw.xfoil_set_paneling(xdg, geom_opts)
    xiw.xfoil_set_buffer_airfoil(xdg, x, z, npoint)
    if xiw.xfoil_smooth_paneling(xdg) != 0:
        print("Error smoothing paneling: xfoil_set_buffer_airfoil must be " +
              "called first.")
        sys.exit(1)

    # Set operating point
    xiw.xfoil_set_reynolds_number(xdg, re)
    xiw.xfoil_set_mach_number(xdg, mach)
    noper = 61
    alpha = np.linspace(-10., 10., noper)
    res = {'re': re, 'mach': mach, 'alpha': np.linspace(-10., 10., noper),
           'cl': np.zeros((noper)), 'cd': np.zeros((noper)),
           'cm': np.zeros((noper))}
    for i in range(noper):
        alpha_out, res['cl'][i], res['cd'][i], res['cm'][i], converged, stat = \
        xiw.xfoil_specal(xdg, res['alpha'][i])
        if stat == 0:
            if not converged:
                print("Warning: convergence failed for alpha {:.1f}." \
                      .format(res['alpha'][i]))
        else:
            print("Error: xfoil_smooth_paneling must be called first.")
            sys.exit(1)

    return res

def plot(results):
    fig, axarr = plt.subplots(2,2)

    axarr[0,0].set_xlabel("Angle of attack")
    axarr[0,0].set_ylabel("Lift coefficient")
    lines = []
    labels = []
    for result in results:
        line, = axarr[0,0].plot(result['alpha'], result['cl'])
        lines.append(line)
        labels.append("Re = {:.1e}".format(result['re']))

    axarr[0,1].set_xlabel("Drag coefficient")
    axarr[0,1].set_ylabel("Lift coefficient")
    for result in results:
        axarr[0,1].plot(result['cd'], result['cl'])

    axarr[1,0].set_xlabel("Angle of attack")
    axarr[1,0].set_ylabel("Pitching moment coefficient")
    for result in results:
        axarr[1,0].plot(result['alpha'], result['cm'])

    # Use the last frame for a legend
    axarr[1,1].legend(handles=lines, labels=labels, loc='center', ncol=3)

    plt.show()

# Analyze an airfoil over a range of Reynolds numbers and angles of attack
if __name__ == "__main__":

    nre = 33
    re = np.linspace(2e+5, 6e+5, nre)
    mach = 0.05
    digits = '0010'

    # Execute in serial
    start = time.time()
    results = analyze_all_re(digits, re, mach)
    end = time.time()
    print("Serial execution time: {:.4f} sec".format(end-start))

    # Execute in parallel
    # Note: currently, xdg cannot be passed as an argument to a worker process,
    #   which means each process must set up a new xdg when run in parallel.
    #   This incurs some overhead compared to the serial case, but testing
    #   indicates this overhead is not very significant, and the parallel
    #   version still performs considerably better than the serial version for
    #   this test (greater than 2x speed-up on a quad-core Intel i5 laptop).
    pool = mp.Pool(mp.cpu_count())
    resultslist = []
    start = time.time()
    for i in range(nre):
        resultslist.append(pool.apply_async(analyze_one_re,
                                            args=(digits, re[i], mach)))
    results = []
    for result in resultslist:
        results.append(result.get())
    end = time.time()
    print("Parallel execution time: {:.4f} sec".format(end-start))

    plot(results)
