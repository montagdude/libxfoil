#!/usr/bin/env python
#
# This file is part of libxfoil.
# 
# libxfoil is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# libxfoil is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.
# 
# Copyright (C) 2017 Daniel Prosser
# 
# See xfoil_interface.f90 for descriptions of inputs and outputs.
#
#===============================================================================
#
# Abstracts swig datatypes so the functions can be called with native types.
#

import xfoil_interface as xi

def naca_4_digit(des, npointside):

  npointside_p = xi.copy_intp(npointside)
  xout_a = xi.new_doublea(2*npointside)
  zout_a = xi.new_doublea(2*npointside)
  nout_p = xi.new_intp()

  xi.naca_4_digit(des, npointside_p, xout_a, zout_a, nout_p)

  nout = xi.intp_value(nout_p)
  xout = nout*[0]
  zout = nout*[0]
  for i in range(nout):
    xout[i] = xi.doublea_getitem(xout_a, i)
    zout[i] = xi.doublea_getitem(zout_a, i)

  xi.delete_intp(npointside_p)
  xi.delete_doublea(xout_a)
  xi.delete_doublea(zout_a)
  xi.delete_intp(nout_p)

  return xout, zout, nout

def naca_5_digit(des, npointside):

  npointside_p = xi.copy_intp(npointside)
  xout_a = xi.new_doublea(2*npointside)
  zout_a = xi.new_doublea(2*npointside)
  nout_p = xi.new_intp()
  stat_p = xi.new_intp()

  xi.naca_5_digit(des, npointside_p, xout_a, zout_a, nout_p, stat_p)

  nout = xi.intp_value(nout_p)
  stat = xi.intp_value(stat_p)
  xout = nout*[0]
  zout = nout*[0]
  for i in range(nout):
    xout[i] = xi.doublea_getitem(xout_a, i)
    zout[i] = xi.doublea_getitem(zout_a, i)

  xi.delete_intp(npointside_p)
  xi.delete_doublea(xout_a)
  xi.delete_doublea(zout_a)
  xi.delete_intp(nout_p)
  xi.delete_intp(stat_p)

  return xout, zout, nout, stat

def smooth_paneling(xin, zin, npointin, npointout):

  npointin_p = xi.copy_intp(npointin)
  npointout_p = xi.copy_intp(npointout)
  xin_a = xi.new_doublea(npointin)
  zin_a = xi.new_doublea(npointin)
  xout_a = xi.new_doublea(npointout)
  zout_a = xi.new_doublea(npointout)
  for i in range(npointin):
    xi.doublea_setitem(xin_a, i, xin[i])
    xi.doublea_setitem(zin_a, i, zin[i])

  xi.smooth_paneling(xin_a, zin_a, npointin_p, npointout_p, xout_a, zout_a)

  xout = npointout*[0]
  zout = npointout*[0]
  for i in range(npointout):
    xout[i] = xi.doublea_getitem(xout_a, i)
    zout[i] = xi.doublea_getitem(zout_a, i)

  xi.delete_intp(npointin_p)
  xi.delete_intp(npointout_p)
  xi.delete_doublea(xin_a)
  xi.delete_doublea(zin_a)
  xi.delete_doublea(xout_a)
  xi.delete_doublea(zout_a)

  return xout, zout

def xfoil_init():

  xi.xfoil_init()

def xfoil_set_airfoil(xin, zin, npointin):

  xin_a = xi.new_doublea(npointin)
  zin_a = xi.new_doublea(npointin)
  npointin_p = xi.copy_intp(npointin)
  for i in range(npointin):
    xi.doublea_setitem(xin_a, i, xin[i])
    xi.doublea_setitem(zin_a, i, zin[i])

  xi.xfoil_set_airfoil(xin_a, zin_a, npointin_p)

  xi.delete_doublea(xin_a)
  xi.delete_doublea(zin_a)
  xi.delete_intp(npointin_p)

def xfoil_geometry_info():

  maxt_p = xi.new_doublep()
  xmaxt_p = xi.new_doublep()
  maxc_p = xi.new_doublep()
  xmaxc_p = xi.new_doublep()

  xi.xfoil_geometry_info(maxt_p, xmaxt_p, maxc_p, xmaxc_p)

  maxt = xi.doublep_value(maxt_p)
  xmaxt = xi.doublep_value(xmaxt_p)
  maxc = xi.doublep_value(maxc_p)
  xmaxc = xi.doublep_value(xmaxc_p)

  xi.delete_doublep(maxt_p)
  xi.delete_doublep(xmaxt_p)
  xi.delete_doublep(maxc_p)
  xi.delete_doublep(xmaxc_p)

  return maxt, xmaxt, maxc, xmaxc

def run_xfoil(npointin, xin, zin, geom_options, noppoint, operating_points,
              op_modes, reynolds_numbers, mach_numbers, use_flap, x_flap,
              y_flap, y_flap_spec, flap_degrees, xfoil_options,
              ncrit_per_point=None):

  npointin_p = xi.copy_intp(npointin)
  xin_a = xi.new_doublea(npointin)
  zin_a = xi.new_doublea(npointin)
  noppoint_p = xi.copy_intp(noppoint)
  operating_points_a = xi.new_doublea(noppoint)
  op_modes_a = xi.new_inta(noppoint)
  reynolds_numbers_a = xi.new_doublea(noppoint)
  mach_numbers_a = xi.new_doublea(noppoint)
  use_flap_p = xi.copy_boolp(use_flap)
  if use_flap:
    x_flap_p = xi.copy_doublep(x_flap)
    y_flap_p = xi.copy_doublep(y_flap)
    y_flap_spec_p = xi.copy_intp(y_flap_spec)
  else:
    x_flap_p = xi.new_doublep()
    y_flap_p = xi.new_doublep()
    y_flap_spec_p = xi.new_intp()
  flap_degrees_a = xi.new_doublea(noppoint)
  lift_a = xi.new_doublea(noppoint)
  drag_a = xi.new_doublea(noppoint)
  moment_a = xi.new_doublea(noppoint) 
  viscrms_a = xi.new_doublea(noppoint)
  alpha_a = xi.new_doublea(noppoint)
  xtrt_a = xi.new_doublea(noppoint)
  xtrb_a = xi.new_doublea(noppoint)
  ncrit_per_point_a = xi.new_doublea(noppoint) 

  for i in range(npointin):
    xi.doublea_setitem(xin_a, i, xin[i])
    xi.doublea_setitem(zin_a, i, zin[i])
  for i in range(noppoint):
    xi.doublea_setitem(operating_points_a, i, operating_points[i])
    xi.inta_setitem(op_modes_a, i, op_modes[i])
    xi.doublea_setitem(reynolds_numbers_a, i, reynolds_numbers[i])
    xi.doublea_setitem(mach_numbers_a, i, mach_numbers[i])
    if use_flap:
      xi.doublea_setitem(flap_degrees_a, i, flap_degrees[i])
    if ncrit_per_point is not None:
      xi.double_setitem(ncrit_per_point_a, i, ncrit_per_point[i]) 

  if ncrit_per_point is not None:
    xi.run_xfoil(npointin_p, xin_a, zin_a, geom_options, noppoint_p,
                 operating_points_a, op_modes_a, reynolds_numbers_a,
                 mach_numbers_a, use_flap_p, x_flap_p, y_flap_p, y_flap_spec_p,
                 flap_degrees_a, xfoil_options, lift_a, drag_a, moment_a,
                 viscrms_a, alpha_a, xtrt_a, xtrb_a, ncrit_per_point_a)
  else:
    xi.run_xfoil(npointin_p, xin_a, zin_a, geom_options, noppoint_p,
                 operating_points_a, op_modes_a, reynolds_numbers_a,
                 mach_numbers_a, use_flap_p, x_flap_p, y_flap_p, y_flap_spec_p,
                 flap_degrees_a, xfoil_options, lift_a, drag_a, moment_a,
                 viscrms_a, alpha_a, xtrt_a, xtrb_a)

  lift = noppoint*[0]
  drag = noppoint*[0]
  moment = noppoint*[0]
  viscrms = noppoint*[0]
  alpha = noppoint*[0]
  xtrt = noppoint*[0]
  xtrb = noppoint*[0]
  for i in range(noppoint):
    lift[i] = xi.doublea_getitem(lift_a, i)
    drag[i] = xi.doublea_getitem(drag_a, i)
    moment[i] = xi.doublea_getitem(moment_a, i)
    viscrms[i] = xi.doublea_getitem(viscrms_a, i)
    alpha[i] = xi.doublea_getitem(alpha_a, i)
    xtrt[i] = xi.doublea_getitem(xtrt_a, i)
    xtrb[i] = xi.doublea_getitem(xtrb_a, i)

  xi.delete_intp(npointin_p)
  xi.delete_doublea(xin_a)
  xi.delete_doublea(zin_a)
  xi.delete_intp(noppoint_p)
  xi.delete_doublea(operating_points_a)
  xi.delete_inta(op_modes_a)
  xi.delete_doublea(reynolds_numbers_a)
  xi.delete_doublea(mach_numbers_a)
  xi.delete_boolp(use_flap_p)
  xi.delete_doublep(x_flap_p)
  xi.delete_doublep(y_flap_p)
  xi.delete_intp(y_flap_spec_p)
  xi.delete_doublea(flap_degrees_a)
  xi.delete_doublea(lift_a)
  xi.delete_doublea(drag_a)
  xi.delete_doublea(moment_a)
  xi.delete_doublea(viscrms_a)
  xi.delete_doublea(alpha_a)
  xi.delete_doublea(xtrt_a)
  xi.delete_doublea(xtrb_a)
  xi.delete_doublea(ncrit_per_point_a)

  return lift, drag, moment, viscrms, alpha, xtrt, xtrb

def xfoil_cleanup():

  xi.xfoil_cleanup()
