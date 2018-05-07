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
# Copyright (C) 2018 Daniel Prosser
# 
# See xfoil_interface.f90 for descriptions of inputs and outputs.
#
#===============================================================================
#
# Abstracts swig datatypes so the functions can be called with native types.
#

import xfoil_interface as xi

def xfoil_init():

  xi.xfoil_init()

def xfoil_defaults(xfoil_options):

  xi.xfoil_defaults(xfoil_options)

def xfoil_set_paneling(geom_opts):

  xi.xfoil_set_paneling(geom_opts)

def xfoil_set_airfoil(xin, zin, npointin):

  xin_a = xi.new_doublea(npointin)
  zin_a = xi.new_doublea(npointin)
  npointin_p = xi.copy_intp(npointin)
  stat_p = xi.new_intp()
  for i in range(npointin):
    xi.doublea_setitem(xin_a, i, xin[i])
    xi.doublea_setitem(zin_a, i, zin[i])

  xi.xfoil_set_airfoil(xin_a, zin_a, npointin_p, stat_p)
  stat = xi.intp_value(stat_p)

  xi.delete_doublea(xin_a)
  xi.delete_doublea(zin_a)
  xi.delete_intp(npointin_p)
  xi.delete_intp(stat_p)

  return stat

def xfoil_smooth_paneling():

  stat_p = xi.new_intp()

  xi.xfoil_smooth_paneling(stat_p)
  stat = xi.intp_value(stat_p)

  xi.delete_intp(stat_p)

  return stat

def xfoil_apply_flap_deflection(xflap, yflap, y_flap_spec, degrees):

  xflap_p = xi.copy_doublep(xflap)
  yflap_p = xi.copy_doublep(yflap)
  y_flap_spec_p = xi.copy_intp(y_flap_spec)
  degrees_p = xi.copy_doublep(degrees)
  npointout_p = xi.new_intp()
  stat_p = xi.new_intp()

  xi.xfoil_apply_flap_deflection(xflap_p, yflap_p, y_flap_spec_p, degrees_p,
                                 npointout_p, stat_p)

  npointout = xi.intp_value(npointout_p)
  stat = xi.intp_value(stat_p)

  xi.delete_doublep(xflap_p)
  xi.delete_doublep(yflap_p)
  xi.delete_intp(y_flap_spec_p)
  xi.delete_intp(npointout_p)
  xi.delete_intp(stat_p)

  return npointout, stat

def xfoil_modify_tegap(gap, blendloc):

  stat_p = xi.new_intp()
  gap_p = xi.copy_doublep(gap)
  blendloc_p = xi.copy_doublep(blendloc)
  npointout_p = xi.new_intp()

  xi.xfoil_modify_tegap(gap_p, blendloc_p, npointout_p, stat_p)

  npointout = xi.intp_value(npointout_p)
  stat = xi.intp_value(stat_p)

  xi.delete_doublep(gap_p)
  xi.delete_doublep(blendloc_p)
  xi.delete_intp(npointout_p)
  xi.delete_intp(stat_p)

  return npointout, stat

def xfoil_get_airfoil(npoint):

  xout_a = xi.new_doublea(npoint)
  zout_a = xi.new_doublea(npoint)
  npoint_p = xi.copy_intp(npoint)

  xi.xfoil_get_airfoil(xout_a, zout_a, npoint_p)

  xout = npoint*[0]
  zout = npoint*[0]
  for i in range(npoint):
    xout[i] = xi.doublea_getitem(xout_a, i)
    zout[i] = xi.doublea_getitem(zout_a, i)

  xi.delete_doublea(xout_a)
  xi.delete_doublea(zout_a)
  xi.delete_intp(npoint_p)

  return xout, zout

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

def xfoil_set_reynolds_number(re):

  re_p = xi.copy_doublep(re)

  xi.xfoil_set_reynolds_number(re_p)

  xi.delete_doublep(re_p)

def xfoil_set_mach_number(mach):

  mach_p = xi.copy_doublep(mach)

  xi.xfoil_set_mach_number(mach_p)

  xi.delete_doublep(mach_p)

def xfoil_reinitialize_bl():

  xi.xfoil_reinitialize_bl()

def xfoil_specal(alpha_spec):

  alpha_spec_p = xi.copy_doublep(alpha_spec)
  alpha_p = xi.new_doublep()
  lift_p = xi.new_doublep()
  drag_p = xi.new_doublep() 
  moment_p = xi.new_doublep()
  converged_p = xi.new_boolp()
  stat_p = xi.new_intp()

  xi.xfoil_specal(alpha_spec_p, alpha_p, lift_p, drag_p, moment_p, converged_p,
                  stat_p)

  alpha = xi.doublep_value(alpha_p)
  lift = xi.doublep_value(lift_p)
  drag = xi.doublep_value(drag_p)
  moment = xi.doublep_value(moment_p)
  converged = xi.boolp_value(converged_p)
  stat = xi.intp_value(stat_p)

  xi.delete_doublep(alpha_spec_p)
  xi.delete_doublep(alpha_p)
  xi.delete_doublep(lift_p)
  xi.delete_doublep(drag_p)
  xi.delete_doublep(moment_p)
  xi.delete_boolp(converged_p)
  xi.delete_intp(stat_p)

  return alpha, lift, drag, moment, converged, stat

def xfoil_speccl(cl_spec):

  cl_spec_p = xi.copy_doublep(cl_spec)
  alpha_p = xi.new_doublep()
  lift_p = xi.new_doublep()
  drag_p = xi.new_doublep() 
  moment_p = xi.new_doublep()
  converged_p = xi.new_boolp()
  stat_p = xi.new_intp()

  xi.xfoil_speccl(cl_spec_p, alpha_p, lift_p, drag_p, moment_p, converged_p,
                  stat_p)

  alpha = xi.doublep_value(alpha_p)
  lift = xi.doublep_value(lift_p)
  drag = xi.doublep_value(drag_p)
  moment = xi.doublep_value(moment_p)
  converged = xi.boolp_value(converged_p)
  stat = xi.intp_value(stat_p)

  xi.delete_doublep(cl_spec_p)
  xi.delete_doublep(alpha_p)
  xi.delete_doublep(lift_p)
  xi.delete_doublep(drag_p)
  xi.delete_doublep(moment_p)
  xi.delete_boolp(converged_p)
  xi.delete_intp(stat_p)

  return alpha, lift, drag, moment, converged, stat

def xfoil_cleanup():

  xi.xfoil_cleanup()

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

def xfoil_spline_coordinates(x, z, npt):

  x_a = xi.new_doublea(npt)
  z_a = xi.new_doublea(npt)
  npt_p = xi.copy_intp(npt)
  s_a = xi.new_doublea(npt) 
  xs_a = xi.new_doublea(npt) 
  zs_a = xi.new_doublea(npt) 
  for i in range(npt):
    xi.doublea_setitem(x_a, i, x[i])
    xi.doublea_setitem(z_a, i, z[i])

  xi.xfoil_spline_coordinates(x_a, z_a, npt_p, s_a, xs_a, zs_a)

  s = npt*[0]
  xs = npt*[0]
  zs = npt*[0]
  for i in range(npt):
    s[i] = xi.doublea_getitem(s_a, i)
    xs[i] = xi.doublea_getitem(xs_a, i)
    zs[i] = xi.doublea_getitem(zs_a, i)

  xi.delete_doublea(x_a)
  xi.delete_doublea(z_a)
  xi.delete_intp(npt_p)
  xi.delete_doublea(s_a)
  xi.delete_doublea(xs_a)
  xi.delete_doublea(zs_a)

  return s, xs, zs

def xfoil_eval_spline(x, z, s, xs, zs, npt, sc):

  x_a = xi.new_doublea(npt)
  z_a = xi.new_doublea(npt)
  s_a = xi.new_doublea(npt) 
  xs_a = xi.new_doublea(npt) 
  zs_a = xi.new_doublea(npt) 
  npt_p = xi.copy_intp(npt)
  sc_p = xi.copy_doublep(sc)
  xc_p = xi.new_doublep()
  zc_p = xi.new_doublep()
  for i in range(npt):
    xi.doublea_setitem(x_a, i, x[i])
    xi.doublea_setitem(z_a, i, z[i])
    xi.doublea_setitem(s_a, i, s[i])
    xi.doublea_setitem(xs_a, i, xs[i])
    xi.doublea_setitem(zs_a, i, zs[i])

  xi.xfoil_eval_spline(x_a, z_a, s_a, xs_a, zs_a, npt_p, sc_p, xc_p, zc_p)

  xc = xi.doublep_value(xc_p) 
  zc = xi.doublep_value(zc_p) 

  xi.delete_doublea(x_a)
  xi.delete_doublea(z_a)
  xi.delete_doublea(s_a)
  xi.delete_doublea(xs_a)
  xi.delete_doublea(zs_a)
  xi.delete_doublep(sc_p)
  xi.delete_doublep(xc_p)
  xi.delete_doublep(zc_p)

  return xc, zc

def xfoil_lefind(x, z, s, xs, zs, npt):

  x_a = xi.new_doublea(npt)
  z_a = xi.new_doublea(npt)
  s_a = xi.new_doublea(npt) 
  xs_a = xi.new_doublea(npt) 
  zs_a = xi.new_doublea(npt) 
  npt_p = xi.copy_intp(npt)
  sle_p = xi.new_doublep()
  xle_p = xi.new_doublep()
  zle_p = xi.new_doublep()
  for i in range(npt):
    xi.doublea_setitem(x_a, i, x[i])
    xi.doublea_setitem(z_a, i, z[i])
    xi.doublea_setitem(s_a, i, s[i])
    xi.doublea_setitem(xs_a, i, xs[i])
    xi.doublea_setitem(zs_a, i, zs[i])

  xi.xfoil_lefind(x_a, z_a, s_a, xs_a, zs_a, npt_p, sle_p, xle_p, zle_p)

  sle = xi.doublep_value(sle_p) 
  xle = xi.doublep_value(xle_p) 
  zle = xi.doublep_value(zle_p) 

  xi.delete_doublea(x_a)
  xi.delete_doublea(z_a)
  xi.delete_doublea(s_a)
  xi.delete_doublea(xs_a)
  xi.delete_doublea(zs_a)
  xi.delete_doublep(sle_p)
  xi.delete_doublep(xle_p)
  xi.delete_doublep(zle_p)

  return sle, xle, zle

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
