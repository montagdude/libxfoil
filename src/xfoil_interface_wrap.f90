! This file is part of libxfoil.
! 
! libxfoil is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! libxfoil is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.
! 
! Copyright (C) 2018 Daniel Prosser
! 
! See xfoil_interface.f90 for descriptions of inputs and outputs.
!
!===============================================================================
! This file contains all the signatures of subroutines and types needed to 
! interface with libxfoil in Fortran90, kind of like a header file in C. Just
! `include "xfoil_interface_wrap.f90"` in your code.

type, bind(c) :: xfoil_options_type
  real(c_double) :: ncrit
  real(c_double) :: xtript, xtripb
  logical(c_bool) :: viscous_mode
  logical(c_bool) :: silent_mode
  integer(c_int) :: maxit
  real(c_double) :: vaccel
  logical(c_bool) :: fix_unconverged
  logical(c_bool) :: reinitialize
end type xfoil_options_type

type, bind(c) :: xfoil_geom_options_type
  integer(c_int) :: npan
  real(c_double) :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
end type xfoil_geom_options_type

interface
subroutine naca_4_digit(des, npointside, xout, zout, nout)                     &
           bind(c, name="naca_4_digit")

  use iso_c_binding
  character(c_char), dimension(4), intent(in) :: des
  integer(c_int), intent(in) :: npointside
  real(c_double), dimension(2*npointside), intent(out) :: xout, zout
  integer(c_int), intent(out) :: nout 

end subroutine naca_4_digit 
end interface

interface
subroutine naca_5_digit(des, npointside, xout, zout, nout, stat)               &
           bind(c, name="naca_5_digit")

  use iso_c_binding
  character(c_char), dimension(5), intent(in) :: des
  integer(c_int), intent(in) :: npointside
  real(c_double), dimension(2*npointside), intent(out) :: xout, zout
  integer(c_int), intent(out) :: nout, stat

end subroutine naca_5_digit 
end interface

interface
subroutine smooth_paneling(xin, zin, npointin, npointout, geom_options, xout,  &
                           zout) bind(c, name="smooth_paneling")

  use iso_c_binding

  type, bind(c) :: xfoil_geom_options_type
    integer(c_int) :: npan
    real(c_double) :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  end type xfoil_geom_options_type

  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin, npointout
  type(xfoil_geom_options_type), intent(in) :: geom_options
  real(c_double), dimension(npointout), intent(out) :: xout, zout
  
end subroutine smooth_paneling
end interface

interface
subroutine xfoil_init() bind(c, name="xfoil_init")
end subroutine xfoil_init
end interface

interface
subroutine xfoil_set_airfoil(xin, zin, npointin)                               &
           bind(c, name="xfoil_set_airfoil")

  use iso_c_binding
  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin

end subroutine xfoil_set_airfoil
end interface

interface
subroutine xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)                       &
           bind(c, name="xfoil_geometry_info")

  use iso_c_binding
  real(c_double), intent(out) :: maxt, xmaxt, maxc, xmaxc

end subroutine xfoil_geometry_info
end interface

interface
subroutine xfoil_spline_coordinates(x, z, npt, s, xs, zs)                      &
           bind(c, name="xfoil_spline_coordinates")

  use iso_c_binding
  real(c_double), dimension(npt), intent(in) :: x, z
  integer(c_int), intent(in) :: npt
  real(c_double), dimension(npt), intent(out) :: s, xs, zs

end subroutine xfoil_spline_coordinates
end interface

interface
subroutine xfoil_eval_spline(x, z, s, xs, zs, npt, sc, xc, zc)                 &
           bind(c, name="xfoil_eval_spline")

  use iso_c_binding
  real(c_double), dimension(npt), intent(in) :: x, z, s, xs, zs
  integer, intent(in) :: npt
  real(c_double), intent(in) :: sc
  real(c_double), intent(out) :: xc, zc

end subroutine xfoil_eval_spline
end interface

interface
subroutine xfoil_lefind(x, z, s, xs, zs, npt, sle, xle, zle)                   &
           bind(c, name="xfoil_lefind")

  use iso_c_binding
  real(c_double), dimension(npt), intent(in) :: x, z, s, xs, zs
  integer, intent(in) :: npt
  real(c_double), intent(out) :: sle, xle, zle

end subroutine xfoil_lefind
end interface

interface
subroutine run_xfoil(npointin, xin, zin, geom_options, noppoint,               &
                     operating_points, op_modes, reynolds_numbers,             &
                     mach_numbers, use_flap, x_flap, y_flap, y_flap_spec,      &
                     flap_degrees, xfoil_options, lift, drag, moment, viscrms, &
                     alpha, xtrt, xtrb, ncrit_per_point)                       &
           bind(c, name="run_xfoil")

  use iso_c_binding

  type, bind(c) :: xfoil_options_type
    real(c_double) :: ncrit
    real(c_double) :: xtript, xtripb
    logical(c_bool) :: viscous_mode
    logical(c_bool) :: silent_mode
    integer(c_int) :: maxit
    real(c_double) :: vaccel
    logical(c_bool) :: fix_unconverged
    logical(c_bool) :: reinitialize
  end type xfoil_options_type

  type, bind(c) :: xfoil_geom_options_type
    integer(c_int) :: npan
    real(c_double) :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  end type xfoil_geom_options_type

  integer(c_int), intent(in) :: npointin, noppoint
  real(c_double), dimension(npointin), intent(in) :: xin, zin
  type(xfoil_geom_options_type), intent(in) :: geom_options
  real(c_double), dimension(noppoint), intent(in) :: operating_points,         &
                                    reynolds_numbers, mach_numbers, flap_degrees
  real(c_double), intent(in) :: x_flap, y_flap
  integer(c_int), intent(in) :: y_flap_spec
  logical(c_bool), intent(in) :: use_flap
  integer(c_int), dimension(noppoint), intent(in) :: op_modes
  type(xfoil_options_type), intent(in) :: xfoil_options
  real(c_double), dimension(noppoint), intent(out) :: lift, drag, moment,      &
                                                      viscrms
  real(c_double), dimension(noppoint), intent(out) :: alpha, xtrt, xtrb
  real(c_double), dimension(noppoint), intent(in), optional :: ncrit_per_point

end subroutine run_xfoil
end interface

interface
subroutine xfoil_cleanup() bind(c, name="xfoil_cleanup")
end subroutine xfoil_cleanup
end interface
