!  This file is part of libxfoil.

!  libxfoil is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  libxfoil is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2018 Daniel Prosser

module xfoil_interface

! Contains subroutines to use XFoil to analyze an airfoil

  use iso_c_binding
  use xfoil_data_mod

  implicit none

  type, bind(c) :: xfoil_options_type

    real(c_double) :: ncrit             !Critical ampl. ratio
    real(c_double) :: xtript, xtripb    !Trip locations
    logical(c_bool) :: viscous_mode
    logical(c_bool) :: silent_mode      !Toggle xfoil screen write
    integer(c_int) :: maxit             !Iterations for BL calcs
    real(c_double) :: vaccel            !Xfoil BL convergence accelerator
                                        !  point (recommended for optimization)

  end type xfoil_options_type

  type, bind(c) :: xfoil_geom_options_type

    integer(c_int) :: npan
    real(c_double) :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2

  end type xfoil_geom_options_type

  contains

!=============================================================================80
!
! Allocates memory for structs in xfoil_data_group using C backend
!
!=============================================================================80
subroutine xfoil_init(xdg) bind(c, name="xfoil_init")

  type(xfoil_data_group), intent(inout) :: xdg

  interface
    subroutine allocate_xdg(xdg) bind(c)
      use xfoil_data_mod, only : xfoil_data_group
      type(xfoil_data_group), intent(inout) :: xdg
    end subroutine allocate_xdg
  end interface

  call allocate_xdg(xdg)

end subroutine xfoil_init

!=============================================================================80
!
! Frees memory for structs in xfoil_data_group using C backend
!
!=============================================================================80
subroutine xfoil_cleanup(xdg) bind(c, name="xfoil_cleanup")

  type(xfoil_data_group), intent(inout) :: xdg

  interface
    subroutine free_xdg(xdg) bind(c)
      use xfoil_data_mod, only : xfoil_data_group
      type(xfoil_data_group), intent(inout) :: xdg
    end subroutine free_xdg
  end interface

  call free_xdg(xdg)

end subroutine xfoil_cleanup

!=============================================================================80
!
! Copies xfoil_data_groups using C backend. xfoil_init must have already been
! called for both.
!
!=============================================================================80
subroutine xfoil_copy(xdg_from, xdg_to) bind(c, name="xfoil_copy")

  type(xfoil_data_group), intent(in) :: xdg_from
  type(xfoil_data_group), intent(inout) :: xdg_to

  interface
    subroutine copy_xdg(xdg_from, xdg_to) bind(c)
      use xfoil_data_mod, only : xfoil_data_group
      type(xfoil_data_group), intent(in) :: xdg_from
      type(xfoil_data_group), intent(inout) :: xdg_to
    end subroutine copy_xdg
  end interface

  call copy_xdg(xdg_from, xdg_to)

end subroutine xfoil_copy

!=============================================================================80
!
! Initializes xfoil variables from settings
!
!=============================================================================80
subroutine xfoil_defaults(xdg, xfoil_options) bind(c, name="xfoil_defaults")

  use iso_c_binding
  type(xfoil_data_group), intent(inout) :: xdg
  type(xfoil_options_type), intent(in) :: xfoil_options

  real(c_double), pointer :: SIG(:)
  real(c_double), pointer :: QF0(:)
  real(c_double), pointer :: QF1(:)
  real(c_double), pointer :: QF2(:)
  real(c_double), pointer :: QF3(:)
  real(c_double), pointer :: XSTRIP(:)
  real(c_double), pointer :: GAMU(:,:)
  real(c_double), pointer :: GAM(:)
  real(c_double), pointer :: APANEL(:)

  call c_f_pointer(xdg%xfd%SIG, SIG, [IZX])
  call c_f_pointer(xdg%xfd%QF0, QF0, [IQX])
  call c_f_pointer(xdg%xfd%QF1, QF1, [IQX])
  call c_f_pointer(xdg%xfd%QF2, QF2, [IQX])
  call c_f_pointer(xdg%xfd%QF3, QF3, [IQX])
  call c_f_pointer(xdg%xfd%XSTRIP, XSTRIP, [ISX])
  call c_f_pointer(xdg%xfd%GAMU, GAMU, [IQX,2])
  call c_f_pointer(xdg%xfd%GAM, GAM, [IQX])
  call c_f_pointer(xdg%xfd%APANEL, APANEL, [IZX])

  xdg%xfd%SILENT_MODE = xfoil_options%silent_mode
  xdg%xfd%VISCOUS_MODE = xfoil_options%viscous_mode
  xdg%xfd%MAXIT = xfoil_options%maxit
  xdg%xfd%N = 0
  xdg%xfd%PI = 4.d0*atan(1.d0)
  xdg%xfd%HOPI = 0.5d0/xdg%xfd%PI
  xdg%xfd%QOPI = 0.25d0/xdg%xfd%PI
  xdg%xfd%DTOR = xdg%xfd%PI/180.d0
  xdg%xfd%QINF = 1.d0
  SIG(:) = 0.d0
  QF0(:) = 0.d0
  QF1(:) = 0.d0
  QF2(:) = 0.d0
  QF3(:) = 0.d0
  xdg%xfd%NW = 0
  xdg%xfd%RETYP = 1
  xdg%xfd%MATYP = 1
  xdg%xfd%GAMMA = 1.4d0
  xdg%xfd%GAMM1 = xdg%xfd%GAMMA - 1.d0
  xdg%xfd%XCMREF = 0.25d0
  xdg%xfd%YCMREF = 0.d0
  xdg%xfd%LVISC = xfoil_options%viscous_mode
  xdg%xfd%AWAKE = 0.d0
  xdg%xfd%AVISC = 0.d0
  xdg%xfd%ITMAX = xfoil_options%maxit
  xdg%xfd%LWDIJ = .false.
  xdg%xfd%LIPAN = .false.
  xdg%xfd%LBLINI = .false.
  xdg%xfd%ACRIT = xfoil_options%ncrit
  xdg%xfd%IDAMP = 0
  XSTRIP(1) = xfoil_options%xtript
  XSTRIP(2) = xfoil_options%xtripb
  xdg%xfd%VACCEL = xfoil_options%vaccel
  xdg%xfd%WAKLEN = 1.d0
  xdg%xfd%PSIO = 0.d0
  GAMU(:,:) = 0.d0
  GAM(:) = 0.d0
  xdg%xfd%SIGTE = 0.d0
  xdg%xfd%GAMTE = 0.d0
  xdg%xfd%SIGTE_A = 0.d0
  xdg%xfd%GAMTE_A = 0.d0
  APANEL(:) = 0.d0

! Set boundary layer calibration parameters

  call BLPINI(xdg%bld)

end subroutine xfoil_defaults

!=============================================================================80
!
! Sets xfoil paneling options
!
!=============================================================================80
subroutine xfoil_set_paneling(xdg, geom_options)                               &
           bind(c, name="xfoil_set_paneling")

  type(xfoil_data_group), intent(inout) :: xdg
  type(xfoil_geom_options_type), intent(in) :: geom_options

  xdg%xfd%NPAN = geom_options%npan
  xdg%xfd%CVPAR = geom_options%cvpar
  xdg%xfd%CTERAT = geom_options%cterat
  xdg%xfd%CTRRAT = geom_options%ctrrat
  xdg%xfd%XSREF1 = geom_options%xsref1
  xdg%xfd%XSREF2 = geom_options%xsref2
  xdg%xfd%XPREF1 = geom_options%xpref1
  xdg%xfd%XPREF2 = geom_options%xpref2

end subroutine xfoil_set_paneling

!=============================================================================80
!
! Sets buffer airfoil for xfoil
!
!=============================================================================80
subroutine xfoil_set_buffer_airfoil(xdg, xin, zin, npointin)                   &
           bind(c, name="xfoil_set_buffer_airfoil")

  use iso_c_binding
  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin

  real(c_double), pointer :: XB(:)
  real(c_double), pointer :: YB(:)

  call c_f_pointer(xdg%xfd%XB, XB, [IBX])
  call c_f_pointer(xdg%xfd%YB, YB, [IBX])

  xdg%xfd%NB = npointin
  XB(1:xdg%xfd%NB) = xin
  YB(1:xdg%xfd%NB) = zin

end subroutine xfoil_set_buffer_airfoil

!=============================================================================80
!
! Returns buffer  airfoil coordinates from Xfoil
! stat: 0 for success, 1 if buffer airfoil is not available (call
!   xfoil_set_buffer_airfoil first)
!
!=============================================================================80
subroutine xfoil_get_buffer_airfoil(xdg, xout, zout, npoint, stat)             &
           bind(c, name="xfoil_get_buffer_airfoil")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: xout, zout
  integer(c_int), intent(out) :: stat

  real(c_double), pointer :: XB(:)
  real(c_double), pointer :: YB(:)

  call c_f_pointer(xdg%xfd%XB, XB, [IBX])
  call c_f_pointer(xdg%xfd%YB, YB, [IBX])

! Check that buffer airfoil is available

  stat = 0
  if (xdg%xfd%NB == 0) then
    stat = 1
    return
  end if

  xout(1:npoint) = XB(1:npoint)
  zout(1:npoint) = YB(1:npoint)

end subroutine xfoil_get_buffer_airfoil

!=============================================================================80
!
! Returns current (not buffer) airfoil coordinates from Xfoil
! stat: 0 for success, 1 if current airfoil is not available (call
!   xfoil_smooth_paneling first)
!
!=============================================================================80
subroutine xfoil_get_current_airfoil(xdg, xout, zout, npoint, stat)            &
           bind(c, name="xfoil_get_current_airfoil")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: xout, zout
  integer(c_int), intent(out) :: stat

  real(c_double), pointer :: X(:)
  real(c_double), pointer :: Y(:)

  call c_f_pointer(xdg%xfd%X, X, [IZX])
  call c_f_pointer(xdg%xfd%Y, Y, [IZX])

! Check that airfoil is available

  stat = 0
  if (xdg%xfd%N == 0) then
    stat = 1
    return
  end if

  xout(1:npoint) = X(1:npoint)
  zout(1:npoint) = Y(1:npoint)

end subroutine xfoil_get_current_airfoil

!=============================================================================80
!
! Smooths buffer airfoil using Xfoil's PANGEN subroutine
! stat: 0 for success, 1 for failure (xfoil_set_buffer_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_smooth_paneling(xdg, stat)                                    &
           bind(c, name="xfoil_smooth_paneling")

  type(xfoil_data_group), intent(inout) :: xdg
  integer(c_int), intent(out) :: stat

! Check that buffer airfoil is set

  stat = 0
  if (xdg%xfd%NB == 0) then
    stat = 1
    return
  end if

! Smooth paneling with PANGEN

  call PANGEN(xdg%xfd, .NOT. xdg%xfd%SILENT_MODE)

end subroutine xfoil_smooth_paneling

!=============================================================================80
!
! Subroutine to apply a flap deflection to the buffer airfoil and set it as the
! current airfoil. It is recommended to call this after xfoil_smooth_paneling.
! z_flap_spec = 0: specified as y/c
!             = 1: specified as y/local thickness
! stat: 0 for success, 1 for failure (xfoil_set_buffer_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_apply_flap_deflection(xdg, xflap, zflap, z_flap_spec, degrees,&
                                       npointout, stat)                        &
           bind(c, name="xfoil_apply_flap_deflection")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: xflap, zflap, degrees
  integer(c_int), intent(in) :: z_flap_spec
  integer(c_int), intent(out) :: npointout, stat

! Check that buffer airfoil is set

  stat = 0
  if (xdg%xfd%NB == 0) then
    stat = 1
    return
  end if

! Apply flap deflection

  call FLAP(xdg%xfd, xflap, zflap, z_flap_spec, degrees)

! Get new buffer airfoil points (may have changed)

  npointout = xdg%xfd%NB

end subroutine xfoil_apply_flap_deflection

!=============================================================================80
!
! Subroutine to modify the trailing edge gap of the buffer airfoil and set it as
! the current airfoil.
! gap: the new TE gap
! blendloc: x/c location where the shape is first modified to accomodate the gap
!   0 < blendloc < 1
! stat: 0 for success, 1 for failure (xfoil_set_buffer_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_modify_tegap(xdg, gap, blendloc, npointout, stat)             &
           bind(c, name="xfoil_modify_tegap")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: gap, blendloc
  integer(c_int), intent(out) :: npointout, stat

! Check that buffer airfoil is set

  stat = 0
  if (xdg%xfd%NB == 0) then
    stat = 1
    return
  end if

! Modify trailing edge gap

  call TGAP(xdg%xfd, gap, blendloc)

! Get new buffer airfoil points (may have changed)

  npointout = xdg%xfd%NB

end subroutine xfoil_modify_tegap

!=============================================================================80
!
! Gets thickness and camber information for the current (not buffer) airfoil
! stat: 0 for success, 1 if current airfoil is not available (call
!   xfoil_smooth_paneling first)
!
!=============================================================================80
subroutine xfoil_geometry_info(xdg, maxt, xmaxt, maxc, xmaxc, stat)            &
           bind(c, name="xfoil_geometry_info")

  type(xfoil_data_group), intent(in) :: xdg
  real(c_double), intent(out) :: maxt, xmaxt, maxc, xmaxc
  integer(c_int), intent(out) :: stat

! Check that airfoil is available

  stat = 0
  if (xdg%xfd%N == 0) then
    stat = 1
    return
  end if

  maxt = xdg%xfd%THICKB
  xmaxt = xdg%xfd%XTHICKB
  maxc = xdg%xfd%CAMBR
  xmaxc = xdg%xfd%XCAMBR

end subroutine xfoil_geometry_info

!=============================================================================80
!
! Sets Reynolds number for viscous calculations
!
!=============================================================================80
subroutine xfoil_set_reynolds_number(xdg, re)                                  &
           bind(c, name="xfoil_set_reynolds_number")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: re

  xdg%xfd%REINF1 = re

end subroutine xfoil_set_reynolds_number

!=============================================================================80
!
! Sets Mach number
!
!=============================================================================80
subroutine xfoil_set_mach_number(xdg, mach)                                    &
           bind(c, name="xfoil_set_mach_number")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: mach

  call MINFSET(xdg%xfd, mach)

end subroutine xfoil_set_mach_number

!=============================================================================80
!
! Resets BL initialization flags in xfoil, so BL will be reinitialized at next
! point
!
!=============================================================================80
subroutine xfoil_reinitialize_bl(xdg) bind(c, name="xfoil_reinitialize_bl")

  type(xfoil_data_group), intent(inout) :: xdg

  xdg%xfd%LIPAN = .false.
  xdg%xfd%LBLINI = .false.

end subroutine xfoil_reinitialize_bl

!=============================================================================80
!
! Runs Xfoil at a specified angle of attack
! Assumes airfoil geometry, reynolds number, and mach number have already been
! set in Xfoil.
! stat: 0 for success, 1 if current airfoil is not available (call
!   xfoil_smooth_paneling first)
!
!=============================================================================80
subroutine xfoil_specal(xdg, alpha_spec, alpha, lift, drag, moment, converged, &
                        stat) bind(c, name="xfoil_specal")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: alpha_spec
  real(c_double), intent(out) :: alpha, lift, drag, moment
  logical(c_bool), intent(out) :: converged
  integer(c_int), intent(out) :: stat

! Check that airfoil is available

  stat = 0
  if (xdg%xfd%N == 0) then
    stat = 1
    return
  end if

! Inviscid calculations for specified angle of attack

  converged = .true.
  xdg%xfd%LALFA = .true.
  xdg%xfd%ALFA = alpha_spec*xdg%xfd%DTOR
  call SPECAL(xdg%xfd)
  if (abs(xdg%xfd%ALFA-xdg%xfd%AWAKE) .GT. 1.0D-5) xdg%xfd%LWAKE  = .false.
  if (abs(xdg%xfd%ALFA-xdg%xfd%AVISC) .GT. 1.0D-5) xdg%xfd%LVCONV = .false.
  if (abs(xdg%xfd%MINF-xdg%xfd%MVISC) .GT. 1.0D-5) xdg%xfd%LVCONV = .false.

! Viscous calculations (if requested)

  if (xdg%xfd%VISCOUS_MODE) then
    call VISCAL(xdg%xfd, xdg%bld, xdg%xbd, xdg%xfd%MAXIT)
    converged = xdg%xfd%LVCONV
  end if

! Outputs

  alpha = xdg%xfd%ALFA/xdg%xfd%DTOR
  lift = xdg%xfd%CL
  moment = xdg%xfd%CM
  if (xdg%xfd%VISCOUS_MODE) then
    drag = xdg%xfd%CD
  else
    drag = xdg%xfd%CDP
  end if

end subroutine xfoil_specal

!=============================================================================80
!
! Runs Xfoil at a specified lift coefficient
! Assumes airfoil geometry, reynolds number, and mach number have already been
! set in Xfoil.
! stat: 0 for success, 1 if current airfoil is not available (call
!   xfoil_smooth_paneling first)
!
!=============================================================================80
subroutine xfoil_speccl(xdg, cl_spec, alpha, lift, drag, moment, converged,    &
           stat) bind(c, name="xfoil_speccl")

  type(xfoil_data_group), intent(inout) :: xdg
  real(c_double), intent(in) :: cl_spec
  real(c_double), intent(out) :: alpha, lift, drag, moment
  logical(c_bool), intent(out) :: converged
  integer(c_int), intent(out) :: stat

! Check that airfoil is available

  stat = 0
  if (xdg%xfd%N == 0) then
    stat = 1
    return
  end if

! Inviscid calculations for specified lift coefficient

  converged = .true.
  xdg%xfd%LALFA = .false.
  xdg%xfd%ALFA = 0.d0
  xdg%xfd%CLSPEC = cl_spec
  call SPECCL(xdg%xfd)
  if (abs(xdg%xfd%ALFA-xdg%xfd%AWAKE) .GT. 1.0D-5) xdg%xfd%LWAKE  = .false.
  if (abs(xdg%xfd%ALFA-xdg%xfd%AVISC) .GT. 1.0D-5) xdg%xfd%LVCONV = .false.
  if (abs(xdg%xfd%MINF-xdg%xfd%MVISC) .GT. 1.0D-5) xdg%xfd%LVCONV = .false.

! Viscous calculations (if requested)

  if (xdg%xfd%VISCOUS_MODE) then
    call VISCAL(xdg%xfd, xdg%bld, xdg%xbd, xdg%xfd%MAXIT)
    converged = xdg%xfd%LVCONV
  end if

! Outputs

  alpha = xdg%xfd%ALFA/xdg%xfd%DTOR
  lift = xdg%xfd%CL
  moment = xdg%xfd%CM
  if (xdg%xfd%VISCOUS_MODE) then
    drag = xdg%xfd%CD
  else
    drag = xdg%xfd%CDP
  end if

end subroutine xfoil_speccl

!=============================================================================80
!
! Returns transition locations on top and bottom in x and z
!
!=============================================================================80
subroutine xfoil_get_transloc(xdg, xtranst, ztranst, xtransb, ztransb)         &
           bind(c, name="xfoil_get_transloc")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  real(c_double), intent(out) :: xtranst, ztranst, xtransb, ztransb

  real(c_double), pointer :: XOCTR(:)
  real(c_double), pointer :: YOCTR(:)

  call c_f_pointer(xdg%xfd%XOCTR, XOCTR, [ISX])
  call c_f_pointer(xdg%xfd%YOCTR, YOCTR, [ISX])

  xtranst = XOCTR(1)
  ztranst = YOCTR(1)
  xtransb = XOCTR(2)
  ztransb = YOCTR(2)

end subroutine xfoil_get_transloc

!=============================================================================80
!
! Returns cp on surface
!
!=============================================================================80
subroutine xfoil_get_cp(xdg, npoint, cp) bind(c, name="xfoil_get_cp")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: cp

  real(c_double), pointer :: CPV(:)
  real(c_double), pointer :: CPI(:)

  call c_f_pointer(xdg%xfd%CPV, CPV, [IZX])
  call c_f_pointer(xdg%xfd%CPI, CPI, [IZX])

  if (xdg%xfd%VISCOUS_MODE) then
    cp(1:npoint) = CPV(1:npoint)
  else
    cp(1:npoint) = CPI(1:npoint)
  end if

end subroutine xfoil_get_cp

!=============================================================================80
!
! Returns skin friction coefficient on surface
!
!=============================================================================80
subroutine xfoil_get_cf(xdg, npoint, cf) bind(c, name="xfoil_get_cf")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: cf

  integer(c_int) :: is, ibl, i
  real(c_double) :: que

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: TAU(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%TAU, TAU, [IVX,ISX])

  que = 0.5d0*xdg%xfd%QINF**2.d0

! Populate skin friction array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) cf(i) = TAU(ibl,is) / que

    end do
  end do

end subroutine xfoil_get_cf

!=============================================================================80
!
! Returns BL edge velocity on surface
!
!=============================================================================80
subroutine xfoil_get_uedge(xdg, npoint, uedge) bind(c, name="xfoil_get_uedge")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: uedge

  integer(c_int) :: is, ibl, i
  real(c_double) :: uei

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: UEDG(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%UEDG, UEDG, [IVX,ISX])

! Populate uedge array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        uei = UEDG(ibl,is)
        uedge(i) = uei * (1.d0-xdg%xfd%TKLAM) /                                &
                         (1.d0 - xdg%xfd%TKLAM*(uei/xdg%xfd%QINF)**2.d0)
      end if

    end do
  end do

end subroutine xfoil_get_uedge

!=============================================================================80
!
! Returns BL displacement thickness on surface
!
!=============================================================================80
subroutine xfoil_get_deltastar(xdg, npoint, deltastar)                         &
           bind(c, name="xfoil_get_deltastar")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: deltastar

  integer(c_int) :: is, ibl, i

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: DSTR(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%DSTR, DSTR, [IVX,ISX])

! Populate deltastar array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        deltastar(i) = DSTR(ibl,is)
      end if

    end do
  end do

end subroutine xfoil_get_deltastar

!=============================================================================80
!
! Returns BL dissipation coefficient on surface
!
!=============================================================================80
subroutine xfoil_get_diss(xdg, npoint, diss) bind(c, name="xfoil_get_diss")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: diss

  integer(c_int) :: is, ibl, i

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: DIS(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%DIS, DIS, [IVX,ISX])

! Populate diss array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        diss(i) = DIS(ibl,is) / xdg%xfd%QINF**3.d0
      end if

    end do
  end do

end subroutine xfoil_get_diss

!=============================================================================80
!
! Returns BL kinematic shape parameter on surface
!
!=============================================================================80
subroutine xfoil_get_hk(xdg, npoint, hk) bind(c, name="xfoil_get_hk")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: hk

  integer(c_int) :: is, ibl, i
  real(c_double) :: thi, dsi, uei, uc, amsq, dummy

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: THET(:,:)
  real(c_double), pointer :: DSTR(:,:)
  real(c_double), pointer :: UEDG(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%THET, THET, [IVX,ISX])
  call c_f_pointer(xdg%xfd%DSTR, DSTR, [IVX,ISX])
  call c_f_pointer(xdg%xfd%UEDG, UEDG, [IVX,ISX])

! Populate hk array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        thi = THET(ibl,is)
        dsi = DSTR(ibl,is)
        uei = UEDG(ibl,is)
        uc = uei * (1.d0-xdg%xfd%TKLAM) /                                      &
                   (1.d0 - xdg%xfd%TKLAM*(uei/xdg%xfd%QINF)**2.d0)
        amsq = uc*uc*xdg%xbd%HSTINV /                                          &
               (xdg%xfd%GAMM1*(1.d0 - 0.5d0*uc*uc*xdg%xbd%HSTINV))
        call HKIN(dsi/thi, amsq, hk(i), dummy, dummy)
      end if

    end do
  end do

end subroutine xfoil_get_hk

!=============================================================================80
!
! Returns BL momentum thickness Reynolds number on surface
!
!=============================================================================80
subroutine xfoil_get_retheta(xdg, npoint, retheta)                             &
           bind(c, name="xfoil_get_retheta")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: retheta

  integer(c_int) :: is, ibl, i
  real(c_double) :: uei, ue, herat, rhoe, amue

! Sutherland's constant/To (assumes stagnation conditions are at STP)

  real(c_double), parameter :: hvrat = 0.35d0

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: UEDG(:,:)
  real(c_double), pointer :: THET(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%UEDG, UEDG, [IVX,ISX])
  call c_f_pointer(xdg%xfd%THET, THET, [IVX,ISX])

! Populate ampl array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        uei = UEDG(ibl,is)
        ue = uei * (1.d0-xdg%xfd%TKLAM) /                                      &
                   (1.d0 - xdg%xfd%TKLAM*(uei/xdg%xfd%QINF)**2.d0)
        herat = (1.d0 - 0.5d0*xdg%xbd%HSTINV*uei**2.d0)                        &
              / (1.d0 - 0.5d0*xdg%xbd%HSTINV*xdg%xfd%QINF**2.d0)
        rhoe = herat**(1.d0/xdg%xfd%GAMM1)
        amue = sqrt(herat**3.d0) * (1.d0+hvrat)/(herat+hvrat)
        retheta(i) = xdg%xfd%REINF * rhoe*ue*THET(ibl,is)/amue
      end if

    end do
  end do

end subroutine xfoil_get_retheta

!=============================================================================80
!
! Returns amplification ratio N
!
!=============================================================================80
subroutine xfoil_get_ampl(xdg, npoint, ampl) bind(c, name="xfoil_get_ampl")

  use iso_c_binding
  type(xfoil_data_group), intent(in) :: xdg
  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: ampl

  integer(c_int) :: is, ibl, i

  integer(c_int), pointer :: NBL(:)
  integer(c_int), pointer :: IPAN(:,:)
  real(c_double), pointer :: X(:)
  real(c_double), pointer :: XOCTR(:)
  real(c_double), pointer :: CTAU(:,:)

  call c_f_pointer(xdg%xfd%NBL, NBL, [ISX])
  call c_f_pointer(xdg%xfd%IPAN, IPAN, [IVX,ISX])
  call c_f_pointer(xdg%xfd%X, X, [IZX])
  call c_f_pointer(xdg%xfd%XOCTR, XOCTR, [ISX])
  call c_f_pointer(xdg%xfd%CTAU, CTAU, [IVX,ISX])

! Populate ampl array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
!       In laminar regions, CTAU is log of amplification ratio. In turbulent
!       regions, we'll just set it to ACRIT
        if (X(i) <= XOCTR(is)) then
          ampl(i) = CTAU(ibl,is)
        else
          ampl(i) = xdg%xfd%ACRIT
        end if
      end if

    end do
  end do

end subroutine xfoil_get_ampl

!=============================================================================80
!
! Subroutine to analyze an airfoil in a loop over a number of operating points.
!
! Inputs:
!   npointin: number of points in airfoil coordinates
!   xin, zin: airfoil coordinates
!   geom_opts: Xfoil paneling settings
!   noppoints: number of operating points to analyze
!   operating_points: specified AoA or Cl for each point
!   op_modes: indicates whether each operating_point specifies an AoA or Cl.
!     0 => AoA, 1 => Cl
!   reynolds_numbers, mach_numbers: Re and Mach for each operating point
!   use_flap: T or F; whether flap deflections will be applied
!   x_flap, z_flap: x and z flap hinge coordinates
!   z_flap_spec: 0 => z_flap = z/c, 1 => z_flap = z/local_thickness
!   flap_degrees: flap deflection at each operating point (+ve down)
!   xfoil_opts: Xfoil run settings
!   reinitialize: T or F; whether to always reinitialize the BL at each point
!   fix_unconverged: T or F; if true, for any points that fail to converge,
!     will try to initialize the BL at a point closer to Cl = 0, and then
!     run again at the original operating point
!
! Outputs:
!   alpha, Cl, Cd, Cm, and x/c transition locations at each operating point
!   viscrms: rms for viscous calculations (check for convergence)
!
!=============================================================================80
subroutine run_xfoil(npointin, xin, zin, geom_opts, noppoint, operating_points,&

                     op_modes, reynolds_numbers, mach_numbers, use_flap,       &
                     x_flap, z_flap, z_flap_spec, flap_degrees, xfoil_opts,    &
                     reinitialize, fix_unconverged, lift, drag, moment,        &
                     viscrms, alpha, xtrt, xtrb, stat) bind(c, name="run_xfoil")

  integer(c_int), intent(in) :: npointin, noppoint
  real(c_double), dimension(npointin), intent(in) :: xin, zin
  type(xfoil_geom_options_type), intent(in) :: geom_opts
  real(c_double), dimension(noppoint), intent(in) :: operating_points,         &
                                    reynolds_numbers, mach_numbers, flap_degrees
  real(c_double), intent(in) :: x_flap, z_flap
  integer(c_int), intent(in) :: z_flap_spec
  type(xfoil_options_type), intent(in) :: xfoil_opts
  logical(c_bool), intent(in) :: use_flap, reinitialize, fix_unconverged
  integer(c_int), dimension(noppoint), intent(in) :: op_modes
  real(c_double), dimension(noppoint), intent(out) :: lift, drag, moment,      &
                                                      viscrms
  real(c_double), dimension(noppoint), intent(out) :: alpha, xtrt, xtrb
  integer(c_int), intent(out) :: stat

  type(xfoil_data_group) :: xdg
  integer(c_int) :: i, dummy
  logical(c_bool), dimension(noppoint) :: point_converged, point_fixed
  real(c_double) :: newpoint, ztrt, ztrb
  character(30) :: text
  character(150) :: message

  if (.not. xfoil_opts%silent_mode) then
    write(*,*)
    write(*,*) 'Analyzing aerodynamics using the XFOIL engine ...'
  end if

  point_converged(:) = .true.
  point_fixed(:) = .false.

! Set xfoil defaults and paneling settings

  call xfoil_init(xdg)
  call xfoil_defaults(xdg, xfoil_opts)
  call xfoil_set_paneling(xdg, geom_opts)

! Set airfoil and smooth paneling

  if (.not. use_flap) then
    call xfoil_set_buffer_airfoil(xdg, xin, zin, npointin)
    call xfoil_smooth_paneling(xdg, stat)
    if (stat /= 0) return
  end if

! Run xfoil for requested operating points

  lift(:) = 0.d0
  drag(:) = 0.d0
  moment(:) = 0.d0
  viscrms(:) = 0.d0

! Run xfoil for requested operating points

  run_oppoints: do i = 1, noppoint

!   Reset airfoil, smooth paneling, and apply flap deflection

    if (use_flap) then
      call xfoil_set_buffer_airfoil(xdg, xin, zin, npointin)
      call xfoil_smooth_paneling(xdg, stat)
      call xfoil_apply_flap_deflection(xdg, x_flap, z_flap, z_flap_spec,       &
                                       flap_degrees(i), dummy, stat)
    end if

    call xfoil_set_reynolds_number(xdg, reynolds_numbers(i))
    call xfoil_set_mach_number(xdg, mach_numbers(i))

    if (reinitialize) call xfoil_reinitialize_bl(xdg)

    if (op_modes(i) == 0) then

      call xfoil_specal(xdg, operating_points(i), alpha(i), lift(i), drag(i),  &
                        moment(i), point_converged(i), stat)

    elseif (op_modes(i) == 1) then

      call xfoil_speccl(xdg, operating_points(i), alpha(i), lift(i), drag(i),  &
                        moment(i), point_converged(i), stat)

    else

      write(*,*)
      write(*,*) "Error in xfoil_interface: op_mode must be 0 or 1."
      write(*,*)
      stop

    end if

!   Additional outputs

    call xfoil_get_transloc(xdg, xtrt(i), ztrt, xtrb(i), ztrb)

!   Handling of unconverged points

    if (xfoil_opts%viscous_mode .and. .not. point_converged(i)) then

      if (fix_unconverged) then

!       Try to initialize BL at new point (in the direction away from stall)

        newpoint = operating_points(i) - 0.5d0*abs(operating_points(i))*sign(  &
                                                   1.d0, operating_points(i))
        if (newpoint == 0.d0) newpoint = 0.1d0

        call xfoil_reinitialize_bl(xdg)
        if (op_modes(i) == 0) then
          call xfoil_specal(xdg, newpoint, alpha(i), lift(i), drag(i),         &
                            moment(i), point_converged(i), stat)
        else
          call xfoil_speccl(xdg, newpoint, alpha(i), lift(i), drag(i),         &
                            moment(i), point_converged(i), stat)
        end if

!       Now try to run again at the old operating point

        if (op_modes(i) == 0) then
          call xfoil_specal(xdg, operating_points(i), alpha(i), lift(i),       &
                            drag(i), moment(i), point_converged(i), stat)
        else
          call xfoil_speccl(xdg, operating_points(i), alpha(i), lift(i),       &
                            drag(i), moment(i), point_converged(i), stat)
        end if

        if (point_converged(i)) point_fixed(i) = .true.

        call xfoil_get_transloc(xdg, xtrt(i), ztrt, xtrb(i), ztrb)

      end if
  end if

!   Convergence check

    viscrms(i) = xdg%xfd%RMSBL

  end do run_oppoints

! Final check for NaNs

  do i = 1, noppoint
    if (isnan(lift(i))) then
      lift(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(drag(i))) then
      drag(i) = 1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(moment(i))) then
      moment(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(viscrms(i))) then
      viscrms(i) = 1.D+08
    end if
  end do

! Print warnings about unconverged points

  if (.not. xfoil_opts%silent_mode) then

    write(*,*)

    do i = 1, noppoint

      write(text,*) i
      text = adjustl(text)

      if (point_converged(i)) then
        message = 'Operating point '//trim(text)//' converged.'
      elseif (.not. point_converged(i) .and. point_fixed(i)) then
        message = 'Operating point '//trim(text)//' initially did not '//      &
                  'converge but was fixed.'
      elseif (.not. point_converged(i) .and. .not. point_fixed(i)) then
        message = 'Operating point '//trim(text)//' initially did not '//      &
                  'converge and was not fixed.'
      end if

      write(*,*) trim(message)

    end do
  end if

  call xfoil_cleanup(xdg)

end subroutine run_xfoil

!=============================================================================80
!
! Subroutine to generate a 4-digit NACA airfoil
! Inputs:
! 	des: 4-digit designation
!	npointside: number of points per side
! Outputs:
!	xout: x coordinates
!	zout: z coordinates
!	nout: total number of points (2*npointside)
!
!=============================================================================80
subroutine naca_4_digit(des, npointside, xout, zout, nout)                     &

           bind(c, name="naca_4_digit")

  character(c_char), dimension(4), intent(in) :: des
  integer(c_int), intent(in) :: npointside
  real(c_double), dimension(2*npointside), intent(out) :: xout, zout
  integer(c_int), intent(out) :: nout

  integer(c_int) :: ides, i
  real(c_double), dimension(npointside) :: xx, yt, yc
  character(4, kind=c_char) :: desfor
  character(9) :: foilname

  do i = 1, 4
    desfor(i:i) = des(i)
  end do
  read(desfor,'(I4)') ides

  call NACA4(ides, xx, yt, yc, npointside, xout, zout, nout, foilname)

end subroutine naca_4_digit

!=============================================================================80
!
! Subroutine to generate a 5-digit NACA airfoil
! Inputs:
! 	des: 5-digit designation
!	npointside: number of points per side
! Outputs:
!	xout: x coordinates
!	zout: z coordinates
!	nout: total number of points (2*npointside)
!       stat: 0 for success, 1 for failure
!
!=============================================================================80
subroutine naca_5_digit(des, npointside, xout, zout, nout, stat)               &

           bind(c, name="naca_5_digit")

  character(c_char), dimension(5), intent(in) :: des
  integer(c_int), intent(in) :: npointside
  real(c_double), dimension(2*npointside), intent(out) :: xout, zout
  integer(c_int), intent(out) :: nout, stat

  integer(c_int) :: ides, i
  real(c_double), dimension(npointside) :: xx, yt, yc
  character(5, kind=c_char) :: desfor
  character(10) :: foilname

  do i = 1, 5
    desfor(i:i) = des(i)
  end do
  read(desfor,'(I5)') ides

  call NACA5(ides, xx, yt, yc, npointside, xout, zout, nout, foilname)
  if (ides == 0) then
    stat = 1
  else
    stat = 0
  end if

end subroutine naca_5_digit

!=============================================================================80
!
! Fits a spline to airfoil coordinates
!
!=============================================================================80
subroutine xfoil_spline_coordinates(x, z, npt, s, xs, zs)                      &
           bind(c, name="xfoil_spline_coordinates")

  real(c_double), dimension(npt), intent(in) :: x, z
  integer(c_int), intent(in) :: npt
  real(c_double), dimension(npt), intent(out) :: s, xs, zs

  call SCALC(x, z, s, npt)
  call SEGSPL(x, xs, s, npt)
  call SEGSPL(z, zs, s, npt)

end subroutine xfoil_spline_coordinates

!=============================================================================80
!
! Computes x and z at a given spline coordinate
!
! Inputs
! x, z: buffer airfoil coordinates
! s: arc length array (from xfoil_spline_coordinates)
! xs, zs: splined airfoil coordinates (from xfoil_spline_coordinates)
! npt: number of points in buffer coordinates
! sc: arc length value to calculate xc and zc
!
! Outputs
! xc, zc: coordinates at sc
!
!=============================================================================80
subroutine xfoil_eval_spline(x, z, s, xs, zs, npt, sc, xc, zc)                 &
           bind(c, name="xfoil_eval_spline")

  real(c_double), dimension(npt), intent(in) :: x, z, s, xs, zs
  integer(c_int), intent(in) :: npt
  real(c_double), intent(in) :: sc
  real(c_double), intent(out) :: xc, zc

  interface
    double precision function SEVAL(SS, X, XS, S, N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface

  xc = SEVAL(sc, x, xs, s, npt)
  zc = SEVAL(sc, z, zs, s, npt)

end subroutine xfoil_eval_spline

!=============================================================================80
!
! Computes leading edge arc length, x, and z
!
! Inputs
! x, z: buffer airfoil coordinates
! s: arc length array (from xfoil_spline_coordinates)
! xs, zs: splined airfoil coordinates (from xfoil_spline_coordinates)
! npt: number of points in buffer coordinates
!
! Outputs
! sle: leading edge arc length
! xle, zle: leading edge coordinates
!
!=============================================================================80
subroutine xfoil_lefind(x, z, s, xs, zs, npt, sle, xle, zle)                   &
           bind(c, name="xfoil_lefind")

  real(c_double), dimension(npt), intent(in) :: x, z, s, xs, zs
  integer(c_int), intent(in) :: npt
  real(c_double), intent(out) :: sle, xle, zle

  call LEFIND(sle, x, xs, z, zs, s, npt, .true.)
  call xfoil_eval_spline(x, z, s, xs, zs, npt, sle, xle, zle)

end subroutine xfoil_lefind

end module xfoil_interface
