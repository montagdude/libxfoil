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
! Allocates xfoil variables that may be too big for the stack in OpenMP
!
!=============================================================================80
subroutine xfoil_init() bind(c, name="xfoil_init")

  use xfoil_inc

! Allocate variables that may be too big for the stack in OpenMP

  if (allocated(AIJ)) return

  allocate(AIJ(IQX,IQX))
  allocate(BIJ(IQX,IZX))
  allocate(DIJ(IZX,IZX))
  allocate(CIJ(IWX,IQX))
  allocate(IPAN(IVX,ISX))
  allocate(ISYS(IVX,ISX))
  allocate(W1(6*IQX))
  allocate(W2(6*IQX))
  allocate(W3(6*IQX))
  allocate(W4(6*IQX))
  allocate(W5(6*IQX))
  allocate(W6(6*IQX))
  allocate(VTI(IVX,ISX))
  allocate(XSSI(IVX,ISX))
  allocate(UINV(IVX,ISX))
  allocate(UINV_A(IVX,ISX))
  allocate(UEDG(IVX,ISX))
  allocate(THET(IVX,ISX))
  allocate(DSTR(IVX,ISX))
  allocate(CTAU(IVX,ISX))
  allocate(MASS(IVX,ISX))
  allocate(TAU(IVX,ISX))
  allocate(DIS(IVX,ISX))
  allocate(CTQ(IVX,ISX))
  allocate(DELT(IVX,ISX))
  allocate(TSTR(IVX,ISX))
  allocate(USLP(IVX,ISX))
  allocate(VM(3,IZX,IZX))
  allocate(VA(3,2,IZX))
  allocate(VB(3,2,IZX))
  allocate(VDEL(3,2,IZX))

end subroutine xfoil_init

!=============================================================================80
!
! Initializes xfoil variables from settings
!
!=============================================================================80
subroutine xfoil_defaults(xfoil_options) bind(c, name="xfoil_defaults")

  use xfoil_inc

  type(xfoil_options_type), intent(in) :: xfoil_options

  SILENT_MODE = xfoil_options%silent_mode
  VISCOUS_MODE = xfoil_options%viscous_mode
  MAXIT = xfoil_options%maxit
  N = 0
  PI = 4.d0*atan(1.d0)
  HOPI = 0.5d0/PI
  QOPI = 0.25d0/PI
  DTOR = PI/180.d0
  QINF = 1.d0
  SIG(:) = 0.d0
  QF0(:) = 0.d0
  QF1(:) = 0.d0
  QF2(:) = 0.d0
  QF3(:) = 0.d0
  NW = 0
  RETYP = 1
  MATYP = 1
  GAMMA = 1.4d0
  GAMM1 = GAMMA - 1.d0
  XCMREF = 0.25d0
  YCMREF = 0.d0
  LVISC = xfoil_options%viscous_mode
  AWAKE = 0.d0
  AVISC = 0.d0
  ITMAX = xfoil_options%maxit
  LWDIJ = .false.
  LIPAN = .false.
  LBLINI = .false.
  ACRIT = xfoil_options%ncrit
  IDAMP = 0
  XSTRIP(1) = xfoil_options%xtript
  XSTRIP(2) = xfoil_options%xtripb
  VACCEL = xfoil_options%vaccel
  WAKLEN = 1.d0
  PSIO = 0.d0
  GAMU(:,:) = 0.d0
  GAM(:) = 0.d0
  SIGTE = 0.d0
  GAMTE = 0.d0
  SIGTE_A = 0.d0
  GAMTE_A = 0.d0
  APANEL(:) = 0.d0

! Set boundary layer calibration parameters

  call BLPINI

end subroutine xfoil_defaults

!=============================================================================80
!
! Sets xfoil paneling options
!
!=============================================================================80
subroutine xfoil_set_paneling(geom_options) bind(c, name="xfoil_set_paneling")

  use xfoil_inc, only : NPAN, CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1,   &
                        XPREF2

  type(xfoil_geom_options_type), intent(in) :: geom_options

  NPAN = geom_options%npan
  CVPAR = geom_options%cvpar
  CTERAT = geom_options%cterat
  CTRRAT = geom_options%ctrrat
  XSREF1 = geom_options%xsref1
  XSREF2 = geom_options%xsref2
  XPREF1 = geom_options%xpref1
  XPREF2 = geom_options%xpref2
  
end subroutine xfoil_set_paneling

!=============================================================================80
!
! Sets buffer airfoil for xfoil.
! stat: 0 for success, 1 for failure (xfoil_init not called yet)
!
!=============================================================================80
subroutine xfoil_set_airfoil(xin, zin, npointin, stat)                         &
           bind(c, name="xfoil_set_airfoil")

  use xfoil_inc, only : AIJ, XB, YB, NB

  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin
  integer(c_int), intent(out) :: stat

! Check to make sure xfoil is initialized

  stat = 0
  if (.not. allocated(AIJ)) then
    stat = 1
    return
  end if

  NB = npointin
  XB(1:NB) = xin
  YB(1:NB) = zin

end subroutine xfoil_set_airfoil

!=============================================================================80
!
! Smooths buffer airfoil using Xfoil's PANGEN subroutine
! stat: 0 for success, 1 for failure (xfoil_set_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_smooth_paneling(stat) bind(c, name="xfoil_smooth_paneling")

  use xfoil_inc, only : NB, SILENT_MODE

  integer(c_int), intent(out) :: stat

! Check that buffer airfoil is set

  stat = 0
  if (NB == 0) then
    stat = 1
    return
  end if

! Smooth paneling with PANGEN

  call PANGEN(.NOT. SILENT_MODE)

! Overwrite buffer airfoil

  call GSET

end subroutine xfoil_smooth_paneling

!=============================================================================80
!
! Subroutine to apply a flap deflection to the buffer airfoil and set it as the
! current airfoil. It is recommended to call this after xfoil_smooth_paneling.
! y_flap_spec = 0: specified as y/c
!             = 1: specified as y/local thickness
! stat: 0 for success, 1 for failure (xfoil_set_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_apply_flap_deflection(xflap, yflap, y_flap_spec, degrees,     &
                                       npointout, stat)                        &
           bind(c, name="xfoil_apply_flap_deflection")

  use xfoil_inc, only : NB

  real(c_double), intent(in) :: xflap, yflap, degrees
  integer(c_int), intent(in) :: y_flap_spec
  integer(c_int), intent(out) :: npointout, stat

! Check that buffer airfoil is set

  stat = 0
  if (NB == 0) then
    stat = 1
    return
  end if

! Apply flap deflection

  call FLAP(xflap, yflap, y_flap_spec, degrees)

! Get new buffer airfoil points (may have changed)

  npointout = NB

end subroutine xfoil_apply_flap_deflection

!=============================================================================80
!
! Subroutine to modify the trailing edge gap of the buffer airfoil and set it as
! the current airfoil.
! gap: the new TE gap
! blendloc: x/c location where the shape is first modified to accomodate the gap
!   0 < blendloc < 1
! stat: 0 for success, 1 for failure (xfoil_set_airfoil not called yet)
!
!=============================================================================80
subroutine xfoil_modify_tegap(gap, blendloc, npointout, stat)                  &
           bind(c, name="xfoil_modify_tegap")

  use xfoil_inc, only : NB

  real(c_double), intent(in) :: gap, blendloc
  integer(c_int), intent(out) :: npointout, stat

! Check that buffer airfoil is set

  stat = 0
  if (NB == 0) then
    stat = 1
    return
  end if

! Modify trailing edge gap

  call TGAP(gap, blendloc)

! Get new buffer airfoil points (may have changed)

  npointout = NB

end subroutine xfoil_modify_tegap

!=============================================================================80
!
! Returns current (not buffer) airfoil coordinates from Xfoil
!
!=============================================================================80
subroutine xfoil_get_airfoil(xout, zout, npoint)                               &
           bind(c, name="xfoil_get_airfoil")

  use xfoil_inc, only : X, Y

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: xout, zout

  xout(1:npoint) = X(1:npoint)
  zout(1:npoint) = Y(1:npoint)
   
end subroutine xfoil_get_airfoil

!=============================================================================80
!
! Gets thickness and camber information for the current (not buffer) airfoil
!
!=============================================================================80
subroutine xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)                       &
           bind(c, name="xfoil_geometry_info")

  use xfoil_inc, only : THICKB, XTHICKB, CAMBR, XCAMBR

  real(c_double), intent(out) :: maxt, xmaxt, maxc, xmaxc

  maxt = THICKB
  xmaxt = XTHICKB
  maxc = CAMBR
  xmaxc = XCAMBR 

end subroutine xfoil_geometry_info

!=============================================================================80
!
! Sets Reynolds number for viscous calculations
!
!=============================================================================80
subroutine xfoil_set_reynolds_number(re)                                       &
           bind(c, name="xfoil_set_reynolds_number")

  use xfoil_inc, only : REINF1

  real(c_double), intent(in) :: re

  REINF1 = re

end subroutine xfoil_set_reynolds_number

!=============================================================================80
!
! Sets Mach number
!
!=============================================================================80
subroutine xfoil_set_mach_number(mach) bind(c, name="xfoil_set_mach_number")

  real(c_double), intent(in) :: mach

  call MINFSET(mach)

end subroutine xfoil_set_mach_number

!=============================================================================80
!
! Resets BL initialization flags in xfoil, so BL will be reinitialized at next
! point
!
!=============================================================================80
subroutine xfoil_reinitialize_bl() bind(c, name="xfoil_reinitialize_bl")

  use xfoil_inc, only : LIPAN, LBLINI

  LIPAN = .false.
  LBLINI = .false.

end subroutine xfoil_reinitialize_bl

!=============================================================================80
!
! Runs Xfoil at a specified angle of attack
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
! stat: 0 for success, 1 for failure (xfoil_init not called yet)
!
!=============================================================================80
subroutine xfoil_specal(alpha_spec, alpha, lift, drag, moment, converged, stat)&
           bind(c, name="xfoil_specal")

  use xfoil_inc

  real(c_double), intent(in) :: alpha_spec
  real(c_double), intent(out) :: alpha, lift, drag, moment
  logical(c_bool), intent(out) :: converged
  integer(c_int), intent(out) :: stat

! Check to make sure xfoil is initialized

  stat = 0
  if (.not. allocated(AIJ)) then
    stat = 1
    return
  end if

! Inviscid calculations for specified angle of attack

  converged = .true.
  LALFA = .true.
  ALFA = alpha_spec*DTOR
  call SPECAL
  if (abs(ALFA-AWAKE) .GT. 1.0D-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0D-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0D-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (VISCOUS_MODE) then
    call VISCAL(MAXIT)
    converged = LVCONV
  end if

! Outputs

  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  if (VISCOUS_MODE) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_specal

!=============================================================================80
!
! Runs Xfoil at a specified lift coefficient
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
! stat: 0 for success, 1 for failure (xfoil_init not called yet)
!
!=============================================================================80
subroutine xfoil_speccl(cl_spec, alpha, lift, drag, moment, converged, stat)   &
           bind(c, name="xfoil_speccl")

  use xfoil_inc

  real(c_double), intent(in) :: cl_spec
  real(c_double), intent(out) :: alpha, lift, drag, moment
  logical(c_bool), intent(out) :: converged
  integer(c_int), intent(out) :: stat

! Check to make sure xfoil is initialized

  stat = 0
  if (.not. allocated(AIJ)) then
    stat = 1
    return
  end if

! Inviscid calculations for specified lift coefficient

  converged = .true.
  LALFA = .false.
  ALFA = 0.d0
  CLSPEC = cl_spec
  call SPECCL
  if (abs(ALFA-AWAKE) .GT. 1.0D-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0D-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0D-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (VISCOUS_MODE) then
    call VISCAL(MAXIT)
    converged = LVCONV
  end if

! Outputs

  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  if (VISCOUS_MODE) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_speccl

!=============================================================================80
!
! Returns transition locations on top and bottom in x and z
!
!=============================================================================80
subroutine xfoil_get_transloc(xtranst, ztranst, xtransb, ztransb)              &
           bind(c, name="xfoil_get_transloc")

  use xfoil_inc, only : XOCTR, YOCTR

  real(c_double), intent(out) :: xtranst, ztranst, xtransb, ztransb

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
subroutine xfoil_get_cp(npoint, cp) bind(c, name="xfoil_get_cp")

  use xfoil_inc, only : CPI, CPV, VISCOUS_MODE

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: cp

  if (VISCOUS_MODE) then
    cp(1:npoint) = CPV(1:NPOINT)
  else
    cp(1:npoint) = CPI(1:NPOINT)
  end if

end subroutine xfoil_get_cp

!=============================================================================80
!
! Returns skin friction coefficient on surface
!
!=============================================================================80
subroutine xfoil_get_cf(npoint, cf) bind(c, name="xfoil_get_cf")

  use xfoil_inc, only : TAU, QINF, IPAN, NBL

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: cf

  integer(c_int) :: is, ibl, i
  real(c_double) :: que

  que = 0.5d0*QINF**2.d0
  
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
subroutine xfoil_get_uedge(npoint, uedge) bind(c, name="xfoil_get_uedge")

  use xfoil_inc, only : UEDG, TKLAM, QINF, IPAN, NBL

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: uedge

  integer(c_int) :: is, ibl, i
  real(c_double) :: uei

! Populate uedge array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        uei = UEDG(ibl,is)
        uedge(i) = uei * (1.d0-TKLAM) / (1.d0 - TKLAM*(uei/QINF)**2.d0)
      end if

    end do
  end do

end subroutine xfoil_get_uedge

!=============================================================================80
!
! Returns BL displacement thickness on surface
!
!=============================================================================80
subroutine xfoil_get_deltastar(npoint, deltastar)                              &
           bind(c, name="xfoil_get_deltastar")

  use xfoil_inc, only : DSTR, IPAN, NBL

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: deltastar

  integer(c_int) :: is, ibl, i

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
subroutine xfoil_get_diss(npoint, diss) bind(c, name="xfoil_get_diss")

  use xfoil_inc, only : DIS, QINF, IPAN, NBL

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: diss

  integer(c_int) :: is, ibl, i

! Populate diss array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        diss(i) = DIS(ibl,is) / QINF**3.d0
      end if

    end do
  end do

end subroutine xfoil_get_diss

!=============================================================================80
!
! Returns BL kinematic shape parameter on surface
!
!=============================================================================80
subroutine xfoil_get_hk(npoint, hk) bind(c, name="xfoil_get_hk")

  use xfoil_inc, only : THET, DSTR, UEDG, TKLAM, QINF, GAMM1, IPAN, NBL
  use xbl_inc,   only : HSTINV

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: hk

  integer(c_int) :: is, ibl, i
  real(c_double) :: thi, dsi, uei, uc, amsq, dummy 

! Populate hk array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        thi = THET(ibl,is)
        dsi = DSTR(ibl,is)
        uei = UEDG(ibl,is)
        uc = uei * (1.d0-TKLAM) / (1.d0 - TKLAM*(uei/QINF)**2.d0) 
        amsq = uc*uc*HSTINV / (GAMM1*(1.d0 - 0.5d0*uc*uc*HSTINV))
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
subroutine xfoil_get_retheta(npoint, retheta) bind(c, name="xfoil_get_retheta")

  use xfoil_inc, only : UEDG, QINF, TKLAM, GAMM1, REINF, THET, IPAN, NBL, IVX
  use xbl_inc,   only : HSTINV

  integer(c_int), intent(in) :: npoint
  real(c_double), dimension(npoint), intent(out) :: retheta

  integer(c_int) :: is, ibl, i
  real(c_double) :: uei, ue, herat, rhoe, amue 

! Sutherland's constant/To (assumes stagnation conditions are at STP)

  real(c_double), parameter :: hvrat = 0.35d0
  
! Populate ampl array, going over upper surface and then lower surface

  do is = 1, 2
    do ibl = 2, NBL(is)
      i = IPAN(ibl,is)

!     Xfoil BL arrays include wake; only accept surface points here

      if (i <= npoint) then
        uei = UEDG(ibl,is)
        ue = uei * (1.d0-TKLAM) / (1.d0 - TKLAM*(uei/QINF)**2.d0) 
        herat = (1.d0 - 0.5d0*HSTINV*uei**2.d0)                                &
              / (1.d0 - 0.5d0*HSTINV*QINF**2.d0)
        rhoe = herat**(1.d0/GAMM1)
        amue = sqrt(herat**3.d0) * (1.d0+hvrat)/(herat+hvrat)
        retheta(i) = REINF * rhoe*ue*THET(ibl,is)/amue
      end if

    end do
  end do

end subroutine xfoil_get_retheta

!=============================================================================80
!
! Deallocates memory in xfoil
!
!=============================================================================80
subroutine xfoil_cleanup() bind(c, name="xfoil_cleanup")

  use xfoil_inc

! Don't bother if xfoil is not initialized

  if (.not. allocated(AIJ)) return

! Deallocate variables

  deallocate(AIJ)
  deallocate(BIJ)
  deallocate(DIJ)
  deallocate(CIJ)
  deallocate(IPAN)
  deallocate(ISYS)
  deallocate(W1)
  deallocate(W2)
  deallocate(W3)
  deallocate(W4)
  deallocate(W5)
  deallocate(W6)
  deallocate(VTI)
  deallocate(XSSI)
  deallocate(UINV)
  deallocate(UINV_A)
  deallocate(UEDG)
  deallocate(THET)
  deallocate(DSTR)
  deallocate(CTAU)
  deallocate(MASS)
  deallocate(TAU)
  deallocate(DIS)
  deallocate(CTQ)
  deallocate(DELT)
  deallocate(TSTR)
  deallocate(USLP)
  deallocate(VM)
  deallocate(VA)
  deallocate(VB)
  deallocate(VDEL)

end subroutine xfoil_cleanup

!=============================================================================80
!
! Subroutine to get Cl, Cd, Cm for an airfoil from Xfoil at given operating
! conditions.  Reynolds numbers and mach numbers should be specified for each
! operating point.  Additionally, op_mode determines whether each point is run
! at a constant alpha or cl - use 0 for specified alpha and 1 for specified cl.
!
! This is a convenience method to run xfoil at a bunch of different operating
! points, optionally with changing flap deflections and ncrit values. Still
! requires xfoil_init, xfoil_defaults, and xfoil_set_paneling to be called
! first.
! 
! Outputs:
!   alpha, Cl, Cd, Cm each operating point
!   viscrms: rms for viscous calculations (check for convergence)
!
!=============================================================================80
subroutine run_xfoil(npointin, xin, zin, noppoint, operating_points, op_modes, &
                     reynolds_numbers, mach_numbers, use_flap, x_flap, y_flap, &
                     y_flap_spec, flap_degrees, reinitialize, fix_unconverged, &
                     lift, drag, moment, viscrms, alpha, xtrt, xtrb, stat,     &
                     ncrit_per_point) bind(c, name="run_xfoil")

  use xfoil_inc

  integer(c_int), intent(in) :: npointin, noppoint
  real(c_double), dimension(npointin), intent(in) :: xin, zin
  real(c_double), dimension(noppoint), intent(in) :: operating_points,         &
                                    reynolds_numbers, mach_numbers, flap_degrees
  real(c_double), intent(in) :: x_flap, y_flap
  integer(c_int), intent(in) :: y_flap_spec
  logical(c_bool), intent(in) :: use_flap, reinitialize, fix_unconverged
  integer(c_int), dimension(noppoint), intent(in) :: op_modes
  real(c_double), dimension(noppoint), intent(out) :: lift, drag, moment,      &
                                                      viscrms
  real(c_double), dimension(noppoint), intent(out) :: alpha, xtrt, xtrb
  integer(c_int), intent(out) :: stat
  real(c_double), dimension(noppoint), intent(in), optional :: ncrit_per_point

  integer(c_int) :: i, dummy 
  logical(c_bool), dimension(noppoint) :: point_converged, point_fixed 
  real(c_double) :: newpoint, ztrt, ztrb
  character(30) :: text
  character(150) :: message

  if (.not. SILENT_MODE) then
    write(*,*) 
    write(*,*) 'Analyzing aerodynamics using the XFOIL engine ...'
  end if 

  point_converged(:) = .true.
  point_fixed(:) = .false.

! Set airfoil and smooth paneling

  if (.not. use_flap) then
    call xfoil_set_airfoil(xin, zin, npointin, stat)
    if (stat /= 0) return
    call xfoil_smooth_paneling(stat)
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
      call xfoil_set_airfoil(xin, zin, npointin, stat)
      call xfoil_smooth_paneling(stat)
      call xfoil_apply_flap_deflection(x_flap, y_flap, y_flap_spec,            &
                                       flap_degrees(i), dummy, stat)
    end if

    call xfoil_set_reynolds_number(reynolds_numbers(i))
    call xfoil_set_mach_number(mach_numbers(i))

    if (reinitialize) call xfoil_reinitialize_bl()

!   Set ncrit per point

    if (present(ncrit_per_point)) ACRIT = ncrit_per_point(i)

    if (op_modes(i) == 0) then

      call xfoil_specal(operating_points(i), alpha(i), lift(i), drag(i),       &
                        moment(i), point_converged(i), stat)

    elseif (op_modes(i) == 1) then

      call xfoil_speccl(operating_points(i), alpha(i), lift(i), drag(i),       &
                        moment(i), point_converged(i), stat)

    else

      write(*,*)
      write(*,*) "Error in xfoil_interface: op_mode must be 0 or 1."
      write(*,*)
      stop

    end if

!   Additional outputs

    call xfoil_get_transloc(xtrt(i), ztrt, xtrb(i), ztrb)

!   Handling of unconverged points

    if (VISCOUS_MODE .and. .not. point_converged(i)) then

      if (fix_unconverged) then

!       Try to initialize BL at new point (in the direction away from stall)

        newpoint = operating_points(i) - 0.5d0*abs(operating_points(i))*sign(  &
                                                   1.d0, operating_points(i))
        if (newpoint == 0.d0) newpoint = 0.1d0

        call xfoil_reinitialize_bl()
        if (op_modes(i) == 0) then
          call xfoil_specal(newpoint, alpha(i), lift(i), drag(i), moment(i),   &
                            point_converged(i), stat)
        else
          call xfoil_speccl(newpoint, alpha(i), lift(i), drag(i), moment(i),   &
                            point_converged(i), stat)
        end if

!       Now try to run again at the old operating point

        if (op_modes(i) == 0) then
          call xfoil_specal(operating_points(i), alpha(i), lift(i), drag(i),   &
                            moment(i), point_converged(i), stat)
        else
          call xfoil_speccl(operating_points(i), alpha(i), lift(i), drag(i),   &
                            moment(i), point_converged(i), stat)
        end if

        if (point_converged(i)) point_fixed(i) = .true.

        call xfoil_get_transloc(xtrt(i), ztrt, xtrb(i), ztrb)

      end if
  end if

!   Convergence check

    viscrms(i) = RMSBL

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

  if (.not. SILENT_MODE) then

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
