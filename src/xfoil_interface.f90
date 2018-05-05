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
    logical(c_bool) :: fix_unconverged  !Reinitialize to fix unconverged pts.
    logical(c_bool) :: reinitialize     !Reinitialize BLs at every operating
                                        !  point (recommended for optimization)

  end type xfoil_options_type

  type, bind(c) :: xfoil_geom_options_type

    integer(c_int) :: npan
    real(c_double) :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2

  end type xfoil_geom_options_type

  contains

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
! Subroutine to smooth an airfoil using Xfoil's PANGEN subroutine. Note that
! geom_options%npan is ignored in favor of npointout.
!
!=============================================================================80
subroutine smooth_paneling(xin, zin, npointin, npointout, geom_options, xout,  &
                           zout) bind(c, name="smooth_paneling")

  use xfoil_inc

  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin, npointout
  type(xfoil_geom_options_type), intent(in) :: geom_options
  real(c_double), dimension(npointout), intent(out) :: xout, zout
  
  integer(c_int) :: i

! Check to make sure xfoil is initialized

  if (.not. allocated(AIJ)) then
    write(*,*) "Error: xfoil is not initialized!  Call xfoil_init() first."
    stop
  end if

! Set xfoil airfoil and paneling options

  call xfoil_set_airfoil(xin, zin, npointin)
  call xfoil_set_paneling(geom_options)

! Smooth paneling with PANGEN

  call PANGEN(.NOT. SILENT_MODE)

! Save smoothed airfoil coordinates

  do i = 1, npointout
    xout(i) = X(i)
    zout(i) = Y(i)
  end do

end subroutine smooth_paneling

!=============================================================================80
!
! Subroutine to apply a flap deflection to the buffer airfoil and set it as the
! current airfoil.  For best results, this should be called after PANGEN.
! y_flap_spec = 0: specified as y/c
!             = 1: specified as y/local thickness
!
!=============================================================================80
subroutine xfoil_apply_flap_deflection(xflap, yflap, y_flap_spec, degrees)     &
           bind(c, name="xfoil_apply_flap_deflection")

  real(c_double), intent(in) :: xflap, yflap, degrees
  integer(c_int), intent(in) :: y_flap_spec

! Apply flap deflection

  call FLAP(xflap, yflap, y_flap_spec, degrees)

end subroutine xfoil_apply_flap_deflection

!=============================================================================80
!
! Subroutine to modify the trailing edge gap of the buffer airfoil and set it as
! the current airfoil.
! gap: the new TE gap
! blendloc: x/c location where the shape is first modified to accomodate the gap
!   0 < blendloc < 1
!
!=============================================================================80
subroutine xfoil_modify_tegap(gap, blendloc) bind(c, name="xfoil_modify_tegap")

  real(c_double), intent(in) :: gap, blendloc

! Modify trailing edge gap

  call TGAP(gap, blendloc)

end subroutine xfoil_modify_tegap

!=============================================================================80
!
! Subroutine to get Cl, Cd, Cm for an airfoil from Xfoil at given operating
! conditions.  Reynolds numbers and mach numbers should be specified for each
! operating point.  Additionally, op_mode determines whether each point is run
! at a constant alpha or cl - use 0 for specified alpha and 1 for specified cl.
! 
! Outputs:
!   alpha, Cl, Cd, Cm each operating point
!   viscrms: rms for viscous calculations (check for convergence)
!
!=============================================================================80
subroutine run_xfoil(npointin, xin, zin, geom_options, noppoint,               &
                     operating_points, op_modes, reynolds_numbers,             &
                     mach_numbers, use_flap, x_flap, y_flap, y_flap_spec,      &
                     flap_degrees, xfoil_options, lift, drag, moment, viscrms, &
                     alpha, xtrt, xtrb, ncrit_per_point)                       &
           bind(c, name="run_xfoil")

  use xfoil_inc

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

  integer(c_int) :: i
  logical(c_bool), dimension(noppoint) :: point_converged, point_fixed 
  real(c_double) :: newpoint
  character(30) :: text
  character(150) :: message

  if (.not. xfoil_options%silent_mode) then
    write(*,*) 
    write(*,*) 'Analyzing aerodynamics using the XFOIL engine ...'
  end if 

! Check to make sure xfoil is initialized

  if (.not. allocated(AIJ)) then
    write(*,*) "Error: xfoil is not initialized!  Call xfoil_init() first."
    stop
  end if

! Set default Xfoil parameters

  call xfoil_defaults(xfoil_options)

  point_converged(:) = .true.
  point_fixed(:) = .false.

! Set paneling options

  call xfoil_set_paneling(geom_options)

! Set airfoil and smooth paneling

  if (.not. use_flap) then
    call xfoil_set_airfoil(xin, zin, npointin)
    call PANGEN(.not. SILENT_MODE)
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
      call xfoil_set_airfoil(xin, zin, npointin)
      call PANGEN(.not. SILENT_MODE)
      call xfoil_apply_flap_deflection(x_flap, y_flap, y_flap_spec,            &
                                       flap_degrees(i))
    end if

    REINF1 = reynolds_numbers(i)
    call MINFSET(mach_numbers(i))

    if (xfoil_options%reinitialize) then
      LIPAN = .false.
      LBLINI = .false.
    end if

!   Set compressibility parameters from MINF

    CALL COMSET

!   Set ncrit per point

    if (present(ncrit_per_point)) ACRIT = ncrit_per_point(i)

    if (op_modes(i) == 0) then

      call xfoil_specal(operating_points(i), xfoil_options%viscous_mode,       &
                        xfoil_options%maxit, lift(i), drag(i), moment(i))

    elseif (op_modes(i) == 1) then

      call xfoil_speccl(operating_points(i), xfoil_options%viscous_mode,       &
                        xfoil_options%maxit, lift(i), drag(i), moment(i))

    else

      write(*,*)
      write(*,*) "Error in xfoil_interface: op_mode must be 0 or 1."
      write(*,*)
      stop

    end if

!   Additional outputs

    alpha(i) = ALFA/DTOR
    xtrt(i) = XOCTR(1)
    xtrb(i) = XOCTR(2)

!   Handling of unconverged points

    if (xfoil_options%viscous_mode .and. .not. LVCONV) then

      point_converged(i) = .false.

      if (xfoil_options%fix_unconverged) then

!       Try to initialize BL at new point (in the direction away from stall)

        newpoint = operating_points(i) - 0.5d0*abs(operating_points(i))*sign(  &
                                                   1.d0, operating_points(i))
        if (newpoint == 0.d0) newpoint = 0.1d0

        LIPAN = .false.
        LBLINI = .false.
        if (op_modes(i) == 0) then
          call xfoil_specal(newpoint, xfoil_options%viscous_mode,              & 
                            xfoil_options%maxit, lift(i), drag(i), moment(i))
        else
          call xfoil_speccl(newpoint, xfoil_options%viscous_mode,              & 
                            xfoil_options%maxit, lift(i), drag(i), moment(i))
        end if

!       Now try to run again at the old operating point

        if (op_modes(i) == 0) then
          call xfoil_specal(operating_points(i), xfoil_options%viscous_mode,   &
                            xfoil_options%maxit, lift(i), drag(i), moment(i))
        else
          call xfoil_speccl(operating_points(i), xfoil_options%viscous_mode,   &
                            xfoil_options%maxit, lift(i), drag(i), moment(i))
        end if

        if (LVCONV) point_fixed(i) = .true.

        alpha(i) = ALFA/DTOR
        xtrt(i) = XOCTR(1)
        xtrb(i) = XOCTR(2)

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

  if (.not. xfoil_options%silent_mode) then

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
! Runs Xfoil at a specified angle of attack
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!=============================================================================80
subroutine xfoil_specal(angle_of_attack, viscous_mode, maxit, lift, drag,      &
                        moment) bind(c, name="xfoil_specal")

  use xfoil_inc

  real(c_double), intent(in) :: angle_of_attack
  logical(c_bool), intent(in) :: viscous_mode
  integer(c_int), intent(in) :: maxit
  real(c_double), intent(out) :: lift, drag, moment

! Inviscid calculations for specified angle of attack

  LALFA = .TRUE.
  ALFA = angle_of_attack*DTOR
  call SPECAL
  if (abs(ALFA-AWAKE) .GT. 1.0D-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0D-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0D-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
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
!
!=============================================================================80
subroutine xfoil_speccl(cl_spec, viscous_mode, maxit, lift, drag, moment)      &
           bind(c, name="xfoil_speccl")

  use xfoil_inc

  real(c_double), intent(in) :: cl_spec
  logical(c_bool), intent(in) :: viscous_mode
  integer(c_int), intent(in) :: maxit
  real(c_double), intent(out) :: lift, drag, moment

! Inviscid calculations for specified lift coefficient

  LALFA = .FALSE.
  ALFA = 0.d0
  CLSPEC = cl_spec
  call SPECCL
  if (abs(ALFA-AWAKE) .GT. 1.0D-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0D-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0D-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_speccl

!=============================================================================80
!
! Allocates xfoil variables that may be too big for the stack in OpenMP
!
!=============================================================================80
subroutine xfoil_init() bind(c, name="xfoil_init")

  use xfoil_inc

! Allocate variables that may be too big for the stack in OpenMP

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
! Initializes xfoil variables
!
!=============================================================================80
subroutine xfoil_defaults(xfoil_options) bind(c, name="xfoil_defaults")

  use xfoil_inc

  type(xfoil_options_type), intent(in) :: xfoil_options

  N = 0
  SILENT_MODE = xfoil_options%silent_mode
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
! Sets airfoil for xfoil
!
!=============================================================================80
subroutine xfoil_set_airfoil(xin, zin, npointin)                               &
           bind(c, name="xfoil_set_airfoil")

  use xfoil_inc, only : XB, YB, NB

  real(c_double), dimension(npointin), intent(in) :: xin, zin
  integer(c_int), intent(in) :: npointin

  NB = npointin
  XB(1:NB) = xin
  YB(1:NB) = zin

end subroutine xfoil_set_airfoil

!=============================================================================80
!
! Returns current airfoil coordinates from Xfoil
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
! Deallocates memory in xfoil
!
!=============================================================================80
subroutine xfoil_cleanup() bind(c, name="xfoil_cleanup")

  use xfoil_inc

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
! Gets thickness and camber information for the current airfoil
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
