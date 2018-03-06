program main

  use iso_c_binding

  implicit none

  ! For compiling purposes, it is necessary to add interfaces for all the
  ! subroutines in the xfoil library. These are provided for you in the
  ! following file, which is essentially the Fortran equivalent of a C header
  ! file. Be sure that it is in your include path when compiling.
  include "xfoil_interface_wrap.f90"

  type(xfoil_options_type) :: opts
  type(xfoil_geom_options_type) :: geom_opts
  integer :: npoint
  integer, parameter :: noppoint = 8
  double precision, dimension(:), allocatable :: x, z
  double precision :: xle, zle
  double precision, dimension(noppoint) :: &
    oppoints = (/-6.d0, -3.d0, 0.d0, 3.d0, 6.d0, 9.d0, 12.d0, 15.d0/)
  integer, dimension(noppoint) :: opmodes = (/0, 0, 0, 0, 0, 0, 0, 0/)
  double precision, dimension(noppoint) :: re = 3.D+05
  double precision, dimension(noppoint) :: mach = 0.3d0
  double precision, dimension(noppoint) :: flapang = 0.d0
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms, alpha, &
                                           xtrt, xtrb
  integer i 
  logical(c_bool) :: use_flap = .false.

  opts%ncrit = 9.d0
  opts%xtript = 1.d0
  opts%xtripb = 1.d0
  opts%viscous_mode = .true.
  opts%silent_mode = .false.
  opts%maxit = 100
  opts%vaccel = 0.01d0
  opts%fix_unconverged = .false.
  opts%reinitialize = .false.

  geom_opts%npan = 200
  geom_opts%cvpar = 1.d0
  geom_opts%cterat = 0.15d0
  geom_opts%ctrrat = 0.2d0
  geom_opts%xsref1 = 1.d0
  geom_opts%xsref2 = 1.d0
  geom_opts%xpref1 = 1.d0
  geom_opts%xpref2 = 1.d0

  call read_airfoil_points("clarky.dat", npoint)
  allocate(x(npoint))
  allocate(z(npoint))
  call read_airfoil("clarky.dat", npoint, x, z)

  call xfoil_init()
  call xfoil_lefind(x, z, npoint, xle, zle)
  write(*,'(A14,F8.5,A2,F8.5)') "Leading edge: ", xle, ", ", zle
  call run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,     &
                 mach, use_flap, 0.d0, 0.d0, 0, flapang, opts, lift, drag,     &
                 moment, viscrms, alpha, xtrt, xtrb)
  call xfoil_cleanup()

  write(*,*)
  do i = 1, noppoint
    write(*,'(A6,I1,A8,F7.4,A7,F7.4,A7,F7.4,A7,F7.4)') &
      "Point ", i, ": AoA = ", alpha(i), ", Cl = ", lift(i), ", Cd = ",        &
      drag(i), ", Cm = ", moment(i)
  end do

  deallocate(x)
  deallocate(z)

end program main

subroutine read_airfoil_points(filename, npoint)

  character(*), intent(in) :: filename
  integer, intent(out) :: npoint

  integer :: iunit, stat

  iunit = 12
  open(iunit, file=trim(filename), status="old", iostat=stat)
  if (stat /= 0) then
    write(*,*) "Error opening file "//trim(filename)//"."
    stop
  end if 

  npoint = 0
  do while(.true.)
    read(iunit,*,end=500)
    npoint = npoint + 1
  end do

500 close(iunit)

end subroutine read_airfoil_points

subroutine read_airfoil(filename, npoint, x, z)

  character(*), intent(in) :: filename
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(out) :: x, z

  integer :: i, iunit, stat

  iunit = 12
  open(iunit, file=trim(filename), status="old", iostat=stat)
  if (stat /= 0) then
    write(*,*) "Error opening file "//trim(filename)//"."
    stop
  end if 

  do i = 1, npoint
    read(iunit,*) x(i), z(i)
  end do

  close(iunit)

end subroutine read_airfoil
