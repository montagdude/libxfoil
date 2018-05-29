program main

  use iso_c_binding
  use xfoil_interface

  implicit none

  type(xfoil_options_type) :: opts
  type(xfoil_geom_options_type) :: geom_opts
  integer :: npoint
  integer, parameter :: noppoint = 8
  double precision, dimension(:), allocatable :: x, z, s, xs, zs
  double precision :: xle, zle, sle
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
  logical(c_bool) :: fix_unconverged = .false.
  logical(c_bool) :: reinitialize = .false.
  integer(c_int) :: stat
  logical(c_bool) converged
  type(xfoil_data_group) :: xdg

  opts%ncrit = 9.d0
  opts%xtript = 1.d0
  opts%xtripb = 1.d0
  opts%viscous_mode = .true.
  opts%silent_mode = .true.
  opts%maxit = 100
  opts%vaccel = 0.01d0

  geom_opts%npan = 200
  geom_opts%cvpar = 1.d0
  geom_opts%cterat = 0.15d0
  geom_opts%ctrrat = 0.2d0
  geom_opts%xsref1 = 1.d0
  geom_opts%xsref2 = 1.d0
  geom_opts%xpref1 = 1.d0
  geom_opts%xpref2 = 1.d0

  ! Read airfoil from file
  call read_airfoil_points("clarky.dat", npoint)
  allocate(x(npoint))
  allocate(z(npoint))
  allocate(s(npoint))
  allocate(xs(npoint))
  allocate(zs(npoint))
  call read_airfoil("clarky.dat", npoint, x, z)

  ! Get leading edge location
  call xfoil_spline_coordinates(x, z, npoint, s, xs, zs)
  call xfoil_lefind(x, z, s, xs, zs, npoint, sle, xle, zle)
  write(*,'(A14,F8.5,A2,F8.5)') "Leading edge: ", xle, ", ", zle

  ! Set up inputs
  call xfoil_init(xdg)
  call xfoil_defaults(xdg, opts)
  call xfoil_set_paneling(xdg, geom_opts)
  call xfoil_set_buffer_airfoil(xdg, x, z, npoint)
  call xfoil_smooth_paneling(xdg, stat)
  if (stat /= 0) then
    write(*,*) "Error smoothing paneling."
    stop
  end if

  ! Run at all operating points
  write(*,*) "Running Xfoil ..."
  do i = 1, noppoint
    call xfoil_set_reynolds_number(xdg, re(i))
    call xfoil_set_mach_number(xdg, mach(i))
    call xfoil_specal(xdg, oppoints(i), alpha(i), lift(i), drag(i), moment(i), &
                      converged, stat)
    if (stat /= 0) then
      write(*,*) "Error running Xfoil."
      stop
    end if
    if (.not. converged) then
      write(*,"(A7,I1,A18)") " Point ", i, " did not converge."
    else
      write(*,"(A7,I1,A11)") " Point ", i, " converged."
      write(*,'(A7,I1,A8,F7.4,A7,F7.4,A7,F7.4,A7,F7.4)')                       &
        " Point ", i, ": AoA = ", alpha(i), ", Cl = ", lift(i), ", Cd = ",     &
        drag(i), ", Cm = ", moment(i)
    end if
  end do
  call xfoil_cleanup(xdg)

  ! Now, do the same thing but with the run_xfoil method
  write(*,*)
  write(*,*) "Running Xfoil using the run_xfoil method ..."
  call run_xfoil(npoint, x, z, geom_opts, noppoint, oppoints, opmodes, re,     &
                 mach, use_flap, 0.d0, 0.d0, 0, flapang, opts, reinitialize,   &
                 fix_unconverged, lift, drag, moment, viscrms, alpha, xtrt,    &
                 xtrb, stat)
  if (stat /= 0) then
    write(*,*) "Error running Xfoil."
    stop
  end if

  do i = 1, noppoint
    write(*,'(A6,I1,A8,F7.4,A7,F7.4,A7,F7.4,A7,F7.4)')                         &
      "Point ", i, ": AoA = ", alpha(i), ", Cl = ", lift(i), ", Cd = ",        &
      drag(i), ", Cm = ", moment(i)
  end do

  deallocate(x)
  deallocate(z)
  deallocate(s)
  deallocate(xs)
  deallocate(zs)

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
