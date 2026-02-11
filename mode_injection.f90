program mode_injection
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric. This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$

use multifield_globals
use multifield_utils

implicit none

real, parameter    :: N_bound = 100.0             ! Upper bound in N
real, parameter    :: dN = 0.01                   ! Data flushing period
real, parameter    :: dt = 100.0                  ! Time-step
real, dimension(6) :: y                           ! State array
real, dimension(2) :: phi, phidot                 ! Field multiplet
real               :: H, Hdot, H_end, N, slowroll ! Variables
real               :: N_flush, k_phys             ! Variables
integer            :: i, j, k                     ! Indices
character(len=32)  :: arg                         ! Command-line argument
character(len=100) :: filename                    ! Output file name
logical            :: condition                   ! Loop condition

! Parse argument
if (field_space_geometry == "RPM") then
	call get_command_argument(1, arg)
	read (arg, *) energyscale
end if

if (field_space_geometry == "GBM") then
	left_GBM = [ &
20.0, -0.85, 20.0, -0.85, -0.85, &
20.0, 20.0, -0.85, 20.0, 20.0, &
-0.85, 20.0, -0.85, 20.0, -0.85, &
-0.85, 20.0, 20.0, -0.85, 20.0, &
-0.85, 20.0, 20.0, -0.85, -0.85, &
20.0, -0.85, 20.0, 20.0, -0.85 ]
	amp_GBM  = [ left_GBM, left_GBM(n_GBM/2:1:-1) ]
end if

! Background initial conditions
phi      = initialize_phi()    ! $\phi^{A}$
phidot   = 0.0                 ! $\dot{\phi}^{A}(t_{0})$
H        = Hubble(phi, phidot) ! $H(t_{0})$
N        = 0.0                 ! $N(t_{0})$
slowroll = 0.0                 ! Initialize slow-roll parameter $\epsilon(t_{0})$
N_flush  = 0.0                 ! Data writing trigger

! Initialize background arrays
call pack_state_background(y, phi, phidot, H, N)

condition = .true.

!do while (slowroll <= 1.0)
do while (condition)
	! Update functions of time
	call unpack_state_background(y, phi, phidot, H, N)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	call gl8_background(y, dt)

	if (slowroll >= 1.0 .and. convergence(y(6), N)) condition = .false.
end do

H_end  = H             ! Hubble parameter at $\epsilon = 1$
k_phys = 1.0d5 * H_end ! Surface of constant $k_{\mathrm{phys}}$

! Background initial conditions
phi      = initialize_phi()    ! $\phi^{A}$
phidot   = 0.0                 ! $\dot{\phi}^{A}(t_{0})$
H        = Hubble(phi, phidot) ! $H(t_{0})$
N        = 0.0                 ! $N(t_{0})$
slowroll = 0.0                 ! Initialize slow-roll parameter $\epsilon(t_{0})$
N_flush  = 0.0                 ! Data writing trigger

condition = .true.

! Initialize background arrays
call pack_state_background(y, phi, phidot, H, N)

do while (condition) ! (N <= N_bound)
	! Update functions of time
	call unpack_state_background(y, phi, phidot, H, N)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	if (N >= N_flush) then
		write (*, '(7(6e25.10e3))') N, log(1.0 / H), log(2.0 * pi / k_phys), slowroll
		N_flush = N_flush + dN
	end if
!	write (*, '(7(6e25.10e3))') phi, phidot, H, N, slowroll
!	write (*, '(7(6e25.10e3))') N, log(1.0 / H), log(2.0 * pi / k_phys), slowroll
	
	call gl8_background(y, dt)
	
	if (slowroll >= 1.0 .and. convergence(y(6), N)) condition = .false.
!	if (slowroll >= 1.0) condition = .false.
end do

end program mode_injection
