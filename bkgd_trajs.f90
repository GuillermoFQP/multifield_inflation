program bkgd_trajs
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric. This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$

use multifield_globals
use multifield_utils

implicit none

real, parameter    :: N_bound = 200.0                     ! Upper bound in N
real, parameter    :: dN = 0.01                           ! Data flushing period
real, parameter    :: dt = 500.0                          ! Time-step
integer, parameter :: ntrajsx = 6                         ! Size of the grid of initial conditions
integer, parameter :: ntrajsy = 6                         ! Size of the grid of initial conditions
real, dimension(6) :: y                                   ! State array
real, dimension(2) :: phi, phidot, phi_min, phi_max, dphi ! Colormap grid parameters
real               :: H, Hdot, N, slowroll, N_flush       ! Variables
integer            :: i, j, k, u                          ! Indices
character(len=32)  :: arg                                 ! Command-line argument
character(len=100) :: filename                            ! Output file name
logical            :: condition                           ! Loop condition

! Parse argument
if (field_space_geometry == "RPM") then
	call get_command_argument(1, arg)
	read (arg, *) energyscale
end if

! Create destination directory
call execute_command_line('rm -rf trajs; mkdir trajs')

! Define the range for initial conditions
select case (potential_shape)
	case ("ELP")
		phi_min = [-20.0, -20.0]
		phi_max = [ 20.0,  20.0]
	case ("NLP")
		phi_min = [-20.0, -20.0]
		phi_max = [ 20.0,  20.0]
	case ("HYP")
		phi_min = [-20.0, -20.0]
		phi_max = [ 20.0,  20.0]
	case ("MVP")
		phi_min = [-20.0, -20.0]
		phi_max = [ 20.0,  20.0]
end select

! Define the step size (NGRID-1 intervals between NGRID points)
dphi = [(phi_max(1) - phi_min(1)) / real(ntrajsx-1), (phi_max(2) - phi_min(2)) / real(ntrajsy-1)]

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(dynamic) PRIVATE(i, j, phi, phidot, H, N, slowroll, N_flush, k, u, filename, y) SHARED(phi_min, phi_max, dphi)
do j = 1, ntrajsy
	do i = 1, ntrajsx
		! Background initial conditions
		phi    = [phi_min(1) + real(i-1) * dphi(1), phi_min(2) + real(j-1) * dphi(2)] ! $\phi^{A}(t_{0})$
		phidot = [0.0, 0.0]                                                           ! $\dot{\phi}^{A}(t_{0})$
		H      = Hubble(phi, phidot)                                                  ! $H(t_{0})$
		N      = 0.0                                                                  ! $N(t_{0})$
		
		! Initialize background arrays
		call pack_state_background(y, phi, phidot, H, N)
		
		slowroll = 0.0                 ! Initialize slow-roll parameter $\epsilon(t_{0})$
		N_flush  = 0.0                 ! Data writing trigger
		k        = (j-1) * ntrajsx + i ! Output file number
		u        = 10 + k              ! Unit number
		
		write (filename, '("trajs/traj_", I0, ".txt")') k
		open (unit=u, file=filename, status='replace', action='write')
		
		condition = .true.
		
!		do while (slowroll <= 1.0 .or. abs(abs(phi(2))-2.5) >= 1.0d-3 .or. abs(phi(1)) >= 1.0d-3) ! Condition for hybrid potential
!		do while (slowroll <= 1.0 .or. abs(phi(2)) >= 1.0d-3 .or. abs(phi(1)) >= 1.0d-3)          ! Condition for Elliptic potential
!		do while (slowroll <= 1.0 .or. abs(phi(1)) >= 1.0d-3)                                     ! Condition for non-linear potential
!		do while (N <= N_bound)
		do while (condition)
			! Update functions of time
			call unpack_state_background(y, phi, phidot, H, N)
			
			! Update slow-roll parameter
			slowroll = - Hubbledot(phi, phidot) / H**2
			
			if (N >= N_flush) then
				write (u,'(7(6e25.10e3))') phi, potential(phi), N, slowroll
				N_flush = N_flush + dN
			end if
			
			call gl8_background(y, dt)
			
			! Update condition
			select case (potential_shape)
				case ("ELP")
					condition = slowroll <= 1.0 .or. abs(phi(2)) >= 1.0d-3 .or. abs(phi(1)) >= 1.0d-3
				case ("NLP")
					condition = slowroll <= 1.0 .or. abs(phi(1)) >= 1.0d-3
				case ("HYP")
					condition = slowroll <= 1.0 .or. abs(abs(phi(2))-2.5) >= 1.0d-3 .or. abs(phi(1)) >= 1.0d-3
				case ("MVP")
					condition = slowroll <= 1.0 .or. abs(phi(1)) >= 5.0d-1 .or. abs(phi(2)) >= 5.0d-1
			end select
		end do
		
		close (u)
	end do
end do
!$OMP END PARALLEL DO

end program bkgd_trajs
