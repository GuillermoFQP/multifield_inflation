program ps_full
! Solves the set of coupled background and perturbation equations for an inflation model with two scalar fields with a non-trivial field space metric and calculates the curvature power spectrum of the perturbation fields $\mathcal{P}_{\mathcal{R}}$ as a function of $\log(\frac{k}{aH})$. This script evolves 28 functions of time stored in the entries of the array y(1:28).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$
! y(7:10)  stores $\Lambda^{A}_{j}$
! y(11:12) stores $\mathrm{Re} \, \mathbf{r}_{1}$
! y(13:14) stores $\mathrm{Im} \, \mathbf{r}_{1}$
! y(15:16) stores $\mathrm{Re} \, \mathbf{r}_{2}$
! y(17:18) stores $\mathrm{Im} \, \mathbf{r}_{2}$
! y(19:20) stores $\mathrm{Re} \, \mathbf{r}_{1}^{\prime}$
! y(21:22) stores $\mathrm{Im} \, \mathbf{r}_{1}^{\prime}$
! y(23:24) stores $\mathrm{Re} \, \mathbf{r}_{2}^{\prime}$
! y(25:26) stores $\mathrm{Im} \, \mathbf{r}_{2}^{\prime}$
! y(27)    stores $\theta_{1}^{\prime}$
! y(28)    stores $\theta_{2}^{\prime}$

use multifield_globals
use multifield_utils

implicit none

real, parameter              :: N_start = 0.5, N_stop = 100.0                    ! Lower and upper bounds for mode injection
real, parameter              :: N_step = 10.0                                    ! E-fold interval between consecutive mode injections
real, parameter              :: dt_back = 20.0, dt_pert = 20.0                   ! Time steps
real, dimension(6)           :: y_back                                           ! Background state array
real, dimension(28)          :: y_pert                                           ! Background + perturbation state array
real, dimension(2)           :: phi, phidot                                      ! Field multiplet and its time derivative
real, dimension(2)           :: Re_r1, Im_r1, Re_r2, Im_r2                       ! Complex amplitude vector
real, dimension(2)           :: Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p               ! Time derivative of complex amplitude vector
real                         :: theta1_p, theta2_p                               ! Time derivative of real phases
integer, parameter           :: ngrid_max = int((N_stop - N_start) / N_step) + 1 ! E-fold interval between consecutive mode injections
real, dimension(2,ngrid_max) :: phi_grid, phidot_grid                            ! Initial conditions for mode-injection
real, dimension(ngrid_max)   :: H_grid, N_grid, PR_values, x_grid                ! Initial conditions for mode-injection
real, dimension(2,2)         :: phidotphidot, vbein_PT, W2_ij                    ! Variables
real                         :: H, H_end, Hdot, N, slowroll, k_mode, k_phys      ! Variables
real                         :: N_trigger, W_11, W_22                            ! Variables
integer                      :: i, ngrid                                         ! Indices
character(len=32)            :: arg                                              ! Command-line argument

! Parse argument
call get_command_argument(1, arg)
read (arg, *) energyscale

!========================================================================================================
! Background
!========================================================================================================
! Background initial conditions
select case (potential_shape)
	case ("elliptic")
		phi = [12.0, 12.0]
	case ("non-linear")
		phi = [20.0, 20.0]
end select

phidot = [0.0,  0.0]         ! $\dot{\phi}^{A}(t_{0})$
H      = Hubble(phi, phidot) ! $H(t_{0})$
N      = 0.0                 ! $N(t_{0})$

! Initialize background arrays
call pack_state_background(y_back, phi, phidot, H, N)

i         = 0       ! Counter
slowroll  = 0.0     ! Initialize slow-roll parameter $\epsilon(t_{0})$
N_trigger = N_start ! Initialize mode injection trigger

do while (slowroll <= 1.0)
	! Update functions of time
	call unpack_state_background(y_back, phi, phidot, H, N)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	if (N >= N_trigger .and. N <= N_stop) then
		i                = i + 1
		phi_grid(:,i)    = phi
		phidot_grid(:,i) = phidot
		H_grid(i)        = H
		N_grid(i)        = N
		N_trigger        = N_trigger + N_step
	end if
	
	call gl8_background(y_back, dt_back)
end do

H_end = H ! Hubble parameter at $\epsilon = 1$
ngrid = i ! Number of grid points

!========================================================================================================
! Perturbations
!========================================================================================================
k_phys = 1.0d5 * H_end               ! Physical wave vector magnitude

!$OMP PARALLEL DO PRIVATE(i, phi, phidot, H, N, slowroll, vbein_PT, k_mode, W2_ij, W_11, W_22, Re_r1, Re_r2, Im_r1, Im_r2, Re_r1_p, Re_r2_p, Im_r1_p, Im_r2_p, theta1_p, theta2_p, y_pert, N_trigger) SHARED(phi_grid, phidot_grid, H_grid, N_grid, k_phys, x_grid, PR_values) SCHEDULE(dynamic)
do i = 1, ngrid
	! Background initial conditions
	phi      = phi_grid(:,i)                   ! $\phi^{A}(t_{0})$
	phidot   = phidot_grid(:,i)                ! $\dot{\phi}^{A}(t_{0})$
	H        = H_grid(i)                       ! $H(t_{0})$
	N        = N_grid(i)                       ! $N(t_{0})$
	slowroll = - Hubbledot(phi, phidot) / H**2 ! $\epsilon(t_{0})$
	vbein_PT = vielbein_ad_is(phi, phidot)     ! $\Lambda^{A}_{i}(t_{0})$
	
	! Comoving wave vector magnitude
	k_mode = exp(N) * k_phys
	
	! Initial conditions for perturbations
	W2_ij    = freq_matrix(phi, phidot, H, N, vbein_PT, k_mode) ! $\Omega^{2}_{ij}(t_{0})$
	W_11     = sqrt(W2_ij(1,1))                                 ! $\sqrt{\Omega^{2}_{11}}$
	W_22     = sqrt(W2_ij(2,2))                                 ! $\sqrt{\Omega^{2}_{22}}$
	Re_r1    = [sqrt(2.0 * W_11)**(-1) , 0.0]
	Re_r2    = [0.0, sqrt(2.0 * W_22)**(-1)]
	Im_r1    = 0.0
	Im_r2    = 0.0
	Re_r1_p  = 0.0
	Re_r2_p  = 0.0
	Im_r1_p  = 0.0
	Im_r2_p  = 0.0
	theta1_p = W_11
	theta2_p = W_22

	! Initialize perturbation functions
	call pack_state_perturbations(y_pert, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, theta1_p, theta2_p)
	
	! Data writing trigger
	N_trigger = N

	do while (slowroll <= 1.0)
		! Update functions of time
		call unpack_state_perturbations(y_pert, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, theta1_p, theta2_p)
		
		! Update slow-roll parameter
		slowroll = - Hubbledot(phi, phidot) / H**2
		
		call gl8_perturbations(y_pert, dt_pert, k_mode)
	end do
	
	x_grid(i)    = log(k_mode/(exp(N)*H))
	PR_values(i) = powerspectrum(phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, k_mode)
end do
!$OMP END PARALLEL DO

do i = 1, ngrid
	write (*,'(7(6e25.10e3))') x_grid(i), PR_values(i)
end do

end program ps_full
