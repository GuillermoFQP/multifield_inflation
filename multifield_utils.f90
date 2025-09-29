module multifield_utils

use multifield_globals

implicit none

contains

! Get background state array from background functions
subroutine pack_state_background(y, phi, phidot, H, N)
	real, intent(out) :: y(6)
	real, intent(in)  :: phi(2), phidot(2), H, N

	y(1:2) = phi    ! $\varphi^{A}$
	y(3:4) = phidot ! $\dot{\varphi}^{A}$
	y(5)   = H      ! $H(t)$
	y(6)   = N      ! $N(t)$
	
end subroutine pack_state_background

! Get background + perturbation state array from background + perturbation functions
subroutine pack_state_perturbations(y, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, theta1_p, theta2_p)
	real, intent(out) :: y(28)
	real, intent(in)  :: phi(2), phidot(2), H, N
	real, intent(in)  :: vbein_PT(2,2), Re_r1(2), Im_r1(2), Re_r2(2), Im_r2(2), Re_r1_p(2), Im_r1_p(2), Re_r2_p(2), Im_r2_p(2), theta1_p, theta2_p

	y(1:2)   = phi                    ! $\varphi^{A}$
	y(3:4)   = phidot                 ! $\dot{\varphi}^{A}$
	y(5)     = H                      ! $H(t)$
	y(6)     = N                      ! $N(t)$
	y(7:10)  = reshape(vbein_PT, [4]) ! $\Lambda^{A}_{i}$
	y(11:12) = Re_r1                  ! $\mathrm{Re} \, \mathbf{r}_{1}$
	y(13:14) = Im_r1                  ! $\mathrm{Im} \, \mathbf{r}_{1}$
	y(15:16) = Re_r2                  ! $\mathrm{Re} \, \mathbf{r}_{2}$
	y(17:18) = Im_r2                  ! $\mathrm{Im} \, \mathbf{r}_{2}$
	y(19:20) = Re_r1_p                ! $\mathrm{Re} \, \mathbf{r}_{1}^{\prime}$
	y(21:22) = Im_r1_p                ! $\mathrm{Im} \, \mathbf{r}_{1}^{\prime}$
	y(23:24) = Re_r2_p                ! $\mathrm{Re} \, \mathbf{r}_{2}^{\prime}$
	y(25:26) = Im_r2_p                ! $\mathrm{Im} \, \mathbf{r}_{2}^{\prime}$
	y(27)    = theta1_p               ! $\theta_{1}^{\prime}$
	y(28)    = theta2_p               ! $\theta_{2}^{\prime}$
	
end subroutine pack_state_perturbations

! Get background functions from background state array
subroutine unpack_state_background(y, phi, phidot, H, N)
	real, intent(in)  :: y(6)
	real, intent(out) :: phi(2), phidot(2), H, N

	phi    = y(1:2) ! $\varphi^{A}$
	phidot = y(3:4) ! $\dot{\varphi}^{A}$
	H      = y(5)   ! $H(t)$
	N      = y(6)   ! $N(t)$
	
end subroutine unpack_state_background

! Get background + perturbation functions from background + perturbation state array
subroutine unpack_state_perturbations(y, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, theta1_p, theta2_p)
	real, intent(in)  :: y(28)
	real, intent(out) :: phi(2), phidot(2), H, N
	real, intent(out) :: vbein_PT(2,2), Re_r1(2), Im_r1(2), Re_r2(2), Im_r2(2), Re_r1_p(2), Im_r1_p(2), Re_r2_p(2), Im_r2_p(2), theta1_p, theta2_p
	
	phi      = y(1:2)                  ! $\varphi^{A}$
	phidot   = y(3:4)                  ! $\dot{\varphi}^{A}$
	H        = y(5)                    ! $H(t)$
	N        = y(6)                    ! $N(t)$
	vbein_PT = reshape(y(7:10), [2,2]) ! $\Lambda^{A}_{i}$
	Re_r1    = y(11:12)                ! $\mathrm{Re} \, \mathbf{r}_{1}$
	Im_r1    = y(13:14)                ! $\mathrm{Im} \, \mathbf{r}_{1}$
	Re_r2    = y(15:16)                ! $\mathrm{Re} \, \mathbf{r}_{2}$
	Im_r2    = y(17:18)                ! $\mathrm{Im} \, \mathbf{r}_{2}$
	Re_r1_p  = y(19:20)                ! $\mathrm{Re} \, \mathbf{r}_{1}^{\prime}$
	Im_r1_p  = y(21:22)                ! $\mathrm{Im} \, \mathbf{r}_{1}^{\prime}$
	Re_r2_p  = y(23:24)                ! $\mathrm{Re} \, \mathbf{r}_{2}^{\prime}$
	Im_r2_p  = y(25:26)                ! $\mathrm{Im} \, \mathbf{r}_{2}^{\prime}$
	theta1_p = y(27)                   ! $\theta_{1}^{\prime}$
	theta2_p = y(28)                   ! $\theta_{2}^{\prime}$
	
end subroutine unpack_state_perturbations

subroutine evalf_background(y, dydx)
	real, dimension(6)  :: y, dydx
	real, dimension(2)   :: phi, phidot, phidotdot
	real                 :: H, N
	
	! Update functions
	call unpack_state_background(y, phi, phidot, H, N)
	
	! Update arrays
	dydx(1:2) = phidot                      ! $\dot{\varphi}^{A}$
	dydx(3:4) = phi_dot_dot(phi, phidot, H) ! $\ddot{\varphi}^{A}$
	dydx(5)   = Hubbledot(phi, phidot)      ! $\dot{H}$
	dydx(6)   = H                           ! $\dot{N}$
	
end subroutine evalf_background

subroutine evalf_perturbations(y, dydx, k_mode)
	real, dimension(28)  :: y, dydx
	real, dimension(2)   :: phi, phidot, phidotdot, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, Re_r1_pp, Im_r1_pp, Re_r2_pp, Im_r2_pp
	real, dimension(2,2) :: vbein_PT, vbein_PT_dot, Gammaphidot, W2_ij, Id, W2_eff_r1, W2_eff_r2
	real                 :: k_mode, H, N, theta1_p, theta2_p, theta1_pp, theta2_pp, r1_mod2, r2_mod2, r1_mod2_p, r2_mod2_p
	
	! Update functions
	call unpack_state_perturbations(y, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p, theta1_p, theta2_p)
	
	! Calculate time-evolution functions
	phidotdot    = phi_dot_dot(phi, phidot, H)                                               ! $\ddot{\phi}^{A}$
	W2_ij        = freq_matrix(phi, phidot, H, N, vbein_PT, k_mode)                          ! $\Omega^{2}_{ij}$
	Gammaphidot  = Gamma_phidot(phi, phidot)                                                 ! $\Gamma^{A}_{BC} \dot{\varphi}^{C}$
	Id           = reshape([1.0, 0.0, 0.0, 1.0], [2,2])                                      ! $\mathbf{I}$
	vbein_PT_dot = - matmul(Gammaphidot, vbein_PT)                                           ! $\dot{\Lambda}^{A}_{i}$
	r1_mod2      = dot_product(Re_r1, Re_r1) + dot_product(Im_r1, Im_r1)                     ! $r_{1}^{2}$
	r2_mod2      = dot_product(Re_r2, Re_r2) + dot_product(Im_r2, Im_r2)                     ! $r_{2}^{2}$
	r1_mod2_p    = 2.0 * (dot_product(Re_r1_p, Re_r1) + dot_product(Im_r1_p, Im_r1))         ! $(r_{1}^{2})^{\prime}$
	r2_mod2_p    = 2.0 * (dot_product(Re_r2_p, Re_r2) + dot_product(Im_r2_p, Im_r2))         ! $(r_{2}^{2})^{\prime}$
	W2_eff_r1    = W2_ij - theta1_p**2 * Id                                                  ! $\mathbf{\Omega}^{2} - \theta_{1}^{\prime 2} \mathbf{I}$
	W2_eff_r2    = W2_ij - theta2_p**2 * Id                                                  ! $\mathbf{\Omega}^{2} - \theta_{2}^{\prime 2} \mathbf{I}$
	theta1_pp    = - (r1_mod2_p / r1_mod2) * theta1_p                                        ! $\theta_{1}^{\prime\prime}$
	theta2_pp    = - (r2_mod2_p / r2_mod2) * theta2_p                                        ! $\theta_{2}^{\prime\prime}$
	Re_r1_pp     = - matmul(W2_eff_r1, Re_r1) + theta1_pp * Im_r1 + 2.0 * theta1_p * Im_r1_p ! $\mathrm{Re} \, \mathbf{r}_{1}^{\prime\prime}$
	Re_r2_pp     = - matmul(W2_eff_r2, Re_r2) + theta2_pp * Im_r2 + 2.0 * theta2_p * Im_r2_p ! $\mathrm{Re} \, \mathbf{r}_{2}^{\prime\prime}$
	Im_r1_pp     = - matmul(W2_eff_r1, Im_r1) - theta1_pp * Re_r1 - 2.0 * theta1_p * Re_r1_p ! $\mathrm{Im} \, \mathbf{r}_{1}^{\prime\prime}$
	Im_r2_pp     = - matmul(W2_eff_r2, Im_r2) - theta2_pp * Re_r2 - 2.0 * theta2_p * Re_r2_p ! $\mathrm{Im} \, \mathbf{r}_{2}^{\prime\prime}$
	
	! Update arrays
	dydx(1:2)   = phidot                     ! $\dot{\varphi}^{A}$
	dydx(3:4)   = phidotdot                  ! $\ddot{\varphi}^{A}$
	dydx(5)     = Hubbledot(phi, phidot)     ! $\dot{H}$
    dydx(6)     = H                          ! $\dot{N}$
	dydx(7:10)  = reshape(vbein_PT_dot, [4]) ! $\dot{\Lambda}^{A}_{i}$
	dydx(11:12) = exp(-N) * Re_r1_p
	dydx(13:14) = exp(-N) * Im_r1_p
	dydx(15:16) = exp(-N) * Re_r2_p
	dydx(17:18) = exp(-N) * Im_r2_p
	dydx(19:20) = exp(-N) * Re_r1_pp
	dydx(21:22) = exp(-N) * Im_r1_pp
	dydx(23:24) = exp(-N) * Re_r2_pp
	dydx(25:26) = exp(-N) * Im_r2_pp
	dydx(27)    = exp(-N) * theta1_pp
	dydx(28)    = exp(-N) * theta2_pp
	
end subroutine evalf_perturbations

! 8th order implicit Gauss-Legendre integrator
subroutine gl8_background(y, dt)
	integer, parameter  :: s = 4, n = 6
	real, intent(in)    :: dt
	real, intent(inout) :: y(n)
	integer             :: i, j
	real                :: g(n,s)

	! Butcher tableau for 8th order Gauss-Legendre method
	real, parameter :: a(s,s) = reshape((/ &
	0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
	0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
	0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
	-0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
	0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
	0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
	0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
	0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
	real, parameter ::   b(s) = (/ &
	0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
	0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

	! Iterate trial steps
	g = 0.0
	do j = 1, 16
		g = matmul(g,a)
		do i = 1, s
			call evalf_background(y + g(:,i) * dt, g(:,i))
		end do
	end do

	! Update the solution
	y = y + matmul(g,b)*dt
	
end subroutine gl8_background

! 8th order implicit Gauss-Legendre integrator
subroutine gl8_perturbations(y, dt, k_mode)
	integer, parameter  :: s = 4, n = 28
	real, intent(in)    :: dt, k_mode
	real, intent(inout) :: y(n)
	integer             :: i, j
	real                :: g(n,s)

	! Butcher tableau for 8th order Gauss-Legendre method
	real, parameter :: a(s,s) = reshape((/ &
	0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
	0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
	0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
	-0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
	0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
	0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
	0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
	0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
	real, parameter ::   b(s) = (/ &
	0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
	0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

	! Iterate trial steps
	g = 0.0
	do j = 1, 16
		g = matmul(g,a)
		do i = 1, s
			call evalf_perturbations(y + g(:,i) * dt, g(:,i), k_mode)
		end do
	end do

	! Update the solution
	y = y + matmul(g,b)*dt
	
end subroutine gl8_perturbations

! Determinant of a square matrix or order 2
pure function det(matrix) result(det_matrix)
	real, intent(in) :: matrix(2,2)
	real             :: det_matrix
	
	det_matrix = matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1)

end function det

! Field space metric $h_{AB}$
pure function metric(phi) result(h_AB)
	real, intent(in) :: phi(2)
	real             :: h_AB(2,2)
	integer          :: i
	
	select case (field_space_geometry)
		case ("EUM")
			h_AB = reshape([1.0, 0.0, 0.0, 1.0], [2,2])
		case ("RPM")
			h_AB(1,1) = 1.0 + 2.0 * phi(2)**2 / energyscale**2
			h_AB(1,2) = 0.0
			h_AB(2,1) = 0.0
			h_AB(2,2) = 1.0
		case ("AAM")
			h_AB(1,1) = 1.0 / (1.0 - phi(1)**2 / (6.0 * alpha_aa))**2
			h_AB(1,2) = 0.0
			h_AB(2,1) = 0.0
			h_AB(2,2) = 1.0 / (1.0 - phi(2)**2 / (6.0 * beta_aa))**2
		case ("OWM")
			h_AB(1,1) = 1.0
			h_AB(1,2) = 0.0
			h_AB(2,1) = 0.0
			h_AB(2,2) = amp_OWM * cos(frq_OWM * phi(1)) + 1.0
		case ("GBM")
			h_AB(1,1) = 1.0
			h_AB(1,2) = 0.0
			h_AB(2,1) = 0.0
			h_AB(2,2) = 1.0
			do i = 1, n_GBM
				if (abs(phi(1) - x_GBM(i)) <= spp_GBM * std_GBM(i)) then
					h_AB(2,2) = 1.0 + amp_GBM(i) * exp(-0.5 * (phi(1) - x_GBM(i))**2 / std_GBM(i)**2)
				end if
			end do
	end select

end function metric

! Field space inverse metric $h^{AB}$
pure function inv_metric(phi) result(hAB)
	real, intent(in) :: phi(2)
	real             :: hAB(2,2)
	integer          :: i
	
	select case (field_space_geometry)
		case ("EUM")
			hAB = reshape([1.0, 0.0, 0.0, 1.0], [2,2])
		case ("RPM")
			hAB(1,1) = 1.0 * energyscale**2 / (1.0 * energyscale**2 + 2.0 * phi(2)**2)
			hAB(1,2) = 0.0
			hAB(2,1) = 0.0
			hAB(2,2) = 1.0
		case ("AAM")
			hAB(1,1) = (1.0 - phi(1)**2 / (6.0 * alpha_aa))**2
			hAB(1,2) = 0.0
			hAB(2,1) = 0.0
			hAB(2,2) = (1.0 - phi(2)**2 / (6.0 * beta_aa))**2
		case ("OWM")
			hAB(1,1) = 1.0
			hAB(1,2) = 0.0
			hAB(2,1) = 0.0
			hAB(2,2) = 1.0 / (amp_OWM * cos(frq_OWM * phi(1)) + 1.0)
		case ("GBM")
			hAB(1,1) = 1.0
			hAB(1,2) = 0.0
			hAB(2,1) = 0.0
			hAB(2,2) = 1.0
			do i = 1, n_GBM
				if (abs(phi(1) - x_GBM(i)) <= spp_GBM * std_GBM(i)) then
					hAB(2,2) = 1.0 / (amp_GBM(i) * exp(-0.5 * (phi(1) - x_GBM(i))**2 / std_GBM(i)**2) + 1.0)
				end if
			end do
	end select

end function inv_metric

! Christoffel symbols $\Gamma^{A}_{BC}$
pure function Christoffel(phi) result(Gamma)
	real, intent(in) :: phi(2)
	real             :: Gamma(2,2,2)
	integer          :: i
	real             :: dx, sd2, dx2, gbmp ! GBM metric variables
	
	Gamma = 0.0
	
	select case (field_space_geometry)
		case ("EUM")
			! Nothing to do
		case ("RPM")
			Gamma(1,1,2) = 2.0 * phi(2) / (energyscale**2 + 2.0 * phi(2)**2)
			Gamma(1,2,1) = Gamma(1,1,2)
			Gamma(2,1,1) = -2.0 * phi(2) / energyscale**2
		case ("AAM")
			Gamma(1,1,1) = phi(1) / (3.0 * alpha_aa * (1.0 - phi(1)**2 / (6.0 * alpha_aa)))
			Gamma(2,2,2) = phi(2) / (3.0 *  beta_aa * (1.0 - phi(2)**2 / (6.0 *  beta_aa)))
		case ("OWM")
			Gamma(1,2,2) = 0.5 * amp_OWM * frq_OWM * sin(frq_OWM * phi(1))
			Gamma(2,1,2) = -0.5 * amp_OWM * frq_OWM * sin(frq_OWM * phi(1)) / (amp_OWM * cos(frq_OWM * phi(1)) + 1.0)
			Gamma(2,2,1) = Gamma(2,1,2)
		case ("GBM")
			do i = 1, n_GBM
				dx = phi(1) - x_GBM(i)
				if (abs(dx) <= spp_GBM * std_GBM(i)) then
					sd2  = std_GBM(i)**2
					dx2  = dx**2
					gbmp = exp(-0.5 * dx2 / sd2)
					Gamma(1,2,2) =  0.5 * amp_GBM(i) * dx * gbmp / sd2
					Gamma(2,1,2) = -0.5 * amp_GBM(i) * dx * gbmp / (sd2 * (amp_GBM(i) * gbmp + 1.0))
					Gamma(2,2,1) = Gamma(2,1,2)
				end if
			end do
	end select
	
end function Christoffel

! Riemann curvature tensor $R_{ABCD}$
pure function Riemann(phi) result(R)
	real, intent(in) :: phi(2)
	real             :: R(2,2,2,2)
	integer          :: i
	real             :: dx, sd2, dx2, gbmp ! GBM metric variables
	
	R = 0.0
	
	select case (field_space_geometry)
		case ("EUM")
			! Nothing to do
		case ("RPM")
			R(1,2,1,2) = -0.5 / (0.25 * energyscale**2 + 0.5 * phi(2)**2)
		case ("AAM")
			! Nothing to do
		case ("OWM")
			R(1,2,1,2) = amp_OWM * frq_OWM**2 * (-0.25 * amp_OWM * sin(frq_OWM * phi(1))**2 + 0.5 * amp_OWM + 0.5 * cos(frq_OWM * phi(1))) / (amp_OWM * cos(frq_OWM * phi(1)) + 1.0)
		case ("GBM")
			do i = 1, n_GBM
				dx = phi(1) - x_GBM(i)
				if (abs(dx) <= spp_GBM * std_GBM(i)) then
					sd2  = std_GBM(i)**2
					dx2  = dx**2
					gbmp = exp(-0.5 * dx2 / sd2)
					R(1,2,1,2) = amp_GBM(i) * (0.25 * amp_GBM(i) * dx2 * gbmp**2 + 0.5 * sd2 * (amp_GBM(i) * gbmp + 1.0) * gbmp - 0.5 * dx2 * (amp_GBM(i) * gbmp + 1.0) * gbmp) / (sd2**2 * (amp_GBM(i) * gbmp + 1.0))
				end if
			end do
	end select
	
	! Symmetries
	R(1,2,2,1) = -R(1,2,1,2)
	R(2,1,1,2) = -R(1,2,1,2)
	R(2,1,2,1) =  R(1,2,1,2)

end function Riemann

! Field potential
pure function potential(phi) result(V)
	real, intent(in) :: phi(2)
	real             :: V
	real             :: ds2, ds, cs, sn, ang
	
	select case (potential_shape)
		case ("ELP")
			V = 0.5 * m1_ELP**2 * phi(1)**2 + 0.5 * m2_ELP**2 * phi(2)**2
		case ("NLP")
			V = 0.5 * g_NLP * phi(1)**2 * phi(2)**2 + 0.25 * lambda_NLP * phi(1)**4 + 0.25 * mu_NLP * phi(2)**4
		case ("HYP")
			V = 0.25 * m2_HYP**2 * (phi(2)**2 - chi0_HYP**2)**2 / chi0_HYP**2 + 0.5 * m1_HYP**2 * phi(1)**2 + 0.5 * g_HYP**2 * phi(1)**2 * phi(2)**2 + mu_HYP**3 * phi(2)
		case ("MVP")
			! Precompute factors
			ds2 = phi(1)**2 + phi(2)**2
			ds  = sqrt(ds2)
			cs  = cos(f1_mvp * ds)
			sn  = sin(f1_mvp * ds)
			ang = atan2(phi(2), phi(1))
			
			! Final expression
			V = lambda1_mvp**4 * cos(f2_mvp * ang - f3_mvp * sn) + 0.5 * m_mvp**2 * ds2
	end select

end function potential

! Gradient of  the potential in the field space $\mathcal{D}_{A} V$
pure function NablaV(phi) result(DV)
	real, intent(in) :: phi(2)
	real             :: DV(2)
	real             :: ds2, ds, cs, sn, ang
	
	select case (potential_shape)
		case ("ELP")
			DV(1) = m1_ELP**2 * phi(1)
			DV(2) = m2_ELP**2 * phi(2)
		case ("NLP")
			DV(1) = g_NLP * phi(1) * phi(2)**2 + lambda_NLP * phi(1)**3
			DV(2) = g_NLP * phi(1)**2 * phi(2) + mu_NLP * phi(2)**3
		case ("HYP")
			DV(1) = g_HYP**2 * phi(1) * phi(2)**2 + m1_HYP**2 * phi(1)
			DV(2) = m2_HYP**2 * phi(2) * (-chi0_HYP**2 + phi(2)**2) / chi0_HYP**2 + g_HYP**2 * phi(1)**2 * phi(2) + mu_HYP**3
		case ("MVP")
			! Precompute factors
			ds2 = phi(1)**2 + phi(2)**2
			ds  = sqrt(ds2)
			cs  = cos(f1_mvp * ds)
			sn  = sin(f1_mvp * ds)
			ang = atan2(phi(2), phi(1))
			
			! Final expression
			DV(1) = -lambda1_mvp**4 * (-f1_mvp * f3_mvp * phi(1) * cs / ds - f2_mvp * phi(2) / ds2) * sin(f2_mvp * ang - f3_mvp * sn) + m_mvp**2 * phi(1)
			DV(2) = -lambda1_mvp**4 * (-f1_mvp * f3_mvp * phi(2) * cs / ds + f2_mvp * phi(1) / ds2) * sin(f2_mvp * ang - f3_mvp * sn) + m_mvp**2 * phi(2)
	end select
	
end function NablaV

! Hessian of the potential in the field space $\partial_{B} \partial_{A} V$
pure function HessianV(phi) result(D2V)
	real, intent(in) :: phi(2)
	real             :: D2V(2,2)
	real             :: ds2, ds, ds3, cs, sn, css, sns, ang
	
	select case (potential_shape)
		case ("ELP")
			D2V(1,1) = m1_ELP**2
			D2V(1,2) = 0.0
			D2V(2,2) = m2_ELP**2
		case ("NLP")
			D2V(1,1) = g_NLP * phi(2)**2 + 3.0 * lambda_NLP * phi(1)**2
			D2V(1,2) = 2.0 * g_NLP*phi(1) * phi(2)
			D2V(2,2) = g_NLP * phi(1)**2 + 3.0 * mu_NLP * phi(2)**2
		case ("HYP")
			D2V(1,1) = g_HYP**2 * phi(2)**2 + m1_HYP**2
			D2V(1,2) = 2.0 * g_HYP**2 * phi(1) * phi(2)
			D2V(2,2) = m2_HYP**2 * (phi(2)**2 - chi0_HYP**2) / chi0_HYP**2 + 2.0 * m2_HYP**2 * phi(2)**2 / chi0_HYP**2 + g_HYP**2 * phi(1)**2
		case ("MVP")
			! Precompute factors
			ds2 = phi(1)**2 + phi(2)**2
			ds  = sqrt(ds2)
			ds3 = ds**3
			cs  = cos(f1_mvp * ds)
			sn  = sin(f1_mvp * ds)
			css = cos(f2_mvp * ang - f3_mvp * sn)
			sns = sin(f2_mvp * ang - f3_mvp * sn)
			ang = atan2(phi(2), phi(1))
			
			! Final expression
			D2V(1,1) = -lambda1_mvp**4 * (f1_mvp * f3_mvp * phi(1) * cs / ds + f2_mvp * phi(2) / ds2)**2 * css - lambda1_mvp**4 * (f1_mvp**2 * f3_mvp * phi(1)**2 * sn / ds2 + f1_mvp * f3_mvp * phi(1)**2 * cs / ds3 - f1_mvp * f3_mvp * cs / ds + 2.0 * f2_mvp * phi(1) * phi(2) / ds2**2) * sns + m_mvp**2
			D2V(1,2) = -lambda1_mvp**4 * ((f1_mvp * f3_mvp * phi(1) * cs / ds + f2_mvp * phi(2) / ds2) * (f1_mvp * f3_mvp * phi(2) * cs / ds - f2_mvp * phi(1) / ds2) * css + (f1_mvp**2 * f3_mvp * phi(1) * phi(2) * sn / ds2 + f1_mvp * f3_mvp * phi(1) * phi(2) * cs / ds3 + 2.0 * f2_mvp * phi(2)**2 / ds2**2 - f2_mvp / ds2) * sns)
			D2V(2,2) = -lambda1_mvp**4 * (f1_mvp * f3_mvp * phi(2) * cs / ds - f2_mvp * phi(1) / ds2)**2 * css - lambda1_mvp**4 * (f1_mvp**2 * f3_mvp * phi(2)**2 * sn / ds2 + f1_mvp * f3_mvp * phi(2)**2 * cs / ds3 - f1_mvp * f3_mvp * cs / ds - 2.0 * f2_mvp * phi(1) * phi(2) / ds2**2) * sns + m_mvp**2
	end select
	
	! Symmetry
	D2V(2,1) = D2V(1,2)
	
end function HessianV

! Second covariant derivative of the potential in the field space with the first index raised $\mathcal{D}_{B} \mathcal{D}_{A} V$
pure function NablaNablaV(phi) result(DDV)
	real, intent(in) :: phi(2)
	real             :: DDV(2,2), DV(2), D2V(2,2), Gamma(2,2,2), GammaDV(2,2)
	integer          :: i, j
	
	DV    = NablaV(phi)      ! $\partial_{A} V$
	D2V   = HessianV(phi)    ! $\partial_{B} \partial_{A} V$
	Gamma = Christoffel(phi) ! $\Gamma^{A}_{BC}$
	
	! $\Gamma^{C}_{AB} \partial_{C} V$
	do j = 1, 2
		do i = 1, 2
			GammaDV(i,j) = sum(Gamma(:,i,j) * DV(:))
		end do
	end do
	
	DDV = D2V - GammaDV
	
end function NablaNablaV

! Hubble parameter $H$
pure function Hubble(phi, phidot) result(H)
	real, intent(in) :: phi(2), phidot(2)
	real             :: H
	
	H = sqrt((sum(metric(phi) * outer_product(phidot, phidot)) / 2.0 + potential(phi)) / 3.0)

end function Hubble

! Time derivative of Hubble parameter $\dot{H}$
pure function Hubbledot(phi, phidot) result(Hdot)
	real, intent(in) :: phi(2), phidot(2)
	real             :: Hdot
	
	Hdot = - sum(metric(phi) * outer_product(phidot, phidot)) / 2.0

end function Hubbledot

! Second time derivative of field multiplet $\ddot{\phi}^{A}$
pure function phi_dot_dot(phi, phidot, H) result(phidotdot)
	real, intent(in)       :: phi(2), phidot(2), H
	real, dimension(2)     :: phidotdot, DV, Gammaphidotphidot
	real, dimension(2,2)   :: hAB, phidotphidot, Gammaphidot
	real, dimension(2,2,2) :: Gamma
	integer                :: i
	
	hAB               = inv_metric(phi)               ! $h^{AB}$
	DV                = NablaV(phi)                   ! $\mathcal{D}_{A} V$
	phidotphidot      = outer_product(phidot, phidot) ! $\dot{\varphi}^{A} \dot{\varphi}^{B}$
	Gammaphidot       = Gamma_phidot(phi, phidot)     ! $\Gamma^{A}_{BC} \dot{\varphi}^{C}$
	Gammaphidotphidot = matmul(Gammaphidot, phidot)   ! $\Gamma^{A}_{BC} \dot{\varphi}^{B} \dot{\varphi}^{C}$
	
	phidotdot = - Gammaphidotphidot - 3.0 * H * phidot - matmul(hAB, DV) ! $\ddot{\phi}^{A}$

end function phi_dot_dot

! Field squared-mass matrix $\mathcal{M}_{AB}$
pure function mass_matrix(phi, phidot, H) result (M)
	real, intent(in)         :: phi(2), phidot(2), H
	real                     :: C1, C2, Hdot
	real, dimension(2)       :: phidot_lower, Dtphidot, Dtphidot_lower, phidotdot, Gammaphidotphidot
	real, dimension(2,2)     :: h_AB, DDV, Rphidotphidot, term1, term2, M, phidotphidot, Gammaphidot
	real, dimension(2,2,2)   :: Gamma
	real, dimension(2,2,2,2) :: R
	integer                  :: i, j
	
	h_AB         = metric(phi)                   ! $h_{AB}$
	R            = Riemann(phi)                  ! $R_{ABCD}$
	phidotphidot = outer_product(phidot, phidot) ! $\dot{\phi}^{A} \dot{\phi}^{B}$
	phidotdot    = phi_dot_dot(phi, phidot, H)   ! $\ddot{\phi}^{A}$
	Gamma        = Christoffel(phi)              ! $\Gamma^{A}_{BC}$
	Hdot         = Hubbledot(phi, phidot)        ! $\dot{H}$
	DDV          = NablaNablaV(phi)              ! $\mathcal{D}_{B} \mathcal{D}_{A} V$

	! $R_{ACDB} \dot{\varphi}^{C} \dot{\varphi}^{D}$
	do j = 1, 2
		do i = 1, 2
			Rphidotphidot(i,j) = sum(R(i,:,:,j) * phidotphidot)
		end do
	end do
	
	Gammaphidot       = Gamma_phidot(phi, phidot)     ! $\Gamma^{A}_{BC} \dot{\varphi}^{C}$
	Gammaphidotphidot = matmul(Gammaphidot, phidot)   ! $\Gamma^{A}_{BC} \dot{\varphi}^{B} \dot{\varphi}^{C}$
	Dtphidot          = phidotdot + Gammaphidotphidot ! $\mathcal{D}_{t} \varphi^{A}$
	phidot_lower      = matmul(h_AB, phidot)          ! $\dot{\varphi}_{A}$
	Dtphidot_lower    = matmul(h_AB, Dtphidot)        ! $\mathcal{D}_{t} \varphi_{A}$
	
	C1 = 3.0 - Hdot / H**2 ! $\frac{1}{a^{3}} \mathcal{D}_{t} \left( \frac{a^{3}}{H} \right)$
	C2 = H**(-1)           ! $\frac{1}{a^{3}} \left( \frac{a^{3}}{H} \right)$
	
	term1 = outer_product(phidot_lower, phidot_lower)        ! $\dot{\varphi}_{A} \dot{\varphi}_{B}$
	term2 = Sym(outer_product(Dtphidot_lower, phidot_lower)) ! $\mathcal{D}_{t} \dot{\varphi}_{A} \dot{\varphi}_{B} + \dot{\varphi}_{A} \mathcal{D}_{t} \dot{\varphi}_{B}$
	
	M = DDV - Rphidotphidot - (C1*term1 + C2*term2) ! $\mathcal{M}_{AB}$

end function mass_matrix

! Vielbein adiabatic-isocurvature $e^{A}_{i}$
pure function vielbein_ad_is(phi, phidot) result(vbein_AI)
	real, intent(in) :: phi(2), phidot(2)
	real             :: h_AB(2,2), vbein_AI(2,2), phidotphidot(2,2), epsilon(2,2), varepsilon(2,2), phidot_mag

	h_AB         = metric(phi)                           ! $h_{AB}$
	epsilon      = reshape([0.0, -1.0, 1.0, 0.0], [2,2]) ! Levi-Civita symbol $\epsilon^{AB}$
	varepsilon   = sqrt(det(h_AB))**(-1) * epsilon       ! Levi-Civita tensor $\varepsilon^{AB} = \epsilon^{AB} / \sqrt{h}$
	phidotphidot = outer_product(phidot, phidot)         ! $\dot{\varphi}^{A} \dot{\varphi}^{B}$
	phidot_mag   = sqrt(sum(h_AB * phidotphidot))        ! $|\dot{\varphi}| = \sqrt{h_{AB} \dot{\varphi}^{A} \dot{\varphi}^{B}}$
	
	! Adiabatic component
	vbein_AI(:,1) = phidot / phidot_mag                                   ! $e^{A}_{1} = \frac{\dot{\varphi}^{A}}{|\dot{\varphi}|}$
	! Isocurvature component
	vbein_AI(:,2) = matmul(matmul(varepsilon, h_AB), phidot) / phidot_mag ! $e^{A}_{2} = \frac{\varepsilon^{AB} h_{BC} \dot{\varphi}^{C}}{|\dot{\varphi}|}$

end function vielbein_ad_is

! Tensor product
pure function outer_product(v1, v2) result(product)
	real, intent(in) :: v1(2), v2(2)
	real             :: product(2,2)
	
	product(1,1) = v1(1) * v2(1)
	product(2,1) = v1(2) * v2(1)
	product(1,2) = v1(1) * v2(2)
	product(2,2) = v1(2) * v2(2)

end function outer_product

! $\Gamma^{A}_{BC} \dot{\varphi}^{B}$
pure function Gamma_phidot(phi, phidot) result(Gammaphidot)
	real, intent(in) :: phi(2), phidot(2)
	real             :: Gammaphidot(2,2), Gamma(2,2,2)
	integer          :: i, j
	
	Gamma = Christoffel(phi)

	! Christoffel coefficients contracted with the time derivative of the field multiplet $\Gamma^{A}_{BC} \dot{\varphi}^{C}$
	do j = 1, 2
		do i = 1, 2
			Gammaphidot(i,j) = sum(Gamma(i,j,:) * phidot)
		end do
	end do
	
end function Gamma_phidot

! Frequency-squared matrix $\Omega_{ij}$
pure function freq_matrix(phi, phidot, H, N, vbein_PT, k_mode) result(W2_ij)
	real, intent(in) :: phi(2), phidot(2), H, N, vbein_PT(2,2), k_mode
	real             :: Id(2,2), W2_ij(2,2), M_AB(2,2), M_ij(2,2), a2, Hdot
	
	Id   = reshape([1.0, 0.0, 0.0, 1.0], [2,2])
	a2   = exp(2.0 * N)
	Hdot = Hubbledot(phi, phidot)
	M_AB = mass_matrix(phi, phidot, H)
	M_ij = matmul(matmul(transpose(vbein_PT), M_AB), vbein_PT)
	
	! $\Omega^{2}_{ij} \equiv \left( k^{2} - \frac{a^{\prime\prime}}{a} \right) \! \delta_{ij} + a^{2} \mathcal{M}_{ij}$
	! where $\frac{a^{\prime\prime}}{a} = a^{2} (2 H^{2} + \dot{H})$
	W2_ij = (k_mode**2 - a2 * (2.0 * H**2 + Hdot)) * Id + a2 * M_ij
	
end function freq_matrix

! Calculate the power spectrum
pure function powerspectrum(phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, k_mode) result(PR)
	real, intent(in) :: phi(2), phidot(2), H, N, vbein_PT(2,2), Re_r1(2), Im_r1(2), Re_r2(2), Im_r2(2), k_mode
	real             :: phidotphidot(2,2), h_AB(2,2), phidot_mag2, C1, C2, HA(2,2), UU(2,2), HAUUAH(2,2), PR, a
	
	! Field-field correlator (modulo $\delta^{3}(\mathbf{k}+\mathbf{k}^{\prime})$)
	UU = outer_product(Re_r1, Re_r1) + outer_product(Re_r2, Re_r2) + outer_product(Im_r1, Im_r1) + outer_product(Im_r2, Im_r2)
	
	! $\langle \mathcal{R}(\mathbf{k},t) \mathcal{R}(\mathbf{k}^{\prime},t) \rangle$
	a            = exp(N)
	h_AB         = metric(phi)                           ! $h_{AB}$
	phidotphidot = outer_product(phidot, phidot)         ! $\dot{\varphi}^{A} \dot{\varphi}^{B}$
	phidot_mag2  = sum(h_AB * phidotphidot)              ! $|\dot{\varphi}|^{2} = h_{AB} \dot{\varphi}^{A} \dot{\varphi}^{B}$
	HA           = matmul(h_AB, vbein_PT)                ! $\mathbf{H \Lambda}$
	HAUUAH       = matmul(matmul(HA, UU), transpose(HA)) ! $(\mathbf{H \Lambda L}) (\mathbf{H \Lambda L})^{\intercal}$
	C1           = k_mode**3 / (2.0 * pi**2)             ! $\frac{k^{3}}{2 \pi^{2}}$
	C2           = H / (a * phidot_mag2)                 ! $\frac{H}{a |\dot{\varphi}|^{2}}$
	
	! $\mathcal{P}_{\mathcal{R}} (k) = \frac{k^{3}}{2\pi^{2}} \left( \frac{H}{a |\dot{\varphi}|} \right)^{2} \mathbf{e}_{\sigma}^{\intercal} [(\mathbf{H \Lambda L}) (\mathbf{H \Lambda L})^{\intercal}] \mathbf{e}_{\sigma}$
	PR = C1 * C2**2 * sum(HAUUAH * phidotphidot)

end function powerspectrum

! Symmetrizer operator for a 2x2 real matrix
pure function Sym(A) result(SymA)
	real, intent(in) :: A(2,2)
	real             :: SymA(2,2)
	
	SymA = A + transpose(A)

end function Sym

! Compute the orthonormal matrix $\mathbf{U}$ such that $\mathbf{U^{\intercal} M U}$ is diagonal, where $\mathbf{M}$ is symmetric
subroutine diagonalization(M, U, eigenvalues)
	real, intent(in)  :: M(2,2)
	real, intent(out) :: U(2,2), eigenvalues(2)
	real              :: discriminant, v1(2), v2(2)
	
	! Calculate discriminant
	discriminant = (M(1,1) - M(2,2))**2 + 4.0 * M(1,2) * M(2,1)

	! Check for complex eigenvalues
	 if (discriminant < 0.0) then
		 write(*,*) 'ERROR: Complex eigenvalues detected in mass matrix (discriminant = ', discriminant, ')'
		 error stop 'Complex eigenvalues not allowed in s - terminating program'
	 end if

	! Compute eigenvalues
	eigenvalues(1) = (M(1,1) + M(2,2) + sqrt(discriminant)) / 2.0
	eigenvalues(2) = (M(1,1) + M(2,2) - sqrt(discriminant)) / 2.0

	if (abs(M(1,2)) < epsilon(1.0)) then
		! Handle case where $\mathbf{M}$ is already diagonal
		U = reshape([1.0, 0.0, 0.0, 1.0], [2,2])
	else if (abs(M(1,1) - eigenvalues(1)) > epsilon(1.0)) then
		! Calculate first eigenvector: $(\mathbf{M} - \lambda_{1} \mathbf{I}) \mathbf{v}_{1} = 0$
		! We solve the first equation $(M_{11} - \lambda_{1}) x + M_{12} y = 0$
		! Choose $y = 1$, then $x = - M_{12} / (M_{11} - \lambda_{1})$
		v1(1) = - M(1,2) / (M(1,1) - eigenvalues(1))
		v1(2) = 1.0
	else
		! If $M_{11} - \lambda_{1} = 0$, then from $M_{21} x + (M_{22} - \lambda_{1})y = 0$
		! Since M_{12} \ne 0$ (checked above), set $y = 1$ then $x = - (M_{22} - \lambda_{1}) / M_{21}$
		v1(1) = - (M(2,2) - eigenvalues(1)) / M(2,1)
		v1(2) = 1.0
	end if

	! Normalize V1
	v1 = v1 / norm2(v1)

	! Second eigenvector is orthogonal to the first (since $M$ is symmetric)
	v2 = [-v1(2), v1(1)]

	! Construct $U$ matrix
	U = reshape([v1, v2], [2,2])

end subroutine diagonalization

end module multifield_utils
