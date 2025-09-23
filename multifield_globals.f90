module multifield_globals

implicit none

real, parameter             :: pi = 4.0 * atan(1.0)         ! Constant
character(len=*), parameter :: field_space_geometry = "OWM" ! Choose between "EUM" for Euclidean, "RPM" for Renaux-Petel, or "AAM" for alpha-attractor
character(len=*), parameter :: potential_shape = "NLP"      ! Choose between "ELP" for elliptic, "NLP" for non-linear, "HYP" for hybrid or "MVP" for modulated valley potential
real, parameter             :: m1_ELP = 1.4d-6              ! Elliptic potential parameters
real, parameter             :: m2_ELP = 7.0 * m1_ELP        ! Elliptic potential parameters
real, parameter             :: lambda_NLP = 1.d-14          ! Non-linear potential parameters
real, parameter             :: g_NLP = 2.d-14               ! Non-linear potential parameters
real, parameter             :: mu_NLP = 0.0                 ! Non-linear potential parameters
real, parameter             :: chi0_HYP = 2.5               ! Hybrid inflation potential parameters
real, parameter             :: m2_HYP = 1.47d-5             ! Hybrid inflation potential parameters
real, parameter             :: m1_HYP = 0.3 * m2_HYP        ! Hybrid inflation potential parameters
real, parameter             :: g_HYP = 0.8 * m2_HYP         ! Hybrid inflation potential parameters
real, parameter             :: mu_HYP = 0.0                 ! Hybrid inflation potential parameters
!real, parameter             :: mu_HYP = -(5.0d-6 * m2_HYP**2)**(1.0/3.0) ! Hybrid inflation potential parameters
real, parameter             :: m_mvp = 1.4d-6               ! Modulated valley potential parameters
real, parameter             :: lambda1_mvp = 2.0d-3         ! Modulated valley potential parameters
real, parameter             :: f1_mvp = 1.0                 ! Modulated valley potential parameters
real, parameter             :: f2_mvp = 5.0                 ! Modulated valley potential parameters
real, parameter             :: f3_mvp = 2.0                 ! Modulated valley potential parameters
real, parameter             :: alpha_aa = 1.0               ! $\alpha$-attractor metric parameters
real, parameter             :: beta_aa = 1.0d10             ! $\alpha$-attractor metric parameters
real, parameter             :: amp_OWM = 0.75               ! Oscillatory warp metric parameters
real, parameter             :: frq_OWM = 10.0               ! Oscillatory warp metric parameters
real                        :: energyscale                  ! Renaux-Petel metric parameter (set at runtime)

end module multifield_globals
