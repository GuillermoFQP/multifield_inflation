module multifield_globals

implicit none

real, parameter             :: pi = 4.0 * atan(1.0)                     ! Constant
character(len=*), parameter :: field_space_geometry = "alpha-attractor" ! Choose between "Euclidean", "Renaux-Petel" or "alpha-attractor"
character(len=*), parameter :: potential_shape = "hybrid"               ! Choose between "elliptic", "non-linear" or "hybrid"
real, parameter             :: massfield1 = 1.4d-6                      ! Elliptic potential parameters
real, parameter             :: massfield2 = 7.0 * massfield1            ! Elliptic potential parameters
real, parameter             :: lambda_cons = 1.d-14                     ! Non-linear potential parameters
real, parameter             :: g_cons = 2.d-14                          ! Non-linear potential parameters
real, parameter             :: mu_cons = 0.0                            ! Non-linear potential parameters
real, parameter             :: chi0_hp = 2.5                            ! Hybrid inflation potential parameters
real, parameter             :: m2_hp = 1.47d-5                          ! Hybrid inflation potential parameters
real, parameter             :: m1_hp = 0.3 * m2_hp                      ! Hybrid inflation potential parameters
real, parameter             :: g_hp = 0.8 * m2_hp                       ! Hybrid inflation potential parameters
real, parameter             :: mu_hp = 0.0                              ! Hybrid inflation potential parameters
!real, parameter             :: mu_hp = -(5.0d-6 * m2_hp**2)**(1.0/3.0) ! Hybrid inflation potential parameters
real, parameter             :: alpha_aa = 1.0                           ! $\alpha$-attractor metric parameters
real, parameter             :: beta_aa = 1.0d10                         ! $\alpha$-attractor metric parameters
real                        :: energyscale                              ! Renaux-Petel metric parameter (set at runtime)

end module multifield_globals
