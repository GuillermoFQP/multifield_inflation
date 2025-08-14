module multifield_globals

implicit none

real, parameter             :: pi = 4.0 * atan(1.0)                                 ! Constant
character(len=*), parameter :: field_space_geometry = "Euclidean"                   ! Choose between "Euclidean", "Renaux-Petel" or "alpha-attractor"
character(len=*), parameter :: potential_shape = "hybrid"                           ! Choose between "elliptic", "non-linear" or "hybrid"
real, parameter             :: massfield1 = 1.4d-6, massfield2 = 7.0 * massfield1   ! Elliptic potential parameters
real, parameter             :: lambda_cons = 1.d-14, g_cons = 2.d-14, mu_cons = 0.0 ! Non-linear potential parameters
real, parameter             :: lambda_hp = 1.0d-12, m1_hp = 1.95d-6                 ! Hybrid inflation potential parameters
real, parameter             :: M2_hp = 4.0d-3 * lambda_hp**0.25                     ! Hybrid inflation potential parameters
real, parameter             :: g_hp = 1.8d-3 * lambda_hp**0.25                      ! Hybrid inflation potential parameters
real, parameter             :: alpha_aa = 1.0, beta_aa = 1.0                        ! $\alpha$-attractor metric parameters
real                        :: energyscale                                          ! Renaux-Petel metric parameter (set at runtime)

end module multifield_globals
