module multifield_globals

implicit none

real, parameter             :: pi = 4.0 * atan(1.0)                                 ! Constant
character(len=*), parameter :: field_space_geometry = "Renaux-Petel"                ! Choose between "Euclidean" or "Renaux-Petel"
character(len=*), parameter :: potential_shape = "elliptic"                         ! Choose between "elliptic" or "non-linear"
real, parameter             :: massfield1 = 1.4d-6, massfield2 = 7.0 * massfield1   ! Elliptic potential parameters
real, parameter             :: lambda_cons = 1.d-14, g_cons = 2.d-14, mu_cons = 0.0 ! Non-linear potential parameters
real                        :: energyscale                                          ! Renaux-Petel metri parameter (set at runtime)

end module multifield_globals
