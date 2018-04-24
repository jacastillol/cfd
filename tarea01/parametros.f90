module parametros
  implicit none
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: Long = selected_real_kind(33, 4931)
  integer :: dim
  real(Long), parameter :: pi = acos(-1.0_Long)
  real(Long) :: dx, k1, k2  
end module parametros

