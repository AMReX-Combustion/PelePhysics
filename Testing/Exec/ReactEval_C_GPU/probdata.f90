module probdata_module

  use amrex_fort_module, only : amrex_real
  implicit none

  real(amrex_real), save :: thermal_conductivity, diff_coeff
  real(amrex_real), save :: T1, T2, rho0
  real(amrex_real), save :: t_0

end module probdata_module
