module chemistry_module

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, parameter :: nelements = 1  ! number of elements
  integer, parameter :: nspecies = 1   ! number of species
  integer, parameter :: nreactions = 0 ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), save :: elem_names(nelements)

  integer, parameter :: L_spec_name = 8 ! Each species name has at most 8 characters
  character*(L_spec_name), save :: spec_names(nspecies)

  real(amrex_real), save :: molecular_weight(nspecies), inv_mwt(nspecies)

  real(amrex_real), save :: Ru, Ruc, Patm

contains

end module chemistry_module
