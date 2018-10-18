module chemistry_module

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, save :: nelements   ! number of elements
  integer, save :: nspecies    ! number of species
  integer, save :: nreactions  ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), allocatable, save :: elem_names(:)

  integer, parameter :: L_spec_name = 8 ! Each species name has at most 8 characters
  character*(L_spec_name), allocatable, save :: spec_names(:)

  real(amrex_real), allocatable, save :: molecular_weight(:), inv_mwt(:)

  real(amrex_real), save :: Ru, Ruc, Patm

contains

end module chemistry_module
