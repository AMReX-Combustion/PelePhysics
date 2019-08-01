module react_type_module

  use network, only: nspecies

  use amrex_fort_module, only : amrex_real
  implicit none

  type :: reaction_stat_t
     real(amrex_real) :: cost_value
     logical          :: reactions_succesful = .false.
  end type reaction_stat_t

  type :: react_t
     real(amrex_real) :: rhoedot_ext
     real(amrex_real) :: rhohdot_ext
     real(amrex_real) :: rhoY(nspecies)
     real(amrex_real) :: rhoYdot_ext(nspecies)
     real(amrex_real) :: rho
     real(amrex_real) :: T
     real(amrex_real) :: e
     real(amrex_real) :: h
     real(amrex_real) :: p
     integer :: i, j, k
  end type react_t
  
  interface build
     module procedure react_build
  end interface build

  interface destroy
     module procedure react_destroy
  end interface destroy


contains

  subroutine react_build(r)
    type(react_t), intent(inout) :: r
  end subroutine react_build

  subroutine react_destroy(r)
    type(react_t), intent(inout) :: r
  end subroutine react_destroy
  
end module react_type_module

