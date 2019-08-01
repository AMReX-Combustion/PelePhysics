module actual_network

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, parameter :: nelements   = 1 ! number of elements
  integer, parameter :: nspecies    = 1 ! number of species
  integer, parameter :: nreactions  = 0 ! number of reactions
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspecies) 
  character (len=16), save :: aux_names(naux)

  real(amrex_real), save   :: molec_wt(nspecies)

contains
  
  subroutine actual_network_init

    spec_names(1) = "X"

    molec_wt(1) = 1.0

  end subroutine actual_network_init

  subroutine actual_network_close

  end subroutine actual_network_close

end module actual_network
