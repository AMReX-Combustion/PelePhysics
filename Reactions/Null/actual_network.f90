! An "null" network.  This provides the properties of a set of non-reacting species.
!
! nspec            -- the number of species
! naux             -- the number of auxiliary variables
!
! spec_names       -- the name of the species
!
! aux_names        -- the name of the auxiliary variable
!
! molec_wt         -- molecular weight of species

module actual_network

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, parameter :: nspec = 1
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec) 
  character (len=16), save :: aux_names(naux)

  real(amrex_real), save   :: molec_wt(nspec)

contains
  
  subroutine actual_network_init

    spec_names(1) = "X"

    molec_wt(1) = 1.0

  end subroutine actual_network_init

  subroutine actual_network_close

  end subroutine actual_network_close

end module actual_network
