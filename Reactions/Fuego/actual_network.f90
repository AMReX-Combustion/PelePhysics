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

  use chemistry_module, only : nspecies, chemistry_init, chemistry_close, chemistry_initialized, spec_names, elem_names, rwrk, iwrk

  implicit none

  integer :: nspec, nelem, nreac, nfit, naux
  character (len=16), save, allocatable :: aux_names(:)

contains
  
  subroutine actual_network_init

    if (.not. chemistry_initialized)  call chemistry_init()
    call ckindx(iwrk,rwrk,nelem,nspec,nreac,nfit)
    nspecies = nspec
    naux = 0
    allocate(aux_names(naux))

  end subroutine actual_network_init

  subroutine actual_network_close

    call chemistry_close()
    deallocate(aux_names)
    nspecies = 0
    naux = 0

  end subroutine actual_network_close

end module actual_network
