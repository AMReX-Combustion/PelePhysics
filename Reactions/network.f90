! the network module provides the information about the species we are
! advecting: 
!
! nspecies      -- the number of species
!
! nreactons      -- the number of reactions
!
! spec_names -- the name of the chemical species
!
! This module contains the following routines:
!
!  network_init          -- initialize the chemical properties
!
!  network_species_index -- return the index of the species given its name

module network

  use actual_network
  
  implicit none

  logical :: network_initialized = .false.

contains
  
  subroutine network_init

    implicit none
    
    ! First, we call the specific network initialization.
    ! This should set number of aux variables (nspecies, nreactions should be parameters
    ! defined there), and the components of the species.
    
    call actual_network_init

    ! Check to make sure, and if not, throw an error.

    if ( nspecies .le. 1 ) then
       call bl_error("Network cannot have a nonpositive number of species.")
    endif

    if ( nreactions .le. 0 ) then
       call bl_error("Network cannot have a negative number of reactions.")
    endif

    if ( naux .lt. 0 ) then
       call bl_error("Network cannot have a negative number of auxiliary variables.")
    endif

    network_initialized = .true.

  end subroutine network_init


  subroutine network_close

    implicit none
    
    call actual_network_close

    network_initialized = .false.

  end subroutine network_close


  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspecies
       if (name == spec_names(n)) then
          r = n
          return
       endif
    enddo

  end function network_species_index

end module network
