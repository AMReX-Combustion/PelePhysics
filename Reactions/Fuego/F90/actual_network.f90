module actual_network

  use fuego_chemistry 
  use chemistry_module, only : nspecies, nelements, nreactions, chemistry_init, &
       chemistry_close, chemistry_initialized, spec_names, elem_names, &
       L_spec_name, L_elem_name

  implicit none

  integer :: naux
  character (len=16), save, allocatable :: aux_names(:)

contains
  
  subroutine actual_network_init

    use extern_probin_module, only: numaux, auxnamesin
    use amrex_error_module

    implicit none

    integer :: iaux,nnames

    if (.not. chemistry_initialized)  call chemistry_init()

    ! Get auxiliary variables (if any)
    naux = numaux
    allocate(aux_names(naux))
    nnames = count(transfer(auxnamesin, 'a', len(auxnamesin)) == ",") +1
    if (naux .gt. 0) then 
       if (naux .ne. nnames) then
          call amrex_error('simple ::actual_network_init wrong number of aux variable names')
       end if
       read(auxnamesin,*) aux_names
       do iaux = 1,naux
          aux_names(iaux) = trim(adjustl(aux_names(iaux)))
       end do
    end if

  end subroutine actual_network_init

  subroutine actual_network_close

    call chemistry_close()
    deallocate(aux_names)

  end subroutine actual_network_close

end module actual_network
