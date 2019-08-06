module actual_network

  use fuego_chemistry 
  use network, only : nspecies
  use chemistry_module, only : nreactions, chemistry_init, &
       chemistry_close, chemistry_initialized, spec_names, elem_names, &
       L_spec_name, L_elem_name

  implicit none

  integer :: naux
  character (len=16), save, allocatable :: aux_names(:)

contains
  
  subroutine actual_network_init

    if (.not. chemistry_initialized)  call chemistry_init()
    naux = 0
    allocate(aux_names(naux))

  end subroutine actual_network_init

  subroutine actual_network_close

    call chemistry_close()
    deallocate(aux_names)

  end subroutine actual_network_close

end module actual_network
