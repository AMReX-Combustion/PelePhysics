module transport_coeff_module

  implicit none

contains

  subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use network
    use eos_module
    use transport_module

    integer :: namlen
    integer :: name(namlen)

    real (kind=dp_t) :: small_temp = 1.d-200
    real (kind=dp_t) :: small_dens = 1.d-200

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

  end subroutine extern_init

  function get_num_species() bind(C, name="get_num_species") result (n) 
    use network, only: nspec
    integer :: n
    n = nspec
  end function get_num_species


end module transport_coeff_module
