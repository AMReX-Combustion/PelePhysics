module main_module

  use amrex_fort_module, only : amrex_real

#include "mechanism.h"

  implicit none

  integer ::  iE_main, fuel_ID, oxy_ID, bath_ID

contains

    subroutine extern_init(name,namlen,fuel_ID_in,oxy_ID_in,bath_ID_in,cvode_iE_in) bind(C, name="extern_init")

    use, intrinsic :: iso_c_binding
    use eos_module
    use transport_module

    implicit none
    integer :: namlen
    integer :: name(namlen)

    integer(c_int), intent(in) :: cvode_iE_in, fuel_ID_in, oxy_ID_in, bath_ID_in

    real (kind=amrex_real) :: small_temp = 1.d-200
    real (kind=amrex_real) :: small_dens = 1.d-200

    iE_main = cvode_iE_in
    fuel_ID = fuel_ID_in + 1
    oxy_ID  = oxy_ID_in + 1
    bath_ID = bath_ID_in + 1

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call ckinit()

    call eos_init(small_temp, small_dens)

    call transport_init_F()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close_F()

  end subroutine extern_close

end module main_module
