module actual_transport_module

  use eos_type_module
  use transport_type_module

  implicit none

  character (len=64) :: transport_name = "constant"

contains

  subroutine actual_transport_init
    implicit none
  end subroutine actual_transport_init


  subroutine actual_transport_close
    implicit none
  end subroutine actual_transport_close

  
  subroutine actual_transport(which, coeff)

    use extern_probin_module, only: const_conductivity, const_viscosity,&
         const_bulk_viscosity, const_diffusivity, mks_unit

    implicit none

    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff

!   By default the units are in CGS
!   Below mks_unit is a flag to let the user to provide MKS data, and it makes the conversion to CGS

    if (which % wtr_get_lam) then
       coeff % lam = const_conductivity
       if (mks_unit) coeff % lam = coeff % lam * 1.0d5
    endif

    if (which % wtr_get_mu) then
       coeff % mu = const_viscosity
       if (mks_unit) coeff % mu = coeff % mu * 10.0d0
    endif

    if (which % wtr_get_xi) then
       coeff % xi = const_bulk_viscosity
       if (mks_unit) coeff % xi = coeff % xi * 10.0d0
    endif

    if (which % wtr_get_Ddiag) then
       coeff % Ddiag = const_diffusivity
       if (mks_unit) coeff % Ddiag = coeff % Ddiag * 10.0d0
    endif

  end subroutine actual_transport
  
end module actual_transport_module

