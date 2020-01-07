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

    if (which % wtr_get_lam) then
       coeff % lam = const_conductivity
    endif

    if (which % wtr_get_mu) then
       coeff % mu = const_viscosity
    endif

    if (which % wtr_get_xi) then
       coeff % xi = const_bulk_viscosity
    endif

    if (which % wtr_get_Ddiag) then
       coeff % Ddiag = const_diffusivity
    endif

! By default the units are in CGS
! Below is a flag to let the user to provide MKS data, and it makes the conversion to CGS
    if (mks_unit) then
      if (which % wtr_get_lam) coeff % lam = coeff % lam * 1.0d5
      if (which % wtr_get_mu) coeff % mu = coeff % mu * 10.0d0
      if (which % wtr_get_xi) coeff % xi = coeff % xi * 10.0d0
      if (which % wtr_get_Ddiag) coeff % Ddiag = coeff % Ddiag * 10.0d0
    endif

  end subroutine actual_transport
  
end module actual_transport_module

