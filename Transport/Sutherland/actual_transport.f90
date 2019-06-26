module actual_transport_module

  use amrex_error_module
  use amrex_constants_module
  use eos_type_module
  use transport_type_module

  use amrex_fort_module, only : amrex_real

  implicit none

  character (len=64) :: transport_name = "sutherland"
  logical, save, private :: suth_initialized = .false.
  integer, save, private :: npts_suth = 0
  real(amrex_real), save, allocatable :: Tp(:), Cpp(:)

  !$omp threadprivate(suth_initialized,npts_suth,Tp,Cpp)

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine actual_transport_init

    implicit none

    suth_initialized = .true.

  end subroutine actual_transport_init


  subroutine actual_transport_close

    implicit none

    suth_initialized = .false.

  end subroutine actual_transport_close


  subroutine build_internal(npts)
    integer, intent(in) :: npts

    if (npts_suth .ne. npts .and. npts.gt.0) then
       if (npts_suth .ne. 0) then
          call destroy_internal()
       endif
       allocate(Tp(npts))
       allocate(Cpp(npts))
       npts_suth = npts
    endif

  end subroutine build_internal


  subroutine destroy_internal

    deallocate(Tp)
    deallocate(Cpp)
    npts_suth = 0

  end subroutine destroy_internal



  subroutine actual_transport(which, coeff)

    use amrex_error_module
    use extern_probin_module, only: Prandtl_number, viscosity_mu_ref, viscosity_T_ref, viscosity_S,&
         const_bulk_viscosity, const_diffusivity

    implicit none
    
    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff
    integer :: i

    if (.not. suth_initialized) then
       call amrex_error('Sutherland::actual_transport called before initialized')
    endif

    if (npts_suth .ne. coeff%npts) then
       call build_internal(coeff%npts)
    endif

    ! load input data into local vectors
    do i=1,coeff%npts
       call eos_cp(coeff % eos_state(i))
       Cpp(i) = coeff % eos_state(i) % cp
    enddo

    if (which % wtr_get_mu) then
       coeff % mu(:) = viscosity_mu_ref &
            * (coeff % eos_state(:) % T / viscosity_T_ref)**1.5 &
            * (    viscosity_T_ref      + viscosity_S) &
            / (coeff % eos_state(:) % T + viscosity_S)
    endif


    if (which % wtr_get_lam) then
       if (.not. which % wtr_get_mu) then
          coeff % mu(:) = viscosity_mu_ref &
               * (coeff % eos_state(:) % T / viscosity_T_ref)**1.5 &
               * (    viscosity_T_ref      + viscosity_S) &
               / (coeff % eos_state(:) % T + viscosity_S)
       endif
       coeff % lam = coeff % mu * Cpp / Prandtl_number
    endif

    if (which % wtr_get_xi) then
       coeff % xi = const_bulk_viscosity
    endif

    if (which % wtr_get_Ddiag) then
       coeff % Ddiag = const_diffusivity
    endif


  end subroutine actual_transport

end module actual_transport_module

