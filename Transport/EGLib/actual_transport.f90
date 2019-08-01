module actual_transport_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use eos_type_module
  use transport_type_module
  use egz_module, only : egz_init, egz_close, egzpar, egze3, egzk3,&
       egzl1, egzvr1, egzini, egzk1

  implicit none

  character (len=64) :: transport_name = "egz"
  logical, save, private :: egz_initialized = .false.
  integer, save, private :: iflag = 4
  integer, save, private :: npts_egz = 0
  real(amrex_real), save, allocatable :: Tp(:), Yp(:,:), Xp(:,:), Cpp(:,:), L1(:), L2(:)
  logical, parameter :: use_bulk_viscosity = .true.

  !$omp threadprivate(iflag,npts_egz,Tp,Yp,Xp,Cpp,L1,L2)

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine actual_transport_init

    implicit none

    call egz_init(use_bulk_viscosity)

    egz_initialized = .true.

  end subroutine actual_transport_init


  subroutine actual_transport_close

    implicit none

    call egz_close()

    egz_initialized = .false.

  end subroutine actual_transport_close


  subroutine build_internal(npts)
    integer, intent(in) :: npts

    call egzini(npts)
    if (npts_egz .ne. npts .and. npts.gt.0) then
       if (npts_egz .ne. 0) then
          call destroy_internal()
       endif
       allocate(Tp(npts))
       allocate(L1(npts))
       allocate(L2(npts))
       allocate(Yp(npts,nspecies))
       allocate(Xp(npts,nspecies))
       allocate(Cpp(npts,nspecies))
       npts_egz = npts
    endif

  end subroutine build_internal


  subroutine destroy_internal

    deallocate(Tp)
    deallocate(L1)
    deallocate(L2)
    deallocate(Yp)
    deallocate(Xp)
    deallocate(Cpp)
    npts_egz = 0

  end subroutine destroy_internal



  subroutine actual_transport(which, coeff)

    use amrex_error_module
    
    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff
    integer :: i

    if (.not. egz_initialized) then
       call amrex_error('EGLib::actual_transport called before initialized')
    endif

    if (npts_egz .ne. coeff%npts) then
       call build_internal(coeff%npts)
    endif

    ! load input data into local vectors
    if (iflag > 3) then    
       do i=1,coeff%npts
          call eos_cpi(coeff % eos_state(i))
          Cpp(i,:) = coeff % eos_state(i) % cpi(:)
       enddo
    else
       Cpp = 0.d0
    endif
    do i=1,coeff%npts
       call eos_ytx(coeff % eos_state(i))
       Xp(i,:)  = coeff % eos_state(i) % molefrac(:)    
       Yp(i,:)  = coeff % eos_state(i) % massfrac(:)    
       Tp(i)    = coeff % eos_state(i) % T
    enddo

    call egzpar(Tp, Xp, Cpp)
    
    if (which % wtr_get_mu) then
       call egze3(Tp, coeff % mu(:))
    endif

    if (which % wtr_get_xi) then
       call egzk3(Tp, coeff % xi(:))
!      call egzk1(0.75d0,Xp, coeff % xi(:))
    endif

    if (which % wtr_get_lam) then
       call egzl1( .25d0, Xp, L1)
 !     call egzl1( 1.d0, Xp, L1)
 !     call egzl1(-1.d0, Xp, L2)
 !     coeff % lam = 0.5d0 * (L1 + L2)
       coeff % lam = L1
    endif

    if (which % wtr_get_Ddiag) then
       call egzvr1(Tp, coeff % Ddiag)
    endif

  end subroutine actual_transport

end module actual_transport_module

