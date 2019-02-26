module transport_type_module

  use amrex_fort_module, only : amrex_real
  use network, only: nspec
  use eos_module

  implicit none

  type :: wtr_t
     logical :: wtr_get_Ddiag = .false.
     logical :: wtr_get_Dmat = .false.
     logical :: wtr_get_mu = .false.
     logical :: wtr_get_xi = .false.
     logical :: wtr_get_lam = .false.
  end type wtr_t
  
  type :: trv_t
     integer                  :: npts
     type(eos_t), allocatable :: eos_state(:)
     real(amrex_real),  allocatable :: mu(:)
     real(amrex_real),  allocatable :: xi(:)
     real(amrex_real),  allocatable :: lam(:)
     real(amrex_real),  allocatable :: Ddiag(:,:)
     real(amrex_real),  allocatable :: Dmat(:,:,:)
  end type trv_t
  
  interface build
     module procedure trv_build
  end interface build

  interface destroy
     module procedure trv_destroy
  end interface destroy


contains

  subroutine trv_build(trv,n)    
    type(trv_t), intent(inout) :: trv
    integer, intent(in) :: n
    integer :: i
    if (n.gt.0) then
       allocate(trv%eos_state(n))
       allocate(trv%mu(n))
       allocate(trv%xi(n))
       allocate(trv%lam(n))
!   need to set this up to build different D depending on flag
       allocate(trv%Ddiag(n,nspec))
       do i=1,n
          call build(trv%eos_state(i))
       enddo
    endif
    trv%npts = n
  end subroutine trv_build
  
  subroutine trv_destroy(trv)
    type(trv_t), intent(inout) :: trv
    integer :: i
    if (trv%npts .gt. 0) then
       do i=1,trv%npts
          call destroy(trv%eos_state(i))
       enddo
       deallocate(trv%eos_state)
       deallocate(trv%mu)
       deallocate(trv%xi)
       deallocate(trv%lam)
       deallocate(trv%Ddiag)
       trv%npts = 0
    endif
  end subroutine trv_destroy
  
end module transport_type_module
