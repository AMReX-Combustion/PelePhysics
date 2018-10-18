module vode_module

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, save :: verbose, itol, neq, order, maxstep
  integer, save :: MF

  logical, save :: use_ajac, save_ajac, always_new_j, stiff

  real(amrex_real), save :: rtol(1), atol(1)

  real(amrex_real), allocatable, save :: voderwork(:), voderpar(:)
  integer, allocatable, save :: vodeiwork(:), vodeipar(:)
  integer, save :: lvoderwork, lvodeiwork

!$omp threadprivate(voderwork,vodeiwork,voderpar,vodeipar)

contains

  subroutine vode_init(neq_in,vode_verbose,vode_itol,vode_rtol,vode_atol,vode_order,&
       vode_maxstep,vode_use_ajac,vode_save_ajac,vode_always_new_j,vode_stiff)
    implicit none

    integer, intent(in) :: neq_in, vode_verbose, vode_itol, vode_order, vode_maxstep
    logical, intent(in) :: vode_use_ajac,vode_save_ajac,vode_always_new_j,vode_stiff
    real(amrex_real), intent(in) :: vode_rtol,vode_atol

    neq          = neq_in
    verbose      = vode_verbose
    itol         = vode_itol
    rtol(1)      = vode_rtol
    atol(1)      = vode_atol
    order        = vode_order
    maxstep      = vode_maxstep
    use_ajac     = vode_use_ajac
    save_ajac    = vode_save_ajac
    always_new_j = vode_always_new_j
    stiff        = vode_stiff

    if (.not. stiff) then

       MF = 10
       lvoderwork = 20+16*NEQ
       lvodeiwork = 30

    else 

       if (use_ajac) then
          if (save_ajac) then
             MF = 21
             lvoderwork = 22 + 9*NEQ + 2*NEQ**2
          else
             MF = -21
             lvoderwork = 22 + 9*NEQ +   NEQ**2
          end if
       else
          MF = 22
          lvoderwork = 22 + 9*NEQ + 2*NEQ**2
       end if

       lvodeiwork = 30 + NEQ

    end if

    !$omp parallel
    allocate(voderwork(lvoderwork))
    allocate(vodeiwork(lvodeiwork))

    voderwork = 0.d0
    vodeiwork = 0
    vodeiwork(5) = order
    vodeiwork(6) = maxstep

    allocate(voderpar(2))
    allocate(vodeipar(1))

    call setfirst(.true.)
    !$omp end parallel

  end subroutine vode_init


  subroutine vode_close()
    neq = -1
    !$omp parallel
    if (allocated(voderwork)) deallocate(voderwork)
    if (allocated(vodeiwork)) deallocate(vodeiwork)
    if (allocated(voderpar))  deallocate(voderpar)
    if (allocated(vodeipar))  deallocate(vodeipar)
    !$omp end parallel
  end subroutine vode_close

end module vode_module

