module main_module

  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use network
    use eos_module
    use transport_module

    implicit none
    integer :: namlen
    integer :: name(namlen)

    real(amrex_real) :: small_temp = 1.d-200
    real(amrex_real) :: small_dens = 1.d-200

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close()

  end subroutine extern_close


  subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

    use network, only : nspec

    implicit none

    integer, intent(out) :: nspec_out

    nspec_out = nspec

  end subroutine get_num_spec


  subroutine initialize_data( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspec

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::    plo(3),   phi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspec)
    real(amrex_real), intent(inout) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! local variables
    integer          :: i, j, k, n
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y, z, Y_lo(nspec), Y_hi(nspec)
    real(amrex_real) :: sum

    Temp_lo = 500.d0
    Temp_hi = 500.d0
    dTemp = 5.d0

    Y_lo(:) = ZERO
    Y_lo(1) = ONE
    Y_hi(:) = ONE / nspec
    
    L(:) = phi(:) - plo(:)
    P(:) = L(:) / 4
    
    do k = lo(3),hi(3)
       z = plo(3) + (k+HALF)*dx(3)
       do j = lo(2),hi(2)
          y = plo(2) + (j+HALF)*dx(2)
          do i = lo(1),hi(1)
             x = plo(1) + (i+HALF)*dx(1)

             temperature(i,j,k) = Temp_lo + (Temp_hi-Temp_lo)*y/L(2) + dTemp*SIN(TWO*M_PI*y/P(2))
             massfrac(i,j,k,:) = Y_lo(:) + (Y_hi(:)-Y_lo(:))*x/L(1)

             sum = ZERO
             do n=1,nspec-1
                sum = sum + massfrac(i,j,k,n)
             enddo
             massfrac(i,j,k,nspec) = ONE - sum

          end do
       end do
    end do
  end subroutine initialize_data

  subroutine get_cp( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       cp,          cp_lo, cp_hi) &
       bind(C, name="get_cp")

    use network, only: nspec
    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::  cp_lo(3), cp_hi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspec)
    real(amrex_real), intent(inout) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(  out) :: cp(cp_lo(1):cp_hi(1),cp_lo(2):cp_hi(2),cp_lo(3):cp_hi(3))

    ! local variables
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state % T        = temperature(i,j,k)
             eos_state % massfrac = massfrac(i,j,k,1:nspec)

             call eos_cp(eos_state)

             cp(i,j,k) = eos_state % cp

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine get_cp


end module main_module
