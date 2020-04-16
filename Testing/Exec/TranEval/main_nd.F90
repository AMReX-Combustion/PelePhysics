module main_module

  use amrex_fort_module, only : amrex_real

#include "mechanism.h"

  implicit none

contains

    subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use eos_module
    use transport_module
    use fuego_chemistry, only : network_init 

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

    call transport_init_F()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close_F()

  end subroutine extern_close


  subroutine get_num_spec(nspecies_out) bind(C, name="get_num_spec")

    implicit none

    integer, intent(out) :: nspecies_out

    nspecies_out = NUM_SPECIES

  end subroutine get_num_spec

  subroutine initialize_data( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::    plo(3),   phi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),NUM_SPECIES)
    real(amrex_real), intent(inout) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(inout) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))

    ! local variables
    integer          :: i, j, k, n
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y, z, Y_lo(NUM_SPECIES), Y_hi(NUM_SPECIES)
    real(amrex_real) :: Rho_lo, Rho_hi, dRho
    real(amrex_real) :: sum

    Temp_lo = 500.d0
    Temp_hi = 500.d0
    dTemp = 5.d0

    Rho_lo = 1.d-2
    Rho_hi = 1.d-2
    dRho = 5.d-3

    Y_lo(:) = ZERO
    Y_lo(1) = ONE
    Y_hi(:) = ONE / NUM_SPECIES
    
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
             density(i,j,k) = Rho_lo + (Rho_hi-Rho_lo)*y/L(2) + dRho*SIN(TWO*M_PI*y/P(2))

             sum = ZERO
             do n=1,NUM_SPECIES-1
                sum = sum + massfrac(i,j,k,n)
             enddo
             massfrac(i,j,k,NUM_SPECIES) = ONE - sum

          end do
       end do
    end do

  end subroutine initialize_data


end module main_module
