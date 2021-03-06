module main_module

  use amrex_fort_module, only : amrex_real

#include "mechanism.H"

  implicit none

contains

    subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use eos_module, only : eos_init
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

  subroutine initialize_data( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       energy,      e_lo,  e_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use eos_type_module
    use eos_module, only: eos_rt

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    integer         , intent(in   ) ::   e_lo(3),  e_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::    plo(3),   phi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),NUM_SPECIES)
    real(amrex_real), intent(inout) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(inout) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(amrex_real), intent(inout) :: energy(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))

    ! local variables
    integer          :: i, j, k, n
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y, z, Y_lo(NUM_SPECIES), Y_hi(NUM_SPECIES)
    real(amrex_real) :: Rho_lo, Rho_hi, dRho
    real(amrex_real) :: energy_lo, energy_hi, denergy
    real(amrex_real) :: sum

    type (eos_t) :: eos_state

    call build(eos_state)

    Temp_lo = 500.d0
    Temp_hi = 500.d0
    dTemp = 5.d0

    Rho_lo = 1.d-2
    Rho_hi = 1.d-2
    dRho = 5.d-3

    Y_lo(:) = ZERO
    Y_lo(1) = ONE
    Y_hi(:) = ONE / NUM_SPECIES

    energy_lo   = 0.0d0
    energy_hi   = 0.0d0 
    denergy     = 100.0d0 
    
    L(:) = phi(:) - plo(:)
    P(:) = L(:) / 4
    
    do k = lo(3),hi(3)
       z = plo(3) + (k+HALF)*dx(3)
       do j = lo(2),hi(2)
          y = plo(2) + (j+HALF)*dx(2)
          do i = lo(1),hi(1)
             x = plo(1) + (i+HALF)*dx(1)

#if ( AMREX_SPACEDIM == 1 )
             temperature(i,j,k) = Temp_lo
             density(i,j,k)     = Rho_lo
#else
             temperature(i,j,k) = Temp_lo + (Temp_hi-Temp_lo)*y/L(2) + dTemp*SIN(TWO*M_PI*y/P(2))
             density(i,j,k)     = Rho_lo + (Rho_hi-Rho_lo)*y/L(2) + dRho*SIN(TWO*M_PI*y/P(2))
#endif
             massfrac(i,j,k,:) = Y_lo(:) + (Y_hi(:)-Y_lo(:))*x/L(1)

             sum = ZERO
             do n=1,NUM_SPECIES-1
                sum = sum + massfrac(i,j,k,n)
             enddo
             massfrac(i,j,k,NUM_SPECIES) = ONE - sum

             ! compute energy
             eos_state % T        = temperature(i,j,k)
             eos_state % massfrac = massfrac(i,j,k,1:NUM_SPECIES)
             eos_state % rho      = density(i,j,k)

             call eos_rt(eos_state)

#if ( AMREX_SPACEDIM == 1 )
             energy(i,j,k) = eos_state % e
#else
             energy(i,j,k) = eos_state % e +  (energy_hi-energy_lo)*y/L(2) + denergy*SIN(TWO*M_PI*y/P(2))
#endif

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine initialize_data


  subroutine get_cp( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       cp,          cp_lo, cp_hi) &
       bind(C, name="get_cp")

    use eos_type_module
    use eos_module, only : eos_cp

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::  cp_lo(3), cp_hi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),NUM_SPECIES)
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
             eos_state % massfrac = massfrac(i,j,k,1:NUM_SPECIES)

             call eos_cp(eos_state)

             cp(i,j,k) = eos_state % cp

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine get_cp


  subroutine get_cv( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       cv,          cv_lo, cv_hi) &
       bind(C, name="get_cv")

    use eos_type_module
    use eos_module, only : eos_cv

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::  cv_lo(3), cv_hi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),NUM_SPECIES)
    real(amrex_real), intent(inout) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(  out) :: cv(cv_lo(1):cv_hi(1),cv_lo(2):cv_hi(2),cv_lo(3):cv_hi(3))

    ! local variables
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state % T        = temperature(i,j,k)
             eos_state % massfrac = massfrac(i,j,k,1:NUM_SPECIES)

             call eos_cv(eos_state)

             cv(i,j,k) = eos_state % cv

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine get_cv


  subroutine get_T_from_EY( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       energy,       e_lo, e_hi) &
       bind(C, name="get_T_from_EY")

    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    integer         , intent(in   ) ::   e_lo(3),  e_hi(3)
    real(amrex_real), intent(inout) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),NUM_SPECIES)
    real(amrex_real), intent(  out) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(amrex_real), intent(inout) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(amrex_real), intent(inout) :: energy(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))

    ! local variables
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state % massfrac = massfrac(i,j,k,1:NUM_SPECIES)
             eos_state % e        = energy(i,j,k)
             eos_state % rho      = density(i,j,k)

             call eos_re(eos_state)

             temperature(i,j,k)   = eos_state % T

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine get_T_from_EY

end module main_module
