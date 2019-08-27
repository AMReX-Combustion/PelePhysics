module main_module

  use amrex_fort_module, only : amrex_real

  implicit none

  integer ::  iE_main, fuel_ID, oxy_ID, bath_ID

contains

    subroutine extern_init(name,namlen,fuel_ID_in,oxy_ID_in,bath_ID_in,cvode_iE_in) bind(C, name="extern_init")

    use, intrinsic :: iso_c_binding
    use network
    use eos_module
    use transport_module

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

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close()

  end subroutine extern_close


  subroutine get_num_spec(nspecies_out) bind(C, name="get_num_spec")

    use network, only : nspecies

    implicit none

    integer, intent(out) :: nspecies_out

    nspecies_out = nspecies

  end subroutine get_num_spec


  subroutine initialize_data( &
       lo,hi, &
       rhoY,         rY_lo, rY_hi, &
       rhoY_src,     rY_src_lo, rY_src_hi, &
       rhoE,         rE_lo, rE_hi, &
       rhoEs,        rEs_lo, rEs_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspecies
    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  rY_lo(3), rY_hi(3)
    integer         , intent(in   ) ::  rY_src_lo(3), rY_src_hi(3)
    integer         , intent(in   ) ::  rE_lo(3), rE_hi(3)
    integer         , intent(in   ) ::  rEs_lo(3), rEs_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::    plo(3),   phi(3)
    real(amrex_real), intent(inout) ::  rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies+1)
    real(amrex_real), intent(inout) ::  rhoY_src(rY_src_lo(1):rY_src_hi(1),rY_src_lo(2):rY_src_hi(2),rY_src_lo(3):rY_src_hi(3),nspecies)
    real(amrex_real), intent(inout) ::  rhoE(rE_lo(1):rE_hi(1),rE_lo(2):rE_hi(2),rE_lo(3):rE_hi(3),1)
    real(amrex_real), intent(inout) ::  rhoEs(rEs_lo(1):rEs_hi(1),rEs_lo(2):rEs_hi(2),rEs_lo(3):rEs_hi(3),1)

    ! local variables
    integer          :: i, j, k
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y, z, pressure
    type(eos_t)      :: eos_state

    call build(eos_state)

    Temp_lo = 1500.d0
    Temp_hi = 2000.d0
    dTemp = 5.d0

    if (nspecies.lt.3) then
       stop 'This step assumes that there are at least 3 species'
    endif
    eos_state%molefrac = 0.d0
    eos_state%molefrac(oxy_ID)  = 0.2d0
    eos_state%molefrac(fuel_ID) = 0.1d0
    eos_state%molefrac(bath_ID) = 1.d0 - eos_state%molefrac(fuel_ID) - eos_state%molefrac(oxy_ID)
    call eos_xty(eos_state)
    
    L(:) = phi(:) - plo(:)
    P(:) = L(:) / 4

    pressure = 1013250.d0
    
    do k = lo(3),hi(3)
       z = plo(3) + (k+HALF)*dx(3)
       do j = lo(2),hi(2)
          y = plo(2) + (j+HALF)*dx(2)
          do i = lo(1),hi(1)
             x = plo(1) + (i+HALF)*dx(1)

             eos_state % p        = pressure
             eos_state % T        = Temp_lo + (Temp_hi-Temp_lo)*y/L(2) + dTemp*SIN(TWO*M_PI*y/P(2)) !+ (Temp_hi-Temp_lo)*x/L(1) + (Temp_hi-Temp_lo)*z/L(3) 

             call eos_tp(eos_state)

             ! rhoY(:nspecies) = rhoY, rhoY(nspecies+1) = T
             rhoY(i,j,k,1:nspecies) = eos_state % massfrac * eos_state % rho
             rhoY(i,j,k,nspecies+1) = eos_state % T
             ! rhoY_src(:nspecies) = rhoForcingSpecs
             rhoY_src(i,j,k,1:nspecies) = 0.0d0
             if (iE_main == 1) then
                 ! all in e
#ifdef AMREX_USE_SUNDIALS_3x4x
                 rhoE(i,j,k,1) = eos_state % e * eos_state % rho
#else
                 rhoE(i,j,k,1) = eos_state % e
#endif
             else
                 ! all in h
                 rhoE(i,j,k,1) = eos_state % h * eos_state % rho
             end if
             ! all in h
             !rhoE src ext
             rhoEs(i,j,k,1) = 0.0d0

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine initialize_data

end module main_module
