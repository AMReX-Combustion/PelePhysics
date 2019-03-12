module main_module

  use amrex_fort_module, only : amrex_real

  implicit none

  integer ::  iE_main

contains

    subroutine extern_init(name,namlen,cvode_iE_in) bind(C, name="extern_init")

    use, intrinsic :: iso_c_binding
    use network
    use eos_module
    use transport_module

    integer :: namlen
    integer :: name(namlen)

    integer(c_int), intent(in) :: cvode_iE_in

    real (kind=amrex_real) :: small_temp = 1.d-200
    real (kind=amrex_real) :: small_dens = 1.d-200

    iE_main = cvode_iE_in

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module

    call transport_close()

  end subroutine extern_close


  subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

    use network, only : nspec

    implicit none

    integer, intent(out) :: nspec_out

    nspec_out = nspec

  end subroutine get_num_spec


  subroutine initialize_data( &
       lo,           hi, &
       rhoY,         rY_lo, rY_hi, &
       rhoY_src,     rY_src_lo, rY_src_hi, &
       rhoE,         rE_lo, rE_hi, &
       rhoEs,        rEs_lo, rEs_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspec
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
    real(amrex_real), intent(inout) ::  rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspec+1)
    real(amrex_real), intent(inout) ::  rhoY_src(rY_src_lo(1):rY_src_hi(1),rY_src_lo(2):rY_src_hi(2),rY_src_lo(3):rY_src_hi(3),nspec)
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

    if (nspec.lt.3) then
       stop 'This step assumes that there are at least 3 species'
    endif
    eos_state%molefrac = 0.d0
    eos_state%molefrac(1) = 0.2d0
    eos_state%molefrac(2) = 0.1d0
    eos_state%molefrac(nspec) = 1.d0 - eos_state%molefrac(1) - eos_state%molefrac(2)
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
             eos_state % T        = Temp_lo + (Temp_hi-Temp_lo)*y/L(2) + dTemp*SIN(TWO*M_PI*y/P(2))

             eos_state % massfrac(nspec) = ONE - sum(eos_state % massfrac(1:nspec-1))

             call eos_tp(eos_state)

             ! rhoY(:nspec) = rhoY, rhoY(nspec+1) = T
             rhoY(i,j,k,1:nspec) = eos_state % massfrac * eos_state % rho
             rhoY(i,j,k,nspec+1) = eos_state % T
             ! rhoY_src(:nspec) = rhoForcingSpecs
             rhoY_src(i,j,k,1:nspec) = 0.0d0
             if (iE_main == 1) then
                 ! all in e
                 rhoE(i,j,k,1) = eos_state % e * eos_state % rho
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
