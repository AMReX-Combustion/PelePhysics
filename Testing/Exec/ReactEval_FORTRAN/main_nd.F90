module main_module

  use amrex_fort_module, only : amrex_real

  implicit none

contains

    subroutine extern_init(name,namlen) bind(C, name="extern_init")

    use fuego_chemistry
    use eos_module
    use transport_module

    implicit none
    integer :: namlen
    integer :: name(namlen)

    real (kind=amrex_real) :: small_temp = 1.d-200
    real (kind=amrex_real) :: small_dens = 1.d-200

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

  end subroutine extern_init

  subroutine extern_init_reactor() bind(C, name="extern_init_reactor")

#ifdef USE_SUNDIALS_PP
    use cvode_module      , only : reactor_init

    implicit none

    integer :: ncells(1), iiE(1)

    iiE(1)    = 1
    ncells(1) = 1

    call reactor_init(iiE(1),ncells(1))
#else
    use reactor_module

    implicit none

    call reactor_init(1)
#endif


  end subroutine extern_init_reactor


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    implicit none

    call transport_close()

  end subroutine extern_close


  subroutine get_num_spec(nspecies_out) bind(C, name="get_num_spec")

    use fuego_chemistry, only : nspecies

    implicit none

    integer, intent(out) :: nspecies_out

    nspecies_out = nspecies

  end subroutine get_num_spec

  subroutine initialize_data( &
       lo,hi, &
       rhoY,         rY_lo, rY_hi, &
       temperature,  t_lo,  t_hi, &
       eint,         e_lo,  e_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use fuego_chemistry, only: nspecies
    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  rY_lo(3), rY_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   e_lo(3),  e_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::    plo(3),   phi(3)
    real(amrex_real), intent(inout) ::        rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspecies)
    real(amrex_real), intent(inout) :: temperature( t_lo(1): t_hi(1), t_lo(2): t_hi(2), t_lo(3): t_hi(3))
    real(amrex_real), intent(inout) ::        eint( e_lo(1): e_hi(1), e_lo(2): e_hi(2), e_lo(3): e_hi(3))

    ! local variables
    integer          :: i, j, k
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y, z, pressure
    type(eos_t) :: eos_state

    call build(eos_state)

    Temp_lo = 1500.d0
    Temp_hi = 2000.d0
    dTemp = 5.d0

    if (nspecies.lt.3) then
       stop 'This step assumes that there are at least 3 species'
    endif
    eos_state%molefrac = 0.d0
    eos_state%molefrac(1) = 0.2d0
    eos_state%molefrac(2) = 0.1d0
    eos_state%molefrac(nspecies) = 1.d0 - eos_state%molefrac(1) - eos_state%molefrac(2)
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

             eos_state % massfrac(nspecies) = ONE - sum(eos_state % massfrac(1:nspecies-1))

             call eos_tp(eos_state)

#ifdef USE_SUNDIALS_PP
             eint(i,j,k) = eos_state % e * eos_state % rho
#else
             eint(i,j,k) = eos_state % e
#endif
             rhoY(i,j,k,1:nspecies) = eos_state % massfrac * eos_state % rho
             temperature(i,j,k) = eos_state % T

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine initialize_data


  subroutine react_state(lo,hi, &
                         mold,mo_lo,mo_hi, &
                         eold,eo_lo,eo_hi, &
                         Told,To_lo,To_hi, &
                         mnew,mn_lo,mn_hi, &
                         enew,en_lo,en_hi, &
                         Tnew,Tn_lo,Tn_hi, &
                         ysrc,ys_lo,ys_hi, &
                         esrc,es_lo,es_hi, &
                         mask,m_lo,m_hi, &
                         cost,c_lo,c_hi, &
                         time,dt_react) bind(C, name="react_state")

    use fuego_chemistry           , only : nspecies
#ifdef USE_SUNDIALS_PP
    use cvode_module      , only : react
#else
    use reactor_module    , only : react
#endif

    implicit none

    integer          ::    lo(3),    hi(3)
    integer          :: mo_lo(3), mo_hi(3)
    integer          :: eo_lo(3), eo_hi(3)
    integer          :: To_lo(3), To_hi(3)
    integer          :: mn_lo(3), mn_hi(3)
    integer          :: en_lo(3), en_hi(3)
    integer          :: Tn_lo(3), Tn_hi(3)
    integer          :: ys_lo(3), ys_hi(3)
    integer          :: es_lo(3), es_hi(3)
    integer          ::  m_lo(3),  m_hi(3)
    integer          ::  c_lo(3),  c_hi(3)
    real(amrex_real) :: mold(mo_lo(1):mo_hi(1),mo_lo(2):mo_hi(2),mo_lo(3):mo_hi(3),nspecies)
    real(amrex_real) :: eold(eo_lo(1):eo_hi(1),eo_lo(2):eo_hi(2),eo_lo(3):eo_hi(3))
    real(amrex_real) :: Told(To_lo(1):To_hi(1),To_lo(2):To_hi(2),To_lo(3):To_hi(3))
    real(amrex_real) :: mnew(mn_lo(1):mn_hi(1),mn_lo(2):mn_hi(2),mn_lo(3):mn_hi(3),nspecies)
    real(amrex_real) :: enew(en_lo(1):en_hi(1),en_lo(2):en_hi(2),en_lo(3):en_hi(3))
    real(amrex_real) :: Tnew(Tn_lo(1):Tn_hi(1),Tn_lo(2):Tn_hi(2),Tn_lo(3):Tn_hi(3))
    real(amrex_real) :: ysrc(ys_lo(1):ys_hi(1),ys_lo(2):ys_hi(2),ys_lo(3):ys_hi(3),nspecies)
    real(amrex_real) :: esrc(es_lo(1):es_hi(1),es_lo(2):es_hi(2),es_lo(3):es_hi(3))
    integer          :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(amrex_real) :: cost(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
#ifdef USE_SUNDIALS_PP
    real(amrex_real) :: time, dt_react
#else
    real(amrex_real) :: time, dt_react, pressure
#endif

    integer          :: i, j, k

    real(amrex_real) ::    rY(nspecies+1), rY_src(nspecies)
#ifdef USE_SUNDIALS_PP
    real(amrex_real) ::    energy(1), energy_src(1)
#else
    real(amrex_real) ::    energy, energy_src
#endif


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (mask(i,j,k) .eq. 1) then
                rY(1:nspecies)      = mold(i,j,k,1:nspecies)
                rY_src(1:nspecies)  = ysrc(i,j,k,1:nspecies)
                rY(nspecies+1)      = Told(i,j,k)
#ifdef USE_SUNDIALS_PP
                energy(1)           = eold(i,j,k)
                energy_src(1)       = esrc(i,j,k)
#else
                energy           = eold(i,j,k)
                energy_src       = esrc(i,j,k)

                pressure = 1013250.d0
#endif

                cost(i,j,k) = react(rY, rY_src,&
#ifdef USE_SUNDIALS_PP
                                    energy(1), energy_src(1),&
                                    dt_react,time)
#else
                                    energy, energy_src,&
                                    pressure,&
                                    dt_react,time)
#endif

#ifdef USE_SUNDIALS_PP
                enew(i,j,k)            = energy(1) 
#else
                enew(i,j,k)            = energy 
#endif
                Tnew(i,j,k)            = rY(nspecies+1)
                mnew(i,j,k,1:nspecies) = rY(1:nspecies)
             end if

          end do
       enddo
    enddo

  end subroutine react_state

end module main_module
