module main_module

  use amrex_fort_module, only : amrex_real

  implicit none

  integer ::  fuel_ID, oxy_ID, bath_ID

contains

    subroutine extern_init(name,namlen,fuel_ID_in,oxy_ID_in,bath_ID_in) bind(C, name="extern_init")

    use, intrinsic :: iso_c_binding
    use network
    use eos_module
    use transport_module

    integer :: namlen
    integer :: name(namlen)

    real (kind=amrex_real) :: small_temp = 1.d-200
    real (kind=amrex_real) :: small_dens = 1.d-200

    integer(c_int), intent(in) :: fuel_ID_in, oxy_ID_in, bath_ID_in

    fuel_ID = fuel_ID_in + 1
    oxy_ID  = oxy_ID_in + 1
    bath_ID = bath_ID_in + 1

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    call transport_init()

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module

    call transport_close()

  end subroutine extern_close


  !subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

  !  use network, only : nspec

  !  implicit none

  !  integer, intent(out) :: nspec_out

  !  nspec_out = nspec

  !end subroutine get_num_spec

  subroutine initialize_data( &
       lo,hi, &
       Y,            Y_lo,  Y_hi, &
       temperature,  t_lo,  t_hi, &
       dens,         d_lo,  d_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspec
    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::   Y_lo(3),  Y_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   d_lo(3),  d_hi(3)
    real(amrex_real), intent(in   ) ::   dx(3)
    real(amrex_real), intent(in   ) ::   plo(3),   phi(3)
    real(amrex_real), intent(inout) ::   Y(Y_lo(1):Y_hi(1),Y_lo(2):Y_hi(2),Y_lo(3):Y_hi(3),nspec)
    real(amrex_real), intent(inout) ::   temperature( t_lo(1): t_hi(1), t_lo(2): t_hi(2), t_lo(3): t_hi(3))
    real(amrex_real), intent(inout) ::   dens( d_lo(1): d_hi(1), d_lo(2): d_hi(2), d_lo(3): d_hi(3))

    ! local variables
    integer          :: i, j, k
    real(amrex_real) :: Temp_lo, Temp_hi, dTemp, P(3), L(3), x, y_var, z, pressure
    type(eos_t) :: eos_state

    call build(eos_state)

    Temp_lo = 1500.d0
    Temp_hi = 2000.d0
    dTemp = 5.d0
    
    if (nspec.lt.3) then
       stop 'This step assumes that there are at least 3 species'
    endif
    eos_state%molefrac          = 0.d0
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
          y_var = plo(2) + (j+HALF)*dx(2)
          do i = lo(1),hi(1)
             x = plo(1) + (i+HALF)*dx(1)

             eos_state % p        = pressure
             eos_state % T        = Temp_lo + (Temp_hi-Temp_lo)*y_var/L(2) + dTemp*SIN(TWO*M_PI*y_var/P(2))

             !eos_state % massfrac(nspec) = ONE - sum(eos_state % massfrac(1:nspec-1))

             call eos_tp(eos_state)

             dens(i,j,k)         = eos_state % rho
             Y(i,j,k,1:nspec)    = eos_state % massfrac !* eos_state % rho
             temperature(i,j,k)  = eos_state % T

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine initialize_data


  !subroutine react_state(lo,hi, &
  !                       mold,mo_lo,mo_hi, &
  !                       eold,eo_lo,eo_hi, &
  !                       Told,To_lo,To_hi, &
  !                       mnew,mn_lo,mn_hi, &
  !                       enew,en_lo,en_hi, &
  !                       Tnew,Tn_lo,Tn_hi, &
  !                       ysrc,ys_lo,ys_hi, &
  !                       esrc,es_lo,es_hi, &
  !                       mask,m_lo,m_hi, &
  !                       cost,c_lo,c_hi, &
  !                       time,dt_react) bind(C, name="react_state")

  !  use network           , only : nspec
  !  use react_type_module
  !  use reactor_module, only : react
  !  use react_type_module

  !  implicit none

  !  integer          ::    lo(3),    hi(3)
  !  integer          :: mo_lo(3), mo_hi(3)
  !  integer          :: eo_lo(3), eo_hi(3)
  !  integer          :: To_lo(3), To_hi(3)
  !  integer          :: mn_lo(3), mn_hi(3)
  !  integer          :: en_lo(3), en_hi(3)
  !  integer          :: Tn_lo(3), Tn_hi(3)
  !  integer          :: ys_lo(3), ys_hi(3)
  !  integer          :: es_lo(3), es_hi(3)
  !  integer          ::  m_lo(3),  m_hi(3)
  !  integer          ::  c_lo(3),  c_hi(3)
  !  real(amrex_real) :: mold(mo_lo(1):mo_hi(1),mo_lo(2):mo_hi(2),mo_lo(3):mo_hi(3),nspec)
  !  real(amrex_real) :: eold(eo_lo(1):eo_hi(1),eo_lo(2):eo_hi(2),eo_lo(3):eo_hi(3))
  !  real(amrex_real) :: Told(To_lo(1):To_hi(1),To_lo(2):To_hi(2),To_lo(3):To_hi(3))
  !  real(amrex_real) :: mnew(mn_lo(1):mn_hi(1),mn_lo(2):mn_hi(2),mn_lo(3):mn_hi(3),nspec)
  !  real(amrex_real) :: enew(en_lo(1):en_hi(1),en_lo(2):en_hi(2),en_lo(3):en_hi(3))
  !  real(amrex_real) :: Tnew(Tn_lo(1):Tn_hi(1),Tn_lo(2):Tn_hi(2),Tn_lo(3):Tn_hi(3))
  !  real(amrex_real) :: ysrc(ys_lo(1):ys_hi(1),ys_lo(2):ys_hi(2),ys_lo(3):ys_hi(3),nspec)
  !  real(amrex_real) :: esrc(es_lo(1):es_hi(1),es_lo(2):es_hi(2),es_lo(3):es_hi(3))
  !  integer          :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
  !  real(amrex_real) :: cost(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
  !  real(amrex_real) :: time, dt_react, pressure

  !  integer          :: i, j, k

  !  real(amrex_real) ::    rY(nspec+1), rY_src(nspec)
  !  real(amrex_real) ::    energy, energy_src

  !  pressure = 1013250.d0

  !  do k = lo(3), hi(3)
  !     do j = lo(2), hi(2)
  !        do i = lo(1), hi(1)

  !           if (mask(i,j,k) .eq. 1) then
  !              rY(1:nspec)      = mold(i,j,k,1:nspec)
  !              rY_src(1:nspec)  = ysrc(i,j,k,1:nspec)
  !              rY(nspec+1)      = Told(i,j,k)
  !              energy           = eold(i,j,k)
  !              energy_src       = esrc(i,j,k)

  !              cost(i,j,k) = react(rY, rY_src,&
  !                                  energy, energy_src,&
  !                                  pressure,dt_react,time,0)

  !              enew(i,j,k)         = energy 
  !              Tnew(i,j,k)         = rY(nspec+1)
  !              mnew(i,j,k,1:nspec) = rY(1:nspec)
  !           end if

  !        end do
  !     enddo
  !  enddo

  !end subroutine react_state

end module main_module
