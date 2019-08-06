! This is a constant gamma equation of state
!
! This a simplified version of the more general eos_gamma_general.
!

module eos_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use eos_type_module
  use fuego_chemistry
  use chemistry_module, only : nspecies, Ru, inv_mwt, chemistry_init, chemistry_initialized, spec_names, elem_names

  implicit none
  character (len=64) :: eos_name = "fuego"
  logical, save, private :: initialized = .false.

  real(amrex_real), save, public :: smallT = 1.d-50

  public :: eos_init, eos_xty, eos_ytx, eos_ytx2, eos_ytx_vec, &
          eos_cpi, eos_hi, eos_hi_vec, eos_cv, eos_cp, eos_p_wb, eos_wb,&
          eos_get_activity, eos_rt, eos_tp, eos_rp, eos_re, eos_ps,&
          eos_ph, eos_th, eos_rh, eos_get_transport, eos_h, eos_deriv, &
          eos_mui, eos_get_activity_h
  private :: nspecies, Ru, inv_mwt

  interface
     subroutine amrex_array_init_snan (p, nelem) bind(C,name="amrex_array_init_snan")
       use iso_c_binding, only : c_double, c_size_t
       implicit none
       real(c_double),intent(inout) :: p
       integer (kind=c_size_t),intent(in),value :: nelem
     end subroutine amrex_array_init_snan
  end interface

contains

  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module
    use iso_c_binding, only : c_double, c_size_t

    implicit none

    real(amrex_real), optional :: small_temp
    real(amrex_real), optional :: small_dens

    integer (kind=c_size_t) :: nelem

    nelem = 1
    call amrex_array_init_snan(mintemp,nelem)
    call amrex_array_init_snan(maxtemp,nelem)
    call amrex_array_init_snan(mindens,nelem)
    call amrex_array_init_snan(maxdens,nelem)
    call amrex_array_init_snan(minmassfrac,nelem)
    call amrex_array_init_snan(maxmassfrac,nelem)
    call amrex_array_init_snan(mine,nelem)
    call amrex_array_init_snan(maxe,nelem)
    call amrex_array_init_snan(minp,nelem)
    call amrex_array_init_snan(maxp,nelem)
    call amrex_array_init_snan(mins,nelem)
    call amrex_array_init_snan(maxs,nelem)
    call amrex_array_init_snan(minh,nelem)
    call amrex_array_init_snan(maxh,nelem)

    mintemp     = 1.d-200
    maxtemp     = 1.d200
    mindens     = 1.d-200
    maxdens     = 1.d200
    minmassfrac = 1.d-200
    maxmassfrac = 1.d0
    mine        = -1.d200
    maxe        = +1.d200
    minp        = 1.d-200
    maxp        = +1.d200
    mins        = -1.d200
    maxs        = +1.d200
    minh        = -1.d200
    maxh        = +1.d200

    if (.not. chemistry_initialized)  call chemistry_init()

    if (present(small_temp)) then
       if (small_temp < mintemp) then
          small_temp = mintemp
       else
          mintemp = small_temp
       endif
    endif

    if (present(small_dens)) then
       if (small_dens < mindens) then
          small_dens = mindens
       else
          mindens = small_dens
       endif
    endif

    initialized = .true.

  end subroutine eos_init

  subroutine eos_bottom(state)

    use amrex_constants_module
    use amrex_error_module
    implicit none

    type (eos_t), intent(inout) :: state
    real(amrex_real) :: Cvx

    call ckcvms(state % T, state % cvi)  ! erg/gi.K
    call ckcpms(state % T, state % cpi)  ! erg/gi.K
    call ckhms (state % T, state % hi)    ! erg/gi

    state % cv = sum(state % massfrac(:) * state % cvi(:)) ! erg/g.K
    state % cp = sum(state % massfrac(:) * state % cpi(:)) ! erg/g.K
    state % h  = sum(state % massfrac(:) * state %  hi(:)) ! erg/g

    Cvx = state % wbar  *  state % cv ! erg/mole.K

    state % gam1 = (Cvx + Ru) / Cvx ! -
    state % cs = sqrt(state % gam1 * state % p / state % rho) ! cm/s

    state % dpdr_e = state % p / state % rho
    state % dpde = (state % gam1 - ONE) * state % rho

    ! Try to avoid the expensive log function.  Since we don't need entropy
    ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
    state % s = ONE

    ! Actually not sure what this is...used in composition derivatives
    ! in general system (not yet supported)
    state % dPdr = ZERO

  end subroutine eos_bottom

  subroutine eos_bottom_h(state)

    use amrex_constants_module
    use amrex_error_module
    implicit none

    type (eos_t), intent(inout) :: state
    real(amrex_real) :: Cvx

    call ckcvms(state % T, state % cvi)  ! erg/gi.K
    call ckcpms(state % T, state % cpi)  ! erg/gi.K
    call ckums (state % T, state % ei)    ! erg/gi

    state % cv = sum(state % massfrac(:) * state % cvi(:)) ! erg/g.K
    state % cp = sum(state % massfrac(:) * state % cpi(:)) ! erg/g.K
    state % e  = sum(state % massfrac(:) * state %  ei(:)) ! erg/g

    Cvx = state % wbar  *  state % cv ! erg/mole.K

    state % gam1 = (Cvx + Ru) / Cvx ! -
    state % cs = sqrt(state % gam1 * state % p / state % rho) ! cm/s

    state % dpdr_e = state % p / state % rho
    state % dpde = (state % gam1 - ONE) * state % rho

    ! Try to avoid the expensive log function.  Since we don't need entropy
    ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
    state % s = ONE

    ! Actually not sure what this is...used in composition derivatives
    ! in general system (not yet supported)
    state % dPdr = ZERO

  end subroutine eos_bottom_h

  subroutine eos_xty(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckxty(state % molefrac,state % massfrac)

  end subroutine eos_xty

  subroutine eos_ytx(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckytx (state % massfrac,state % molefrac)

  end subroutine eos_ytx

  subroutine eos_ytx2(Y, X, Nsp)

    implicit none

    double precision, intent(in), dimension(1:Nsp) :: Y
    double precision, intent(out), dimension(1:Nsp) :: X
    integer, intent(in) :: Nsp

    call ckytx(Y(:),X(:))

  end subroutine eos_ytx2

  subroutine eos_ytx_vec(Y, ylo, yhi, X, xlo, xhi, lo, hi, Nsp)

    implicit none

    integer, intent(in) :: ylo(3), yhi(3)
    integer, intent(in) :: xlo(3), xhi(3)    
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: Nsp
    double precision, intent(in), dimension(ylo(1)-1:yhi(1)+1, ylo(2)-1:yhi(2)+1, ylo(3)-1:yhi(3)+1, 1:Nsp) :: Y
    double precision, intent(out), dimension(xlo(1)-1:xhi(1)+1, xlo(2)-1:xhi(2)+1, xlo(3)-1:xhi(3)+1, 1:Nsp) :: X

    integer :: j, k
    integer :: npts

    npts = (hi(1)+1)-(lo(1)-1)+1
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
         call VCKYTX( npts, Y(lo(1)-1:hi(1)+1, j, k, :), X( lo(1)-1:hi(1)+1, j, k, :) )
       enddo
    enddo

  end subroutine eos_ytx_vec

  subroutine eos_cpi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcpms(state % T, state % cpi)

  end subroutine eos_cpi

  subroutine eos_hi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckhms(state % T,state % hi)

  end subroutine eos_hi

  subroutine eos_hi2(T, hi, Nsp)

    implicit none

    double precision, intent(in) :: T
    integer, intent(in)          :: Nsp
    double precision, intent(inout), dimension(1:Nsp) :: hi

    call ckhms(T,hi(:))

  end subroutine eos_hi2

  subroutine eos_hi_vec(mass, masslo, masshi, T, Tlo, Thi, hi, hilo, hihi, low, high, Nsp)

    implicit none

    integer, intent(in) :: masslo(3), masshi(3)
    integer, intent(in) :: Tlo(3), Thi(3)
    integer, intent(in) :: hilo(3), hihi(3)    
    integer, intent(in) :: low(3), high(3)    
    integer, intent(in) :: Nsp

    double precision, intent(in), dimension(masslo(1)-1:masshi(1)+1, masslo(2)-1:masshi(2)+1, masslo(3)-1:masshi(3)+1, 1:Nsp) :: mass
    double precision, intent(in), dimension(Tlo(1)-1:Thi(1)+1, Tlo(2)-1:Thi(2)+1, Tlo(3)-1:Thi(3)+1 ) :: T
    double precision, intent(out), dimension(hilo(1)-1:hihi(1)+1, hilo(2)-1:hihi(2)+1, hilo(3)-1:hihi(3)+1, 1:Nsp) :: hi
    
    integer :: j, k
    integer :: npts

    npts = (high(1)+1)-(low(1)-1)+1
    do k = low(3)-1, high(3)+1
       do j = low(2)-1, high(2)+1
          call VCKHMS( npts, T(low(1)-1:high(1)+1, j, k), hi( low(1)-1:high(1)+1, j, k, :) )
       enddo
    enddo

  end subroutine eos_hi_vec

  subroutine eos_cv(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcvbs(state % T, state % massfrac, state % cv)

  end subroutine eos_cv

  subroutine eos_cp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call ckcpbs(state % T, state % massfrac, state % cp)

  end subroutine eos_cp

  subroutine eos_p_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)
    call ckpy(state % rho, state % T, state % massfrac, state % p)

  end subroutine eos_p_wb

  subroutine eos_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  end subroutine eos_wb

  subroutine eos_get_activity(state)

    implicit none

    type (eos_t), intent(inout) :: state

    double precision :: Cvx

    call ckytcr(state%rho, state % T, state % massfrac, state % Acti)
          
    call eos_wb(state)

    call eos_bottom(state)

  end subroutine eos_get_activity

  subroutine eos_get_activity_h(state)

    implicit none

    type (eos_t), intent(inout) :: state

    double precision :: Cvx

    call ckytcr(state%rho, state % T, state % massfrac, state % Acti)
          
    call eos_wb(state)

    call eos_bottom_h(state)

  end subroutine eos_get_activity_h

  subroutine eos_rt(state)

    implicit none
    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    call ckpy(state % rho, state % T, state % massfrac, state % p)
    call ckums(state % T, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_rt

  subroutine eos_tp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    call ckrhoy(state % p,state % T,state % massfrac,state % rho)
    call ckums(state % T, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_tp

  subroutine eos_rp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_wb(state)

    state % T = state % p * state % wbar / (state % rho * Ru)
    call ckums(state % T, state % ei)
    state % e = sum(state % massfrac(:) * state % ei(:))

    call eos_bottom(state)

  end subroutine eos_rp

  subroutine eos_re(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb(state)

    call get_T_given_eY(state % e, state % massfrac, state % T, lierr)
    if (lierr .ne. 0) then
       print *, 'EOS: get_T_given_eY failed, T, e, Y = ', &
            state % T, state % e, state % massfrac
    end if
    state % T = max(state % T, smallT)
    call ckums(state % T, state % ei)
    call ckpy(state % rho, state % T, state % massfrac, state % p)

    call eos_bottom(state)

  end subroutine eos_re

  subroutine eos_ps(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_ps is not supported in this EOS.')

  end subroutine eos_ps

  subroutine eos_ph(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb(state)

    call get_T_given_hY(state % h, state % massfrac, state % T, lierr)
    if (lierr .ne. 0) then
            print *, 'EOS: get_T_given_hY failed, T, h, Y = ', &
                    state % T, state % h, state % massfrac
    end if
    state % T = max(state % T, smallT)
    call ckhms(state % T, state % hi)
    call ckrhoy(state % p, state % T, state % massfrac, state % rho)

    call eos_bottom_h(state)

  end subroutine eos_ph

  subroutine eos_th(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_th is not supported in this EOS.')

  end subroutine eos_th

  subroutine eos_rh(state)

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: lierr

    call eos_wb(state)

    call get_T_given_hY(state % h, state % massfrac, state % T, lierr)
    if (lierr .ne. 0) then
            print *, 'EOS: get_T_given_hY failed, T, h, Y = ', &
                    state % T, state % h, state % massfrac
    end if
    state % T = max(state % T, smallT)
    call ckhms(state % T, state % hi)
    call ckpy(state % rho, state % T, state % massfrac, state % p)

    call eos_bottom_h(state)

  end subroutine eos_rh

  subroutine eos_get_transport(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_get_transport is not supported in this EOS.')

  end subroutine eos_get_transport

  subroutine eos_h(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_h is not supported in this EOS.')

  end subroutine eos_h

  subroutine eos_deriv(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_deriv is not supported in this EOS.')

  end subroutine eos_deriv

  subroutine eos_mui(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_mui is not supported in this EOS.')

  end subroutine eos_mui

end module eos_module

