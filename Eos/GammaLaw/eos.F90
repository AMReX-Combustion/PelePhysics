! This is a constant gamma equation of state
!
! This a simplified version of the more general eos_gamma_general.
!

module eos_module

  use amrex_fort_module, only : amrex_real
  use eos_type_module

  implicit none

  character (len=64) :: eos_name = "gamma_law"

  real(amrex_real), save :: gamma_const
  double precision, parameter :: R = 8.314462145468952d7

  integer, parameter :: nelements = 1
  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), save :: elem_names(nelements)

  !integer, parameter :: nspecies = 1 ! already in eos_type_module
  integer, parameter :: L_spec_name = 16 ! Each species name has at most 8 characters
  character*(L_spec_name), save :: spec_names(nspecies)

  public :: eos_init, eos_xty, eos_ytx, eos_ytx_vec, eos_cpi, eos_hi, eos_hi_vec, eos_cv, eos_cp, eos_p_wb, eos_wb, eos_get_activity, eos_rt, eos_tp, eos_rp, eos_re, eos_ps, eos_ph, eos_th, eos_rh, eos_get_transport, eos_h, eos_deriv, eos_mui

contains

  subroutine eos_init(small_temp, small_dens)

    !use extern_probin_module
    !use parallel
    use extern_probin_module, only: eos_gamma
    use iso_c_binding, only : c_double, c_size_t

    implicit none

    real(amrex_real), optional :: small_temp
    real(amrex_real), optional :: small_dens

    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       call bl_error("gamma_const cannot be < 0")
    end if

    elem_names(1) = "X"
    spec_names(1) = "X"

    mintemp     = 1.d-200
    maxtemp     = 1.d200
    mindens     = 1.d-200
    maxdens     = 1.d200
    minmassfrac = 1.d-200
    maxmassfrac = 1.d0 + 1.d-12
    mine        = 1.d-200
    maxe        = 1.d200
    minp        = 1.d-200
    maxp        = 1.d200
    mins        = 1.d-200
    maxs        = 1.d200
    minh        = 1.d-200
    maxh        = 1.d200

  end subroutine eos_init

  subroutine eos_wb(state)

    use network, only: molec_wt
    implicit none

    type (eos_t), intent(inout) :: state

    state % wbar = 1.d0 / sum(state % massfrac(:) / molec_wt(:))

  end subroutine eos_wb

  subroutine eos_top(state)

    use amrex_constants_module
    use network, only: molec_wt

    implicit none

    type (eos_t), intent(inout) :: state

    ! Calculate wbar
    call eos_wb(state)
    state % gam1 = gamma_const
    state % cv = R / (state % wbar * (state % gam1 - ONE))
    state % cp = state % gam1 * state % cv
    state % cpi = state % cp

  end subroutine eos_top

  subroutine eos_bottom(state)

    use amrex_constants_module

    implicit none

    type (eos_t), intent(inout) :: state

    state % cs = sqrt(state % gam1 * state % p / state % rho)
    state % dpdr_e = state % p / state % rho
    state % dpde = (state % gam1 - ONE) * state % rho
    state % dPdr = ZERO

    state % hi = state % cv * state % T * state % gam1
    state % ei = state % e

  end subroutine eos_bottom

  subroutine eos_p_wb(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % p = (state % gam1 - 1.d0) * state % rho * state % e

  end subroutine eos_p_wb

  subroutine eos_rt(state)

    use amrex_constants_module

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)

    state % e = state % cv * state % T
    state % p = (state % gam1 - ONE) * state % rho * state % e
    state % s = ONE

    call eos_bottom(state)

  end subroutine eos_rt

  subroutine eos_rh(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_rh is not supported in this EOS.')

  end subroutine eos_rh

  subroutine eos_tp(state)

    use amrex_constants_module

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)

    state % rho = state % p * state % wbar / (R * state % T)
    state % e   = state % p / (state % rho * (state % gam1 - ONE))
    state % s   = ONE

    call eos_bottom(state)

  end subroutine eos_tp

  subroutine eos_rp(state)

    use amrex_constants_module

    implicit none

    type (eos_t), intent(inout) :: state

    real(amrex_real) :: poverrho

    call eos_top(state)

    poverrho = state % p / state % rho
    state % T = poverrho * state % wbar / R
    state % e = poverrho / (state % gam1 - ONE)
    state % s = ONE

    call eos_bottom(state)

  end subroutine eos_rp

  subroutine eos_re(state)

    use amrex_constants_module

    implicit none

    type (eos_t), intent(inout) :: state

    double precision :: poverrho

    call eos_top(state)

    poverrho = (state % gam1 - ONE) * state % e
    state % p = poverrho * state % rho
    state % T = poverrho * state % wbar / R
    state % s = ONE

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

    call bl_error('EOS: eos_ph is not supported in this EOS.')

  end subroutine eos_ph

  subroutine eos_th(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_th is not a valid input for the gamma law EOS.')

  end subroutine eos_th

  subroutine eos_xty(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % massfrac = state % molefrac

  end subroutine eos_xty

  subroutine eos_ytx_vec(Y, ylo, yhi, X, xlo, xhi, lo, hi, Nsp)

    implicit none

    integer, intent(in) :: ylo(3), yhi(3)
    integer, intent(in) :: xlo(3), xhi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: Nsp

    double precision, intent(in), dimension(ylo(1)-1:yhi(1)+1, ylo(2)-1:yhi(2)+1, ylo(3)-1:yhi(3)+1, 1:Nsp) :: Y
    double precision, intent(out), dimension(xlo(1)-1:xhi(1)+1, xlo(2)-1:xhi(2)+1, xlo(3)-1:xhi(3)+1, 1:Nsp) :: X

    integer :: i, j, k, n

    do n = 1, Nsp
      do k = lo(3)-1, hi(3) + 1
        do j = lo(2)-1, hi(2) + 1
          do i = lo(1)-1, hi(1) + 1
            X(i,j,k,n)=Y(i,j,k,n)
          enddo
        enddo
      enddo
    enddo

  end subroutine eos_ytx_vec

  subroutine eos_ytx(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % molefrac = state % massfrac

  end subroutine eos_ytx

  subroutine eos_cpi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)

  end subroutine eos_cpi

  subroutine eos_cv(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)

  end subroutine eos_cv

  subroutine eos_cp(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)

  end subroutine eos_cp

  subroutine eos_get_activity(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_get_activity is not supported in this EOS.')

  end subroutine eos_get_activity

  subroutine eos_get_transport(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call bl_error('EOS: eos_get_transport is not supported in this EOS.')

  end subroutine eos_get_transport

  subroutine eos_hi_vec(mass, masslo, masshi, T, Tlo, Thi, hi, hilo, hihi, low, high, Nsp)

    use amrex_constants_module
    use network, only: molec_wt

    implicit none

    integer, intent(in) :: masslo(3), masshi(3)
    integer, intent(in) :: Tlo(3), Thi(3)
    integer, intent(in) :: hilo(3), hihi(3)
    integer, intent(in) :: low(3), high(3)
    integer, intent(in) :: Nsp

    double precision, intent(in), dimension(masslo(1)-1:masshi(1)+1, masslo(2)-1:masshi(2)+1, masslo(3)-1:masshi(3)+1, 1:Nsp) :: mass
    double precision, intent(in), dimension(Tlo(1)-1:Thi(1)+1, Tlo(2)-1:Thi(2)+1, Tlo(3)-1:Thi(3)+1 ) :: T
    double precision, intent(out), dimension(hilo(1)-1:hihi(1)+1, hilo(2)-1:hihi(2)+1, hilo(3)-1:hihi(3)+1, 1:Nsp) :: hi

    integer :: i, j, k, n

    do n = 1, Nsp
       do k = low(3)-1, high(3) + 1
          do j = low(2)-1, high(2) + 1
             do i = low(1)-1, high(1) + 1
                hi(i,j,k,n) = R / ((1.d0 / sum(mass(i,j,k,:) / molec_wt(:))) * (gamma_const - ONE)) * T(i,j,k) * gamma_const
             enddo
          enddo
       enddo
    enddo

  end subroutine eos_hi_vec

  subroutine eos_hi(state)

    implicit none

    type (eos_t), intent(inout) :: state

    call eos_top(state)
    call eos_bottom(state)    

  end subroutine eos_hi

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
