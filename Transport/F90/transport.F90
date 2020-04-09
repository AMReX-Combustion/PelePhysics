module transport_module

  use eos_type_module
  use transport_type_module
  use actual_transport_module

  implicit none

contains

  ! Transport initialization routine: read in general transport parameters, then 
  ! call any specific initialization used by the transport

  ! This subroutine should be called outside OMP PARALLEL
  subroutine transport_init() bind(C, name="transport_init")

    use extern_probin_module

    implicit none

    ! Set up any specific parameters or initialization steps required by the transport we are using.
    call actual_transport_init()

  end subroutine transport_init


  subroutine transport_close()

    use extern_probin_module

    implicit none

    ! Clean up any specific parameters or initialization steps required by the transport we are using.
    call actual_transport_close

  end subroutine transport_close


  subroutine transport(which, coeff)

    !$acc routine seq

    implicit none

    ! Input arguments

    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff

    ! TODO: May want to check/process input state
    
    ! Call the transport evaluator
    call actual_transport(which, coeff)

  end subroutine transport

  subroutine get_transport_coeffs( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       D,            D_lo,  D_hi, &
       mu,          mu_lo, mu_hi, &
       xi,          xi_lo, xi_hi, &
       lam,        lam_lo,lam_hi) &
       bind(C, name="get_transport_coeffs")

    use fuego_chemistry, only: nspecies, naux
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    integer         , intent(in   ) ::   D_lo(3),  D_hi(3)
    integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
    integer         , intent(in   ) ::  xi_lo(3), xi_hi(3)
    integer         , intent(in   ) :: lam_lo(3),lam_hi(3)
    real (amrex_real), intent(in   ) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspecies)
    real (amrex_real), intent(in   ) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real (amrex_real), intent(in   ) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real (amrex_real), intent(inout) :: D(D_lo(1):D_hi(1),D_lo(2):D_hi(2),D_lo(3):D_hi(3),nspecies)
    real (amrex_real), intent(inout) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
    real (amrex_real), intent(inout) :: xi(xi_lo(1):xi_hi(1),xi_lo(2):xi_hi(2),xi_lo(3):xi_hi(3))
    real (amrex_real), intent(inout) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3))

    ! local variables
    integer      :: i, j, k, n, np
    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    np = hi(1)-lo(1)+1
    call build(coeff,np)

    which_trans % wtr_get_xi    = .true.
    which_trans % wtr_get_mu    = .true.
    which_trans % wtr_get_lam   = .true.
    which_trans % wtr_get_Ddiag = .true.

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)

          do i=1,np
             coeff%eos_state(i)%massfrac(1:nspecies) = massfrac(lo(1)+i-1,j,k,1:nspecies)
          end do

          do i=lo(1),hi(1)
             coeff%eos_state(i-lo(1)+1)%T = temperature(i,j,k)
             coeff%eos_state(i-lo(1)+1)%rho = density(i,j,k)
          end do

          call transport(which_trans, coeff)

          do i=lo(1),hi(1)
             mu(i,j,k)  = coeff %  mu(i-lo(1)+1)
             xi(i,j,k)  = coeff %  xi(i-lo(1)+1)
             lam(i,j,k) = coeff % lam(i-lo(1)+1)
          end do
          do n=1,nspecies
             do i=lo(1),hi(1)
                D(i,j,k,n) = coeff % Ddiag(i-lo(1)+1,n)
             end do
          end do

       end do
    end do

    call destroy(coeff)

  end subroutine get_transport_coeffs

  subroutine get_visco_coeffs( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       mu,          mu_lo, mu_hi) &
       bind(C, name="get_visco_coeffs")

    use fuego_chemistry, only: nspecies
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::  mu_lo(3), mu_hi(3)
    real (amrex_real), intent(in   ) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspecies)
    real (amrex_real), intent(in   ) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real (amrex_real), intent(inout) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))

    ! local variables
    integer      :: i, j, k, n, np
    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    np = hi(1)-lo(1)+1
    call build(coeff,np)

    which_trans % wtr_get_xi    = .false.
    which_trans % wtr_get_mu    = .true.
    which_trans % wtr_get_lam   = .false.
    which_trans % wtr_get_Ddiag = .false.

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)

          do i=1,np
             coeff%eos_state(i)%massfrac(1:nspecies) = massfrac(lo(1)+i-1,j,k,1:nspecies)
          end do

          do i=lo(1),hi(1)
             coeff%eos_state(i-lo(1)+1)%T = temperature(i,j,k)
          end do

          call transport(which_trans, coeff)

          do i=lo(1),hi(1)
             mu(i,j,k)  = coeff %  mu(i-lo(1)+1)
          end do

       end do
    end do

    call destroy(coeff)

  end subroutine get_visco_coeffs


  ! For now, we will assume that auxiliary variables diffuse with unity Le (if they diffuse at all)
  ! Therefore, we look up only Lambda, calculate Cp from the EoS, and then rhoDaux = lambda/Cp.
  ! This could be customized later by adding a which_trans option for auxiliaries and adding
  ! the desired transport proprty calculation to each actual_transport file. For now, this is the least
  ! obtrusive was to add diffusion for auxiliaries.
  subroutine get_transport_coeffs_aux( &
       lo,hi, &
       massfrac,    mf_lo, mf_hi, &
       temperature,  t_lo,  t_hi, &
       density,      r_lo,  r_hi, &
       Daux,       Dau_lo,Dau_hi) &
       bind(C, name="get_transport_coeffs_aux")

    use fuego_chemistry, only: nspecies, naux
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  mf_lo(3), mf_hi(3)
    integer         , intent(in   ) ::   t_lo(3),  t_hi(3)
    integer         , intent(in   ) ::   r_lo(3),  r_hi(3)
    integer         , intent(in   ) :: Dau_lo(3),Dau_hi(3)
    real (amrex_real), intent(in   ) :: massfrac(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3),nspecies)
    real (amrex_real), intent(in   ) :: temperature(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real (amrex_real), intent(in   ) :: density(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real (amrex_real), intent(inout) :: Daux(Dau_lo(1):Dau_hi(1),Dau_lo(2):Dau_hi(2),Dau_lo(3):Dau_hi(3),naux)

    ! local variables
    integer      :: i, j, k, n, np
    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    np = hi(1)-lo(1)+1
    call build(coeff,np)

    which_trans % wtr_get_xi    = .false.
    which_trans % wtr_get_mu    = .false.
    which_trans % wtr_get_lam   = .true.
    which_trans % wtr_get_Ddiag = .false.

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)

          do i=1,np
             coeff%eos_state(i)%massfrac(1:nspecies) = massfrac(lo(1)+i-1,j,k,1:nspecies)
          end do

          do i=lo(1),hi(1)
             coeff%eos_state(i-lo(1)+1)%T = temperature(i,j,k)
             coeff%eos_state(i-lo(1)+1)%rho = density(i,j,k)
             call eos_cp(coeff%eos_state(i-lo(1)+1))
          end do

          call transport(which_trans, coeff)

          do n=1,naux
             do i=lo(1),hi(1)
                Daux(i,j,k,n) = coeff % lam(i-lo(1)+1) / coeff % eos_state(i-lo(1)+1) % cp
             end do
          end do

       end do
    end do

    call destroy(coeff)

  end subroutine get_transport_coeffs_aux


  
end module transport_module
