!!$  The following is an implementation of the Peng-Ronbinson
!!$  Cubic equation of state for mixture of species

!!$  Relies on Species massfractions, Chemkin thermodynamic
!!$  correlations and curve fits. Additionally needs critical parameters
!!$  Tc, Pc, Zc, Vc and accentric factor omega
!!$
!!$  Internal energy, enthalpy, Cp, Cv are calculated using departure functions
!!$  from ideal gas values

module eos_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use eos_type_module
  use network, only : nspecies
  use chemistry_module, only : Ru, inv_mwt, chemistry_init, chemistry_initialized, spec_names, elem_names

  implicit none

  character (len=64) :: eos_name = "PREOS"

  logical, save, private :: initialized = .false.

  real(amrex_real), save, public :: smallT = 1.d-50

  private :: nspecies, Ru, inv_mwt

  real(amrex_real), dimension(:), allocatable :: Pc, Tc, Vc,Zc, omega,Bi, Ai, Fomega

  real(amrex_real), dimension(:,:), allocatable :: AMij

  real(amrex_real), parameter :: f0 = 0.379642
  real(amrex_real), parameter :: f1 = 1.48503
  real(amrex_real), parameter :: f2 = -0.164423
  real(amrex_real), parameter :: f3 = 0.016666

  real(amrex_real), parameter :: afactor = 0.457236
  real(amrex_real), parameter :: bfactor = 0.077796

contains

subroutine actual_eos_init

  implicit none
  integer :: ix, iy

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

 initialized = .true.

 allocate(Pc(nspecies), Tc(nspecies), omega(nspecies), Vc(nspecies),Zc(nspecies),Bi(nspecies),Ai(nspecies),Fomega(nspecies))

 allocate(AMij(nspecies,nspecies))

!!$ call READCRITICALPARAMETERS

! Precompute Bi Ai (without inclusion of the accentric factor based correction)
!!$ call preComputeAiBi
 
end subroutine actual_eos_init

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

  ! Set up any specific parameters or initialization steps required by the EOS we are using.

  call actual_eos_init

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

subroutine eos_xty(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Remains unchanged from Ideal EOS to Peng-Robinson EOS
 call ckxty (state % molefrac,state % massfrac)

 call PR_EOS_GetSpeedOfSound(state)
    
end subroutine eos_xty

subroutine eos_ytx(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Remains unchanged from Ideal EOS to Peng-Robinson EOS
 call ckytx (state % massfrac,state % molefrac)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_ytx


subroutine eos_cpi(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call PR_EOS_GetSpeciesCp(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_cpi

subroutine eos_hi(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate species enthalpy accounting for non-ideal effects
 call PR_Eos_GetSpeciesH(state)

 call PR_EOS_GetSpeedOfSound(state)
   
end subroutine eos_hi

subroutine eos_cv(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate mixture Cv accounting for non-ideal effects
 call PR_EOS_GetMixtureCv(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_cv

subroutine eos_cp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate mixture Cp accounting for non-ideal effects
 call PR_EOS_GetMixtureCp(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_cp

subroutine eos_mui(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & calculate species chemical potential accounting for non-ideal effects
 call PR_EOS_GetSpecies_ChemicalPotential(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_mui

subroutine eos_rt(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! (rho, T, massfrac) are inputs, get (p, e)
 call PR_EOS_GetP_GivenRhoT(state)
 call PR_EOS_GetE_GivenRhoT(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_rt

subroutine eos_rp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! (rho, p, massfrac) are inputs, get (T, e)
 call PR_EOS_Get_TE_GivenRhoP(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_rp

subroutine eos_tp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! (temp, press, massfrac) are inputs, get (rho, e)
 call PR_EOS_Get_rhoE_GivenTP(state)

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_tp

subroutine eos_re(state)

 implicit none

 type (eos_t), intent(inout) :: state

 integer :: lierr

 ! (rho, e, massfrac) are inputs, get (T, p)
 call PR_EOS_Get_TP_GivenRhoE(state,lierr)

 if (lierr .ne. 0) then
    print *, 'PR EOS: get_T_given_e Y failed, T, e, Y = ', &
         state % T, state % e, state % massfrac
 end if

 call PR_EOS_GetSpeedOfSound(state)

end subroutine eos_re

subroutine eos_ps(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! (p, s, massfrac) are inputs unsupported
 call bl_error('EOS: eos_ps is not supported in this EOS.')

end subroutine eos_ps

subroutine eos_ph(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call bl_error('EOS: eos_input_ph is not supported in this EOS.')

end subroutine eos_ph

subroutine eos_th(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call bl_error('EOS: eos_th is not supported in this EOS.')

end subroutine eos_th

subroutine eos_rh(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call bl_error('EOS: eos_input_rh is not supported in this EOS.')

end subroutine eos_rh


!===================================================!
! Compute the EOS species attractive and repulsive  !
! EOS parameters - dependent on critical properties !
!===================================================!
subroutine preComputeAiBi
  implicit none
  integer :: ix
  real(amrex_real) :: rTm1, rTm2, invPc

  do ix = 1,nspecies

     rTm1 = Ru * Tc(ix) * inv_mwt(ix)

     rTm2 = rTm1* Ru * Tc(ix) * inv_mwt(ix)

     invPc = 1.0/Pc(ix)
     
     Ai(ix) = afactor*rTm2*invPc

     Bi(ix) = bfactor*rTm1*invPc
     
  end do
  
end subroutine preComputeAiBi

!=======================================!
! Calculate the Accentric factor based  !
!   Temperature correction              !
!=======================================!
subroutine Calc_Fomega(T)
  implicit none
  real(amrex_real),intent(in) :: T
  integer :: ix
  real(amrex_real) :: OmegaFunc, O1, O2, O3
  
  do ix = 1,nspecies

     O1 = omega(ix)
     O2 = O1 * omega(ix)
     O3 = O2 * omega(ix)
     Fomega(ix) = f0 + f1*O1 + f2*O2 + f3*O3          

!!$     Smoothing of Fomega - Giovangigli proposes tanh function
!!$     if( OmegaFunc*(1.0-sqrt(Tr)) .ge. 0.0 ) then
!!$        Fomega(ix) = OmegaFunc*(1.0-sqrt(Tr))
!!$     else
!!$        Fomega(ix) = tanh(OmegaFunc*(1.0-sqrt(Tr))
!!$     end if
  end do

end subroutine Calc_Fomega

!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleAmBm(T,moleFrac,am,bm)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: moleFrac(:)
  real(amrex_real), intent(out) :: am,bm
  integer :: ix, iy
  real(amrex_real) :: O4, Tr, alphaTr

  call Calc_Fomega(T)
  
  do iy = 1, nspecies
     do ix = 1, nspecies

        Amij(ix,iy) = 1.0 

        ! Reduced Temperature for species, i
        Tr = T/Tc(ix)

        ! Temperature and accentric factor based correction
        O4 = (1.0 + Fomega(ix)*(1.0-sqrt(Tr)))
        alphaTr = O4 * O4

        ! AMij
        Amij(ix,iy) = Amij(ix,iy) * sqrt(Ai(ix) * alphaTr)

        ! Reduced Temperature for species, j
        Tr = T/Tc(iy)

        ! Temperature and accentric factor based correction
        O4 = (1.0 + Fomega(iy)*(1.0-sqrt(Tr)))
        alphaTr = O4 * O4

        ! AMij
        Amij(ix,iy) = Amij(ix,iy) * sqrt(Ai(iy) * alphaTr)
        
     end do
  end do

  ! Molefraction based mixing rule to compute am
  am = 0.0
  do iy = 1, nspecies
     do ix = 1, nspecies
        am = am + moleFrac(ix)*moleFrac(iy)*Amij(ix,iy)
     end do
  end do

  ! Molefraction based mixing rule to compute bm
  bm = 0.0
  do ix = 1, nspecies
     bm = bm + moleFrac(ix)*Bi(ix)
  end do


end subroutine MixingRuleAmBm

!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleAm(T,moleFrac,am)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: moleFrac(:)
  real(amrex_real), intent(out) :: am
  integer :: ix, iy
  real(amrex_real) :: O4,Tr, alphaTr

  call Calc_Fomega(T)

  do iy = 1, nspecies
     do ix = 1, nspecies

        Amij(ix,iy) = 1.0 

        ! Reduced Temperature for species, i
        Tr = T/Tc(ix)

        ! Temperature and accentric factor based correction
        O4 = (1.0 + Fomega(ix)*(1.0-sqrt(Tr)))
        alphaTr = O4 * O4

        ! AMij
        Amij(ix,iy) = Amij(ix,iy) * sqrt(Ai(ix) * alphaTr)

        ! Reduced Temperature for species, j
        Tr = T/Tc(iy)

        ! Temperature and accentric factor based correction
        O4 = (1.0 + Fomega(iy)*(1.0-sqrt(Tr)))
        alphaTr = O4 * O4

        ! AMij
        Amij(ix,iy) = Amij(ix,iy) * sqrt(Ai(iy) * alphaTr)
        
     end do
  end do

  ! Molefraction based mixing rule to compute am
  am = 0.0
  do iy = 1, nspecies
     do ix = 1, nspecies
        am = am + moleFrac(ix)*moleFrac(iy)*Amij(ix,iy)
     end do
  end do

end subroutine MixingRuleAm

!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleBm(moleFrac,bm)
  implicit none
  real(amrex_real),intent(in) :: moleFrac(:)
  real(amrex_real), intent(out) :: bm
  integer :: ix

  ! Molefraction based mixing rule to compute bm
  bm = 0.0
  do ix = 1, nspecies
     bm = bm + moleFrac(ix)*Bi(ix)
  end do

end subroutine MixingRuleBm

!==================================!
!  Compute dAm/dT for PR EOS      !
!==================================!
subroutine Calc_dAmdT(T,moleFrac,am,dAmdT)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: moleFrac(:)
  real(amrex_real), intent(in) :: am
  real(amrex_real), intent(out) :: dAmdT
  real(amrex_real) :: Tr_i, Tr_j, alphaTr_i,alphaTr_j
  real(amrex_real) :: A_j, A_i
  integer :: ix, iy
  
  dAmdT = 0.0

  do iy = 1, nspecies
     do ix = 1, nspecies

        ! Reduced temperature 
        Tr_i = T/Tc(ix)
        Tr_j = T/Tc(iy)
        
        alphaTr_i = (1.0 + Fomega(ix)*(1.0-sqrt(Tr_i)))
        alphaTr_j = (1.0 + Fomega(iy)*(1.0-sqrt(Tr_j)))

        ! Multiplying the reduced temperature based correction factor
        A_j = Ai(iy)*alphaTr_j*alphaTr_j
        A_i = Ai(ix)*alphaTr_i*alphaTr_i
        
        dAmdT =  dAmdT - moleFrac(ix)*moleFrac(iy)*(sqrt(Ai(ix)*A_j)*Fomega(ix)/(2.0*sqrt(T * Tc(ix))) + &
                                                    sqrt(Ai(iy)*A_i)*Fomega(iy)/(2.0*sqrt(T * Tc(iy))) )

     end do
  end do

end subroutine Calc_dAmdT

!==================================!
!  Compute d2Am/dT2 for PR EOS      !
!==================================!
subroutine Calc_d2AmdT2(T,moleFrac,am,dAmdT,d2AmdT2)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: moleFrac(:)
  real(amrex_real), intent(in) :: am
  real(amrex_real), intent(in) :: dAmdT
  real(amrex_real), intent(in) :: d2AmdT2
  real(amrex_real) :: Tr_i, Tr_j, alphaTr_i,alphaTr_j
  real(amrex_real) :: A_j, A_i
  integer :: ix, iy

  d2AmdT2 = 0.0
  do iy = 1, nspecies
     do ix = 1, nspecies
        d2AmdT2 =  d2AmdT2 + moleFrac(ix)*moleFrac(iy)*(sqrt(Ai(ix)*Ai(iy))*Fomega(ix)*Fomega(iy))/sqrt(Tci(ix) * Tc(iy))
     end do
  end do

  d2AmdT2 = (d2AmdT2 - dAmdT)/(2.0*T)
  
end subroutine Calc_d2AmdT2

!=============================================!
! Calculate the roots of the cubic equation   !
!  to compute compressibility factor          !
!=============================================!
subroutine Calc_CompressFactor_Z(Z,am,bm,P,T,Wbar)
  implicit none
  real(amrex_real),intent(out) :: Z
  real(amrex_real),intent(in)  :: am,bm,P,T,Wbar
  real(amrex_real) :: alpha,beta,gamma,Q,R
  real(amrex_real) :: RmT, B1, B2, B3
  real(amrex_real) :: R1, R2, R3

  RmT = (Ru/Wbar)*T

  B1 = (bm*P/RmT)
  B2 = B1 * (bm*P/RmT)
  B3 = B2 * (bm*P/RmT)

  R1 = RmT
  R2 = R1*RmT
  R3 = R2*RmT

  alpha = -1.0 + B1
  
  beta = (am*P - 3.0*(bm*P)*(bm*P))/R2 - 2.0*B1

  gamma = -(am*bm*P*P)/(RmT*RmT*RmT) + B2 + B3
  
  Q = (alpha*alpha - 3.0*beta)/9.0
  
  R = (2.0*alpha*alpha*alpha - 9.0*alpha*beta + 27.0*gamma)/54.0
  
  ! Multiple roots of the cubic equation
  if((Q*Q*Q - R*R)  .gt. 0.0) then
     theta = acos(R/(Q**1.5))
     
     Z1 = -2.0*sqrt(Q)*cos(theta/3.0) - alpha/3.0
     
     Z2 = -2.0*sqrt(Q)*cos((theta+2.0*Pi)/3.0) - alpha/3.0
     
     Z3 = -2.0*sqrt(Q)*cos((theta+4.0*Pi)/3.0) - alpha/3.0
     
     Z = max(Z1,Z2,Z3)
  else
     Z = -sign(1.0,R)*( (sqrt(R*R - Q*Q*Q) + abs(R) )**(1.0/3.0) + Q/((sqrt(R*R - Q*Q*Q) + abs(R) )**(1.0/3.0)) ) - alpha/3.0
  end if
  
end subroutine Calc_CompressFactor_Z

!==========================================================!
! Given a mixture composition calculate mixture density &  !
! internal energy given Pressure and Temperature           !
!      using Peng-Robinson EOS                             !
!==========================================================!
subroutine PR_EOS_Get_rhoE_givenTP(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau,K1
  real(amrex_real) :: Eig
  real(amrex_real) :: sqrtOf2
  
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  ! Calculate the roots of cubic EOS to calculate Compressibility factor
  call Calc_CompressFactor_Z(state%Z,state%am,state%bm,state%P,Tn,state%wbar)
  
  state%rho = state%P * state%Wbar/(state%Z*Ru*state%T)

  ! Specific volume
  tau = 1.0/state%rho

  ! Calculate the first derivative of am w.r.t Temperature,T
  call Calc_dAmdT(state%T,state%moleFrac,state%am,state%dAmdT)

  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

  ! Ideal gas internal energy
  call ckums(state%T, Eig)

  ! Add departure function to the ideal gas mixture internal energy 
  state%e = Eig + (state%am - state%T*state%dAmdT)*K1

end subroutine PR_EOS_Get_rhoE_givenTP

!==========================================================!
! Given a mixture composition calculate Pressure           !
!      given density and Temperature                       !
!      using Peng-Robinson EOS                             !
!==========================================================!
subroutine PR_EOS_GetP_GivenRhoT(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau,Eig 
  real(amrex_real) :: K1
  real(amrex_real) :: sqrtOf2
  
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  tau = 1.0/state%rho

  ! Calculate pressure using Rho, T and Peng-Robinson EOS
  state%p = (Ru/state%wbar) * state%T/(tau - state%bm) - state%am/(tau*tau + 2.0*tau*state%bm - state%bm*state%bm)

end subroutine PR_EOS_GetP_GivenRhoT

!==========================================================!
! Given a mixture composition calculate internal energy    !
!    given density and Temperature                         !
!      using Peng-Robinson EOS                             !
!==========================================================!
subroutine PR_EOS_GetE_GivenRhoT(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau,Eig 
  real(amrex_real) :: K1
  real(amrex_real) :: sqrtOf2
  
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  tau = 1.0/state%rho

  ! Calculate the first derivative of am w.r.t Temperature,T
  call Calc_dAmdT(state%T,state%moleFrac,state%am,state%dAmdT)

  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

  ! Ideal gas internal energy
  call ckums(state%T, Eig)

  ! Add departure function to the ideal gas mixture internal energy 
  state%e = Eig + (state%am - state%T*state%dAmdT)*K1

end subroutine PR_EOS_GetE_GivenRhoT

!==========================================================!
!       Generate mixture Peng Robinson EOS
!  calculate T and internal energy given rho & P
!==========================================================!
subroutine PR_EOS_Get_TE_givenRhoP(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: K1
  real(amrex_real) :: dfdT,dT
  real(amrex_real) :: Eig, tau
  real(amrex_real) :: Tnp1, Tn, Pnp1
  real(amrex_real), parameter :: convCrit = 1e-4
  real(amrex_real), parameter :: accelFac = 1e-5
  integer, parameter :: maxIter = 250
  integer :: nIter
  real(amrex_real) :: sqrtOf2
  sqrtOf2 = sqrt(2.0)
  
  ! Calculate the mixture averaged BM
  call MixingRuleBm(state%moleFrac,state%bm)

  ! Specific volume
  tau = 1.0/state%rho

  ! Ideal gas value as the first guess 
  Tnp1 = state%p*state%wbar/(state%rho*Ru)

  nIter = 0

  ! Start the NewtonSolver Iteration 
  do while (abs(state%p-Pnp1).gt.convCrit .and. nIter .lt. maxIter)

     Tn = Tnp1
     nIter = nIter + 1

     call MixingRuleAm(Tn,state%moleFrac,state%am)

     ! Calculate the roots of cubic EOS to calculate Compressibility factor
     call Calc_CompressFactor_Z(state%Z,state%am,state%bm,state%P,Tn,state%wbar)

     ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_dAmdT(Tn,state%moleFrac,state%am,state%dAmdT)
     
     ! Derivative of the non-linear equation
     dfdT = (Ru/state%wbar)/(tau - Bm) - state%dAmdT/(tau*tau + 2.0*tau*state%bm - state%bm*state%bm)

     ! Temperature correction for the next guess 
     dT = -(state%Z*state%rho*Ru*Tn/state%wbar - state%P)/dfdT

     ! Update the next iterate of Temperature 
     Tnp1 = Tn + dT/(1.0 + abs(dT*accelFac))

     ! Use pressure as the iteration convergence criterion
     Pnp1 = state%rho * state%Z * Ru *Tnp1 / state%wbar
     
  end do

  ! Update temperature in the state 
  state%T = Tnp1

  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

  ! Compute ideal gas internal energy 
  call ckums(state % T, Eig)
  
  ! Update the real gas internal energy using departure functions
  state%e = Eig + (state%am - state%T*state%dAmdT)*K1
  
end subroutine PR_EOS_Get_TE_GivenRhoP

!==========================================================!
!       Generate mixture Peng Robinson EOS
!        calculate T and P given rho & e
!==========================================================!
subroutine PR_EOS_Get_TP_GivenRhoE(state,lierr)
  implicit none
  type (eos_t), intent(in) :: state
  integer, intent(out) :: lierr
  integer :: ix, iy
  real(amrex_real) :: Wmix
  real(amrex_real) :: alphaTr,OmegaFunc, Tr, Pr
  real(amrex_real) :: rutfac
  real(amrex_real) :: K1,Vol
  real(amrex_real) :: dT, O1,O2,O3
  real(amrex_real) :: Eig
  real(amrex_real), parameter :: convCrit = 1e-4
  integer, parameter :: maxIter = 250
  integer :: nIter
  real(amrex_real) :: tau,Tn,Tnp1,fzero
  real(amrex_real) :: sqrtOf2
  sqrtOf2 = sqrt(2.0)

  ! Mixture averaged EOS repulsive parameter 
  call MixingRuleBm(state%moleFrac,state%bm)

  tau = 1.0/state%rho

  ! Use ideal gas as the first guess
  Tnp1 = state%P*state%Wbar/(state%rho*Ru)

  nIter = 0
  fzero = 1.0
  lierr = 0
  
  do while (abs(fzero).gt.convCrit .and. nIter .lt. maxIter)
     
     Tn = Tnp1

     ! Update the iteration counter
     nIter = nIter + 1

     call MixingRuleAm(Tn,state%moleFrac,state%am)

     ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_dAmdT(Tn,state%moleFrac,state%am,state%dAmdT)
     
     ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_d2AmdT2(Tn,state%moleFrac,state%am,state%dAmdT,state%d2AmdT2)

     K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

     ! Ideal gas internal energy and specific heat at constant volume, Cv
     call ckums(Tn, Eig)
     call ckcvbs(Tn, state % massfrac,state % cv)

     ! Calculate real gas Cv
     state%cv = state%cv - Tn*state%d2AmdT2*K1

     ! Compute the non-linear equation
     fzero = state%e - Eig - (state%am - Tn*state%dAmdT)*K1
     dT = fzero/state%cv

     ! Update the temperature
     Tnp1 = Tn + dT

  end do

  if(abs(fzero).gt.convCrit .and. nIter.eq.maxIter) then
     lierr = 3
  end if

  ! Update the state structure with the updated Temperature and Pressure
  state%T = Tnp1
  state%p = (Ru/state%wbar)*state%T/(tau-state%bm) - state%am/(tau*tau + 2.0*tau*state%bm - state%bm*state%bm)

end subroutine PR_EOS_Get_TP_GivenRhoE

!========================================================!
! Given a mixture composition calculate species Cp using !
!                Peng-Robinson EOS                       !
!========================================================!
subroutine PR_EOS_GetSpeciesCp(state)
  implicit none
  type (eos_t), intent(in) :: state




end subroutine PR_EOS_GetSpeciesCp

!=================================================================!
! Given a mixture composition calculate species enthalpy, h using !
!                Peng-Robinson EOS                                !
!=================================================================!
subroutine PR_EOS_GetSpeciesH(state)
  implicit none
  type (eos_t), intent(in) :: state

  
end subroutine PR_EOS_GetSpeciesH

!========================================================================!
! Given a mixture composition calculate species chemical potential, mu_k !
!           using Peng-Robinson EOS                                      !
!========================================================================!
subroutine PR_EOS_GetSpecies_ChemicalPotential(state)
  implicit none
  type (eos_t), intent(in) :: state

  
end subroutine PR_EOS_GetSpecies_ChemicalPotential

!======================================================================!
! Given a mixture composition calculate mixture specific heat Cv using !
!                Peng-Robinson EOS                                     !
!======================================================================!
subroutine PR_EOS_GetMixtureCv(state)
  implicit none
  type (eos_t), intent(in) :: state
  real(amrex_real) :: sqrtOf2
  real(amrex_real) :: tau, K1
  
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  tau = 1.0/state%rho

  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%moleFrac,state%am,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%moleFrac,state%am,state%dAmdT,state%d2AmdT2)
  
  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

  ! Ideal gas specific heat at constant volume
  call ckcvbs(state%T, state % massfrac, state % cv)

  ! Real gas specific heat at constant volume
  state%cv = state%cv - state%T*state%d2AmdT2*K1
  
end subroutine PR_EOS_GetMixtureCv

!======================================================================!
! Given a mixture composition calculate mixture specific heat Cp using !
!                Peng-Robinson EOS                                     !
!======================================================================!
subroutine PR_EOS_GetMixtureCp(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1
  real(amrex_real) :: CviG
  real(amrex_real) :: eosT1Denom, eosT2Denom
  real(amrex_real) :: Rm
  real(amrex_real) :: sqrtOf2
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  tau = 1.0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%moleFrac,state%am,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%moleFrac,state%am,state%dAmdT,state%d2AmdT2)
  
  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )

  eosT1Denom = tau-state%bm
  eosT2Denom = tau*tau + 2.0*tau*state%bm - state%bm*state%bm

  Rm = (Ru/state%wbar)
  
  ! Derivative of Pressure w.r.t to Temperature
  state%dPdT = Rm/eosT1Denom - state%dAmdT/eosT2Denom

  ! Derivative of Pressure w.r.t to tau (specific volume)
  state%dpdtau = -Rm*state%T/(eosT1Denom*eosT1Denom) + 2.0*state%am*(tau+state%bm)/(eosT2Denom*eosT2Denom) 

  ! Ideal gas specific heat at constant volume
  call ckcvbs(state%T, state % massfrac, CviG)

  ! Real gas specific heat at constant pressure
  state%cp = CviG - state%T*state%d2AmdT2*K1 - state%T*state%dPdT*state%dPdT/state%dpdtau
  
end subroutine PR_EOS_GetMixtureCp

!=========================================!
! Calculate speed of sound using PR EOS   !
!=========================================!
subroutine PR_EOS_GetSpeedOfSound(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1
  real(amrex_real) :: CviG
  real(amrex_real) :: eosT1Denom, eosT2Denom
  real(amrex_real) :: Rm
  real(amrex_real) :: sqrtOf2
  sqrtOf2 = sqrt(2.0)

  call MixingRuleAmBm(state%T,state%moleFrac,state%am,state%bm)

  tau = 1.0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%moleFrac,state%am,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%moleFrac,state%am,state%dAmdT,state%d2AmdT2)

  K1 = (1.0/(2.0*sqrtOf2*state%bm))*log( (tau + (1.0-sqrtOf2)*state%bm)/(tau + (1.0+sqrtOf2)*state%bm ) )
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*tau + 2.0*tau*state%bm - state%bm*state%bm

  Rm = (Ru/state%wbar)
  
  ! Derivative of Pressure w.r.t to Temperature
  state%dPdT = Rm/eosT1Denom - state%dAmdT/eosT2Denom
  
  ! Derivative of Pressure w.r.t to tau (specific volume)
  state%dpdtau = -Rm*state%T/(eosT1Denom*eosT1Denom) + 2.0*state%am*(tau+state%bm)/(eosT2Denom*eosT2Denom) 
  
  ! Ideal gas specific heat at constant volume
  call ckcvbs(state%T, state % massfrac, CviG)
  
  ! Real gas specific heat at constant pressure
  state%cp = CviG - state%T*state%d2AmdT2*K1 - state%T*state%dPdT*state%dPdT/state%dpdtau

   ! Real gas specific heat at constant volume
  state%cv = state%cv - state%T*state%d2AmdT2*K1

  ! Real gas speed of sound
  state%cs = tau*sqrt(-state%cp*state%dpdtau/state%cv)

end subroutine PR_EOS_GetSpeedOfSound

end module eos_module
