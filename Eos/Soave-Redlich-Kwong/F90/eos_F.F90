!!$  The following is an implementation of the SRK
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
  use fuego_chemistry
  use eos_type_module
  ! use chemistry_module, only : nspecies, Ru, inv_mwt, chemistry_init, chemistry_initialized, spec_names, elem_names, molecular_weight
  ! With restructured chemistry above all comes from fuego_chemistry module

  implicit none

  character (len=64) :: eos_name = "SRKEOS"

  logical, save, private :: initialized = .false.

  real(amrex_real), save, public :: smallT = 1.d-50

  public :: eos_init, eos_xty, eos_ytx, eos_cpi, eos_hi, eos_cv, eos_cp, eos_p_wb, eos_wb, eos_get_activity, eos_rt, eos_tp, eos_rp, eos_re, eos_ps, eos_ph, eos_th, eos_rh, eos_get_transport, eos_h, eos_deriv, eos_mui
  private :: nspecies, Ru, inv_mwt

  !    real(amrex_real), dimension(:), allocatable :: Pc, Tc, Vc,Zc, omega,Bi, Astari, Fomega
  real(amrex_real), dimension(:), allocatable :: oneOverTc, Tc, Bi, Astari, sqrtAsti, omega, Fomega

  real(amrex_real), parameter :: f0 = 0.48508d0
  real(amrex_real), parameter :: f1 = 1.5517d0
  real(amrex_real), parameter :: f2 = -0.151613d0
 
  !  not needed
  ! real(amrex_real), parameter :: afactor = 0.42748d0
  ! real(amrex_real), parameter :: bfactor = 0.08664d0

  interface
     subroutine amrex_array_init_snan (p, nelem) bind(C,name="amrex_array_init_snan")
       use iso_c_binding, only : c_double, c_size_t
       implicit none
       real(c_double),intent(inout) :: p
       integer (kind=c_size_t),intent(in),value :: nelem
     end subroutine amrex_array_init_snan
  end interface

contains

subroutine actual_eos_init

  implicit none
  integer :: i,j,k

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

!allocate(Pc(nspecies), Tc(nspecies), omega(nspecies), Vc(nspecies),Zc(nspecies),Bi(nspecies),Astari(nspecies),Fomega(nspecies))

 allocate( oneOverTc(nspecies), Tc(nspecies), omega(nspecies), Bi(nspecies), Astari(nspecies), sqrtAsti(nspecies), Fomega(nspecies))

 call get_critparams(Tc,Astari,Bi,omega)
  
  do i = 1,nspecies

     oneOverTc(i) = 1.0d0/Tc(i)

     Fomega(i) = f0 + omega(i)*(f1 + f2*omega(i))          

     sqrtAsti(i) = sqrt(Astari(i))

  end do
  
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

 ! Remains unchanged from Ideal EOS to SRK EOS
 call ckxty (state % molefrac,state % massfrac)
    
end subroutine eos_xty

subroutine eos_ytx(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Remains unchanged from Ideal EOS to SRK EOS
 call ckytx (state % massfrac,state % molefrac)

end subroutine eos_ytx

subroutine eos_cpi(state)

 implicit none

 type (eos_t), intent(inout) :: state

 !  for the momennt we can use the ideal ones .  
 !  these are only used to esimate internal degrees of freedom 
 !  for computing bulk viscosity

 call ckcpms(state % T, state % cpi)

 !  call bl_error('EOS: eos_cpi is not supported in this EOS.')
 ! Construct the EOS for mixture & Calculate species Cp accounting for non-ideal effects
 !    call SRK_Eos_GetSpeciesCp(state)

end subroutine eos_cpi

subroutine eos_hi(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate species enthalpy accounting for non-ideal effects
 call SRK_Eos_GetSpeciesH(state)
   
end subroutine eos_hi

subroutine eos_cv(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate mixture Cv accounting for non-ideal effects
 call SRK_Eos_GetMixtureCv(state)

end subroutine eos_cv

subroutine eos_cp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & Calculate mixture Cp accounting for non-ideal effects
 call SRK_Eos_GetMixtureCp(state)

end subroutine eos_cp

subroutine eos_mui(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Used to test chemical potential expressions
 call SRK_EOS_GetMixture_FreeEnergy(state)

 ! Construct the EOS for mixture & calculate species chemical potential accounting for non-ideal effects
 call SRK_EOS_GetSpecies_ChemicalPotential(state)
    
end subroutine eos_mui

subroutine eos_h(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for mixture & calculate mixture enthalpy accounting for non-ideal effects
 call SRK_EOS_GetMixture_H(state)

end subroutine eos_h

subroutine eos_p_wb(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call eos_top(state)

 ! Calculate am & bm using non-Ideal EOS
 call MixingRuleAmBm(state%T, state%massFrac, state%am, state%bm)

 ! Calculate Pressure given rho and T 
 call SRK_EOS_GetP_givenRhoT(state)
    
end subroutine eos_p_wb

subroutine eos_wb(state)

 implicit none

 type (eos_t), intent(inout) :: state

 state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
    
end subroutine eos_wb

subroutine eos_get_activity(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call SRK_EOS_GetSpecies_Activity(state)

end subroutine eos_get_activity

subroutine eos_get_activity_h(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call SRK_EOS_GetSpecies_Activity(state)

end subroutine eos_get_activity_h

subroutine eos_get_transport(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call SRK_EOS_Calc_Deriv_Chem_Potential(state)

end subroutine eos_get_transport

subroutine eos_deriv(state)

 implicit none

 type (eos_t), intent(inout) :: state

 ! Construct the EOS for the mixture and calculate all the derivatives
 call SRK_EOS_GetDerivative(state)

 ! Calculates dtau/dYk at constant P and T 
 call SRK_EOS_GetSpeciesTau(state)
    
end subroutine eos_deriv

subroutine eos_top(state)

 implicit none

 type (eos_t), intent(inout) :: state

 state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

end subroutine eos_top

subroutine eos_bottom(state)

 implicit none

 type (eos_t), intent(inout) :: state
 real(amrex_real) :: dedtau,tau

 call SRK_Eos_GetMixtureCv(state)
 call SRK_Eos_GetMixtureCp(state)
 call SRK_EOS_GetGamma1(state)

 tau = 1.0d0/state%rho
    
 state % cs = sqrt(state % gam1 * state % p / state % rho) ! cm/s
 state % dpdr = -tau**2*state%dpdtau
 state % dpde = state%dpdt / state%cv

 dedtau = ( 1.0d0/(tau*(tau+state%bm)))*(state%am - state%T * state%dAmdT)
 state%dpdr_e = -tau*tau*(state%dpdtau- dedtau*state%dpdT/state%cv)

end subroutine eos_bottom

subroutine eos_rt(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call eos_top(state)

 ! (rho, T, massfrac) are inputs, get (p, e)
 !  pulled this out so no hidden side effects

 call MixingRuleAmBm(state%T, state%massFrac, state%am, state%bm)
 call Calc_dAmdT(state%T, state%massFrac, state%damdT)
 call SRK_EOS_GetP_GivenRhoT(state)
 call SRK_EOS_GetE_GivenRhoT(state)

 ! Populate the eo_state%ei vector
 call SRK_EOS_GetSpeciesE(state)

 call eos_bottom(state)

end subroutine eos_rt

subroutine eos_tp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call eos_top(state)

 ! (temp, press, massfrac) are inputs, get (rho, e)
 call SRK_EOS_Get_rho_GivenTP(state)
 call SRK_EOS_Get_E_GivenTP(state)

 ! Populate the eo_state%ei vector
 call SRK_EOS_GetSpeciesE(state)

 call eos_bottom(state)

end subroutine eos_tp

subroutine eos_rp(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call eos_top(state)

 ! (rho, p, massfrac) are inputs, get (T, e)
 call SRK_EOS_Get_TE_GivenRhoP(state)

 ! Populate the eo_state%ei vector
 call SRK_EOS_GetSpeciesE(state)

 call eos_bottom(state)

end subroutine eos_rp

subroutine eos_re(state)

 implicit none

 type (eos_t), intent(inout) :: state
 real(amrex_real) :: dedtau,tau,K1
 real(amrex_real) :: Cpig
 real(amrex_real) :: eosT1Denom,eosT2Denom,eosT3Denom 
 real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
 real(amrex_real) :: Rm,Rmk,dhmdtau,dhmdT
 real(amrex_real) :: Temp1,Temp2
 real(amrex_real) :: ek(nspecies)
 integer :: lierr
 integer :: i,j

 !begin eos_top(state)
 state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
 !end eos_top(state)

 call SRK_EOS_Get_T_GivenRhoE(state,lierr)
 !if (lierr .ne. 0) then
 !   print *, 'SRK EOS: get_T_given_e Y failed, T, e, Y = ', &
 !        state % T, state % e, state % massfrac
 !end if

 call MixingRuleAmBm(state%T, state%massFrac, state%am, state%bm)
 call Calc_dAmdT(state%T, state%massFrac, state%damdT)
 tau = 1.d0/state%rho
 state%p = (Ru/state%wbar) * state%T/(tau - state%bm) - &
             state%am/(tau*(tau + state%bm))
 call Calc_dAmdY(state%T, state%massFrac,state%dAmdYk)
 call Calc_d2AmdTY(state%T, state%massFrac,state%d2amdYkdT)
 K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
 eosT1Denom = tau-state%bm
 eosT2Denom = tau*(tau+state%bm)
 eosT3Denom = tau+state%bm
 InvEosT1Denom = 1.0d0/eosT1Denom
 InvEosT2Denom = 1.0d0/eosT2Denom
 InvEosT3Denom = 1.0d0/eosT3Denom
 Rm = (Ru/state%wbar)
 call ckums(state % T, ek)
 do i = 1,nspecies
    Rmk = Ru*inv_mwt(i)
    Temp1 = (state%T * state%d2amdYkdT(i) - state%dAmdYk(i))*K1
    Temp2 = (state%T * state%dAmdT - state%am)*(-Bi(i)*K1 + InvEosT3Denom*Bi(i))/state%bm
    state%ei(i) = ek(i) + Temp1 + Temp2
 end do
 ! Non-ideal Cv is already calculated in SRK_EOS_Get_T_GivenRhoE
 ! state%cv = state%cv + state%T*state%d2AmdT2* (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
 state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom
 state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom
 call ckcpbs(state % T, state % massfrac, Cpig)
 dhmdT = Cpig + state%T*state%d2AmdT2*K1 - state%dAmdT*InvEosT3Denom + Rm*state%bm*InvEosT1Denom
 dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
           Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom
 state%cp = dhmdT - (dhmdtau/state%dpdtau)*state%dPdT
 state%gam1 = -tau*state%cp*state%dpdtau/(state%p*state%cv)
 state % cs = sqrt(state % gam1 * state % p / state % rho) ! cm/s
 state % dpdr = -tau**2*state%dpdtau
 state % dpde = state%dpdt / state%cv
 dedtau = ( 1.0d0/(tau*(tau+state%bm)))*(state%am - state%T * state%dAmdT)
 state%dpdr_e = -tau*tau*(state%dpdtau- dedtau*state%dpdT/state%cv)
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

 call bl_error('EOS: eos_ph is not supported in this EOS.')

end subroutine eos_ph

subroutine eos_th(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call bl_error('EOS: eos_th is not supported in this EOS.')

end subroutine eos_th

subroutine eos_rh(state)

 implicit none

 type (eos_t), intent(inout) :: state

 call bl_error('EOS: eos_rh is not supported in this EOS.')

end subroutine eos_rh

!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleAmBm(T,massFrac,am,bm)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real), intent(out) :: am,bm
  integer :: i, j
  real(amrex_real) ::  Tr 

  real(amrex_real) :: amloc(nspecies)

  am = 0.d0
  bm = 0.d0

  do i = 1, nspecies

  ! may want to later add molifier here for am

     Tr = T*oneOverTc(i)
     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

     bm = bm + massFrac(i)*Bi(i)

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        am = am + massFrac(i)*massFrac(j)*amloc(i)*amloc(j)

     end do
  end do
 
end subroutine MixingRuleAmBm
!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleAm(T,massFrac,am)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real), intent(out) :: am
  integer :: i, j
  real(amrex_real) :: Tr
  real(amrex_real) :: amloc(nspecies)

  am = 0.d0

  do i = 1, nspecies

  ! may want to later add molifier here for am

     Tr = T*oneOverTc(i)
     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        am = am + massFrac(i)*massFrac(j)*amloc(i)*amloc(j)

     end do
  end do

end subroutine MixingRuleAm
!==================================!
!  Compute Am, Bm for the mixture  !
!==================================!
subroutine MixingRuleBm(massFrac,bm)
  implicit none
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real), intent(out) :: bm
  integer :: i

  ! Massfraction based mixing rule to compute bm
  bm = 0.0d0
  do i = 1, nspecies
     bm = bm + massFrac(i)*Bi(i)
  end do

end subroutine MixingRuleBm
!==================================!
!  Compute dAm/dT for SRK EOS      !
!==================================!
subroutine Calc_dAmdT(T,massFrac,dAmdT)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real),intent(out) :: dAmdT
  real(amrex_real) :: Tr, oneOverT
  integer :: i,j

  real(amrex_real) :: amloc(nspecies), amlocder(nspecies)
  
  dAmdT = 0.0d0
  oneOverT = 1.0d0/T

  do i = 1, nspecies

  ! may want to later add molifier here for am

     ! Reduced temperature 
     Tr = T*oneOverTc(i)
     amloc(i)    = (1.0d0+Fomega(i)*(1.0d0-sqrt(Tr)))*sqrtAsti(i)
     amlocder(i) = -0.5d0*Fomega(i)*sqrtAsti(i)*oneOverT*oneOverTc(i)*sqrt(T*Tc(i))

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        dAmdT =  dAmdT +massFrac(i)*massFrac(j)*(amloc(i)*amlocder(j)+amloc(j)*amlocder(i))

     end do
  end do

end subroutine Calc_dAmdT
!==================================!
!  Compute dAm/dY for SRK EOS      !
!==================================!
subroutine Calc_dAmdY(T,massFrac,dAmdY)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real), intent(out) :: dAmdY(:)
  real(amrex_real) :: Tr
  integer :: i,j

  real(amrex_real) :: amloc(nspecies)
  
  dAmdY = 0.0d0

  do i = 1, nspecies

  ! may want to later add molifier here for am

     ! Reduced temperature 
     Tr = T*oneOverTc(i)

     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        dAmdY(i) =  dAmdY(i) + 2.d0*massFrac(j)*amloc(i)*amloc(j)

     end do
  end do

end subroutine Calc_dAmdY
!==================================!
!  Compute dAm/dY for SRK EOS      !
!==================================!
subroutine Calc_d2AmdY2(T,massFrac,d2AmdY2)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real), intent(out) :: d2AmdY2(:,:)
  real(amrex_real) :: Tr
  integer :: i,j
  real(amrex_real) :: amloc(nspecies)
  

  do i = 1, nspecies

  ! may want to later add molifier here for am

     ! Reduced temperature 
     Tr = T*oneOverTc(i)

     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        d2AmdY2(i,j) =  2.d0*amloc(i)*amloc(j)

     end do
  end do

end subroutine Calc_d2AmdY2
!==================================!
!  Compute dAm/dT for SRK EOS      !
!==================================!
subroutine Calc_d2AmdTY(T,massFrac,dAmdTY)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real),intent(out) :: dAmdTY(:)
  real(amrex_real) ::  Tr, oneOverT
  integer :: i,j

  real(amrex_real) :: amloc(nspecies), amlocder(nspecies)
  
  dAmdTY = 0.0d0
  oneOverT = 1.0d0/T

  do i = 1, nspecies

  ! may want to later add molifier here for am

     ! Reduced temperature 
     Tr = T*oneOverTc(i)
     amloc(i)    = (1.0d0+Fomega(i)*(1.0d0-sqrt(Tr)))*sqrtAsti(i)
     amlocder(i) = -0.5d0*Fomega(i)*sqrtAsti(i)*oneOverT*oneOverTc(i)*sqrt(T*Tc(i))

  enddo

  do j = 1, nspecies
     do i = 1, nspecies
        
        dAmdTY(i) =  dAmdTY(i) + 2.d0*massFrac(j)*(amloc(i)*amlocder(j)+amloc(j)*amlocder(i))

     end do
  end do

end subroutine Calc_d2AmdTY
!==================================!
!  Compute d2Am/dT2 for SRK EOS    !
!==================================!
subroutine Calc_d2AmdT2(T,massFrac,d2AmdT2)
  implicit none
  real(amrex_real),intent(in) :: T
  real(amrex_real),intent(in) :: massFrac(:)
  real(amrex_real),intent(out) :: d2AmdT2
  real(amrex_real) :: Tr, tmp1, oneOverT
  integer :: i, j
  real(amrex_real) :: amloc(nspecies), amlocder(nspecies)

  d2AmdT2 = 0.0d0
  oneOverT = 1.0d0/T
  tmp1 = -0.5d0*oneOverT

  do i = 1, nspecies
     Tr = T*oneOverTc(i)
     amloc(i)    = (1.0d0+Fomega(i)*(1.0d0-sqrt(Tr)))*sqrtAsti(i)
     amlocder(i) = -0.5d0*Fomega(i)*sqrtAsti(i)*oneOverT*oneOverTc(i)*sqrt(T*Tc(i))
  enddo

  do j = 1, nspecies-1
     do i = j+1, nspecies
        d2AmdT2 = d2AmdT2+tmp1*massFrac(i)*massFrac(j)*(-4.0d0*T*amlocder(i)*amlocder(j)+amloc(i)*amlocder(j)+amloc(j)*amlocder(i))
     end do
  end do

  d2AmdT2 = 2.0d0*d2AmdT2

  do i = 1, nspecies
     d2AmdT2 = d2AmdT2+tmp1*massFrac(i)*massFrac(i)*(-4.0d0*T*amlocder(i)*amlocder(i)+amloc(i)*amlocder(i)+amloc(i)*amlocder(i))
  end do

end subroutine Calc_d2AmdT2
!=============================================!
! Calculate the roots of the cubic equation   !
!  to compute compressibility factor          !
!=============================================!
subroutine Calc_CompressFactor_Z(Z,am,bm,P,T,Wbar)
  implicit none
  real(amrex_real),intent(out) :: Z
  real(amrex_real),intent(in)  :: am,bm,P,T,Wbar
  real(amrex_real) :: alpha,beta,gamma,Q,R, theta
  real(amrex_real) :: Z1, Z2, Z3
  real(amrex_real) :: RmT, B1
  real(amrex_real) :: R1, R2, R3
  real(amrex_real), parameter :: Pi = 4.0d0 * atan (1.0d0)

!  JBB not checked

  RmT = (Ru/Wbar)*T

  B1 = (bm*P/RmT)

  R1 = RmT
  R2 = R1*RmT
  R3 = R2*RmT

  alpha = -1.0d0
  
  beta = (am*P - bm*P*bm*P)/R2 - B1

  gamma = -(am*bm*P*P)/R3 
  
  Q = (alpha*alpha - 3.0d0*beta)/9.0d0
  
  R = (2.0d0*alpha*alpha*alpha - 9.0d0*alpha*beta + 27.0d0*gamma)/54.0d0
  
  ! Multiple roots of the cubic equation
  if((Q*Q*Q - R*R)  .gt. 0.0d0) then
     theta = acos(R/(Q**1.5d0))
     
     Z1 = -2.0d0*sqrt(Q)*cos(theta/3.0d0) - alpha/3.0d0
     
     Z2 = -2.0d0*sqrt(Q)*cos((theta+2.0d0*Pi)/3.0d0) - alpha/3.0d0
     
     Z3 = -2.0d0*sqrt(Q)*cos((theta+4.0d0*Pi)/3.0d0) - alpha/3.0d0
     
     Z = max(Z1,Z2,Z3)
  else
     Z = -sign(1.0d0,R)*( (sqrt(R*R - Q*Q*Q) + abs(R) )**(1.0d0/3.0d0) + Q/((sqrt(R*R - Q*Q*Q) + abs(R) )**(1.0d0/3.0d0)) ) - alpha/3.0d0
  end if
  
end subroutine Calc_CompressFactor_Z
!==========================================================!
! Given a mixture composition calculate mixture density    !
!         given Pressure and Temperature                   !
!           using SRK EOS                                  !
!==========================================================!
subroutine SRK_EOS_Get_rho_givenTP(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau,K1
  real(amrex_real) :: Eig

  ! not checked

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Calculate the roots of cubic EOS to calculate Compressibility factor
  call Calc_CompressFactor_Z(state%Z,state%am,state%bm,state%P,state%T,state%wbar)

  state%rho = state%P * state%Wbar/(state%Z*Ru*state%T)

  ! Specific volume
  tau = 1.0d0/state%rho

end subroutine SRK_EOS_Get_rho_givenTP
!================================================================!
! Given a mixture composition calculate mixture internal energy  !
!  given Pressure and Temperature using SRK EOS                  !
!================================================================!
subroutine SRK_EOS_Get_E_givenTP(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1
  real(amrex_real) :: Eig
  real(amrex_real) :: ei_ig(nspecies)

!not checked

  ! Specific volume
  tau = 1.0d0/state%rho

  ! Calculate the first derivative of am w.r.t Temperature,T
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  ! Ideal gas internal energy
  call ckums(state%T, ei_ig)
  Eig = sum(state % massfrac(:) * ei_ig(:))

  ! Add departure function to the ideal gas mixture internal energy 
  state%e = Eig + (state%T*state%dAmdT - state%am)*K1

end subroutine SRK_EOS_Get_E_givenTP
!==========================================================!
! Given a mixture composition calculate Pressure           !
!      given density and Temperature                       !
!      using SRK EOS                                       !
!==========================================================!
subroutine SRK_EOS_GetP_givenRhoT(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau

  tau = 1.d0/state%rho

  !  called outside
  !  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Calculate pressure using Rho, T and SRK EOS
  state%p = (Ru/state%wbar) * state%T/(tau - state%bm) - &
              state%am/(tau*(tau + state%bm))

end subroutine SRK_EOS_GetP_GivenRhoT
!==========================================================!
! Given a mixture composition calculate internal energy    !
!    given density and Temperature                         !
!      using SRK EOS                                       !
!==========================================================!
subroutine SRK_EOS_GetE_givenRhoT(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau,Eig 
  real(amrex_real) :: K1
  real(amrex_real) :: ei_ig(nspecies)

  !  mixing rule call outside
  !  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  tau = 1.0d0/state%rho

  ! Calculate the first derivative of am w.r.t Temperature,T
  !  call Calc_dAmdT(state%T,state%massFrac,state%am,state%dAmdT)

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  ! Ideal gas internal energy
  call ckums(state%T, ei_ig)
  Eig = sum(state % massfrac(:) * ei_ig(:))

  ! Add departure function to the ideal gas mixture internal energy
  state%e = Eig + (state%T*state%dAmdT - state%am)*K1

end subroutine SRK_EOS_GetE_GivenRhoT
!==========================================================!
!       Generate mixture SRK EOS
!  calculate T and internal energy given rho & P
!==========================================================!
subroutine SRK_EOS_Get_TE_givenRhoP(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: K1
  real(amrex_real) :: dfdT,dT
  real(amrex_real) :: Eig, tau
  real(amrex_real) :: Tnp1, Tn, Pnp1
  real(amrex_real), parameter :: convCrit = 1e-4
  integer, parameter :: maxIter = 250
  integer :: nIter
  real(amrex_real) :: eosT2Denom,InvEosT2Denom
  real(amrex_real) :: eosT1Denom,InvEosT1Denom
  real(amrex_real) :: Rm
  real(amrex_real) :: ei_ig(nspecies)
 
  !  not checked

  ! Calculate the mixture averaged BM
  call MixingRuleBm(state%massFrac,state%bm)

  ! Specific volume
  tau = 1.0d0/state%rho

  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT1Denom = 1.0d0/eosT1Denom

  Rm = (Ru/state%wbar)

  ! Ideal gas value as the first guess 
  Tnp1 = state%p*tau/Rm
  call MixingRuleAm(Tnp1,state%massFrac,state%am)
  Pnp1 = Rm * state%T*InvEosT1Denom - state%am*InvEosT2Denom

  nIter = 0

  ! Start the NewtonSolver Iteration 
  do while (abs(state%p-Pnp1).gt.convCrit .and. nIter .lt. maxIter)

     Tn = Tnp1
     nIter = nIter + 1

     !  JBB  WHY ISN'T THIS JUST NEWTWON 
     !  Tnp1 = Tn - (p(T,y,tau)-pgiven)/ dpdT

     call MixingRuleAm(Tn,state%massFrac,state%am)

     ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_dAmdT(Tn,state%massFrac,state%dAmdT)

     Pnp1 = Rm * Tn * InvEosT1Denom - state%am * InvEosT2Denom

     ! Derivative of Pressure w.r.t to Temperature
     state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom

     ! Update the next iterate of Temperature 
     Tnp1 = Tn - (Pnp1-state%p)/state%dPdT

     call MixingRuleAm(Tnp1,state%massFrac,state%am)

     ! Use pressure as the iteration convergence criterion
     Pnp1 = Rm * Tnp1 * InvEosT1Denom - state%am * InvEosT2Denom
     
  end do

  ! Update temperature in the state 
  state%T = Tnp1
  ! recompute dAmdT at updated temperature (we've already done this for Am)
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  ! Compute ideal gas internal energy 
  call ckums(state%T, ei_ig)
  Eig = sum(state % massfrac(:) * ei_ig(:))
  
  ! Update the real gas internal energy using departure functions
  state%e = Eig + (state%T*state%dAmdT-state%am)*K1
  
end subroutine SRK_EOS_Get_TE_GivenRhoP
!==========================================================!
!       Generate mixture SRK EOS
!        calculate T and P given rho & e
!==========================================================!
subroutine SRK_EOS_Get_T_GivenRhoE(state,lierr)
  implicit none
  type (eos_t), intent(inout) :: state
  integer, intent(out) :: lierr
  integer :: ix, iy
  real(amrex_real) :: Wmix
  real(amrex_real) :: alphaTr,OmegaFunc, Tr, Pr
  real(amrex_real) :: rutfac
  real(amrex_real) :: K1,Vol
  real(amrex_real) :: dT, O1,O2,O3
  real(amrex_real) :: Eig, tau
  real(amrex_real), parameter :: convCrit = 1d-4
  integer, parameter :: maxIter = 2000
  integer :: nIter
  real(amrex_real) :: Tn,Tnp1,fzero
  real(amrex_real) :: Rm
  real(amrex_real) :: ei_ig(nspecies)
  
  !  not checked

  ! Mixture averaged EOS repulsive parameter 
  call MixingRuleBm(state%massFrac,state%bm)

  ! Specific volume
  tau = 1.0d0/state%rho

  Rm = (Ru/state%wbar)
  
  ! Use ideal gas as the first guess
!  can't use p . . . not define
! Tnp1 = state%P*tau/Rm
   Tnp1 = state % T

  nIter = 0
  fzero = 1.0d0
  lierr = 0
  !write(6,*)" in eos_re", state%rho,state%wbar,
  
  do while (abs(fzero).gt.convCrit .and. nIter .lt. maxIter)
     
     Tn = Tnp1

     ! Update the iteration counter
     nIter = nIter + 1

     call MixingRuleAm(Tn,state%massFrac,state%am)

     ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_dAmdT(Tn,state%massFrac,state%dAmdT)
     
     ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
     call Calc_d2AmdT2(Tn,state%massFrac,state%d2AmdT2)

     K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

     ! Ideal gas internal energy and specific heat at constant volume, Cv
     call ckums(Tn, ei_ig)
     Eig = sum(state % massfrac(:) * ei_ig(:))
     
     call ckcvbs(Tn, state % massfrac, state % cv)

     ! Calculate real gas Cv
     state%cv = state%cv + Tn*state%d2AmdT2*K1

     ! Compute the non-linear equation
     fzero = -state%e + Eig + (Tn*state%dAmdT-state%am)*K1
     
     dT = fzero/state%cv
     ! Update the temperature
     Tnp1 = Tn - dT

  end do

  if(abs(fzero).gt.convCrit .and. nIter.eq.maxIter) then
     lierr = 3
  end if

  ! Update the state structure with the updated Temperature and Pressure
  state%T = Tnp1

end subroutine SRK_EOS_Get_T_GivenRhoE
!========================================================!
! Given a mixture composition calculate species Cp using !
!                SRK EOS                                 !
!========================================================!
subroutine SRK_EOS_GetSpeciesCp(state)
  implicit none
  type (eos_t), intent(in) :: state

        !  not sure if we need this  

end subroutine SRK_EOS_GetSpeciesCp
!=================================================================!
! Given a mixture composition calculate species enthalpy, h using !
!                SRK EOS                                          !
!=================================================================!
subroutine SRK_EOS_GetSpeciesH(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau
  real(amrex_real) :: dAmdTYk(nspecies),dhmdYk(nspecies)
  real(amrex_real) :: Temp1,Temp2, Temp3

  real(amrex_real) :: hmix

!  assumes T, rho and massFrac known

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Derivative of am w.r.t to Yk, species k massfraction

  call Calc_dAmdY(state%T, state%massFrac,state%dAmdYk)
  call Calc_d2AmdTY(state%T, state%massFrac,dAmdTYk)

  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)
 
  ! Derivative of Pressure w.r.t to tau (specific volume)
  state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

  ! Derivative of enthalpy w.r.t to tau (specific volume)
  dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
            Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom

  call ckhms(state % T, state % hi)
  
  ! Derivative of Pressure w.r.t to Yk, species k massfraction
  do i = 1,nspecies

     Rmk = Ru*inv_mwt(i)

     state%dPdYk(i) = Rmk*state%T*InvEosT1Denom &
                     - state%dAmdYk(i)*InvEosT2Denom       &
                     + Bi(i)*(Rm*state%T*InvEosT1Denom*InvEosT1Denom &
                              + state%am*InvEosT2Denom*InvEosT3Denom)

!    dhmdYk(i) = state%hi(i) + (state%T*state%d2amdYkdT(i)-state%dAmdYk(i))*K1  &
     dhmdYk(i) = state%hi(i) + (state%T*dAmdTYk(i)-state%dAmdYk(i))*K1  &
                       - Bi(i)*(state%T*state%dAmdT-state%am)*(K1/state%bm - InvEosT3Denom/state%bm) &
                       + state%am*Bi(i)*InvEosT3Denom*InvEosT3Denom   &
                       - InvEosT3Denom*state%damdYk(i)    &
                       + Rmk*state%T*state%bm*InvEosT1Denom &
                       + Rm*state%T*Bi(i)*(InvEosT1Denom + state%bm*InvEosT1Denom*InvEosT1Denom )

     state%hi(i) = dhmdYk(i) - (dhmdtau/state%dpdtau)*state%dPdYk(i)

  end do
  
end subroutine SRK_EOS_GetSpeciesH
!========================================================================!
! Given a mixture composition calculate species chemical potential, mu_k !
!           using SRK EOS                                                !
!========================================================================!
subroutine SRK_EOS_GetSpecies_ChemicalPotential(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: RmT,RmkT,logTauSt
  real(amrex_real), parameter :: Psat = 1013250.0d0

! assumes rho, T, massFrac known
  
  ! Ideal gas Species enthalpy 
  call ckhms(state % T, state % hi)

  ! Ideal gas Species entropy
  call cksms(state % T, state % si)

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  ! Generate Am, Bm for the mixture
  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Calculate the first derivative of Am w.r.t to Yk, species k massfraction
  call Calc_dAmdy(state%T,state%massFrac,state%damdYk)

  logTauSt =  log(Ru*state%T/Psat)

  tau = 1.0d0/state%rho
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom
  
  RmT = (Ru*state%T/state%wbar)

  ! Add departure function to Ideal gas chemical potential mu_id
  do i = 1,nspecies
     
     RmkT = Ru*inv_mwt(i)*state%T

     state%mui(i) = state%hi(i) - state%T * state%si(i) + RmkT*logTauSt

     state%mui(i) = state%mui(i) + RmT*Bi(i)*InvEosT1Denom &
                                 - K1*state%damdYk(i) + (state%am/state%bm)*K1*Bi(i)   &
                                 - state%am*InvEosT3Denom*Bi(i)/state%bm

     if(state%massFrac(i).gt. 1.0e-16) then
        state%mui(i) = state%mui(i) + RmkT*log(state%massFrac(i)*inv_mwt(i)*InvEosT1Denom)
     end if
        
  end do

end subroutine SRK_EOS_GetSpecies_ChemicalPotential
!========================================================================!
! Given a mixture composition calculate mixture free energy, f           !
!           using SRK EOS                                                !
!========================================================================!
subroutine SRK_EOS_GetMixture_FreeEnergy(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, InvEosT1Denom
  real(amrex_real) :: RmT,RmkT
  real(amrex_real) :: sumTerm1, sumTerm2,invPsat
  real(amrex_real), parameter :: Psat = 1013250.0d0

  ! assumes rho, T, massFrac known
  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  ! Generate Am, Bm for the mixture
  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)
  
  ! Ideal gas Species Internal energy 
  call ckums(state % T, state % ei)

  ! Ideal gas Species entropy
  call cksms(state % T, state % si)

  tau = 1.0d0/state%rho
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  InvEosT1Denom = 1.0d0/eosT1Denom
  
  RmT = (Ru*state%T/state%wbar)

  invPsat = 1.0d0/Psat
  
  sumTerm1 = 0.0d0

  sumTerm2 = 0.0d0
  
  do i = 1,nspecies

     if(state%massFrac(i).gt. 1.0e-16) then 

        RmkT = Ru*inv_mwt(i)*state%T
        
        sumTerm1 = sumTerm1 + state%massFrac(i)*(state%ei(i) - state%T * state%si(i))
        
        sumTerm2 = sumTerm2 + state%massFrac(i)*RmkT* log(state%massFrac(i) * RmkT * InvEosT1Denom * invPsat)

     end if
     
  end do
  
  state%f = sumTerm1 + sumTerm2 - state%am*K1 
  
end subroutine SRK_EOS_GetMixture_FreeEnergy
!=====================================================================!
! Given a mixture composition calculate species activity coefficients
!              using SRK EOS                                          !
!=====================================================================!
subroutine SRK_EOS_GetSpecies_Activity(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau,dhmdYk

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)
  
  call Calc_dAmdY(state%T,state%massFrac,state%damdYk)
  
  tau = 1.0d0/state%rho

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)

! write(6,*)" in activity", state%rho, state%am,state%bm, tau

  do i = 1,nspecies

     state%Acti(i) = exp( (molecular_weight(i)/(Ru*state%T))*( Rm*state%T*InvEosT1Denom*Bi(i) &
                         - K1*state%damdYk(i) + (state%am/state%bm)*K1*Bi(i)   &
                         - state%am*InvEosT3Denom*Bi(i)/state%bm ) ) * state%massFrac(i)*inv_mwt(i)*InvEosT1Denom
!    write(6,*)" in loop", i, state%massFrac(i)*inv_mwt(i)*InvEosT1Denom
  end do

end subroutine SRK_EOS_GetSpecies_Activity
!======================================================================!
! Given a mixture composition calculate d mu_k /dtau and d mu_j dYk    !
!                SRK EOS                                               !
!====================================================================== !
subroutine SRK_EOS_Calc_Deriv_Chem_Potential(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j,k
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm
  real(amrex_real) :: d2AmdY2(nspecies,nspecies)
  real(amrex_real) :: dpdtau0, dpdy0(nspecies)

!  should already be there
    
! state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
  
! call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)
  
! call Calc_dAmdY(state%T,state%massFrac,state%damdYk)

  call Calc_d2AmdY2(state%T,state%massFrac,d2AmdY2)
  
  tau = 1.0d0/state%rho

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)

  dpdtau0 = - Rm * state%T *state%rho **2

  do k= 1, nspecies

        dpdy0(k) = Ru*state%T *state%rho * inv_mwt(k)
  enddo

!  should be there
  ! Derivative of Pressure w.r.t to tau (specific volume)
! state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

  state%dijY = 0.d0

  do k = 1,nspecies

     state%diP(k) = state%massfrac(k) * state%wbar* ((state%dAmdYk(k) - state%am *Bi(k) /state%bm)*InvEosT2Denom  &
               +  state%am *Bi(k) /state%bm *InvEosT1Denom**2 ) / (Rm * state%T) &
               - state%massfrac(k) * (Bi(k)*InvEosT1Denom**2 + state%wbar*inv_mwt(k)*InvEosT1Denom)

    state%diP(k) = state%diP(k) / state%dpdtau


!    state%diP(k) = -state%wbar*state%massfrac(k)*inv_mwt(k)*state%rho

!    state%diP(k) = state%diP(k) / state%dpdtau

!    state%diP(k) = state%diP(k) / dpdtau0

     state%dijY(k,k) =  state%wbar*inv_mwt(k)
     
     
     do j = 1,nspecies


        state%dijY(k,j) =  state%dijY(k,j) -  state%diP(k)* state%dPdYk(j)   &
              + state%wbar*state%massfrac(k)*InvEosT1Denom *( Bi(k)*inv_mwt(j) + Bi(j)*inv_mwt(k))  &
              + state%massfrac(k)*InvEosT1Denom**2*Bi(k)*Bi(j)  & 
              + state%massfrac(k)/(Rm * state%T) * ( (K1 - InvEosT3Denom)*state%damdYk(j)*Bi(k)/state%bm &
              - K1*d2AmdY2(k,j) +(K1*state%damdYk(k) - InvEosT3Denom)*Bi(k)/state%bm  &
              + (-2.0d0*state%am*K1 + state%am*InvEosT3Denom)*Bi(j)*Bi(k)/(state%bm*state%bm) &
              + (state%am*InvEosT3Denom/state%bm + state%am* InvEosT3Denom*InvEosT3Denom)*Bi(j)*Bi(k)/state%bm )

!       state%dijY(k,j) =  state%dijY(k,j) -  state%diP(k)* state%dPdYk(j)

!       state%dijY(k,j) =  state%dijY(k,j) -  state%diP(k)* state%dPdYk(j)

     end do
     
  end do

!  implicit none
!  type (eos_t), intent(inout) :: state
!  integer :: i,j,k
!  real(amrex_real) :: tau, K1
!  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom
!  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
!  real(amrex_real) :: Rm,Rmk, RmJ
!  real(amrex_real) :: d2AmdY2(nspecies,nspecies)
!
!  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
!
!  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)
!
!  call Calc_dAmdY(state%T,state%massFrac,state%am,state%damdYk)
!
!  call Calc_d2AmdY2(state%T,state%massFrac,state%am,d2AmdY2)
!
!  tau = 1.0d0/state%rho
!
!  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
!
!  eosT1Denom = tau-state%bm
!  eosT2Denom = tau*(tau+state%bm)
!  eosT3Denom = tau+state%bm
!
!  InvEosT1Denom = 1.0d0/eosT1Denom
!  InvEosT2Denom = 1.0d0/eosT2Denom
!  InvEosT3Denom = 1.0d0/eosT3Denom
!
!  Rm = (Ru/state%wbar)*state%T
!
!  ! Derivative of Pressure w.r.t to tau (specific volume)
!  state%dpdtau = -Rm*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom
!
!
!  do j = 1,nspecies
!
!     RmJ = Ru * state%T * inv_mwt(j)
!
!
!     state%diP(j) = RmJ*InvEosT1Denom - state%dAmdYk(j)*InvEosT2Denom  &
!                      + Bi(j)*(Rm*InvEosT1Denom*InvEosT1Denom            &
!                      + state%am*InvEosT2Denom*InvEosT3Denom)
!
!
!     state%dijY(j) = -Rm*Bi(j)*InvEosT1Denom*InvEosT1Denom &
!                         -RmJ*InvEosT1Denom + (state%damdYk(j) - state%am*Bi(j)/state%bm)*InvEosT2Denom &
!                         + state%am*Bi(j)*InvEosT1Denom*InvEosT1Denom/state%bm
!
!     ! Multiplying by the diagonal matrix of Massfractions
!     state%dijY(j) = state%dijY(j)*state%massFrac(j)
!
!
!
!     do k = 1,nspecies
!
!        RmK = Ru * state%T * inv_mwt(k)
!
!        state%dijY(k,j) =  RmK * Kdelta(k,j) + state%massFrac(k)*( &
!                                                     InvEosT1Denom*Ru*state%T*(inv_mwt(k)*Bi(j)+ inv_mwt(j)*Bi(k)) &
!                                                   + Rm*InvEosT1Denom*InvEosT1Denom*Bi(k)*Bi(j) &
!                                                   - K1*d2AmdY2(k,j) + (K1 - InvEosT3Denom)*state%damdYk(j)*Bi(k)/state%bm &
!                                                   + (K1-InvEosT3Denom)*state%damdYk(k)*Bi(k)/state%bm  &
!                                                   + (-2.0d0*state%am*K1 + state%am*InvEosT3Denom)*Bi(j)*Bi(k)/(state%bm*state%bm) &
!                                                   + (state%am*InvEosT3Denom/state%bm + state%am* InvEosT3Denom*InvEosT3Denom)*Bi(j)*Bi(k)/state%bm )
!     end do
!
!  end do

end subroutine SRK_EOS_Calc_Deriv_Chem_Potential
!=====================================================================!
! Given a mixture composition calculate mixture Cv                    !
!              using SRK EOS                                          !
!=====================================================================!
subroutine SRK_EOS_GetMixtureCv(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  tau = 1.0d0/state%rho

  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  ! call Calc_dAmdT(state%T,state%massFrac,state%am,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)

  ! Ideal gas specific heat at constant volume
  call ckcvbs(state%T, state % massfrac, state % cv)

  ! Real gas specific heat at constant volume
  state%cv = state%cv + state%T*state%d2AmdT2* (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
end subroutine SRK_EOS_GetMixtureCv
!======================================================================!
! Given a mixture composition calculate mixture specific heat Cp using !
!                SRK EOS                                               !  
!======================================================================!
subroutine SRK_EOS_GetMixtureCp(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1
  real(amrex_real) :: Cpig
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: dhmdT,dhmdtau
  real(amrex_real) :: Rm

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm

  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)
  
  ! Derivative of Pressure w.r.t to Temperature
  state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom

  ! Derivative of Pressure w.r.t to tau (specific volume)
  state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

  ! Ideal gas specific heat at constant pressure
  call ckcpbs(state % T, state % massfrac, Cpig)
  
  ! Derivative of enthalpy w.r.t to Temperature
  dhmdT = Cpig + state%T*state%d2AmdT2*K1 - state%dAmdT*InvEosT3Denom + Rm*state%bm*InvEosT1Denom
  
  ! Derivative of enthalpy w.r.t to tau (specific volume)
  dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
            Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom

  ! Real gas specific heat at constant pressure
  state%cp = dhmdT - (dhmdtau/state%dpdtau)*state%dPdT
  
end subroutine SRK_EOS_GetMixtureCp
!=======================================================================!
! Given a mixture composition calculate mixture enthalpy using SRK EOS  !  
!=======================================================================!
subroutine SRK_EOS_GetMixture_H(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT3Denom
  real(amrex_real) :: Rm
  real(amrex_real) :: hmix

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  ! Mixing rule
  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Specific volume
  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)

  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  eosT1Denom = tau-state%bm
  eosT3Denom = tau+state%bm

  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  ! Specific gas constant
  Rm = (Ru/state%wbar)

  ! Ideal gas mixture enthalpy 
  call ckhms(state % T, state % hi)
  state%h = sum(state % massfrac(:) * state % hi(:))
  call ckhbms(state % T, state%massfrac(:), hmix)

  ! Adding non-ideal departure function 
  state%h = state%h + (state%T*state%dAmdT - state%am)*K1 &
                    + Rm * state%T * state%bm * InvEosT1Denom   &
                    - state%am * InvEosT3Denom
  
end subroutine SRK_EOS_GetMixture_H
!=========================================!
! Calculate speed of sound using PR EOS   !
!=========================================!
subroutine SRK_EOS_GetGamma1(state)
  implicit none
  type (eos_t), intent(inout) :: state
  real(amrex_real) :: tau !, K1
  !real(amrex_real) :: Cpig
  !real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  !real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  !real(amrex_real) :: dhmdT,dhmdtau
  !real(amrex_real) :: Rm

  !  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  ! call Calc_dAmdT(state%T,state%massFrac,state%am,state%dAmdT)
  
  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  !call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)
  !
  !K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)

  !eosT1Denom = tau-state%bm
  !eosT2Denom = tau*(tau+state%bm)
  !eosT3Denom = tau+state%bm

  !InvEosT1Denom = 1.0d0/eosT1Denom
  !InvEosT2Denom = 1.0d0/eosT2Denom
  !InvEosT3Denom = 1.0d0/eosT3Denom

  !Rm = (Ru/state%wbar)

  ! JBB  all of these are known at call
  
  ! Derivative of Pressure w.r.t to Temperature
  ! state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom

  ! Derivative of Pressure w.r.t to tau (specific volume)
  ! state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

  ! Ideal gas specific heat at constant pressure
  ! call ckcpbs(state % T, state % massfrac, Cpig)
  
  ! Derivative of enthalpy w.r.t to Temperature
  ! dhmdT = Cpig + state%T*state%d2AmdT2*K1 - state%dAmdT*InvEosT3Denom + Rm*state%bm*InvEosT1Denom
  
  ! Derivative of enthalpy w.r.t to tau (specific volume)
  ! dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
  !          Rm*T*state%bm*InvEosT1Denom*InvEosT1Denom

  ! Real gas specific heat at constant pressure
  ! state%cp = dhmdT - (dhmdtau/state%dpdtau)*state%dPdT

  ! Real gas specific heat at constant volume
  ! state%cv = Cpig - R + state%T*state%d2AmdT2*K1

  ! Real gas speed of sound
  state%gam1 = -tau*state%cp*state%dpdtau/(state%p*state%cv)

end subroutine SRK_EOS_GetGamma1
!=========================================================================!
! Given a mixture composition calculate species specific volume,tau using !
!                SRK EOS                                                  !
!=========================================================================!
subroutine SRK_EOS_GetSpeciesTau(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau
  real(amrex_real) :: dAmdTYk(nspecies),dhmdYk(nspecies)
  real(amrex_real) :: Temp1,Temp2, Temp3

!  assumes T, rho and massFrac known

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Derivative of am w.r.t to Yk, species k massfraction
  call Calc_dAmdY(state%T, state%massFrac,state%dAmdYk)

  tau = 1.0d0/state%rho
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)

  do i = 1,nspecies

     Rmk = Ru*inv_mwt(i)
     
     Temp1 = Rm*state%T*InvEosT1Denom*InvEosT1Denom - state%am*(2.0d0*tau + state%bm)*InvEosT2Denom*InvEosT2Denom

     Temp2 = Rmk*state%T*InvEosT1Denom + Rm*state%T*Bi(i)*InvEosT1Denom*InvEosT1Denom &
             -state%dAmdYk(i)*InvEosT2Denom + tau*state%am*Bi(i)*InvEosT2Denom*InvEosT2Denom

     state%taui(i) = Temp2/Temp1

  end do

end subroutine SRK_EOS_GetSpeciesTau
!=================================================================!
! Given a mixture composition calculate species enthalpy, h using !
!                SRK EOS                                          !
!=================================================================!
subroutine SRK_EOS_GetSpeciesH_V2(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau
  real(amrex_real) :: dAmdTYk(nspecies)
  real(amrex_real) :: Temp1, Temp2, Temp3, Temp4 
  real(amrex_real) :: hk(nspecies)

!  assumes T, rho and massFrac known

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Derivative of am w.r.t to Yk, species k massfraction
  call Calc_dAmdY(state%T, state%massFrac,state%dAmdYk)
  
  call Calc_d2AmdTY(state%T, state%massFrac,dAmdTYk)

  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)
 
  call ckhms(state % T, hk)
  
  ! Derivative of Pressure w.r.t to Yk, species k massfraction
  do i = 1,nspecies

     Rmk = Ru*inv_mwt(i)

     Temp1 = (state%T * dAmdTYk(i) - state%dAmdYk(i))*K1

     Temp2 = (state%T * state%dAmdT - state%am)*(-Bi(i)*K1/state%bm + InvEosT3Denom*(Bi(i) - state%bm * state%taui(i)/tau)/state%bm)

     Temp3 = Rmk*state%T*state%bm*InvEosT1Denom + Rm*state%T*( Bi(i)*InvEosT1Denom &
          - state%bm*(state%taui(i)-Bi(i))*InvEosT1Denom*InvEosT1Denom )

     Temp4 = -state%dAmdYk(i) * InvEosT3Denom + state%am*(state%taui(i) + Bi(i))* InvEosT3Denom*InvEosT3Denom
     
     state%hi(i) = hk(i) + Temp1 + Temp2 + Temp3 + Temp4 

  end do
  
end subroutine SRK_EOS_GetSpeciesH_V2
!=========================================================================!
! Given a mixture composition calculate species internal energy, ei using !
!                SRK EOS                                                  !
!=========================================================================!
subroutine SRK_EOS_GetSpeciesE(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau
  real(amrex_real) :: Temp1, Temp2, Temp3, Temp4 
  real(amrex_real) :: ek(nspecies)

!  assumes T, rho and massFrac known

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Derivative of am w.r.t to Yk, species k massfraction
  call Calc_dAmdY(state%T, state%massFrac,state%dAmdYk)
  
  call Calc_d2AmdTY(state%T, state%massFrac,state%d2amdYkdT)

  tau = 1.0d0/state%rho
  
  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)
 
  call ckums(state % T, ek)
  
  ! Derivative of Pressure w.r.t to Yk, species k massfraction
  do i = 1,nspecies

     Rmk = Ru*inv_mwt(i)

     Temp1 = (state%T * state%d2amdYkdT(i) - state%dAmdYk(i))*K1

     Temp2 = (state%T * state%dAmdT - state%am)*(-Bi(i)*K1 + InvEosT3Denom*Bi(i))/state%bm

     state%ei(i) = ek(i) + Temp1 + Temp2
     
  end do
  
end subroutine SRK_EOS_GetSpeciesE
!=================================================================!
! Given a mixture composition, T, P calculate using SRK           !
!  the following derivatives for testing, dP/dT, dP/dtau, dH/dtau !
!=================================================================!
subroutine  SRK_EOS_GetDerivative(state)
  implicit none
  type (eos_t), intent(inout) :: state
  integer :: i,j
  real(amrex_real) :: tau, K1
  real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
  real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
  real(amrex_real) :: Rm,Rmk,dhmdtau
  real(amrex_real) :: dAmdTYk(nspecies)
  real(amrex_real) :: Temp1, Temp2, Temp3, Temp4 

  state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))
  
  call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

  ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)

  ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
  call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)

  tau = 1.0d0/state%rho
  
  K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
  eosT1Denom = tau-state%bm
  eosT2Denom = tau*(tau+state%bm)
  eosT3Denom = tau+state%bm
  
  InvEosT1Denom = 1.0d0/eosT1Denom
  InvEosT2Denom = 1.0d0/eosT2Denom
  InvEosT3Denom = 1.0d0/eosT3Denom

  Rm = (Ru/state%wbar)

  ! Derivative of Pressure w.r.t to Temperature
  state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom
  
  ! Derivative of Pressure w.r.t to tau (specific volume)
  state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

  ! Derivative of enthalpy w.r.t to tau (specific volume)
  state%dhdr = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
       Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom
  
end subroutine SRK_EOS_GetDerivative

end module eos_module

