module spray_module

  use amrex_fort_module

  implicit none

  public

contains

  subroutine get_num_fuel_spec(n_fuel_spec) &
       bind(C, name="get_num_fuel_spec")

    use fuel_properties

    integer, intent(out) :: n_fuel_spec

    n_fuel_spec = nspec_f

  end subroutine get_num_fuel_spec

  ! ------------------------------------------------------------------------------
  !  Import fuel properties from input file
  ! ------------------------------------------------------------------------------
  subroutine import_fuel_properties(ccnspec_f, ccfuel_mass_frac, ccfuel_density, &
       ccfuel_crit_temp, ccfuel_latent, ccfuel_boil_temp, &
       ccfuel_cp, ccfuel_molwt, ccfuel_indx) &
       bind(C, name="import_fuel_properties")

    use fuel_properties

    integer :: ccnspec_f
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_mass_frac
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_density
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_crit_temp
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_latent
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_boil_temp
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_cp
    real(amrex_real), dimension(ccnspec_f) :: ccfuel_molwt
    integer, dimension(ccnspec_f) :: ccfuel_indx

    nspec_f = ccnspec_f
    fuel_mass_frac(1:nspec_f)=ccfuel_mass_frac
    fuel_density(1:nspec_f)=ccfuel_density
    fuel_crit_temp(1:nspec_f)=ccfuel_crit_temp
    fuel_latent(1:nspec_f)=ccfuel_latent
    fuel_boil_temp(1:nspec_f)=ccfuel_boil_temp
    fuel_cp(1:nspec_f)=ccfuel_cp
    fuel_molwt(1:nspec_f)=ccfuel_molwt
    fuel_indx(1:nspec_f)=ccfuel_indx

  end subroutine import_fuel_properties

  ! ------------------------------------------------------------------------------
  !  Import controls for spray model
  ! ------------------------------------------------------------------------------
  subroutine import_control_parameters(ccheat_transfer, ccmass_transfer, ccmom_transfer) &
       bind(C, name="import_control_parameters")

    use control_parameters

    integer, intent(in) :: ccheat_transfer, ccmass_transfer, ccmom_transfer

    is_heat_tran = ccheat_transfer
    is_mass_tran = ccmass_transfer
    is_mom_tran  = ccmom_transfer

  end subroutine import_control_parameters

  ! ------------------------------------------------------------------------------
  ! > Calculate specific Cp of fuel
  ! > This actually calculates cp_i i.e. cp of all species, and then discards most
  ! > This may be pushed into Eos package as a general wrapper like get_mixcp
  ! ------------------------------------------------------------------------------

  subroutine calc_spec_mix_cp_spray (eos_state,massfrac,temp,mix_cp,fuel_spec_cp)

    use eos_module
    use chemistry_module, only : nspecies
    use fuel_properties, only: nspec_f, fuel_indx

    double precision, dimension(nspecies), intent(in) :: massfrac
    double precision, intent(in)  :: temp
    double precision, intent(out) :: mix_cp
    double precision, dimension(nspec_f), intent(out) :: fuel_spec_cp

    ! Local variables
    type(eos_t), intent(inout) :: eos_state

    integer L

    eos_state%T = temp

    !Actual call. This returns values in %cpi

    call eos_cpi(eos_state)

    !Loop over fuel indices and 
    !assign to the corresponding index

    do L = 1, nspec_f
       fuel_spec_cp(L) = eos_state%cpi(fuel_indx(L))
    enddo

    !Expression for gas mixture cp. These are extra flops if
    !we're using simple Gammalaw, but we save elsewhere
    !HK: Can I use the fortran intrinsic sum cleverly
    mix_cp = 0.0
    do L = 1, nspecies
       mix_cp = mix_cp + massfrac(L)*eos_state%cpi(L) 
    enddo

  end subroutine calc_spec_mix_cp_spray

  ! ------------------------------------------------------------------------------
  ! > Calculate the Latent Heat of fuel
  ! ------------------------------------------------------------------------------
  subroutine calc_fuel_latent(crit_temp,inv_diff_temp,ref_latent,temp_d,L_f)

    ! Input
    double precision, intent(in) :: crit_temp
    double precision, intent(in) :: inv_diff_temp
    double precision, intent(in) :: ref_latent
    double precision, intent(in) :: temp_d

    ! Output
    double precision, intent(out) :: L_f

    ! Executables
    L_f = (crit_temp-temp_d)*inv_diff_temp
    L_f = max(L_f,0d0)
    L_f = ref_latent*L_f**0.38
    !L_f = 0.97*9.453*(619.-temp_d)**0.38*4.184d7

  end subroutine calc_fuel_latent



  ! ------------------------------------------------------------------------------
  ! > Calculate the spalding mass transfer number
  ! ------------------------------------------------------------------------------
  subroutine calc_spalding_num(L_f,temp_d,gas_press,gas_ysp,gas_molwt,fuel_MW,fuel_invMW,&
       inv_boil,invRu,ref_press,B)

    ! Input
    double precision, intent(in) :: L_f
    double precision, intent(in) :: temp_d
    double precision, intent(in) :: gas_press
    double precision, intent(in) :: gas_ysp
    double precision, intent(in) :: gas_molwt
    double precision, intent(in) :: fuel_MW
    double precision, intent(in) :: fuel_invMW
    double precision, intent(in) :: inv_boil
    double precision, intent(in) :: invRu
    double precision, intent(in) :: ref_press

    ! Output
    double precision, intent(out) :: B

    ! Locals
    double precision, parameter :: eps = 1.0e-15
    double precision :: press_sat
    double precision :: Y_surf

    ! Calculate the saturation pressure
    press_sat = L_f*fuel_MW*invRu*(inv_boil - 1.0d0/temp_d)
    press_sat = ref_press*exp(press_sat)

    ! Calculate Surface fuel mass fraction
    Y_surf = (gas_press/(eps+press_sat)-1.0d0)
    Y_surf = 1.0d0+Y_surf*gas_molwt*fuel_invMW
    Y_surf = 1.0d0/Y_surf

    ! Make sure the surface fuel mass fraction is always a tiny bit smaller than
    ! 1.0.
    ! Otherwise you'll have convergence issues for particles that are in high gas
    ! fuel regions.
    Y_surf = min(Y_surf,1.0d0-eps)

    ! Calculate Spalding number
    B = max((Y_surf-gas_ysp)/(1.0d0-Y_surf),0d0) ! Eq. (12)
    B = min(B,20d0)

  end subroutine calc_spalding_num


  ! ------------------------------------------------------------------------------
  ! > Calculate evaporation rate
  ! ------------------------------------------------------------------------------

  subroutine calc_spec_evap_rate(diameter,B,Re,Sc,diff_coeff,Sh,evap)

    ! Input
    double precision, intent(in) :: diameter
    double precision, intent(in) :: B
    double precision, intent(in) :: Re
    double precision, intent(in) :: Sc
    double precision, intent(in) :: diff_coeff

    ! Output
    double precision, intent(out) :: evap
    double precision, intent(out) :: Sh

    ! Locals
    double precision, parameter :: eps = 1.0e-15
    double precision :: log_B
    double precision :: F_M

    ! Constants
    double precision, parameter :: one_third = 1.0d0/3.0d0
    double precision, parameter :: Pi = 4*atan(1.0d0)

    ! Executables
    ! log_B = log(1+B_m)             
    log_B = log(1.0d0+B)

    F_M = log_B*(1.0d0+B)**(0.7)
    F_M = B/F_M
    F_M = max(F_M,0.0d0) ! The inverse of F_M from Eq. (17)

    !-----------------------------------------------------------------------------------------            
    ! Calculate Sherwood Number                   
    Sh = 1.0d0+max(Re**0.077,1.0d0)*(1.0d0+Re*Sc)**(one_third) ! Eq. (21)
    Sh = 2.0d0+(Sh-2.0d0)*F_M ! Eq. (10)

    !-----------------------------------------------------------------------------------------    
    ! Calculate Evaporation Rate          
    ! Calculate rate of change in mass
    evap = -Pi*max(diff_coeff*diameter*Sh*log_B,0d0) ! Eq. (8)

  end subroutine calc_spec_evap_rate


  !-----------------------------------------------------------------------------------------
  ! SUBROUTINE TO ITERATE AND CALCULATE B_T
  !-----------------------------------------------------------------------------------------

  subroutine calc_thermal_B(Nu, Sh, Pr, Sc,cp_fuel,cp_gas, B_m, B_T)

    ! Input Variables
    double precision, intent(   in) :: Sh
    double precision, intent(   in) :: Pr
    double precision, intent(   in) :: Sc
    double precision, intent(   in) :: cp_fuel
    double precision, intent(   in) :: cp_gas
    double precision, intent(   in) :: B_m

    ! Output variables
    double precision, intent(inout) :: Nu
    double precision, intent(  out) :: B_T

    ! Local
    double precision :: B_T_old
    double precision :: F_T
    double precision :: phi
    double precision :: ratio
    double precision :: error
    double precision :: Nu_tmp, Nu_0
    double precision :: log_BT
    double precision :: eps
    integer, parameter :: imax = 5000
    integer :: ii
    integer :: icount

    ! Executables
    eps = 1.0e-5

    ! Collect all the constants together or else we'll waste a lot of computations
    ratio = (cp_fuel/cp_gas)*(Pr/Sc)*Sh

    ! Copy Nu data. Since Nu is constantly updated I can't guarantee it will be
    ! the correct
    ! answer if Nu is used in the iteration, though it probably should be.
    Nu_0 = Nu
    Nu_tmp = Nu

    ! Take an initial guess for B_T
    phi = ratio*(1.0d0/Nu_tmp)
    B_T = ((1.0d0+B_m)**phi)-1.0d0 ! (Eq. 23)

    icount = 0

    iteration: do ii = 1,imax

       B_T_old = B_T

       log_BT = log(1.0d0+B_T)

       F_T = B_T/(log_BT*(1.0d0+B_T)**0.7) ! inverse of Eq. 17

       ! This is needed incase B_T = 0.0, some compilers will not be able to find
       ! a solution.
       F_T = min(F_T,1.0d0)

       !Nu_tmp = 2.0d0+(Nu_tmp-2.0d0)*F_T ! (Eq. 11)
       Nu_tmp = 2.0d0+(Nu_0-2.0d0)*F_T ! (Eq. 11)

       phi = ratio*(1.0d0/Nu_tmp)

       B_T = (1.0d0+B_m)**(phi)-1.0d0 ! (Eq. 23)

       error = abs(B_T-B_T_old)

       icount = icount+1

       if(error.lt.eps) exit iteration

    end do iteration

    if(error.gt.eps) then
       write(*,*) 'B_T iteration did not converge'
       write(*,*) 'error ', error
       write(*,*) 'B_T ', B_T, B_T
    !  stop
    endif

    B_T = max(B_T,0d0)
    B_T = min(B_T,20d0)  ! Abr & Sir page 1608

    log_BT = log(1.0d0+B_T)

    F_T = B_T/(log_BT*(1.0d0+B_T)**0.7)

    F_T = min(F_T,1.0d0)

    ! Nu = 2.0d0+(Nu-2.0d0)*F_T ! Nu*
    Nu = 2.0d0+(Nu_0-2.0d0)*F_T ! Nu*

    F_T = log_BT/B_T

    F_T = min(F_T,1.0d0)

    Nu = Nu*F_T ! the equivalent of Eq. 15

  end subroutine calc_thermal_B


  !=======================================================================
  ! Calculate the skin enthalpy (Make this easier to call later)
  !=======================================================================
  subroutine calc_skin_enth(n,n_indx,indx,temp_d,temp_g,latent_fuel,enth_skin) &
       bind(C, name="calc_skin_enth")

    use fuel_properties
    ! Input
    integer, intent(in) :: n
    integer, intent(in) :: n_indx
    integer, dimension(n_indx), intent(in) :: indx
    double precision, dimension(n), intent(in) :: temp_d
    double precision, dimension(n), intent(in) :: temp_g
    double precision, dimension(n,n_indx), intent(in) :: latent_fuel

    ! Output
    double precision, dimension(n,n_indx), intent(out) :: enth_skin

    ! Locals
    double precision, dimension(n) :: skin_temp
    double precision, dimension(n,n_indx) :: h_vap
    real, parameter :: one_third = 1.0d0/3.0d0

    ! Executables
    ! Calculate the skin temperature
    skin_temp(1:n) = min(temp_d(1:n) + one_third*(temp_g(1:n)-temp_d(1:n)),fuel_boil_temp(1))

    ! Calculate Enthalpy at skin temp for each species
    !call calc_specEnth_spray(n_indx,indx,n,skin_temp(1:n),enth_skin(1:n,:))
    enth_skin(1:n,1) = (-225.2e7+(skin_temp(1:n)-298.)*224.64e4)/fuel_molwt(1)

    ! Calculate Enthalpy at drop temp for each species
    !call calc_specEnth_spray(n_indx,indx,n,temp_d(1:n),h_vap(1:n,:))
    h_vap(1:n,1) = (-225.2e7+(temp_d(1:n)-298.)*224.64e4)/fuel_molwt(1)

    ! Calculate energy needed to evaporate plus energy needed to raise temperature
    ! of vapor.
    enth_skin(1:n,:) = -(enth_skin(1:n,:)-(h_vap(1:n,:)-latent_fuel(1:n,:)))

  end subroutine calc_skin_enth

end module spray_module
