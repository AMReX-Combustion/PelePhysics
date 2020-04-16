module eos_type_module

  use amrex_fort_module, only : amrex_real
  use fuego_chemistry, only: nspecies, naux

  implicit none

  ! Minimum and maximum thermodynamic quantities permitted by the EOS.
  real(amrex_real), save :: mintemp
  real(amrex_real), save :: maxtemp
  real(amrex_real), save :: mindens
  real(amrex_real), save :: maxdens
  real(amrex_real), save :: minmassfrac
  real(amrex_real), save :: maxmassfrac
  real(amrex_real), save :: mine
  real(amrex_real), save :: maxe
  real(amrex_real), save :: minp
  real(amrex_real), save :: maxp
  real(amrex_real), save :: mins
  real(amrex_real), save :: maxs
  real(amrex_real), save :: minh
  real(amrex_real), save :: maxh

  !$acc declare &
  !$acc create(mintemp, maxtemp, mindens, maxdens, minmassfrac, maxmassfrac, minye, maxye) &
  !$acc create(mine, maxe, minp, maxp, mins, maxs, minh, maxh)

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! massfrac -- the mass fractions of the individual species
  ! molefrac -- the mole fractions of the individual species
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the mixture enthalpy (erg/g)
  ! hi       -- the species enthalpy (erg/g) of species i
  ! e        -- the mixture internal energy (erg/g)
  ! ei       -- the internal energy (erg/g) of species i
  ! s        -- the entropy (erg/g/K)
  ! cv       -- specific heat at constant volume
  ! cvi      -- specific heat at constant volume for each species i
  ! cp       -- specific heat at constant pressure
  ! cpi      -- specific heat at constant pressure for each species i
  ! wbar     -- mean molecular weight
  ! mui      -- Chemical potential of each species i
  ! Acti     -- Activity Coefficient of each species i
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdY     -- d pressure / d massfrac
  ! dhdY     -- d enthalpy / d massfrac at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed

  ! Specific to non-Ideal EOS (Peng-Robinson & SRK)
  ! am       -- Mixture EOS attractive parameter
  ! bm       -- Mixture EOS repulsive parameter
  ! damdT    -- Derivative of am w.r.t Temperature
  ! damdyk   -- Derivative of am w.r.t Species k massfraction
  ! dbmdyk   -- Derivative of bm w.r.t Species k massfraction
  ! d2amdT2  -- Second derivative of am w.r.t Temperature
  ! dpdtau   -- Derivative of Pressure w.r.t specific volume (tau = 1/rho)
  ! Z        -- Compressibility factor 

  
  type :: eos_t

    real(amrex_real) :: rho
    real(amrex_real) :: T
    real(amrex_real) :: p
    real(amrex_real) :: e
    real(amrex_real) :: h
    real(amrex_real) :: s
    real(amrex_real) :: f
    real(amrex_real) :: massfrac(nspecies)
    real(amrex_real) :: molefrac(nspecies)
    real(amrex_real),allocatable :: aux(:)

    real(amrex_real) :: dpdT
    real(amrex_real) :: dpdr
    real(amrex_real) :: dedT
    real(amrex_real) :: dedr
    real(amrex_real) :: dhdT
    real(amrex_real) :: dhdr
    real(amrex_real) :: dsdT
    real(amrex_real) :: dsdr
    real(amrex_real) :: dpde
    real(amrex_real) :: dpdr_e

    real(amrex_real) :: cv
    real(amrex_real) :: cp
    real(amrex_real) :: cpi(nspecies)
    real(amrex_real) :: cvi(nspecies)
    real(amrex_real) :: hi(nspecies)
    real(amrex_real) :: ei(nspecies)
    real(amrex_real) :: si(nspecies)
    real(amrex_real) :: wbar
    real(amrex_real) :: mui(nspecies)
    real(amrex_real) :: Acti(nspecies)
    real(amrex_real) :: dedY(nspecies)
    real(amrex_real) :: dpdY(nspecies)
    real(amrex_real) :: dhdY(nspecies)
    real(amrex_real) :: gam1
    real(amrex_real) :: cs

    ! Quantities used for non-Ideal EOS
    real(amrex_real) :: am
    real(amrex_real) :: bm
    real(amrex_real) :: damdYk(nspecies)
    real(amrex_real) :: d2amdYkdT(nspecies)
    real(amrex_real) :: dPdYk(nspecies)
    real(amrex_real) :: damdT
    real(amrex_real) :: d2amdT2
    real(amrex_real) :: dpdtau
    real(amrex_real) :: Z
    real(amrex_real) :: taui(nspecies)
    real(amrex_real) :: diP(nspecies)
    real(amrex_real) :: dijY(nspecies,nspecies)

  end type eos_t

  interface build
     module procedure eos_build
  end interface build

  interface destroy
     module procedure eos_destroy
  end interface destroy


contains

  subroutine eos_build(eos)
    type(eos_t), intent(inout) :: eos
    if (.not. allocated(eos%aux) .and. naux.gt.0) then
       allocate(eos%aux(naux))
    endif
  end subroutine eos_build
  
  subroutine eos_destroy(eos)
    type(eos_t), intent(inout) :: eos
    if (naux.gt.0) then
       deallocate(eos%aux)
    endif
  end subroutine eos_destroy
  
  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition ... nothing to do here yet
  subroutine composition(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state

  end subroutine composition
  
  ! Compute thermodynamic derivatives with respect to massfrac(:)
  subroutine composition_derivatives(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state

    ! No need for dpdY, dEdY, dhdY yet

  end subroutine composition_derivatives

  ! Normalize the mass fractions: they must be individually positive
  ! and less than one, and they must all sum to unity.
  subroutine normalize_massfracs(state)

    !$acc routine seq

    use amrex_constants_module, only: ONE
    use extern_probin_module, only: small_massfrac

    implicit none

    type (eos_t), intent(inout) :: state

    state % massfrac = max(small_massfrac, min(ONE, state % massfrac))

    state % massfrac = state % massfrac / sum(state % massfrac)

  end subroutine normalize_massfracs
  
  ! Ensure that inputs are within reasonable limits.
  subroutine clean_state(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine clean_state

  ! Print out details of the state.
  subroutine print_state(state)

    implicit none

    type (eos_t), intent(in) :: state

    print *, 'DENS = ', state % rho
    print *, 'TEMP = ', state % T
    print *, 'Y    = ', state % massfrac

  end subroutine print_state


  subroutine eos_get_small_temp(small_temp_out)

    !$acc routine seq

    implicit none

    real(amrex_real), intent(out) :: small_temp_out

    small_temp_out = mintemp

  end subroutine eos_get_small_temp


  subroutine eos_get_small_dens(small_dens_out)

    !$acc routine seq

    implicit none

    real(amrex_real), intent(out) :: small_dens_out

    small_dens_out = mindens

  end subroutine eos_get_small_dens



  subroutine eos_get_max_temp(max_temp_out)

    !$acc routine seq

    implicit none

    real(amrex_real), intent(out) :: max_temp_out

    max_temp_out = maxtemp

  end subroutine eos_get_max_temp



  subroutine eos_get_max_dens(max_dens_out)

    !$acc routine seq

    implicit none

    real(amrex_real), intent(out) :: max_dens_out

    max_dens_out = maxdens

  end subroutine eos_get_max_dens

end module eos_type_module
