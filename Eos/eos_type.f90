module eos_type_module

  use amrex_fort_module, only : amrex_real
  use network, only: nspec, naux

  implicit none

  !integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  !integer, parameter :: eos_input_rh = 2  ! rho, h are inputs
  !integer, parameter :: eos_input_tp = 3  ! T, p are inputs
  !integer, parameter :: eos_input_rp = 4  ! rho, p are inputs
  !integer, parameter :: eos_input_re = 5  ! rho, e are inputs
  !integer, parameter :: eos_input_ps = 6  ! p, s are inputs
  !integer, parameter :: eos_input_ph = 7  ! p, h are inputs
  !integer, parameter :: eos_input_th = 8  ! T, h are inputs

  !integer, parameter :: eos_xty = 9  ! molefrac in, massfrac out only
  !integer, parameter :: eos_ytx = 10 ! massfrac in, molefrac out only
  !integer, parameter :: eos_cpi = 11 ! temp in, cp_i out only
  !integer, parameter :: eos_cp  = 12 ! temp,massfrac in, cvmix out only
  !integer, parameter :: eos_cv  = 13 ! temp,massfrac in, cvmix out only
  !integer, parameter :: eos_hi  = 14 ! temp in, h_i out only
  !integer, parameter :: eos_ei  = 15 ! temp in, e_i out only
  !integer, parameter :: eos_mui = 16 ! temp, Pressure and MassFraction in, chemical potential (mu_i) out
  !integer, parameter :: eos_h   = 17 ! temp in, massfrac in and pressure in, Mixture enthalpy out only 
  !integer, parameter :: eos_f   = 18 ! Mixture gibbs free energy
  !integer, parameter :: eos_deriv   = 19 ! temp in, massfrac in and pressure in, dP/dT, dP/dtau and dh/dtau out
  !integer, parameter :: eos_p_wb= 20 ! rho/T/massfrac in, p/wbar out
  !integer, parameter :: eos_get_activity = 21 !  get activty coefficient
  !integer, parameter :: eos_get_transport = 22 !  get transport terms

  ! these are used to allow for a generic interface to the 
  ! root finding
  integer, parameter :: itemp = 1
  integer, parameter :: idens = 2
  integer, parameter :: iener = 3
  integer, parameter :: ienth = 4
  integer, parameter :: ientr = 5
  integer, parameter :: ipres = 6

  ! error codes
  integer, parameter :: ierr_general         = 1
  integer, parameter :: ierr_input           = 2
  integer, parameter :: ierr_iter_conv       = 3
  integer, parameter :: ierr_neg_e           = 4
  integer, parameter :: ierr_neg_p           = 5
  integer, parameter :: ierr_neg_h           = 6
  integer, parameter :: ierr_neg_s           = 7
  integer, parameter :: ierr_iter_var        = 8
  integer, parameter :: ierr_init            = 9
  integer, parameter :: ierr_init_massfrac   = 10
  integer, parameter :: ierr_out_of_bounds   = 11
  integer, parameter :: ierr_not_implemented = 12

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
    real(amrex_real),allocatable :: massfrac(:)
    real(amrex_real),allocatable :: molefrac(:)
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
    real(amrex_real),allocatable :: cpi(:)
    real(amrex_real),allocatable :: cvi(:)
    real(amrex_real),allocatable :: hi(:)
    real(amrex_real),allocatable :: ei(:)
    real(amrex_real),allocatable :: si(:)
    real(amrex_real) :: wbar
    real(amrex_real),allocatable :: mui(:)
    real(amrex_real),allocatable :: Acti(:)
    real(amrex_real),allocatable :: dedY(:)
    real(amrex_real),allocatable :: dpdY(:)
    real(amrex_real),allocatable :: dhdY(:)
    real(amrex_real) :: gam1
    real(amrex_real) :: cs

    ! Quantities used for non-Ideal EOS
    real(amrex_real) :: am
    real(amrex_real) :: bm
    real(amrex_real),allocatable :: damdYk(:)
    real(amrex_real),allocatable :: d2amdYkdT(:)
    real(amrex_real),allocatable :: dPdYk(:)
    real(amrex_real) :: damdT
    real(amrex_real) :: d2amdT2
    real(amrex_real) :: dpdtau
    real(amrex_real) :: Z
    real(amrex_real),allocatable :: taui(:)
    real(amrex_real),allocatable :: diP(:)
    real(amrex_real),allocatable :: dijY(:,:)

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
    if (.not. allocated(eos%massfrac)) then
       allocate(eos%massfrac(nspec))
    endif
    if (.not. allocated(eos%cpi)) then
       allocate(eos%cpi(nspec))
    endif
    if (.not. allocated(eos%cvi)) then
       allocate(eos%cvi(nspec))
    endif
    if (.not. allocated(eos%hi)) then
       allocate(eos%hi(nspec))
    endif
    if (.not. allocated(eos%ei)) then
       allocate(eos%ei(nspec))
    endif
    if (.not. allocated(eos%molefrac)) then
       allocate(eos%molefrac(nspec))
    endif
    if (.not. allocated(eos%aux)) then
       allocate(eos%aux(naux))
    endif
    if (.not. allocated(eos%dedY)) then
       allocate(eos%dedY(nspec))
    endif
    if (.not. allocated(eos%dpdY)) then
       allocate(eos%dpdY(nspec))
    endif
    if (.not. allocated(eos%dhdY)) then
       allocate(eos%dhdY(nspec))
    endif
    if (.not. allocated(eos%dhdY)) then
       allocate(eos%dhdY(nspec))
    endif
    if (.not. allocated(eos%mui)) then
       allocate(eos%mui(nspec))
    end if
    if (.not. allocated(eos%damdYk)) then
       allocate(eos%damdYk(nspec))
    endif
    if(.not.allocated(eos%d2amdYkdT)) then
       allocate(eos%d2amdYkdT(nspec))
    end if
    if(.not.allocated(eos%dPdYk)) then
       allocate(eos%dPdYk(nspec))
    end if
    if (.not.allocated(eos%Acti)) then
       allocate(eos%Acti(nspec))
    end if
    if(.not.allocated(eos%si)) then
       allocate(eos%si(nspec))
    end if
    if(.not.allocated(eos%taui)) then 
       allocate(eos%taui(nspec))
    end if
    if(.not.allocated(eos%diP)) then
       allocate(eos%diP(nspec))
    end if
    if(.not.allocated(eos%dijY)) then
       allocate(eos%dijY(nspec,nspec))
    end if

  end subroutine eos_build
  
  subroutine eos_destroy(eos)
    type(eos_t), intent(inout) :: eos
    deallocate(eos%massfrac)
    deallocate(eos%cpi)
    deallocate(eos%cvi)
    deallocate(eos%hi)
    deallocate(eos%ei)
    deallocate(eos%molefrac)
    deallocate(eos%dedY)
    deallocate(eos%dpdY)
    deallocate(eos%dhdY)
    deallocate(eos%aux)
    deallocate(eos%damdYk)
    deallocate(eos%mui)
    deallocate(eos%d2amdYkdT)
    deallocate(eos%dPdYk)
    deallocate(eos%Acti)
    deallocate(eos%si)
    deallocate(eos%taui)
    deallocate(eos%diP)
    deallocate(eos%dijY)

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
