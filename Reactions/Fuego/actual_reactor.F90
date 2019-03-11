module actual_reactor_module

  use amrex_fort_module, only : amrex_real
  use network, only: nspec, spec_names
  use react_type_module
  use eos_type_module

  implicit none

  real(amrex_real), private, allocatable :: vodeVec(:),cdot(:),rhoydot_ext(:),ydot_ext(:)
  real(amrex_real), private:: rhoedot_ext, rhoe_init, time_init, time_out, rhohdot_ext, &
                              rhoh_init, time_old, hdot_ext, h_init, pressureInit
  integer,private :: iloc, jloc, kloc, iE, iDense
  type (eos_t) :: eos_state
  !$omp threadprivate(vodeVec,cdot,rhoydot_ext,ydot_ext,rhoedot_ext,rhoe_init,time_init,time_out,rhohdot_ext,rhoh_init,hdot_ext,h_init,time_old,iloc,jloc,kloc,eos_state)

contains

        
!*** INITIALISATION ROUTINES ***!
  !DVODE VERSION
  subroutine actual_reactor_init(iE_in)

    use, intrinsic :: iso_c_binding
    use vode_module, only : vode_init

    integer(c_int),  intent(in   ) :: iE_in
    integer :: neq, verbose, itol, order, maxstep
    real(amrex_real) :: rtol, atol
    logical :: use_ajac, save_ajac, always_new_j, stiff

    neq = nspec + 1
    verbose = 0
    itol = 1
    order = 2
    maxstep = 10000
    use_ajac = .false.
    save_ajac = .false.
    always_new_J = .false.
    stiff = .true.
    rtol = 1.d-10
    atol = 1.d-10

    call vode_init(neq,verbose,itol,rtol,atol,order,&
         maxstep,use_ajac,save_ajac,always_new_j,stiff)

    print *,"Using good ol' dvode"
    print *,"--> DENSE solver without Analytical J"
    iE = iE_in
    if (iE == 1) then
        print *," ->with internal energy (UV cst)"
        allocate(rhoydot_ext(nspec))
    else if (iE == 5) then
        print *," ->with enthalpy (HP cst)"
        allocate(ydot_ext(nspec))
    else
        print *," ->with enthalpy (sort of rhoP cst)"
        allocate(rhoydot_ext(nspec))
    end if 

    allocate(vodeVec(neq))
    allocate(cdot(nspec))

    call build(eos_state)

  end subroutine actual_reactor_init

!*** REACTION ROUTINES ***!
  ! Original DVODE version
  function actual_react(react_state_in, react_state_out, dt_react, time)
    
    use amrex_error_module
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    use chemistry_module, only : molecular_weight
    use eos_module

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: actual_react

    external dvode

    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail, neq
    real(amrex_real) :: vodeTime, vodeEndTime, rhoInv

    eos_state % rho               = sum(react_state_in % rhoY(:))
    eos_state % T                 = react_state_in % T
    rhoInv                         = 1.d0 / eos_state % rho
    eos_state % massfrac(1:nspec) = react_state_in % rhoY(1:nspec) * rhoInv

    if (iE == 1) then
        eos_state % e = react_state_in % e
        call eos_re(eos_state)
    else if (iE == 5) then
        ! for cst HP
        pressureInit  = react_state_in % p
        eos_state % p = react_state_in % p
        eos_state % h = react_state_in % h
        call eos_ph(eos_state)
    else
        eos_state % h = react_state_in % h
        call eos_rh(eos_state)
    end if

    if (always_new_j) call setfirst(.true.)

    MF          = vode_MF
    vodeTime    = time
    vodeEndTime = time + dt_react
    neq         = nspec + 1
    time_old    = time

    vodeVec(neq)     = eos_state % T

    if (iE == 1) then
        rhoe_init            = eos_state % e  *  eos_state % rho
        rhoedot_ext          = react_state_in % rhoedot_ext
        rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
        vodeVec(1:nspec)     = react_state_in % rhoY(:)
    else if (iE == 5) then
        h_init            = eos_state % h  
        hdot_ext          = react_state_in % rhohdot_ext / eos_state % rho
        ydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec) / eos_state % rho
        vodeVec(1:nspec)  = react_state_in % rhoY(:) / eos_state % rho
    else
        rhoh_init            = eos_state % h  *  eos_state % rho
        rhohdot_ext          = react_state_in % rhohdot_ext 
        rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
        vodeVec(1:nspec)     = react_state_in % rhoY(:)
    end if

    time_init = time
    iloc      = react_state_in % i
    jloc      = react_state_in % j
    kloc      = react_state_in % k

    ! Vode: istate
    ! in:  1: init, 2: continue (no change), 3: continue (w/changes)
    ! out: 1: nothing was done, 2: success
    !      -1: excessive work, -2: too much accuracy, -3: bad input
    !      -4: repeated step failure, -5: repeated conv failure
    !      -6: EWT became 0
    istate = 1

    call dvode(f_rhs, neq, vodeVec(:), vodeTime, vodeEndTime,&
         itol, rtol, atol, itask, istate, iopt, voderwork, lvoderwork, &
         vodeiwork, lvodeiwork, f_jac, MF, voderpar, vodeipar)

    if (verbose .ge. 1) then
       write(6,*) '......dvode done:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '    number of fs(RHS EVAL) = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (istate.eq.-4 .or. istate.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspec+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if
    end if

    if (istate > 0) then

       actual_react % reactions_succesful = .true.
       actual_react % cost_value = DBLE(vodeiwork(12)) ! number of f evaluations

       if (iE == 1) then
           eos_state % rho               = sum(vodeVec(1:nspec))
           rhoInv                        = 1.d0 / eos_state % rho
           eos_state % massfrac(1:nspec) = vodeVec(1:nspec) * rhoInv
           eos_state % T                 = vodeVec(neq)
           eos_state % e                 = (rhoe_init  +  dt_react*rhoedot_ext) /eos_state % rho
           call eos_re(eos_state)
           react_state_out % rhoY(:)     = vodeVec(1:nspec) 
           react_state_out % rho         = sum(vodeVec(1:nspec))
           react_state_out % rhoedot_ext = rhoedot_ext
           react_state_out % rhoydot_ext(1:nspec) = rhoydot_ext(1:nspec)
       else if (iE == 5) then
           eos_state % p                 = pressureInit  
           eos_state % massfrac(1:nspec) = vodeVec(1:nspec)
           eos_state % T                 = vodeVec(neq)
           eos_state % h                 = (h_init  +  dt_react*hdot_ext)
           call eos_ph(eos_state)
           react_state_out % rhoY(:)     = vodeVec(1:nspec) * eos_state % rho
           react_state_out % rho         = eos_state % rho 
           react_state_out % rhohdot_ext = hdot_ext * eos_state % rho
           !react_state_out % rhohdot_ext = react_state_in % rhohdot_ext
           react_state_out % rhoydot_ext(1:nspec) = ydot_ext(1:nspec) *  eos_state % rho
           !react_state_out % rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
       else
           eos_state % rho               = sum(vodeVec(1:nspec))
           rhoInv                        = 1.d0 / eos_state % rho
           eos_state % massfrac(1:nspec) = vodeVec(1:nspec) * rhoInv
           eos_state % T                 = vodeVec(neq)
           eos_state % h                 = (rhoh_init  +  dt_react*rhohdot_ext) * rhoInv
           call eos_rh(eos_state)
           react_state_out % rhoY(:)     = vodeVec(1:nspec)
           react_state_out % rho         = sum(vodeVec(1:nspec))
           react_state_out % rhohdot_ext = rhohdot_ext 
           react_state_out % rhoydot_ext(1:nspec) = rhoydot_ext(1:nspec)
       end if

       react_state_out % T = eos_state % T
       react_state_out % e = eos_state % e
       react_state_out % h = eos_state % h
       react_state_out % p = eos_state % p

    else

       actual_react % reactions_succesful = .false.

       print *,'vode failed at',react_state_in % i,react_state_in % j,react_state_in % k
       print *,'input state:'
       print *,'T',react_state_in%T
       if (iE == 1) then
           print *,'e',react_state_in%e
       else
           print *,'h',react_state_in%h
       end if
       print *,'rho',eos_state%rho
       print *,'rhoY',react_state_in%rhoY
       if (iE == 1) then
           print *,'rhoe forcing',react_state_in%rhoedot_ext
       else
           print *,'rhoh forcing',react_state_in%rhohdot_ext
       end if
       print *,'rhoY forcing',react_state_in%rhoydot_ext(1:nspec)

       write(6,*) '......dvode data:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '              number of fs = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (istate.eq.-4 .or. istate.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspec+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if

       print *,'Final T',vodeVec(neq)
       print *,'Final rhoY',vodeVec(1:nspec)

       call amrex_error('vode failed')

    end if
  end function actual_react

  ! Original DVODE version
  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use chemistry_module, only : molecular_weight
    use eos_module

    integer,         intent(in)   :: neq, ipar(*)
    real(amrex_real), intent(in)  :: y(neq), time, rpar(*)
    real(amrex_real), intent(out) :: ydot(neq)
    integer          :: n
    real(amrex_real) :: rhoInv


    if (iE == 1) then
        eos_state % rho               = sum(y(1:nspec))
        rhoInv                        = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspec) = y(1:nspec) * rhoInv
        eos_state % T                 = y(neq) ! guess
        eos_state % e = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        call eos_get_activity(eos_state)
    else if (iE == 5) then
        eos_state % massfrac(1:nspec) = y(1:nspec) 
        eos_state % T                 = y(neq) ! guess
        eos_state % h = h_init + (time - time_init) * hdot_ext
        eos_state % p                 = pressureInit  
        call eos_ph(eos_state)
        call eos_get_activity_h(eos_state)
    else
        eos_state % rho               = sum(y(1:nspec))
        rhoInv                        = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspec) = y(1:nspec) * rhoInv 
        eos_state % T                 = y(neq) ! guess
        eos_state % h = (rhoh_init + (time - time_init) * rhohdot_ext) * rhoInv
        call eos_rh(eos_state)
        call eos_get_activity_h(eos_state)
    end if

    call ckwc(eos_state % T, eos_state % Acti, iwrk, rwrk, cdot)

    if (iE == 1) then
        ydot(neq)    = rhoedot_ext 
        do n=1,nspec
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%ei(n)*ydot(n)
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cv)
    else if (iE == 5) then
        ydot(neq)    = hdot_ext
        do n=1,nspec
           ydot(n)   = cdot(n) * molecular_weight(n) / eos_state%rho + ydot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/eos_state%cp
    else
        ydot(neq)    = rhohdot_ext
        do n=1,nspec
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cp)
    end if

  end subroutine f_rhs

  ! Original DVODE version
  subroutine f_jac(neq, npt, y, t, pd)
    use amrex_error_module

    integer,        intent(in)  :: neq, npt
    real(amrex_real),intent(in)  :: y(neq,npt), t
    real(amrex_real),intent(out) :: pd(neq,neq)

    call amrex_error('DVODE version: Analytic Jacobian not yet implemented')

  end subroutine f_jac

!*** FINALIZE ROUTINES ***!
  subroutine actual_reactor_close()

    if (allocated(vodeVec)) deallocate(vodeVec)
    if (allocated(cdot)) deallocate(cdot)
    if (allocated(rhoydot_ext)) deallocate(rhoydot_ext)
    if (allocated(ydot_ext)) deallocate(ydot_ext)

    call destroy(eos_state)
   
  end subroutine actual_reactor_close


!*** SPECIFIC ROUTINES ***!
  function actual_ok_to_react(state)

    use extern_probin_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: actual_ok_to_react
    real(amrex_real)           :: rho

    actual_ok_to_react = .true.

    rho = sum(state % rhoY)
    if (state % T   < react_T_min   .or. state % T   > react_T_max .or. &
        rho         < react_rho_min .or. rho         > react_rho_max) then

       actual_ok_to_react = .false.

    endif

  end function actual_ok_to_react


  function actual_react_null(react_state_in, react_state_out, dt_react, time)
    
    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: actual_react_null

    react_state_out = react_state_in
    actual_react_null % cost_value = 0.d0
    actual_react_null % reactions_succesful = .true.

  end function actual_react_null


end module actual_reactor_module
