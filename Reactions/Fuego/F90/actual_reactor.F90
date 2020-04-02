module reactor_module

  use amrex_fort_module, only : amrex_real
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use network, only: nspecies, spec_names
  use react_type_module
  use eos_type_module

  implicit none

  integer, parameter :: neq = nspecies + 1
  real(amrex_real), private :: vodeVec(neq),cdot(nspecies),rhoydot_ext(nspecies),ydot_ext(nspecies)
  real(amrex_real), private:: rhoedot_ext, rhoe_init, time_init, rhohdot_ext, &
                              rhoh_init, time_old, hdot_ext, h_init, pressureInit
  integer,private :: iloc, jloc, kloc, iE
  type (eos_t) :: eos_state

  logical, save, private :: reactor_initialized = .false.

  !$omp threadprivate(vodeVec,cdot,rhoydot_ext,ydot_ext,rhoedot_ext,rhoe_init,time_init,rhohdot_ext,rhoh_init,hdot_ext,h_init,pressureInit,time_old,iloc,jloc,kloc,ie,eos_state,reactor_initialized)

contains

        
!*** INITIALISATION ROUTINES ***!
  !DVODE VERSION
  subroutine reactor_init(iE_in, Ncells) bind(C, name="reactor_init")

    use, intrinsic :: iso_c_binding
    use vode_module, only : vode_init
    use extern_probin_module, only : new_Jacobian_each_cell
    use amrex_omp_module

    implicit none
    integer(c_int),  intent(in   ) :: iE_in
    integer(c_int),  intent(in   ), optional :: Ncells
    integer :: verbose, itol, order, maxstep
    real(amrex_real) :: rtol, atol
    logical :: use_ajac, save_ajac, always_new_j_loc, stiff, isio

    verbose = 0
    itol = 1
    order = 2
    maxstep = 10000
    use_ajac = .false.
    save_ajac = .false.
    if (new_Jacobian_each_cell .ne. 0) then
       always_new_j_loc = .true.
    else
       always_new_j_loc = .false.
    endif
    stiff = .true.
    rtol = 1.d-10
    atol = 1.d-10

    call vode_init(neq,verbose,itol,rtol,atol,order,&
         maxstep,use_ajac,save_ajac,always_new_j_loc,stiff)

    isio = parallel_IOProcessor() .and. omp_get_thread_num().eq.0
    if (isio) then
       print *,"Using good ol' dvode"
       print *,"--> DENSE solver without Analytical J"
       print *,"--> Always new J ? ",always_new_j_loc
    endif
    iE = iE_in
    if (iE == 1) then
       if (isio) print *," ->with internal energy (UV cst)"
    else if (iE == 5) then
       if (isio) print *," ->with enthalpy (HP cst)"
    else
       if (isio) print *," ->with enthalpy (sort of rhoP cst)"
    end if 

    call build(eos_state)

    reactor_initialized = .true.

  end subroutine reactor_init

!*** REACTION ROUTINES ***!
  ! Original DVODE version
  !function react(react_state_in, react_state_out, dt_react, time) bind(C, name="react") result(stat)
  function react(rY_in,rY_src_in,rX_in,rX_src_in,P_in,dt_react,time) bind(C, name="react") result(cost_value)
    
    use amrex_error_module
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    use eos_module

    implicit none
    real(amrex_real),   intent(inout) :: rY_in(nspecies+1),rY_src_in(nspecies)
    real(amrex_real),   intent(inout) :: rX_in,rX_src_in,P_in
    real(amrex_real)                  :: dt_react, time
    integer                           :: cost_value
    
    ! For compatibility to remove later
    type(react_t) :: react_state_in

    external dvode

    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail, neq
    real(amrex_real) :: vodeTime, vodeEndTime, rhoInv

    ! For compatibility to remove later
    call build(react_state_in)

    react_state_in %              T = rY_in(nspecies+1)
    react_state_in %        rhoY(:) = rY_in(1:nspecies)
    react_state_in %            rho = sum(react_state_in % rhoY(:))
    react_state_in % rhoYdot_ext(:) = rY_src_in(1:nspecies)

    if (iE == 1) then
        react_state_in %              e = rX_in !/ react_state_in % rho
        react_state_in %    rhoedot_ext = rX_src_in
    else if (iE == 5) then
        react_state_in %              p = P_in
        react_state_in %              h = rX_in !/ react_state_in % rho
        react_state_in %    rhohdot_ext = rX_src_in
    else
        react_state_in %              h = rX_in / react_state_in % rho
        react_state_in %    rhohdot_ext = rX_src_in
    end if
    ! END For compatibility to remove later


    if (.not. reactor_initialized) then
       call amrex_error('reactor::react called before initialized')
    endif

    if ( .not. ok_to_react(react_state_in) ) then
       cost_value = 0.d0
       return
    end if

    eos_state % rho               = sum(react_state_in % rhoY(:))
    eos_state % T                 = react_state_in % T
    rhoInv                         = 1.d0 / eos_state % rho
    eos_state % massfrac(1:nspecies) = react_state_in % rhoY(1:nspecies) * rhoInv

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
    neq         = nspecies + 1
    time_old    = time

    vodeVec(neq)     = eos_state % T

    if (iE == 1) then
        rhoe_init            = eos_state % e  *  eos_state % rho
        rhoedot_ext          = react_state_in % rhoedot_ext
        rhoydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies)
        vodeVec(1:nspecies)     = react_state_in % rhoY(:)
    else if (iE == 5) then
        h_init            = eos_state % h  
        hdot_ext          = react_state_in % rhohdot_ext / eos_state % rho
        ydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies) / eos_state % rho
        vodeVec(1:nspecies)  = react_state_in % rhoY(:) / eos_state % rho
    else
        rhoh_init            = eos_state % h  *  eos_state % rho
        rhohdot_ext          = react_state_in % rhohdot_ext 
        rhoydot_ext(1:nspecies) = react_state_in % rhoydot_ext(1:nspecies)
        vodeVec(1:nspecies)     = react_state_in % rhoY(:)
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

#ifdef MOD_REACTOR
    time = vodeTime
#endif

    if (verbose .ge. 1) then
       write(6,*) '......dvode done:', time
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
          if (ifail .eq. nspecies+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if
    end if

    if (istate > 0) then

       cost_value = DBLE(vodeiwork(12)) ! number of f evaluations

       if (iE == 1) then
           eos_state % rho               = sum(vodeVec(1:nspecies))
           rhoInv                        = 1.d0 / eos_state % rho
           eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies) * rhoInv
           eos_state % T                 = vodeVec(neq)
           eos_state % e                 = (rhoe_init  +  dt_react*rhoedot_ext) /eos_state % rho
           call eos_re(eos_state)
           rY_in(1:nspecies)                = vodeVec(1:nspecies)
           rX_in                         = eos_state % e
           rX_src_in                     = rhoedot_ext
           rY_src_in(1:nspecies)            = rhoydot_ext(1:nspecies)
       else if (iE == 5) then
           eos_state % p                 = pressureInit  
           eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies)
           eos_state % T                 = vodeVec(neq)
           eos_state % h                 = (h_init  +  dt_react*hdot_ext)
           call eos_ph(eos_state)
           rY_in(1:nspecies)                = vodeVec(1:nspecies) * eos_state % rho
           rX_in                         = eos_state % h
           rX_src_in                     = hdot_ext * eos_state % rho
           rY_src_in(1:nspecies)            = ydot_ext(1:nspecies) * eos_state % rho
       else
           eos_state % rho               = sum(vodeVec(1:nspecies))
           rhoInv                        = 1.d0 / eos_state % rho
           eos_state % massfrac(1:nspecies) = vodeVec(1:nspecies) * rhoInv
           eos_state % T                 = vodeVec(neq)
           eos_state % h                 = (rhoh_init  +  dt_react*rhohdot_ext) * rhoInv
           call eos_rh(eos_state)
           rY_in(1:nspecies)                = vodeVec(1:nspecies)
           rX_in                         = eos_state % h * eos_state % rho
           rX_src_in                     = rhohdot_ext
           rY_src_in(1:nspecies)            = rhoydot_ext(1:nspecies)
       end if

       rY_in(1+nspecies)      = eos_state % T

    else


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
       print *,'rhoY forcing',react_state_in%rhoydot_ext(1:nspecies)

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
          if (ifail .eq. nspecies+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if

       print *,'Final T',vodeVec(neq)
       print *,'Final rhoY',vodeVec(1:nspecies)

       call amrex_error('vode failed')

    end if
  end function react

  ! Original DVODE version
  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use chemistry_module, only : molecular_weight
    use eos_module

    implicit none
    integer,         intent(in)   :: neq, ipar(*)
    real(amrex_real), intent(in)  :: y(neq), time, rpar(*)
    real(amrex_real), intent(out) :: ydot(neq)
    integer          :: n
    real(amrex_real) :: rhoInv

    if (iE == 1) then
        eos_state % rho               = sum(y(1:nspecies))
        rhoInv                        = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspecies) = y(1:nspecies) * rhoInv
        eos_state % T                 = y(neq) ! guess
        eos_state % e = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        call eos_get_activity(eos_state)
    else if (iE == 5) then
        eos_state % massfrac(1:nspecies) = y(1:nspecies) 
        eos_state % T                 = y(neq) ! guess
        eos_state % h = h_init + (time - time_init) * hdot_ext
        eos_state % p                 = pressureInit  
        call eos_ph(eos_state)
        call eos_get_activity_h(eos_state)
    else
        eos_state % rho               = sum(y(1:nspecies))
        rhoInv                        = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspecies) = y(1:nspecies) * rhoInv 
        eos_state % T                 = y(neq) ! guess
        eos_state % h = (rhoh_init + (time - time_init) * rhohdot_ext) * rhoInv
        call eos_rh(eos_state)
        call eos_get_activity_h(eos_state)
    end if

    call ckwc(eos_state % T, eos_state % Acti, cdot)

    if (iE == 1) then
        ydot(neq)    = rhoedot_ext 
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%ei(n)*ydot(n)
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cv)
    else if (iE == 5) then
        ydot(neq)    = hdot_ext
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) / eos_state%rho + ydot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/eos_state%cp
    else
        ydot(neq)    = rhohdot_ext
        do n=1,nspecies
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cp)
    end if

  end subroutine f_rhs

  ! Original DVODE version
  subroutine f_jac(neq, npt, y, t, pd)
    use amrex_error_module

    implicit none
    integer,        intent(in)  :: neq, npt
    real(amrex_real),intent(in)  :: y(neq,npt), t
    real(amrex_real),intent(inout) :: pd(neq,neq)

    call amrex_error('DVODE version: Analytic Jacobian not yet implemented')

  end subroutine f_jac

!*** FINALIZE ROUTINES ***!
  subroutine reactor_close() bind(C, name="reactor_close")

    implicit none
    call destroy(eos_state)

    reactor_initialized = .false.
   
  end subroutine reactor_close


!*** SPECIFIC ROUTINES ***!
  function ok_to_react(state)

    use extern_probin_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: ok_to_react
    real(amrex_real)           :: rho

    ok_to_react = .true.

    rho = sum(state % rhoY)
    if (state % T   < react_T_min   .or. state % T   > react_T_max .or. &
        rho         < react_rho_min .or. rho         > react_rho_max) then

       ok_to_react = .false.

    endif

  end function ok_to_react

end module reactor_module
