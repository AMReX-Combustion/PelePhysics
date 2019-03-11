module reactor_module

  use, intrinsic :: iso_c_binding
  use amrex_fort_module, only : amrex_real
  use network
  use eos_module
  use react_type_module
  use actual_reactor_module

  implicit none

  logical, save, private :: reactor_initialized = .false.

contains

  !Original DVODE version
  subroutine reactor_init(iE) bind(C, name="reactor_init")

    implicit none
    integer(c_int),  intent(in   ) :: iE

    !$omp parallel
    call actual_reactor_init(iE)
    !$omp end parallel
    
    reactor_initialized = .true.

  end subroutine reactor_init


  subroutine reactor_close() bind(C, name="reactor_close")

    implicit none

    call actual_reactor_close()
    
    reactor_initialized = .false.

  end subroutine reactor_close


  function ok_to_react(state)

    implicit none
    type (react_t),intent(in) :: state
    logical                   :: ok_to_react

    ok_to_react = actual_ok_to_react(state)

  end function ok_to_react


  function react(react_state_in, react_state_out, dt_react, time)

    use amrex_error_module

    type(react_t),  intent(in    ) :: react_state_in
    type(react_t),  intent(inout ) :: react_state_out
    real(c_double), intent(in    ) :: dt_react, time
    type(reaction_stat_t)          :: react

    if (.not. reactor_initialized) then
       call amrex_error('reactor::react called before initialized')
    endif

    if ( ok_to_react(react_state_in) ) then

       react = actual_react(react_state_in, react_state_out, dt_react, time)

    else

       react = actual_react_null(react_state_in, react_state_out, dt_react, time)

    endif

  end function react

end module reactor_module
