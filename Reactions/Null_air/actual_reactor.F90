module actual_reactor_module

  use amrex_fort_module, only : amrex_real
  use react_type_module

  implicit none

contains

  subroutine actual_reactor_init(iE_in)

    use, intrinsic :: iso_c_binding

    integer(c_int),  intent(in   ) :: iE_in

    !nothing needed here

  end subroutine actual_reactor_init


  subroutine actual_reactor_close()

    ! nothing needed here

  end subroutine actual_reactor_close


  function actual_ok_to_react(state)

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: actual_ok_to_react

    actual_ok_to_react = .true.

  end function actual_ok_to_react


  function actual_react_null(react_state_in, react_state_out, dt_react, time) result(stat)
    
    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: stat

    react_state_out = react_state_in
    stat % cost_value = 0.d0
    stat % reactions_succesful = .true.

  end function actual_react_null


  function actual_react(react_state_in, react_state_out, dt_react, time) result(stat)
    
    use eos_module

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: stat

    stat = actual_react_null(react_state_in, react_state_out, dt_react, time)

  end function actual_react


end module actual_reactor_module
