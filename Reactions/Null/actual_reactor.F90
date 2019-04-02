module reactor_module

  use network, only: nspec
  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine reactor_init(iE_in) bind(C, name="reactor_init") 

    use, intrinsic :: iso_c_binding

    integer(c_int),  intent(in   ) :: iE_in

    !nothing needed here

  end subroutine reactor_init


  subroutine reactor_close() bind(C, name="reactor_close")

    ! nothing needed here

  end subroutine reactor_close


  function react(rY_in,rY_src_in,rX_in,rX_src_in,P_in,dt_react,time,Init) bind(C, name="react") result(cost_value)
    
    use eos_module

    real(amrex_real),   intent(inout) :: rY_in(nspec+1),rY_src_in(nspec)
    real(amrex_real),   intent(inout) :: rX_in,rX_src_in,P_in
    real(amrex_real),   intent(inout) :: dt_react, time
    integer                           :: Init, cost_value

    cost_value = 0

  end function react


end module reactor_module
