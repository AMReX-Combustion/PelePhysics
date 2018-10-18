module actual_burner_module

  use amrex_fort_module, only : amrex_real
  use network
  use burn_type_module

contains

  subroutine actual_burner_init()

    ! Do nothing in this burner.

  end subroutine actual_burner_init


  subroutine actual_burner(state_in, state_out, dt, time)

    implicit none

    type (burn_t), intent(in)    :: state_in
    type (burn_t), intent(inout) :: state_out
    real(amrex_real), intent(in) :: dt, time

    ! Do nothing in this burner.

  end subroutine actual_burner

end module actual_burner_module
