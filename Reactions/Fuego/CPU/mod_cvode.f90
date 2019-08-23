module cvode_module

    interface

        subroutine reactor_init(cvode_iE, Ncells) bind(c,name='reactor_init')
            integer,          intent(in)    :: cvode_iE, Ncells
        end subroutine


        integer(c_int) function react(rY_in,rY_src_in, &
                        rX_in,rX_src_in,&
                        P_in,dt_react,time,Init) &
                        bind(c,name='react')
                use, intrinsic :: iso_c_binding
                real(c_double), dimension(*), intent(inout) :: rY_in
                real(c_double), dimension(*), intent(inout) :: rY_src_in
                real(c_double), dimension(*), intent(inout) :: rX_in
                real(c_double), dimension(*), intent(inout) :: rX_src_in
                real(c_double),  intent(inout) :: P_in
                real(c_double),  intent(inout) :: dt_react
                real(c_double),  intent(inout) :: time
                integer, intent(inout)        :: Init
        end function react

        subroutine reactor_close() bind(c,name='reactor_close')
        end subroutine

    end interface

end module cvode_module


!        integer(c_int) function react(rY_in, rY_src_in, rX_in, rX_src_in, &
!                        P_in, dt_react, time, Init) bind(c,name='react')
!            use amrex_fort_module, only : amrex_real
!            use, intrinsic :: iso_c_binding
!            implicit none
!            real(amrex_real), intent(inout) :: rY_in(*), rX_in(*)
!            real(amrex_real), intent(in   ) :: rY_src_in(*), rX_src_in(*)
!            real(amrex_real), intent(inout) :: P_in
!            real(amrex_real), intent(inout) :: dt_react
!            real(amrex_real), intent(in   ) :: time
!            integer,          intent(in   ) :: Init
!        end function react


