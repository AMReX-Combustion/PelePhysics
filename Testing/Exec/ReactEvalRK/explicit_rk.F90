module explicit_chemistry_module

contains

    subroutine explicit_react(nc,nspecies,        &
            react_vect,      &
            spec_src_vect,   &
            energy_vect,     &
            energy_src_vect, &
            dt_react,nsubsteps_guess,    &
            nsubsteps_min,nsubsteps_max, & 
            nsubsteps_final,errtol) bind(C, name="explicit_react")


        use chemistry_module, only : molecular_weight
        use rk_params_module
        use fuego_chemistry, only : CKWC,CKCVBS,CKUMS,CKYTCR

        implicit none

        integer           ::    nc,   nspecies
        double precision  ::    react_vect(nspecies+1,nc) !species and temperature
        double precision  ::    spec_src_vect(nspecies,nc)
        double precision  ::    energy_vect(nc)
        double precision  ::    energy_src_vect(nc)
        double precision  ::    dt_react
        double precision  ::    dt_guess
        double precision  ::    errtol
        integer           ::    nsubsteps_guess,nsubsteps_min,nsubsteps_max
        integer           ::    nsubsteps_final


        integer           :: cellid,spec,ns,stage,steps
        integer           :: TEMP_INDX
        double precision  :: soln_reg(nspecies+1,nc)
        double precision  :: carryover_reg(nspecies+1,nc)
        double precision  :: error_reg(nspecies+1,nc)
        double precision  :: rhs(nspecies+1,nc)
        double precision  :: rho
        double precision  :: y(nspecies),c(nspecies),eint(nspecies),cv

        double precision :: updt_time
        double precision :: dt_rk,dt_rk_min,dt_rk_max

        !Euler explicit for testing
        !integer,parameter :: nrkstages=1
        !double precision,parameter,dimension(nrkstages)  :: rkcoeffs=(/1.d0/)

        !===============================================================
        dt_rk     = dt_react/nsubsteps_guess
        dt_rk_min = dt_react/nsubsteps_max
        dt_rk_max = dt_react/nsubsteps_min
        !===============================================================

        !write(6,*) "dt_rk_guess:",dt_rk
        !write(6,*) "============================"
        !write(6,*) "============================"
        !flush(6)

        steps=0
        TEMP_INDX     = nspecies+1
        soln_reg      = react_vect

        do while(updt_time .lt. dt_react)

                carryover_reg = soln_reg
                error_reg     = 0.d0
                rhs           = 0.d0
                updt_time     = updt_time+dt_rk
                steps         = steps+1

                !write(6,*) "updt_time,dt_rk:",updt_time,dt_rk,dt_rk_min,dt_rk_max,dt_react
                !flush(6)

                do stage=1,rk64_stages


                        do cellid=1,nc

                                rho = sum(soln_reg(1:nspecies,cellid))
                                y   = soln_reg(1:nspecies,cellid)/rho

                                call CKYTCR(rho,soln_reg(TEMP_INDX,cellid),y,c)
                                call CKWC(soln_reg(TEMP_INDX,cellid),c,rhs(1:nspecies,cellid))
                                call CKUMS(soln_reg(TEMP_INDX,cellid),eint)
                                call CKCVBS(soln_reg(TEMP_INDX,cellid),y,cv)

                                rhs(1:nspecies,cellid) = rhs(1:nspecies,cellid)*molecular_weight(1:nspecies) &
                                + spec_src_vect(1:nspecies,cellid)

                                rhs(TEMP_INDX,cellid)=energy_src_vect(cellid)
                                do ns=1,nspecies
                                        rhs(TEMP_INDX,cellid)=rhs(TEMP_INDX,cellid)-rhs(ns,cellid)*eint(ns)
                                enddo
                                rhs(TEMP_INDX,cellid)=rhs(TEMP_INDX,cellid)/(rho*cv)

                        enddo

                        !time stepping
                        error_reg(:,:) = error_reg(:,:)     + err_rk64(stage)  *dt_rk*rhs(:,:)
                        soln_reg(:,:)  = carryover_reg(:,:) + alpha_rk64(stage)*dt_rk*rhs(:,:)
                        carryover_reg(:,:) = soln_reg(:,:)  + beta_rk64(stage) *dt_rk*rhs(:,:)

                enddo !stage loop

                call adapt_timestep(nc,nspecies,error_reg,dt_rk,dt_rk_min,dt_rk_max,errtol)

                !write(6,*)"dt_rk:",dt_rk,dt_rk_min,dt_rk_max
                !write(6,*)"================================"
                !flush(6)

        enddo !substep loop
        react_vect=soln_reg

        nsubsteps_final=steps
        !write(6,*)"No: of chemistry substeps:",steps
        !write(6,*)"================================"
        !flush(6)


    end subroutine explicit_react

    subroutine adapt_timestep(nc,nspecies,error_reg,dt_rk4,dt_rk4_min,dt_rk4_max,tol)

        implicit none

        integer           :: nc,nspecies
        double precision  :: error_reg(nspecies+1,nc)
        double precision  :: dt_rk4, dt_rk4_max, dt_rk4_min, tol

        integer :: cellid
        double precision :: max_err,change_factor
        double precision,parameter :: safety_fac=1e4
        double precision,parameter :: exp1=0.25
        double precision,parameter :: exp2=0.2
        double precision,parameter :: beta=1.d0

        double precision :: maxerr_cell


        max_err=tiny(max_err)

        do cellid=1,nc

        maxerr_cell=maxval(abs(error_reg(:,nc)))
        if(maxerr_cell .gt. max_err) then
            max_err=maxerr_cell
        endif

        enddo

        !write(6,*)"max_err:",max_err
        !flush(6)


        !chance to increase time step
        if(max_err .lt. tol) then
            !limit max_err,can't be 0
            change_factor=beta*(tol/max_err)**(exp1)
            dt_rk4=min(dt_rk4_max,dt_rk4*change_factor)

            !reduce time step (error is high!)
        else
            change_factor=beta*(tol/max_err)**(exp2)
            dt_rk4=max(dt_rk4_min,dt_rk4*change_factor)
        endif

    end subroutine adapt_timestep

end module explicit_chemistry_module
