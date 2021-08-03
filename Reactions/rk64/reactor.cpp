#include "reactor.H"

int reactor_init(int reactor_type,int Ncells)
{
    amrex::ParmParse pp("ode");
    pp.query("atol", absTol);
    pp.query("rtol", relTol);
    pp.query("rk64_nsubsteps_guess",rk64_nsubsteps_guess);
    pp.query("rk64_nsubsteps_min",    rk64_nsubsteps_min);
    pp.query("rk64_nsubsteps_max",    rk64_nsubsteps_max);
    return (0);
}

int
react(
        amrex::Real* rY_in,
        amrex::Real* rY_src_in,
        amrex::Real* rX_in,
        amrex::Real* rX_src_in,
        amrex::Real& dt_react,
        amrex::Real& time
        int reactor_type,
        int Ncells,
#ifdef AMREX_USE_GPU
        ,
        amrex::gpuStream_t stream)
#endif
{

#include "rk64params.H"

    amrex::Real time_init = time;
    amrex::Real time_out = time + dt_react;
    amrex::Real current_time = time_init;
    amrex::Real *soln_reg, *carryover_reg, *error_reg, *zero_reg, *rhs;
    amrex::Real dt_rk, dt_rk_min, dt_rk_max, change_factor;
    const amrex::Real tinyval = 1e-50;
    int neq=(NUM_SPECIES+1);


    //capture reactor type
    int captured_reactor_type    = reactor_type;
    int captured_nsubsteps_guess = rk64_nsubsteps_guess;
    int captured_nsubsteps_min   = rk64_nsusbteps_min;
    int captured_nsubsteps_max   = rk64_nsusbteps_max;
    
    int *nstepsvec;
    nstepsvec = new int[NCells]();

    amrex::ParallelFor(Ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept 
    {
        amrex::Real soln_reg[NUM_SPECIES+1];
        amrex::Real carryover_reg[NUM_SPECIES+1];
        amrex::Real error_reg[NUM_SPECIES+1];
        amrex::Real zero_reg[NUM_SPECIES+1];
        amrex::Real rhs[NUM_SPECIES+1];
        amrex::Real rYsrc[NUM_SPECIES];

        for(int sp=0;sp<neq;sp++)
        {
            soln_reg[sp]      = rY_in[icell*neq+sp];
            carryover_reg[sp] = soln_reg[sp];
            error_reg[sp]     = 0.0;
            rhs[sp]           = 0.0;
        }

        dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
        dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
        dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);
        
        current_time = time_init;
        amrex::Real rhoe_init[]={rX_in[icell]};
        amrex::Real rhoesrc_ext[]={rX_src_in[icell]};

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            rYsrc[sp]=rY_src_in[icell*NUM_SPECIES+sp];
        }

        int nsteps = 0;
        amrex::Real change_factor;
        while (current_time < time_out) 
        {
            for(int sp=0;sp<neq;sp++)
            {
                    eror_reg[sp]=0.0;
            }
            for (int stage = 0; stage < nstages_rk64; stage++) 
            {
                fKernelSpec(0, current_time-time_init, 
                        captured_reactor_type, soln_reg, rhs, 
                        rhoe_init, rhoesrc_ext, rYsrc);

                for (int i = 0; i < neq; i++) 
                {
                    error_reg[i]    += err_rk64[stage] * dt_rk * rhs[i];
                    soln_reg[i]      = carryover_reg[i] + alpha_rk64[stage] * dt_rk * rhs[i];
                    carryover_reg[i] = soln_reg[i] + beta_rk64[stage] * dt_rk * rhs[i];
                }
            }

            current_time += dt_rk;
            nsteps++;

            amrex::Real max_err = tinyval;
            for (int i = 0; i < neq; i++) 
            {
                if (fabs(error_reg[i]) > max_err) 
                {
                    max_err = fabs(error_reg[i]);
                }
            }

            if (max_err < captured_abstol) 
            {
                change_factor = beta * pow((data->errtol / max_err), exp1_rk64);
                dt_rk = std::min(dt_rk_max, dt_rk * change_factor);
            } 
            else 
            {
                change_factor = beta * pow((data->errtol / max_err), exp2_rk64);
                dt_rk = std::max(dt_rk_min, dt_rk * change_factor);
            }
        }
        nstepsvec[icell]=nsteps;
        //copy data back
        for(int sp=0;sp<neq;sp++)
        {
            rY_in[icell*neq+sp]=soln_reg[sp];
        }
        rX_in[icell] = rhoe_init[0] + dt_react*rhoesrc_init[0];
    });

#ifdef MOD_REACTOR
time = time_out;
#endif

    amrex::Real avgsteps=0.0;
    for(int i=0;i<Ncells;i++)
    {
        avgsteps += nstepsvec[i];
    } 

    return(int(avgsteps/amrex::Real(Ncells)));
}

    int
react(
        const amrex::Box& box,
        amrex::Array4<amrex::Real> const& rY_in,
        amrex::Array4<amrex::Real> const& rY_src_in,
        amrex::Array4<amrex::Real> const& T_in,
        amrex::Array4<amrex::Real> const& rEner_in,
        amrex::Array4<amrex::Real> const& rEner_src_in,
        amrex::Array4<amrex::Real> const& FC_in,
        amrex::Array4<int> const& mask,
        amrex::Real& dt_react,
        amrex::Real& time,
        const int& reactor_type
#ifdef AMREX_USE_GPU
        ,
        amrex::gpuStream_t stream
#endif
     )
{
    amrex::Real time_init = time;
    amrex::Real time_out = time + dt_react;
    amrex::Real current_time = time_init;
    const amrex::Real tinyval = 1e-50;
    int neq=(NUM_SPECIES+1);


    //capture reactor type
    int captured_reactor_type    = reactor_type;
    int captured_nsubsteps_guess = rk64_nsubsteps_guess;
    int captured_nsubsteps_min   = rk64_nsusbteps_min;
    int captured_nsubsteps_max   = rk64_nsusbteps_max;
    
    int *nstepsvec;
    int Ncells=box.numPts();
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    nstepsvec = new int[NCells]();

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept 
    {
        amrex::Real soln_reg[NUM_SPECIES+1];
        amrex::Real carryover_reg[NUM_SPECIES+1];
        amrex::Real error_reg[NUM_SPECIES+1];
        amrex::Real zero_reg[NUM_SPECIES+1];
        amrex::Real rhs[NUM_SPECIES+1];
        amrex::Real rYsrc[NUM_SPECIES];
    
        amrex::Real dt_rk, dt_rk_min, dt_rk_max, change_factor;
        amrex::Real rho;

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            soln_reg[sp]      = rY_in(i,j,k,sp);
            carryover_reg[sp] = soln_reg[sp];
            rho += rY_in(i,j,k,sp);
            error_reg[sp]     = 0.0;
            rhs[sp]           = 0.0;
        }
        rhs[NUM_SPECIES]       = 0.0;
        error_reg[NUM_SPECIES] = 0.0;
        amrex::Real rho_inv = 1.0 / rho;
        amrex::Real temp = T_in(i, j, k, 0);
    
        amrex::Real Enrg_loc = rEner_in(i,j,k,0)/rho;
        auto eos = pele::physics::PhysicsType::eos();
        if (captured_reactor_type == 1) 
        {
            eos.EY2T(Enrg_loc, mass_frac, temp);
        } 
        else 
        {
            eos.HY2T(Enrg_loc, mass_frac, temp);
        }
        soln_reg[NUM_SPECIES] = temp;

        dt_rk = dt_react / amrex::Real(captured_nsubsteps_guess);
        dt_rk_min = dt_react / amrex::Real(captured_nsubsteps_max);
        dt_rk_max = dt_react / amrex::Real(captured_nsubsteps_min);
        
        current_time = time_init;
        amrex::Real rhoe_init[]={rEner_in(i,j,k,0)};
        amrex::Real rhoesrc_ext[]={rEner_src_in(i,j,k,0)};

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            rYsrc[sp]=rY_src_in(i,j,k,sp);
        }

        int nsteps = 0;
        while (current_time < time_out) 
        {
            for(int sp=0;sp<neq;sp++)
            {
                    eror_reg[sp]=0.0;
            }
            for (int stage = 0; stage < nstages_rk64; stage++) 
            {
                fKernelSpec(0, current_time-time_init, 
                        captured_reactor_type, soln_reg, rhs, 
                        rhoe_init, rhoesrc_ext, rYsrc);

                for (int i = 0; i < neq; i++) 
                {
                    error_reg[i]    += err_rk64[stage] * dt_rk * rhs[i];
                    soln_reg[i]      = carryover_reg[i] + alpha_rk64[stage] * dt_rk * rhs[i];
                    carryover_reg[i] = soln_reg[i] + beta_rk64[stage] * dt_rk * rhs[i];
                }
            }

            current_time += dt_rk;
            nsteps++;

            amrex::Real max_err = tinyval;
            for (int i = 0; i < neq; i++) 
            {
                if (fabs(error_reg[i]) > max_err) 
                {
                    max_err = fabs(error_reg[i]);
                }
            }

            if (max_err < captured_abstol) 
            {
                change_factor = beta * pow((data->errtol / max_err), exp1_rk64);
                dt_rk = std::min(dt_rk_max, dt_rk * change_factor);
            } 
            else 
            {
                change_factor = beta * pow((data->errtol / max_err), exp2_rk64);
                dt_rk = std::max(dt_rk_min, dt_rk * change_factor);
            }
        }

        //copy data back
        int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
        nstepsvec[icell]=nsteps;
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            rY_in(i,j,k,sp) = soln_reg[sp];
        }
        T_in(i,j,k,0)     = soln_reg[NUM_SPECIES];
        rEner_in(i,j,k,0) = rhoe_init[0] + dt_react*rhoesrc_init[0];
        FC_in(i,j,k,0)    = nsteps;
    });

#ifdef MOD_REACTOR
time = time_out;
#endif

    amrex::Real avgsteps=0.0;
    for(int i=0;i<Ncells;i++)
    {
        avgsteps += nstepsvec[i];
    }   
    avgsteps = int(avgsteps/amrex::Real(Ncells));

    return(int(avgsteps/amrex::Real(Ncells)));

}
