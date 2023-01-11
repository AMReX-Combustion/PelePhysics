#include "AMReX_Reduce.H"
#include "ReactorBDF.H"
#include "ReactorBDFsolver.H"

namespace pele::physics::reactions {

    int
    ReactorBDF::init(int reactor_type, int /*ncells*/)
    {
        BL_PROFILE("Pele::ReactorBDF::init()");
        m_reactor_type = reactor_type;
        ReactorTypes::check_reactor_type(m_reactor_type);
        amrex::ParmParse pp("ode");
        pp.query("verbose", verbose);
        pp.query("atol", absTol);
        pp.query("nonlinear_iters",bdf_nonlinear_iters);
        pp.query("bdf_nsubsteps", bdf_nsubsteps);
        pp.query("bdf_gmres_restarts",bdf_gmres_restarts);
        pp.query("bdf_gmres_tol",bdf_gmres_tol);
        pp.query("bdf_gmres_precond",bdf_gmres_precond);
        pp.query("clean_init_massfrac", m_clean_init_massfrac);
        return (0);
    }

    int
    ReactorBDF::react(
        amrex::Real* rY_in,
        amrex::Real* rYsrc_in,
        amrex::Real* rX_in,
        amrex::Real* rX_src_in,
        amrex::Real& dt_react,
        amrex::Real& time,
        int ncells
#ifdef AMREX_USE_GPU
        ,
        amrex::gpuStream_t /*stream*/
#endif
        )
    {
        BL_PROFILE("Pele::ReactorBDF::react()");

        amrex::Real time_init = time;
        amrex::Real time_out = time + dt_react;

        // Copy to device
        amrex::Gpu::DeviceVector<amrex::Real> rY(ncells * (NUM_SPECIES + 1), 0);
        amrex::Gpu::DeviceVector<amrex::Real> rYsrc(ncells * NUM_SPECIES, 0);
        amrex::Gpu::DeviceVector<amrex::Real> rX(ncells, 0);
        amrex::Gpu::DeviceVector<amrex::Real> rX_src(ncells, 0);
        amrex::Real* d_rY = rY.data();
        amrex::Real* d_rYsrc = rYsrc.data();
        amrex::Real* d_rX = rX.data();
        amrex::Real* d_rX_src = rX_src.data();
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, rY_in, rY_in + ncells * (NUM_SPECIES + 1), d_rY);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, rYsrc_in, rYsrc_in + ncells * NUM_SPECIES,
            d_rYsrc);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, rX_in, rX_in + ncells, d_rX);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, rX_src_in, rX_src_in + ncells, d_rX_src);

        // capture variables
        const int captured_reactor_type = m_reactor_type;
        const int captured_nsubsteps = bdf_nsubsteps;
        const amrex::Real captured_abstol = absTol;
        const amrex::Real captured_gmres_tol = bdf_gmres_tol;
        const int captured_gmres_restarts=bdf_gmres_restarts;
        const int captured_nonlinear_iters=bdf_nonlinear_iters;
        const int captured_gmres_precond=bdf_gmres_precond;

        amrex::Gpu::DeviceVector<int> v_nsteps(ncells, 0);
        int* d_nsteps = v_nsteps.data();

        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
            amrex::Real soln_n[NUM_SPECIES + 1] = {0.0}; //at time level n
            amrex::Real soln[NUM_SPECIES + 1] = {0.0};
            amrex::Real dsoln[NUM_SPECIES + 1] = {0.0};  //newton_soln_k+1 - newton_soln_k
            amrex::Real dsoln0[NUM_SPECIES + 1] = {0.0}; //initial newton_soln_k+1 -newton_soln_k
            amrex::Real dsoln_n[NUM_SPECIES+1] = {0.0};  //newton_soln_k-newton_soln_n
            amrex::Real rYsrc_ext[NUM_SPECIES] = {0.0};
            amrex::Real ydot[NUM_SPECIES + 1] = {0.0};
            amrex::Real current_time = time_init;
            const int neq = (NUM_SPECIES + 1);

            for (int sp = 0; sp < neq; sp++) {
                soln_n[sp] = d_rY[icell * neq + sp];
                soln[sp] = soln_n[sp];
            }
    
            amrex::Real rho = 0.0;
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                 rho = rho + soln[sp];
            }
            amrex::Real temp = soln[NUM_SPECIES];
    
            amrex::Real massfrac[NUM_SPECIES] = {0.0};
            amrex::Real rhoinv = 1.0 / rho;
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                massfrac[sp] = soln[sp] * rhoinv;
            }   

            amrex::Real dt_bdf = dt_react / amrex::Real(captured_nsubsteps);

            amrex::Real rhoe_init[] = {d_rX[icell]};
            amrex::Real rhoesrc_ext[] = {d_rX_src[icell]};

            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                rYsrc_ext[sp] = d_rYsrc[icell * NUM_SPECIES + sp];
            }

            int nsteps = 0;
            int printflag=0;
            const int consP = (captured_reactor_type 
                               == ReactorTypes::h_reactor_type);
            auto eos = pele::physics::PhysicsType::eos();
            while (current_time < time_out) 
            {
                for(int nlit=0;nlit<bdf_nonlinear_iters;nlit++)
                {
                    for (int sp = 0; sp < neq; sp++)
                    {
                        dsoln0[sp] = dsoln[sp];
                        dsoln_n[sp] = soln[sp]-soln_n[sp];
                    }
                    
                    amrex::Real Jmat1d[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
                    amrex::Real mw[NUM_SPECIES] = {0.0};
                    get_mw(mw);
                    amrex::Real Jmat2d[NUM_SPECIES+1][NUM_SPECIES+1]={0.0};
                    eos.RTY2JAC(rho, temp, massfrac, Jmat1d, consP);
                    for (int i = 0; i < NUM_SPECIES; i++) 
                    {
                        for (int j = 0; j < NUM_SPECIES; j++) 
                        {
                            Jmat2d[i][j]=Jmat1d[j*(NUM_SPECIES+1) + i] * mw[i] / mw[j];
                        }
                        Jmat2d[i][NUM_SPECIES] = Jmat1d[NUM_SPECIES*(NUM_SPECIES+1)+i]*mw[i];
                        Jmat2d[NUM_SPECIES][i] = Jmat1d[i*(NUM_SPECIES+1)+NUM_SPECIES]/mw[i];
                    }
                    Jmat2d[NUM_SPECIES][NUM_SPECIES]=Jmat1d[(NUM_SPECIES+1)*(NUM_SPECIES+1)-1];
                    
                    utils::fKernelSpec<Ordering>(
                        0, 1, current_time - time_init, 
                        captured_reactor_type, soln, ydot,
                        rhoe_init, rhoesrc_ext, rYsrc_ext);

                    performgmres(Jmat2d,ydot,dsoln0,dsoln,dsoln_n,
                                 dt_bdf,captured_gmres_precond,
                                 captured_gmres_restarts,captured_gmres_tol,printflag);

                    for (int sp = 0; sp < neq; sp++) {
                        soln[sp] += dsoln[sp];
                    }
                }
                current_time += dt_bdf;
                nsteps++;

            }

            //ideally should be cost
            d_nsteps[icell] = nsteps;

            // copy data back
            for (int sp = 0; sp < neq; sp++) {
                d_rY[icell * neq + sp] = soln[sp];
            }
            d_rX[icell] = rhoe_init[0] + dt_react * rhoesrc_ext[0];
        });

#ifdef MOD_REACTOR
        time = time_out;
#endif

        const int avgsteps = amrex::Reduce::Sum<int>(
            ncells, [=] AMREX_GPU_DEVICE(int i) noexcept -> int { return d_nsteps[i]; },
            0);

        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, d_rY, d_rY + ncells * (NUM_SPECIES + 1), rY_in);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, d_rYsrc, d_rYsrc + ncells * NUM_SPECIES,
            rYsrc_in);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_rX, d_rX + ncells, rX_in);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, d_rX_src, d_rX_src + ncells, rX_src_in);

        //return cost here
        return (int(avgsteps / amrex::Real(ncells)));
    }

    int
    ReactorBDF::react(
        const amrex::Box& box,
        amrex::Array4<amrex::Real> const& rY_in,
        amrex::Array4<amrex::Real> const& rYsrc_in,
        amrex::Array4<amrex::Real> const& T_in,
        amrex::Array4<amrex::Real> const& rEner_in,
        amrex::Array4<amrex::Real> const& rEner_src_in,
        amrex::Array4<amrex::Real> const& FC_in,
        amrex::Array4<int> const& /*mask*/,
        amrex::Real& dt_react,
        amrex::Real& time
#ifdef AMREX_USE_GPU
        ,
        amrex::gpuStream_t /*stream*/
#endif
        )
    {
        BL_PROFILE("Pele::ReactorBDF::react()");

        amrex::Real time_init = time;
        amrex::Real time_out = time + dt_react;

        // capture variables
        const int captured_reactor_type = m_reactor_type;
        const int captured_nsubsteps = bdf_nsubsteps;
        const amrex::Real captured_abstol = absTol;
        const amrex::Real captured_gmres_tol = bdf_gmres_tol;
        const int captured_gmres_restarts=bdf_gmres_restarts;
        const int captured_nonlinear_iters=bdf_nonlinear_iters;
        const int captured_gmres_precond=bdf_gmres_precond;

        int ncells = static_cast<int>(box.numPts());
        const auto len = amrex::length(box);
        const auto lo = amrex::lbound(box);

        amrex::Gpu::DeviceVector<int> v_nsteps(ncells, 0);
        int* d_nsteps = v_nsteps.data();

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real soln_n[NUM_SPECIES + 1] = {0.0}; //at time level n
            amrex::Real soln[NUM_SPECIES + 1] = {0.0};
            amrex::Real dsoln[NUM_SPECIES + 1] = {0.0};  //newton_soln_k+1 - newton_soln_k
            amrex::Real dsoln0[NUM_SPECIES + 1] = {0.0}; //initial newton_soln_k+1 -newton_soln_k
            amrex::Real dsoln_n[NUM_SPECIES+1] = {0.0};  //newton_soln_k-newton_soln_n
            amrex::Real ydot[NUM_SPECIES + 1] = {0.0};
            amrex::Real rYsrc_ext[NUM_SPECIES] = {0.0};
            amrex::Real current_time = time_init;
            const int neq = (NUM_SPECIES + 1);

            amrex::Real rho = 0.0;
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                soln_n[sp] = rY_in(i, j, k, sp);
                soln[sp] = soln_n[sp];
                rho += rY_in(i, j, k, sp);
            }
            amrex::Real rho_inv = 1.0 / rho;
            amrex::Real massfrac[NUM_SPECIES] = {0.0};
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                massfrac[sp] = rY_in(i, j, k, sp) * rho_inv;
            }
            amrex::Real temp = T_in(i, j, k, 0);

            amrex::Real Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;
            auto eos = pele::physics::PhysicsType::eos();
            if (captured_reactor_type == ReactorTypes::e_reactor_type) {
                eos.REY2T(rho, Enrg_loc, massfrac, temp);
            } else if (captured_reactor_type == ReactorTypes::h_reactor_type) {
                eos.RHY2T(rho, Enrg_loc, massfrac, temp);
            } else {
                amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
            }
            soln_n[NUM_SPECIES] = temp;
            soln[NUM_SPECIES] = temp;

            amrex::Real dt_bdf = dt_react / amrex::Real(captured_nsubsteps);

            amrex::Real rhoe_init[] = {rEner_in(i, j, k, 0)};
            amrex::Real rhoesrc_ext[] = {rEner_src_in(i, j, k, 0)};

            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                rYsrc_ext[sp] = rYsrc_in(i, j, k, sp);
            }

            int nsteps = 0;
            int printflag=1;
            const int consP = (captured_reactor_type 
                               == ReactorTypes::h_reactor_type);
            while (current_time < time_out) 
            {
                for(int nlit=0;nlit<bdf_nonlinear_iters;nlit++)
                {
                    for (int sp = 0; sp < neq; sp++) {
                        dsoln0[sp] = dsoln[sp];
                        dsoln_n[sp] = soln[sp]-soln_n[sp];
                    }

                    amrex::Real Jmat1d[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)] = {0.0};
                    amrex::Real mw[NUM_SPECIES] = {0.0};
                    get_mw(mw);
                    amrex::Real Jmat2d[NUM_SPECIES+1][NUM_SPECIES+1]={0.0};
                    eos.RTY2JAC(rho, temp, massfrac, Jmat1d, consP);
                    for (int i = 0; i < NUM_SPECIES; i++) 
                    {
                        for (int j = 0; j < NUM_SPECIES; j++) 
                        {
                            Jmat2d[i][j]=Jmat1d[j*(NUM_SPECIES+1) + i] * mw[i] / mw[j];
                        }
                        Jmat2d[i][NUM_SPECIES] = Jmat1d[NUM_SPECIES*(NUM_SPECIES+1)+i]*mw[i];
                        Jmat2d[NUM_SPECIES][i] = Jmat1d[i*(NUM_SPECIES+1)+NUM_SPECIES]/mw[i];
                    }
                    Jmat2d[NUM_SPECIES][NUM_SPECIES]=Jmat1d[(NUM_SPECIES+1)*(NUM_SPECIES+1)-1];

                    utils::fKernelSpec<Ordering>(
                        0, 1, current_time - time_init, captured_reactor_type, soln, ydot,
                        rhoe_init, rhoesrc_ext, rYsrc_ext);

                    performgmres(Jmat2d,ydot,dsoln0,dsoln,dsoln_n,
                                 dt_bdf,captured_gmres_precond,
                                 captured_gmres_restarts,captured_gmres_tol,printflag);

                    for (int sp = 0; sp < neq; sp++) {
                        soln[sp] += dsoln[sp];
                    }
                }
                current_time += dt_bdf;
                nsteps++;
            }

            // copy data back
            int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
            d_nsteps[icell] = nsteps;
            rho = 0.0;
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                rY_in(i, j, k, sp) = soln[sp];
                rho += rY_in(i, j, k, sp);
            }
            rho_inv = 1.0 / rho;
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                massfrac[sp] = rY_in(i, j, k, sp) * rho_inv;
            }
            temp = soln[NUM_SPECIES];
            rEner_in(i, j, k, 0) = rhoe_init[0] + dt_react * rhoesrc_ext[0];
            Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;

            if (captured_reactor_type == ReactorTypes::e_reactor_type) {
                eos.REY2T(rho, Enrg_loc, massfrac, temp);
            } else if (captured_reactor_type == ReactorTypes::h_reactor_type) {
                eos.RHY2T(rho, Enrg_loc, massfrac, temp);
            } else {
                amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
            }
            T_in(i, j, k, 0) = temp;
            FC_in(i, j, k, 0) = nsteps;
        });

#ifdef MOD_REACTOR
        time = time_out;
#endif

        const int avgsteps = amrex::Reduce::Sum<int>(
            ncells, [=] AMREX_GPU_DEVICE(int i) noexcept -> int { return d_nsteps[i]; },
            0);
        return (int(avgsteps / amrex::Real(ncells)));
    }

} // namespace pele::physics::reactions
