#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include "mechanism.h"

using namespace amrex;

#include <Transport_F.H>
#include <main_F.H>
#include <PlotFileFromMF.H>
#ifdef USE_SUNDIALS_PP
#if USE_ARKODE_PP
        #include <actual_CARKODE.h>
#else
        #include <actual_Creactor.h>
#endif
#else
#include <actual_reactor.H> 
#endif

/**********************************/
int main (int   argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    Real timer_tot = amrex::second();
    Real timer_init = 0.;
    Real timer_advance = 0.;
    Real timer_print = 0.;
    Real timer_print_tmp = 0.;

    {

        timer_init = amrex::second();

        //set default inputs
        //=================
        Real dt = 1.e-5; 
        int ndt = 1; 
        int cvode_iE = 1;
        std::string fuel_name="H2";
        Real rk64_errtol=1e-13;
        int rk64_nsubsteps_guess=40;
        int rk64_nsubsteps_max=100;
        int rk64_nsubsteps_min=2;
        int rk64_nsubsteps_final;
        int max_grid_size = 16;
        int chem_integrator=1;

        //ARK parameters
        Real ark_rtol=1e-9;
        Real ark_atol=1e-9;
        int implicitflag=1;
        int use_erkode=0;


        std::string probin_file="probin";
        std::string pltfile("plt");
        int print_temp=0;

        /* CVODE inputs */
        int cvode_ncells = 1;
        int fuel_idx = -1;
        int oxy_idx = -1;
        int bath_idx = -1;
        
        //number of cells/points
        std::vector<int> npts(3,1);
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            npts[i] = 1;
        }

        //parse through inputs file
        //===========================
        {
            /* ParmParse from the inputs file */
            ParmParse pp;

            pp.query("dt",dt);
            pp.query("ndt",ndt); 
            pp.query("cvode_iE",cvode_iE);
            pp.query("cvode_ncells",cvode_ncells);
            pp.query("max_grid_size",max_grid_size);
            pp.query("probin_file",probin_file);

            pp.query("chem_integrator",chem_integrator);
            pp.query("rk64_errtol",rk64_errtol);
            pp.query("rk64_nsubsteps_guess",rk64_nsubsteps_guess);
            pp.query("rk64_nsubsteps_min",rk64_nsubsteps_min);
            pp.query("rk64_nsubsteps_max",rk64_nsubsteps_max);
            pp.get("fuel_name", fuel_name);
            pp.query("print_temp",print_temp);
            pp.queryarr("npts", npts, 0, AMREX_SPACEDIM);

            pp.query("implicitflag",implicitflag);
            pp.query("use_erkode",use_erkode);
            pp.query("ark_rtol",ark_rtol);
            pp.query("ark_atol",ark_atol);

        }

        //print some run parameters
        //=========================
        amrex::Print() << "Fuel: ";
        amrex::Print() << fuel_name << ", Oxy: O2";
        amrex::Print() << std::endl;

        amrex::Print() << "Integration method: ";
        if(chem_integrator==1)
        {
#if USE_ARKODE_PP 
            amrex::Print() << "ARK ODE";
#else
            amrex::Print() << "CVODE BDF";
#endif
            amrex::Print() << std::endl;
        }
        else
        {
            amrex::Print()<<"Using RK64 explicit\n"; 
        }

        amrex::Print() << "Type of reactor: ";
        amrex::Print() << cvode_iE;
        amrex::Print() << std::endl;


        /* take care of probin init to initialize problem */
        int probin_file_length = probin_file.length();
        std::vector<int> probin_file_name(probin_file_length);

        for (int i = 0; i < probin_file_length; i++)
            probin_file_name[i] = probin_file[i];

        //get fuel id (by default it is H2, 
        //most chemistries have it)
        //================================
        if (fuel_name == "H2") 
        {
            fuel_idx  = H2_ID;
        } 
        else if (fuel_name == "CH4") 
        {
#ifdef CH4_ID
            fuel_idx  = CH4_ID;
#endif
        } 
        else if (fuel_name == "NC12H26") 
        {
#ifdef NC12H26_ID
            fuel_idx  = NC12H26_ID;
#endif
        }
        oxy_idx   = O2_ID;
        bath_idx  = N2_ID;
        extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx,&cvode_iE);

        /* Initialize D/CVODE reactor */
#ifdef _OPENMP
#pragma omp parallel
#endif

#if USE_ARKODE_PP 
        reactor_init(&cvode_iE, &cvode_ncells,implicitflag,use_erkode,ark_rtol,ark_atol);
        //reactor_init(&cvode_iE, &cvode_ncells);
#else
        reactor_init(&cvode_iE, &cvode_ncells);
#endif

        amrex::Print() << "Integrating "<<npts[0]<< "x"<<npts[1]<< "x"<<npts[2]<< "  box for: ";
        amrex::Print() << dt << " seconds";
        amrex::Print() << std::endl;


        amrex::Print() << std::endl;

        Box domain(IntVect(D_DECL(0,0,0)),
                IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        /* Additional defs to initialize domain */
        std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
        for (int i=0; i<BL_SPACEDIM; ++i) {
            plo[i] = 0.0; //(i+1)*0.35;
            phi[i] = domain.length(i);
            dx[i] = (phi[i] - plo[i])/domain.length(i);
        }

        int Ncomp;
        Ncomp = NUM_SPECIES;

        /* Create MultiFabs with no ghost cells */
        DistributionMapping dm(ba);
        MultiFab mf(ba, dm, Ncomp+1, 0);
        MultiFab rY_source_ext(ba,dm,Ncomp,0);
        MultiFab mfE(ba, dm, 1, 0);
        MultiFab rY_source_energy_ext(ba,dm,1,0);
        MultiFab temperature(ba,dm,1,0);
        MultiFab fctCount(ba,dm,1,0);

        IntVect tilesize(D_DECL(10240,10240,32));

        /* INITIALIZE DATA */
#ifdef _OPENMP
#pragma omp parallel 
#endif
        for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi ){
            const Box& box = mfi.tilebox();
            initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                    BL_TO_FORTRAN_N_3D(mf[mfi],0),
                    BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
                    BL_TO_FORTRAN_N_3D(mfE[mfi],0),
                    BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
                    &(dx[0]), &(plo[0]), &(phi[0]));
        }
       

        timer_init = amrex::second() - timer_init; 

        timer_print = amrex::second();

        ParmParse ppa("amr");
        ppa.query("plot_file",pltfile);
        std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
        // Specs
        PlotFileFromMF(mf,outfile);

        timer_print = amrex::second() - timer_print;

        /* EVALUATE */
        amrex::Print() << " \n STARTING THE ADVANCE \n";

        timer_advance = amrex::second();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi )
        {
            const Box& box = mfi.tilebox();
            //int ncells = box.numPts();

            const auto&      rhoY     = mf.array(mfi);
            const auto&      frcExt   =  rY_source_ext.array(mfi);
            const auto&      rhoE     = mfE.array(mfi);
            const auto&      frcEExt  = rY_source_energy_ext.array(mfi);
            const auto&      fc       = fctCount.array(mfi);

            const auto len     = amrex::length(box);
            const auto lo      = amrex::lbound(box);

            int num_cells_in_box = len.x*len.y*len.z;

            /* Pack the data */
            // rhoY,T
            double state_vect[num_cells_in_box*(Ncomp+1)];
            // rhoY_src_ext
            double state_src_vect[num_cells_in_box*(Ncomp)];
            // rhoE/rhoH
            double state_vect_energy[num_cells_in_box];
            double state_src_vect_energy[num_cells_in_box];

            int indx_i[num_cells_in_box];
            int indx_j[num_cells_in_box];
            int indx_k[num_cells_in_box];

            Print()<<"Num cells in box:"<<num_cells_in_box<<"\t"<<len.x<<"\t"<<len.y<<"\t"<<len.z<<"\n";

            int nc = 0;
            for(int k = 0; k < len.z; ++k)
            {
                for(int j = 0; j < len.y; ++j) 
                {
                    for(int i = 0; i < len.x; ++i) 
                    {
                        /* Fill the vectors */
                        for (int sp=0;sp<Ncomp; sp++)
                        {
                            state_vect[nc*(Ncomp+1) + sp]   = rhoY(i+lo.x,j+lo.y,k+lo.z,sp);
                            state_src_vect[nc*Ncomp + sp]   = frcExt(i+lo.x,j+lo.y,k+lo.z,sp);
                        }
                        state_vect[nc*(Ncomp+1) + Ncomp]    = rhoY(i+lo.x,j+lo.y,k+lo.z,Ncomp);
                        state_vect_energy[nc]               = rhoE(i+lo.x,j+lo.y,k+lo.z,0);
                        state_src_vect_energy[nc]           = frcEExt(i+lo.x,j+lo.y,k+lo.z,0);
                        //
                        indx_i[nc] = i+lo.x;
                        indx_j[nc] = j+lo.y;
                        indx_k[nc] = k+lo.z;
                        //
                        nc = nc+1;
                    }
                }
            }

            if(chem_integrator==1)
            {
                // rhoY,T
                double tmp_vect[cvode_ncells*(Ncomp+1)];
                // rhoY_src_ext
                double tmp_src_vect[cvode_ncells*(Ncomp)];
                // rhoE/rhoH
                double tmp_vect_energy[cvode_ncells];
                double tmp_src_vect_energy[cvode_ncells];

                for(int i = 0; i < num_cells_in_box; i+=cvode_ncells)
                { 
                    /* Fill the vectors */
                    for(int l=0;l<cvode_ncells;l++)
                    {
                        for(int sp=0;sp<Ncomp;sp++)
                        {
                            tmp_vect[l*(Ncomp+1)+sp] = state_vect[(i+l)*(Ncomp+1)+sp]; 
                            tmp_src_vect[l*Ncomp+sp] = state_src_vect[(i+l)*Ncomp+sp];
                        }
                        tmp_vect[l*(Ncomp+1)+Ncomp]      = state_vect[(i+l)*(Ncomp+1)+Ncomp];
                        tmp_vect_energy[l]               = state_vect_energy[i+l];
                        tmp_src_vect_energy[l]           = state_src_vect_energy[i+l];
                    }

                    //
                    Real time = 0.0;
                    Real dt_incr =  dt/ndt;
                    for (int ii = 0; ii < ndt; ++ii) 
                    {
#ifdef USE_SUNDIALS_PP 
                        double cost = react(tmp_vect, tmp_src_vect,
                                tmp_vect_energy, tmp_src_vect_energy,
                                &dt_incr, &time);
#else
                        double pressure = 1013250.0;
                        double cost = react(tmp_vect, tmp_src_vect,
                                tmp_vect_energy, tmp_src_vect_energy,
                                &pressure,&dt_incr, &time);
#endif
                        dt_incr =  dt/ndt;
                        //printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp]);
                        if(print_temp)
                        {
                            PrintToFile("time_vs_temp")<<time<<"\t"<<tmp_vect[Ncomp]<<"\n";
                        }
                    }


                    for (int l = 0; l < cvode_ncells ; ++l)
                    {
                        for (int sp=0;sp<(Ncomp+1); sp++)
                        {
                            state_vect[(i+l)*(Ncomp+1)+sp] = tmp_vect[l*(Ncomp+1)+sp];
                        }
                        state_vect_energy[i+l]             = tmp_vect_energy[l];
                    }

                }

            }
            else
            {
                Real time = 0.0;
                Real dt_incr =  dt/ndt;
                int guess_substeps=rk64_nsubsteps_guess;
                for (int ii = 0; ii < ndt; ++ii) 
                {
                    explicit_react(&nc, &Ncomp, state_vect, state_src_vect,
                            state_vect_energy, state_src_vect_energy,
                            &dt_incr,&guess_substeps,&rk64_nsubsteps_min,
                            &rk64_nsubsteps_max,&rk64_nsubsteps_final,&rk64_errtol);

                    time=time+dt_incr;
                    guess_substeps=rk64_nsubsteps_final;
                    //printf("%14.6e %14.6e \n", time, state_vect[Ncomp]);
                    if(print_temp)
                    {
                        PrintToFile("time_vs_temp")<<time<<"\t"<<state_vect[Ncomp]<<"\n";
                    }
                }
            }

            for (int l = 0; l < num_cells_in_box; ++l)
            {
                for (int sp=0;sp<Ncomp; sp++)
                {
                    rhoY(indx_i[l],indx_j[l],indx_k[l],sp) = state_vect[l*(Ncomp+1) + sp];
                }
                rhoY(indx_i[l],indx_j[l],indx_k[l],Ncomp)  = state_vect[l*(Ncomp+1) + Ncomp];
                rhoE(indx_i[l],indx_j[l],indx_k[l],0)      = state_vect_energy[l];
            }
        }

        timer_advance = amrex::second() - timer_advance;


        timer_print_tmp = amrex::second();

        outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
        // Specs
        PlotFileFromMF(mf,outfile);

        timer_print = amrex::second() - timer_print_tmp + timer_print;

        extern_close();

    }

    timer_tot = amrex::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance, timer_print},
            ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
        << "Run Time init         = " << timer_init    << "\n"
        << "Run Time advance      = " << timer_advance << "\n"
        << "Run Time print plt    = " << timer_print << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
