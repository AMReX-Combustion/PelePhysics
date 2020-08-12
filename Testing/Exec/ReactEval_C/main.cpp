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

#include <main_F.H>
#include <PlotFileFromMF.H>
#include <EOS.H>
#include <main_K.H> 

#if defined(USE_SUNDIALS_PP)
// ARKODE or CVODE
  #include <Transport.H>
  #include <reactor.h>
#else
  #if defined(USE_RK64_PP)
    #include <Transport.H>
// Expl RK solver
    #include <reactor.h>
  #endif
#endif

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main::main()", pmain);

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //INITIAL TIME
    Real timer_init = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_init,IOProc);
    //TOTAL TIME
    Real timer_tot = 0.;
    //INITIALIZATION
    Real timer_initialize_stop = 0.;
    // ADVANCE
    Real timer_adv = 0.;
    Real timer_adv_stop = 0.;
    //Print
    Real timer_print = 0.;
    Real timer_print_stop = 0.;

    {

    int max_grid_size = 16;
    std::string pltfile("plt");
    std::string fuel_name="none";
    /* ODE inputs */
    int ode_ncells = 1;
    int ode_iE     = -1;
    int third_dim  = 1024;
    int ndt        = 1; 
    Real dt        = 1.e-5;
    /* Scaling */
    Real rtol      = 1e-10;
    Real atol      = 1e-10;
    int use_typ_vals = 0;

    {
        /* ParmParse from the inputs file */
        ParmParse pp;

        // domain size
        pp.query("max_grid_size",max_grid_size);

        // third dim
        pp.query("third_dim",third_dim);

        // Get name of fuel
        pp.get("fuel_name", fuel_name);
    }

    {
        /* ParmParse from the inputs file */
        ParmParse ppode("ode");

        // final time
        ppode.query("dt",dt);

        // time stepping
        ppode.query("ndt",ndt); 

        ppode.query("reactor_type",ode_iE);
        /* Select ODE type of energy employed */
        //1 for UV, 2 for HP
        //   1 = Internal energy
        //   anything else = enthalpy (PeleLM restart)

        // nb of cells to integrate
        ppode.query("ode_ncells",ode_ncells);

        /* Additional scaling queries */
        ppode.query("rtol",rtol);
        ppode.query("atol",atol);
        ppode.query("use_typ_vals",use_typ_vals);

    }

    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);

    /* PRINT ODE INFO */
    amrex::Print() << "ODE solver: ";
#ifdef USE_SUNDIALS_PP
#ifdef USE_ARKODE_PP 
    amrex::Print() << "Using ARKODE (impl/expl solver)";
#else
    amrex::Print() << "Using CVODE (implicit solver)";
#endif
#else
#ifdef USE_RK64_PP
    amrex::Print()<<"Using custom RK64 (explicit solver)";
#endif
#endif
    amrex::Print() << std::endl;

    amrex::Print() << "Type of reactor: ";
        amrex::Print() << ode_iE;
    amrex::Print() << std::endl;

    amrex::Print() << "Fuel: ";
    amrex::Print() << fuel_name << ", Oxy: O2";
    amrex::Print() << std::endl;

    /* Mixture info */
    int fuel_idx   = -1;
    int oxy_idx    = -1;
    int bath_idx   = -1;
    if (fuel_name == "H2") {
        fuel_idx  = H2_ID;
#ifdef CH4_ID
    } else if (fuel_name == "CH4") {
        fuel_idx  = CH4_ID;
#endif
#ifdef NC12H26_ID
    } else if (fuel_name == "NC12H26") {
        fuel_idx  = NC12H26_ID;
#endif
    }
    oxy_idx   = O2_ID;
    bath_idx  = N2_ID;

    EOS::init();
    transport_init();

    BL_PROFILE_VAR("main::reactor_info()", reactInfo);

    /* Initialize reactor object */
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
    // Set ODE r/a tolerances
    SetTolFactODE(rtol,atol);
    // Set species-specific abs tolerances
    if (use_typ_vals) {
        amrex::Print() << "Using user-defined typical values for the absolute tolerances of the ode solver.\n";
        amrex::ParmParse pptv("ode");
        int nb_typ_vals = pptv.countval("typ_vals");
        if (nb_typ_vals != (NUM_SPECIES + 1)){
            printf("%d %d\n", nb_typ_vals, (NUM_SPECIES + 1));
            amrex::Abort("Not enough/too many typical values");
        }
        std::vector<double> typ_vals(nb_typ_vals);
        for (int i = 0; i < nb_typ_vals; ++i) {
            pptv.get("typ_vals", typ_vals[i],i);
        }
        SetTypValsODE(typ_vals);
    }
#ifdef USE_CUDA
    reactor_info(ode_iE, ode_ncells);
#else
    reactor_init(ode_iE, ode_ncells);
#endif
}

    BL_PROFILE_VAR_STOP(reactInfo);

    /* make domain and BoxArray */
    std::vector<int> npts(3,1);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        npts[i] = 2;
    }
    npts[1] = third_dim;

    amrex::Print() << "Integrating "<<npts[0]<< "x"<<npts[1]<< "x"<<npts[2]<< "  box for: ";
    amrex::Print() << dt << " seconds";
    amrex::Print() << std::endl;

    Box domain(IntVect(D_DECL(0,0,0)),
    IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    /* Additional defs to initialize domain */
    GpuArray<Real, AMREX_SPACEDIM>  dx;
    GpuArray<Real, AMREX_SPACEDIM>  plo;
    GpuArray<Real, AMREX_SPACEDIM>  phi;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
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
    iMultiFab dummyMask(ba,dm,1,0);

    BL_PROFILE_VAR("main::initialize_data()", InitData);

    /* INITIALIZE DATA */
#ifdef USE_CUDA
    int count_mf = 0;
#endif
    IntVect tilesize(D_DECL(1024,1024,1024));
    FabArrayBase::mfiter_tile_size = tilesize;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& gbox = mfi.tilebox();

        Array4<Real> const& rY_a    = mf.array(mfi);
        Array4<Real> const& rYs_a   = rY_source_ext.array(mfi);
        Array4<Real> const& E_a     = mfE.array(mfi);
        Array4<Real> const& rE_a    = rY_source_energy_ext.array(mfi);

        amrex::ParallelFor(gbox,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                initialize_data(i, j, k, fuel_idx, 
                                rY_a, rYs_a, E_a, rE_a,
                                dx, plo, phi);
        });

#ifdef USE_CUDA
        count_mf = count_mf + 1;
#endif
    }

    BL_PROFILE_VAR_STOP(InitData);

#ifdef USE_CUDA
    amrex::Print() << "That many boxes: " << count_mf<< "\n";
#endif

    timer_initialize_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_initialize_stop,IOProc);

    BL_PROFILE_VAR("main::PlotFileFromMF()", PlotFile);
    std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);
    BL_PROFILE_VAR_STOP(PlotFile);

    /* EVALUATE */
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    BL_PROFILE_VAR("main::Malloc()", Allocs);
    BL_PROFILE_VAR_STOP(Allocs);

    BL_PROFILE_VAR("main::React()", ReactInLoop);
    BL_PROFILE_VAR_STOP(ReactInLoop);

    BL_PROFILE_VAR("main::(un)flatten()", mainflatten);
    BL_PROFILE_VAR_STOP(mainflatten);

    BL_PROFILE_VAR("main::advance()", Advance);

    timer_adv = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_adv,IOProc);

    /* ADVANCE */
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& box  = mfi.tilebox();
        int ncells      = box.numPts();
        int extra_cells = 0;

        const auto len     = amrex::length(box);
        const auto lo      = amrex::lbound(box);

        auto const& rhoY    = mf.array(mfi);
        auto const& T       = mf.array(mfi, NUM_SPECIES);
        auto const& rhoE    = mfE.array(mfi);
        auto const& frcExt  = rY_source_ext.array(mfi);
        auto const& frcEExt = rY_source_energy_ext.array(mfi);
        auto const& fc      = fctCount.array(mfi);
        auto const& mask    = dummyMask.array(mfi);

#ifdef USE_CUDA
        cudaError_t cuda_status = cudaSuccess;
        ode_ncells    = ncells;
#else
        extra_cells = ncells - ncells / ode_ncells * ode_ncells; 
#endif

        amrex::Print() << " Integrating " << ncells << " cells with a "<<ode_ncells<< " ode cell buffer \n";
        amrex::Print() << "("<< extra_cells<<" extra cells) \n";

#ifndef CVODE_BOXINTEG
        /* ALLOCS */
        BL_PROFILE_VAR_START(Allocs);
        // rhoY,T
        amrex::Real *tmp_vect; 
        // rhoY_src_ext
        amrex::Real *tmp_src_vect;
        // rhoE/rhoH
        amrex::Real *tmp_vect_energy;
        amrex::Real *tmp_src_vect_energy;

        tmp_vect            =  new amrex::Real[(ncells+extra_cells)*(NUM_SPECIES+1)];
        tmp_src_vect        =  new amrex::Real[(ncells+extra_cells)*(NUM_SPECIES)];
        tmp_vect_energy     =  new amrex::Real[(ncells+extra_cells)];
        tmp_src_vect_energy =  new amrex::Real[(ncells+extra_cells)];
        BL_PROFILE_VAR_STOP(Allocs);

        /* Packing of data */
        BL_PROFILE_VAR_START(mainflatten);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
            for(int sp=0; sp<NUM_SPECIES; sp++) {
                tmp_vect[icell*(NUM_SPECIES+1)+sp]     = rhoY(i,j,k,sp);
                tmp_src_vect[icell*NUM_SPECIES+sp]     = frcExt(i,j,k,sp);
            }
            tmp_vect[icell*(NUM_SPECIES+1)+NUM_SPECIES] = rhoY(i,j,k,NUM_SPECIES);
            tmp_vect_energy[icell]                      = rhoE(i,j,k,0);
            tmp_src_vect_energy[icell]                  = frcEExt(i,j,k,0);
        });

        for (int icell=ncells; icell<ncells+extra_cells; icell++) {
            for(int sp=0; sp<NUM_SPECIES; sp++) {
                tmp_vect[icell*(NUM_SPECIES+1)+sp]     = rhoY(0,0,0,sp);
                tmp_src_vect[icell*NUM_SPECIES+sp]     = frcExt(0,0,0,sp);
            }
            tmp_vect[icell*(NUM_SPECIES+1)+NUM_SPECIES] = rhoY(0,0,0,NUM_SPECIES);
            tmp_vect_energy[icell]                      = rhoE(0,0,0,0); 
            tmp_src_vect_energy[icell]                  = frcEExt(0,0,0,0);
        }
        BL_PROFILE_VAR_STOP(mainflatten);
#endif


        /* Solve */
        BL_PROFILE_VAR_START(ReactInLoop);
#ifndef CVODE_BOXINTEG
        Real fc_tmp_lcl = 0.0;
        for(int i = 0; i < ncells+extra_cells; i+=ode_ncells) {
#endif
            Real time      = 0.0;
            Real dt_incr   = dt/ndt;
            for (int ii = 0; ii < ndt; ++ii) {
#ifndef CVODE_BOXINTEG
                fc_tmp_lcl = react(tmp_vect + i*(NUM_SPECIES+1), tmp_src_vect + i*NUM_SPECIES,
                                   tmp_vect_energy + i, tmp_src_vect_energy + i,
                                   dt_incr, time);
                //printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp + (NUM_SPECIES+1)]);
#else

#ifdef USE_CUDA
                react(box,
                      rhoY, frcExt, T,
                      rhoE, frcEExt,
                      fc, mask,
                      dt_incr, time,
                      ode_iE, amrex::Gpu::gpuStream());
#else
                react_1(box,
                      rhoY, frcExt, T,
                      rhoE, frcEExt,
                      fc, mask,
                      dt_incr, time);
#endif

#endif
                dt_incr =  dt/ndt;
            }
#ifndef CVODE_BOXINTEG
        }   
#endif
        BL_PROFILE_VAR_STOP(ReactInLoop);

#ifdef USE_CUDA
        cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());  
#endif


#ifndef CVODE_BOXINTEG
        /* Unpacking of data */
        BL_PROFILE_VAR_START(mainflatten);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
               for(int sp=0; sp<NUM_SPECIES; sp++) {
                   rhoY(i,j,k,sp) = tmp_vect[icell*(NUM_SPECIES+1)+sp];
               }
               rhoY(i,j,k,NUM_SPECIES) = tmp_vect[icell*(NUM_SPECIES+1) + NUM_SPECIES];
               rhoE(i,j,k,0)           = tmp_vect_energy[icell];
               fc(i,j,k,0)             = fc_tmp_lcl;
        });
        BL_PROFILE_VAR_STOP(mainflatten);
#endif

        printf("DONE");
#ifndef CVODE_BOXINTEG 
        delete(tmp_vect);
        delete(tmp_src_vect);
        delete(tmp_vect_energy);
        delete(tmp_src_vect_energy);
#endif

    }
    BL_PROFILE_VAR_STOP(Advance);

    timer_adv_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_adv_stop,IOProc);


    timer_print = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_print,IOProc);

    BL_PROFILE_VAR_START(PlotFile);
    outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);
    BL_PROFILE_VAR_STOP(PlotFile);

    timer_print_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_print_stop,IOProc);
    
    EOS::close();
    transport_close();

    }

    timer_tot = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_tot,IOProc);


    amrex::Print() << "Run Time total     (main())             = " << timer_tot             - timer_init    << "\n"
                   << "Run Time init      (initialize_data())  = " << timer_initialize_stop - timer_init    << "\n"
                   << "Run Time advance   (advance())          = " << timer_adv_stop        - timer_adv << "\n"
                   << "Run Time print plt (PlotFileFromMF())   = " << timer_print_stop      - timer_print << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
