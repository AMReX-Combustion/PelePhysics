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
#include <Transport.H>
#include <reactor.h>
#ifdef USE_ARKODE_PP 
static std::string ODE_SOLVER = "ARKODE";
#else
static std::string ODE_SOLVER = "CVODE";
#endif
#else
#if defined(USE_RK64_PP)
#include <Transport.H>
#include <reactor.h>
static std::string ODE_SOLVER = "RK64";
#endif
#endif

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    BL_PROFILE_VAR("main::main()", pmain);

    ParmParse pp;
    std::string fuel_name;
    pp.get("fuel_name", fuel_name);

    int max_grid_size = 16;
    pp.query("max_grid_size",max_grid_size);

    int ncells = 16;
    pp.query("ncells",ncells);
    
    /* ODE inputs */
    ParmParse ppode("ode");
    int ode_ncells = 1;
    ppode.query("ode_ncells",ode_ncells); // number of cells to integrate per call

    Real dt = 1.e-5;
    ppode.query("dt",dt);
    
    int ndt = 1;
    ppode.query("ndt",ndt); // number of solver calls per dt 
    
    int ode_iE = -1;
    ppode.query("reactor_type",ode_iE); // RHS type, 1: e (PeleC), !1: h (PeleLM)  <------ FIXME!
    
    Real rtol = 1e-10;
    ppode.query("rtol",rtol);
    
    Real atol = 1e-10;
    ppode.query("atol",atol);

    int use_typ_vals = 0;
    ppode.query("use_typ_vals",use_typ_vals);

    ParmParse ppa("amr");
    std::string pltfile("plt");
    ppa.query("plot_file",pltfile);

    /* PRINT ODE INFO */
    Print() << "ODE solver: " << ODE_SOLVER << std::endl;
    Print() << "Type of reactor: " << (ode_iE == 1 ? "e (PeleC)" : "h (PeleLM)") << std::endl; // <---- FIXME
    Print() << "Fuel: " << fuel_name << ", Oxy: O2"  << std::endl;

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

    /* Initialize reactor object inside OMP region, including tolerances */
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      // Set ODE tolerances
      SetTolFactODE(rtol,atol);

      // Set species-specific abs tolerances
      if (use_typ_vals) {
        Print() << "Using user-defined typical values for the absolute tolerances of the ode solver.\n";
        ParmParse pptv("ode");
        int nb_typ_vals = pptv.countval("typ_vals");
        if (nb_typ_vals != (NUM_SPECIES + 1)){
          printf("%d %d\n", nb_typ_vals, (NUM_SPECIES + 1));
          Abort("Not enough/too many typical values");
        }
        std::vector<double> typ_vals(nb_typ_vals);
        for (int i = 0; i < nb_typ_vals; ++i) {
          pptv.get("typ_vals", typ_vals[i],i);
        }
        SetTypValsODE(typ_vals);
      }
#ifdef USE_CUDA_SUNDIALS_PP
      reactor_info(ode_iE, ode_ncells);
#else
      reactor_init(ode_iE, ode_ncells);
#endif
    }
    BL_PROFILE_VAR_STOP(reactInfo);

    /* make domain and BoxArray */
    std::array<int,3> npts = {D_DECL(ncells,ncells,ncells)};

    Box domain(IntVect(D_DECL(0,0,0)),
               IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

    Print() << "Integrating "<< domain.numPts() << " cells for: " << dt << " seconds" << std::endl;

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    /* Additional defs to initialize domain */
    GpuArray<Real, AMREX_SPACEDIM>  plo = {D_DECL(0,0,0)};
    GpuArray<Real, AMREX_SPACEDIM>  dx  = {D_DECL(1,1,1)};
    GpuArray<Real, AMREX_SPACEDIM>  phi = {D_DECL(Real(domain.length(0)),
                                                  Real(domain.length(1)),
                                                  Real(domain.length(2)))};

    /* Create MultiFabs with no ghost cells */
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, NUM_SPECIES + 1, 0);
    MultiFab rY_source_ext(ba,dm,NUM_SPECIES,0);
    MultiFab mfE(ba, dm, 1, 0);
    MultiFab rY_source_energy_ext(ba,dm,1,0);
    MultiFab temperature(ba,dm,1,0);
    MultiFab fctCount(ba,dm,1,0);
    iMultiFab dummyMask(ba,dm,1,0);

    int count_mf = mf.local_size();
    FabArrayBase::mfiter_tile_size = IntVect(D_DECL(1024,1024,1024));

    /* INIT */
    BL_PROFILE_VAR("main::initialize_data()", InitData);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& gbox = mfi.tilebox();

      Array4<Real> const& rY_a    = mf.array(mfi);
      Array4<Real> const& rYs_a   = rY_source_ext.array(mfi);
      Array4<Real> const& E_a     = mfE.array(mfi);
      Array4<Real> const& rE_a    = rY_source_energy_ext.array(mfi);

      ParallelFor(gbox,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        initialize_data(i, j, k, fuel_idx, 
                        rY_a, rYs_a, E_a, rE_a,
                        dx, plo, phi);
      });
    }
    BL_PROFILE_VAR_STOP(InitData);

    BL_PROFILE_VAR("main::PlotFileFromMF()", PlotFile);
    std::string outfile = Concatenate(pltfile,0);
    PlotFileFromMF(mf,outfile);
    BL_PROFILE_VAR_STOP(PlotFile);

    Print() << " \n STARTING THE ADVANCE \n";

    /* REACT */
    BL_PROFILE_VAR("Advance",Advance);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box  = mfi.tilebox();
      int nc          = box.numPts();
      int extra_cells = 0;

      const auto len     = length(box);
      const auto lo      = lbound(box);

      auto const& rhoY    = mf.array(mfi);
      auto const& T       = mf.array(mfi, NUM_SPECIES);
      auto const& rhoE    = mfE.array(mfi);
      auto const& frcExt  = rY_source_ext.array(mfi);
      auto const& frcEExt = rY_source_energy_ext.array(mfi);
      auto const& fc      = fctCount.array(mfi);
      auto const& mask    = dummyMask.array(mfi);

#ifdef USE_CUDA_SUNDIALS_PP
      cudaError_t cuda_status = cudaSuccess;
      ode_ncells    = nc;
#else
      extra_cells = nc - (nc / ode_ncells) * ode_ncells; 
#endif

      Print() << " Integrating " << nc << " cells with a "<<ode_ncells<< " ode cell buffer \n";
      Print() << "("<< extra_cells<<" extra cells) \n";

#ifndef CVODE_BOXINTEG
      BL_PROFILE_VAR("Allocs",Allocs);
      Real *tmp_vect            =  new Real[(nc+extra_cells)*(NUM_SPECIES+1)];
      Real *tmp_src_vect        =  new Real[(nc+extra_cells)*(NUM_SPECIES)];
      Real *tmp_vect_energy     =  new Real[(nc+extra_cells)];
      Real *tmp_src_vect_energy =  new Real[(nc+extra_cells)];
      Real *tmp_fc              =  new Real[(nc+extra_cells)];
      BL_PROFILE_VAR_STOP(Allocs);

      BL_PROFILE_VAR("Flatten",mainflatten);
      ParallelFor(box,
      [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        for(int sp=0; sp<NUM_SPECIES; sp++) {
          tmp_vect[icell*(NUM_SPECIES+1)+sp]     = rhoY(i,j,k,sp);
          tmp_src_vect[icell*NUM_SPECIES+sp]     = frcExt(i,j,k,sp);
        }
        tmp_vect[icell*(NUM_SPECIES+1)+NUM_SPECIES] = rhoY(i,j,k,NUM_SPECIES);
        tmp_vect_energy[icell]                      = rhoE(i,j,k,0);
        tmp_src_vect_energy[icell]                  = frcEExt(i,j,k,0);
      });

      for (int icell=nc; icell<nc+extra_cells; icell++) {
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
      BL_PROFILE_VAR("React",ReactInLoop);
#ifndef CVODE_BOXINTEG
      for(int i = 0; i < nc+extra_cells; i+=ode_ncells) {
        tmp_fc[i] = 0.0;
        Real time = 0.0;
        Real dt_incr = dt/ndt;
        for (int ii = 0; ii < ndt; ++ii) {
          tmp_fc[i] += react(tmp_vect + i*(NUM_SPECIES+1), tmp_src_vect + i*NUM_SPECIES,
                             tmp_vect_energy + i, tmp_src_vect_energy + i,
                             dt_incr, time);
          dt_incr =  dt/ndt;
        }
      }
#else
      {      
        Real time = 0.0;
        Real dt_incr = dt/ndt;
        for (int ii = 0; ii < ndt; ++ii)
        {
#ifdef USE_CUDA_SUNDIALS_PP
          react(box,
                rhoY, frcExt, T,
                rhoE, frcEExt,
                fc, mask,
                dt_incr, time,
                ode_iE, Gpu::gpuStream());
#else
          react_1(box,
                  rhoY, frcExt, T,
                  rhoE, frcEExt,
                  fc, mask,
                  dt_incr, time);
#endif
          dt_incr =  dt/ndt;
        }
      }
#endif
      BL_PROFILE_VAR_STOP(ReactInLoop);

#ifndef CVODE_BOXINTEG
      BL_PROFILE_VAR_START(mainflatten);
      ParallelFor(box,
      [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        for(int sp=0; sp<NUM_SPECIES; sp++) {
          rhoY(i,j,k,sp)        = tmp_vect[icell*(NUM_SPECIES+1)+sp];
        }
        rhoY(i,j,k,NUM_SPECIES) = tmp_vect[icell*(NUM_SPECIES+1) + NUM_SPECIES];
        rhoE(i,j,k)             = tmp_vect_energy[icell];
        fc(i,j,k)               = tmp_fc[icell];
      });
      BL_PROFILE_VAR_STOP(mainflatten);

      delete(tmp_vect);
      delete(tmp_src_vect);
      delete(tmp_vect_energy);
      delete(tmp_src_vect_energy);
#endif

    }
    BL_PROFILE_VAR_STOP(Advance);
    
    outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass

    BL_PROFILE_VAR_START(PlotFile);
    MultiFab out(mf.boxArray(),mf.DistributionMap(),mf.nComp()+1,mf.nGrow());
    MultiFab::Copy(out,mf,0,0,mf.nComp(),mf.nGrow());
    AMREX_ALWAYS_ASSERT(fctCount.boxArray()==mf.boxArray() && fctCount.nGrow()>=mf.nGrow());
    MultiFab::Copy(out,fctCount,0,mf.nComp(),1,mf.nGrow());
    PlotFileFromMF(out,outfile);
    BL_PROFILE_VAR_STOP(PlotFile);

    EOS::close();
    transport_close();

    BL_PROFILE_VAR_STOP(pmain);
  }
  Finalize();

  return 0;
}
