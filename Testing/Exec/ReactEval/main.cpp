#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif
#include "mechanism.H"
#include <GPU_misc.H>

#include <PlotFileFromMF.H>
#include <PelePhysics.H>
#include <reactor.H>

#ifndef USE_RK64_PP
#ifdef USE_ARKODE_PP 
static std::string ODE_SOLVER = "ARKODE";
#else
static std::string ODE_SOLVER = "CVODE";
#endif
#else
static std::string ODE_SOLVER = "RK64";
#endif

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
#ifdef AMREX_USE_GPU
  amrex::sundials::MemoryHelper::Initialize(); /* TODO: this ideally (I think) will go in the amrex::Initialize */
#endif
  {

    Real strt_time = ParallelDescriptor::second();

    BL_PROFILE_VAR("main::main()", pmain);

    ParmParse pp;
    std::string fuel_name;
    pp.get("fuel_name", fuel_name);

    std::string pltfile;
    bool do_plt = false;
    if (pp.countval("plotfile")>0) {
      pp.get("plotfile",pltfile);
      do_plt = true;
    }

    /* react() function version */
    // 1 -> Array4 version of react()  (Default)
    // 2 -> 1d raw pointer version of react()
    std::string reactFormat = "Array4";
    int reactFunc;
    pp.query("reactFormat", reactFormat);
    if ( reactFormat == "Array4") {
       reactFunc = 1;
    } else if ( reactFormat == "1dArray" ) {
       reactFunc = 2;
    } else {
       Abort(" --> reactFormat can only be 'Array4' or '1dArray' !");
    }

    /* ODE inputs */
    ParmParse ppode("ode");
    int ode_ncells = 1;
    ppode.query("ode_ncells",ode_ncells); // number of cells to integrate per call, used only if reactFunc = 2

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

    Print() << "ODE solver: " << ODE_SOLVER << std::endl;
    Print() << "Type of reactor: " << (ode_iE == 1 ? "e (PeleC)" : "h (PeleLM)") << std::endl; // <---- FIXME
    Print() << "Fuel: " << fuel_name << ", Oxy: O2"  << std::endl;

    /* Mixture info */
    int fuel_idx   = -1;
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
#ifdef IC8H18_ID
    } else if (fuel_name == "IC8H18") {
      fuel_idx  = IC8H18_ID;
#endif
    }

    pele::physics::transport::InitTransport<
      pele::physics::PhysicsType::eos_type>()();

    /* Initialize reactor object inside OMP region, including tolerances */
    BL_PROFILE_VAR("main::reactor_info()", reactInfo);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      // Set ODE tolerances
#ifndef USE_ARKODE_PP
      //SetTolFactODE(rtol,atol);
#endif

#ifdef AMREX_USE_GPU
      reactor_init(ode_iE, ode_ncells);
#endif
    }
    BL_PROFILE_VAR_STOP(reactInfo);

    // -----------------------------------------------------------------------------
    // Define geometry
    // -----------------------------------------------------------------------------
    std::array<int,3> ncells = {AMREX_D_DECL(1,1,1)};
    if (pp.countval("ncells") == 1) {
      pp.get("ncells",ncells[0]);
      ncells = {AMREX_D_DECL(ncells[0],1,1)};
    }
    else if (pp.countval("ncells") >= AMREX_SPACEDIM) {
      Vector<int> nc(AMREX_SPACEDIM);
      pp.getarr("ncells",nc,0,AMREX_SPACEDIM);
      ncells = {AMREX_D_DECL(nc[0],nc[1],nc[2])};
    }
    else {
      Abort("ncells has to have length 1 or spacedim");
    }
    
    Box domain(IntVect(AMREX_D_DECL(0,0,0)),
               IntVect(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1)));


    RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                     {AMREX_D_DECL( 1.0, 1.0, 1.0)});

    int coord = 0;   

    Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};

    Geometry geom(domain, real_box, coord, is_periodic);

    Print() << "Integrating "<< domain.numPts() << " cells for: " << dt << " seconds with " << ndt << " substeps \n";

    // -----------------------------------------------------------------------------
    // Define BoxArray / Dmap
    // -----------------------------------------------------------------------------
    int max_grid_size = 16;
    pp.query("max_grid_size",max_grid_size);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    // -----------------------------------------------------------------------------
    // Create MultiFabs with no ghost cells
    // -----------------------------------------------------------------------------
    int num_grow = 0;
    MultiFab mf(ba, dm, NUM_SPECIES + 1, num_grow);
    MultiFab rY_source_ext(ba,dm,NUM_SPECIES, num_grow);
    MultiFab mfE(ba, dm, 1, num_grow);
    MultiFab rY_source_energy_ext(ba,dm,1,num_grow);
    MultiFab temperature(ba,dm,1,num_grow);
    MultiFab fctCount(ba,dm,1,num_grow);
    iMultiFab dummyMask(ba,dm,1,num_grow);
    dummyMask.setVal(1);

    FabArrayBase::mfiter_tile_size = IntVect(AMREX_D_DECL(1024,1024,1024));

    // -----------------------------------------------------------------------------
    // Initialize data
    // -----------------------------------------------------------------------------
    BL_PROFILE_VAR("main::initialize_data()", InitData);
    const auto geomdata = geom.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box = mfi.tilebox();

      Array4<Real> const& rY_a    = mf.array(mfi);
      Array4<Real> const& rYs_a   = rY_source_ext.array(mfi);
      Array4<Real> const& E_a     = mfE.array(mfi);
      Array4<Real> const& rE_a    = rY_source_energy_ext.array(mfi);

      ParallelFor(box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        initialize_data(i, j, k, fuel_idx, 
                        rY_a, rYs_a, E_a, rE_a,
                        geomdata);
      });
    }
    BL_PROFILE_VAR_STOP(InitData);

    // Print initial data if needed
    BL_PROFILE_VAR_NS("PlotFile",PlotFile);
    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = Concatenate(pltfile,0);
      PlotFileFromMF(mf,outfile);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    // -----------------------------------------------------------------------------
    // Set typical values
    // -----------------------------------------------------------------------------
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
#ifndef USE_ARKODE_PP
      /*if (use_typ_vals) {
        Print() << "Using user-defined typical values for the absolute tolerances of the ode solver.\n";
        Vector<double> typ_vals(NUM_SPECIES+1);
        ppode.getarr("typ_vals", typ_vals,0,NUM_SPECIES+1);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          typ_vals[i] = std::max(typ_vals[i],1.e-10);
        }
        SetTypValsODE(typ_vals);
      }*/
#endif
    }

    Print() << " \n STARTING THE ADVANCE \n";

    /* REACT */
    BL_PROFILE_VAR("Advance",Advance);
    BL_PROFILE_VAR_NS("React",ReactInLoop);
    BL_PROFILE_VAR_NS("Allocs",Allocs);
    BL_PROFILE_VAR_NS("Flatten",mainflatten);
#ifdef AMREX_USE_OMP
    const auto tiling = MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
    const bool tiling = TilingIfNotGPU();
#endif
    for ( MFIter mfi(mf,tiling); mfi.isValid(); ++mfi) {

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

      // -------------------------------------------------------------
      // Integration with Array4 react function
      if (reactFunc == 1) {
        Real time = 0.0;
        Real dt_incr = dt/ndt;

        Print() << " Integrating " << nc << " cells \n";
        /* Solve */
        BL_PROFILE_VAR_START(ReactInLoop);
        for (int ii = 0; ii < ndt; ++ii)
        {
          react(box, rhoY, frcExt, T,
                rhoE, frcEExt, fc, mask,
                dt_incr, time, ode_iE
#ifdef AMREX_USE_GPU
                , amrex::Gpu::gpuStream()
#endif
               );
          dt_incr =  dt/ndt;
        }
        BL_PROFILE_VAR_STOP(ReactInLoop);

      // -------------------------------------------------------------
      // Integration with 1dArray raw pointer react function
      }
#ifndef AMREX_USE_GPU 
      else if (reactFunc == 2) {

          // On GPU, integrate the entirely box at once
          // othewise use the user-input ode_ncells
#ifdef AMREX_USE_GPU
          ode_ncells    = nc;
#endif
          extra_cells = nc - (nc / ode_ncells) * ode_ncells;
          Print() << " Integrating " << nc << " cells with a "<< ode_ncells << " ode cell buffer ";
          Print() << "("<< extra_cells<<" extra cells) \n";

          BL_PROFILE_VAR_START(Allocs);
          int nCells               =  nc+extra_cells;
          auto tmp_vect            =  new Real[nCells * (NUM_SPECIES+1)];
          auto tmp_src_vect        =  new Real[nCells * NUM_SPECIES];
          auto tmp_vect_energy     =  new Real[nCells];
          auto tmp_src_vect_energy =  new Real[nCells];
          auto tmp_fc              =  new int[nCells];
          auto tmp_mask            =  new int[nCells];
          BL_PROFILE_VAR_STOP(Allocs);

          BL_PROFILE_VAR_START(mainflatten);
          ParallelFor(box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
            int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);

            box_flatten(icell, i, j, k, ode_iE, 
                        rhoY, frcExt, T, rhoE, frcEExt, mask,
                        tmp_vect, tmp_src_vect, tmp_vect_energy, tmp_src_vect_energy, tmp_mask);
          });

          for (int icell=nc; icell<nc+extra_cells; icell++) {
            for(int sp=0; sp<NUM_SPECIES; sp++) {
              tmp_vect[icell*(NUM_SPECIES+1)+sp]     = rhoY(0,0,0,sp);
              tmp_src_vect[icell*NUM_SPECIES+sp]     = frcExt(0,0,0,sp);
            }
            tmp_vect[icell*(NUM_SPECIES+1)+NUM_SPECIES] = T(0,0,0);
            tmp_vect_energy[icell]                      = rhoE(0,0,0); 
            tmp_src_vect_energy[icell]                  = frcEExt(0,0,0);
            tmp_mask[icell]                             = mask(0,0,0);
          }
          BL_PROFILE_VAR_STOP(mainflatten);

          /* Solve */
          BL_PROFILE_VAR_START(ReactInLoop);
          for(int i = 0; i < nCells; i+=ode_ncells) {
             tmp_fc[i] = 0;
             Real time = 0.0;
             Real dt_incr = dt/ndt;
             for (int ii = 0; ii < ndt; ++ii) {
               tmp_fc[i] += react(&tmp_vect[i*(NUM_SPECIES+1)], &tmp_src_vect[i*NUM_SPECIES],
                                  &tmp_vect_energy[i], &tmp_src_vect_energy[i],
                                  dt_incr,time);
               dt_incr =  dt/ndt;
               for (int ic = i+1; ic < i+ode_ncells ; ++ic) {
                  tmp_fc[ic] = tmp_fc[i];
               }
             }
          }
          BL_PROFILE_VAR_STOP(ReactInLoop);
      
          BL_PROFILE_VAR_START(mainflatten);
          ParallelFor(box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
            int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
            box_unflatten(icell, i, j, k, ode_iE,
                          rhoY, T, rhoE, frcEExt, fc,
                          tmp_vect, tmp_vect_energy, tmp_fc, dt);
          });
          BL_PROFILE_VAR_STOP(mainflatten);

          delete[] tmp_vect;
          delete[] tmp_src_vect;
          delete[] tmp_vect_energy;
          delete[] tmp_src_vect_energy;
          delete[] tmp_fc;
          delete[] tmp_mask;
      }
#endif
    }
    BL_PROFILE_VAR_STOP(Advance);

    {
      Vector<double> typ_vals(NUM_SPECIES+1);
      Print() << "ode.typ_vals= ";
      for (int i = 0; i < NUM_SPECIES+1; ++i) {
        Print() << std::max(1.e-10,mf.max(i)) << " ";
      }
      Print() << std::endl;
    }

    if (do_plt) {
        BL_PROFILE_VAR_START(PlotFile);
        MultiFab out(mf.boxArray(),mf.DistributionMap(),mf.nComp()+1,mf.nGrow());
        MultiFab::Copy(out,mf,0,0,mf.nComp(),mf.nGrow());
        AMREX_ALWAYS_ASSERT(fctCount.boxArray()==mf.boxArray() && fctCount.nGrow()>=mf.nGrow());
        MultiFab::Copy(out,fctCount,0,mf.nComp(),1,mf.nGrow());
        PlotFileFromMF(out,Concatenate(pltfile,1));
        BL_PROFILE_VAR_STOP(PlotFile);
    }
    
    pele::physics::transport::CloseTransport<
      pele::physics::PhysicsType::eos_type>()();

    BL_PROFILE_VAR_STOP(pmain);

    Real run_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << " \n >> React::main() " << run_time << "\n\n";
  }
  Finalize();

  return 0;
}
