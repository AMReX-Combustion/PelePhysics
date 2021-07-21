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

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
       constexpr std::streamsize bl_ignore_max{100000};
           is.ignore(bl_ignore_max, '\n');
}

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

    /* initial data */
    // Can either generate single level data
    // or read in data from a checkpoint like file
    int initFromChk = 0;
    pp.query("initFromFile",initFromChk);
    std::string chkfile = "";
    if (initFromChk) {
       pp.query("initFile",chkfile);
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
      SetTolFactODE(rtol,atol);
#endif

      reactor_init(ode_iE, ode_ncells);
    }
    BL_PROFILE_VAR_STOP(reactInfo);

    // -----------------------------------------------------------------------------
    // Initialize geom/data
    // When initFromChk = 0, default is single level data provided by initialize_data
    // with a -1:1 unit length realbox in each dir and cell count provided by the user
    // otherwise number of levels, geom and all are read in from the chkfile
    // -----------------------------------------------------------------------------
    int finest_level = 0;
    Vector<Geometry> geoms;
    Vector<BoxArray> grids;
    Vector<DistributionMapping> dmaps;

    if ( initFromChk == 0 ) {
        // -----------------------------------------------------------------------------
        // Resize vectors
        // -----------------------------------------------------------------------------
        geoms.resize(finest_level+1);
        grids.resize(finest_level+1);
        dmaps.resize(finest_level+1);

        // -----------------------------------------------------------------------------
        // Define the geometry
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

        geoms[0] = Geometry(domain, real_box, coord, is_periodic);

        Print() << "Integrating "<< domain.numPts() << " cells for: " << dt << " seconds with " << ndt << " substeps \n";

        // -----------------------------------------------------------------------------
        // Define BoxArray / Dmap
        // -----------------------------------------------------------------------------
        int max_grid_size = 16;
        pp.query("max_grid_size",max_grid_size);
        grids[0] = BoxArray(domain);
        grids[0].maxSize(max_grid_size);
        dmaps[0] = DistributionMapping(grids[0], ParallelDescriptor::NProcs());
    } else {
        // Read chkfile header to get the geometry/BAs info

        std::string File(chkfile + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);
        std::string line, word;

        //--------------------------------------------
        // General info
        //--------------------------------------------
        Real crse_dt = 0.0;
        std::getline(is, line);  // Dummy title
        is >> finest_level;      // Finest level
        GotoNextLine(is);
        is >> crse_dt;         // Coarse level dt
        GotoNextLine(is);

        Print() << "  Warning: dt and ndt are overwritten when using data from a checkfile ! \n";
        ndt = 1;
        dt = crse_dt;

        //--------------------------------------------
        // Geometry
        //--------------------------------------------
        Real prob_lo[AMREX_SPACEDIM];
        Real prob_hi[AMREX_SPACEDIM];
        // Low coordinates of domain bounding box
        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while(lis >> word)
            {
                prob_lo[i++] = std::stod(word);
            }
        }

        // High coordinates of domain bounding box
        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while(lis >> word)
            {
                prob_hi[i++] = std::stod(word);
            }
        }
        int coord = 0;   
        Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};
        RealBox domainSize(prob_lo, prob_hi);
        Box domain;            // Read domain Box
        is >> domain;
        GotoNextLine(is);

        // -----------------------------------------------------------------------------
        // Resize vectors
        // -----------------------------------------------------------------------------
        geoms.resize(finest_level+1);
        grids.resize(finest_level+1);
        dmaps.resize(finest_level+1);
        
        // -----------------------------------------------------------------------------
        // Define geoms, read BoxArray and define dmap
        // -----------------------------------------------------------------------------
        geoms[0] = Geometry(domain, domainSize, coord, is_periodic);
        for(int lev = 1; lev <= finest_level; ++lev) {
           geoms[lev] = amrex::refine(geoms[lev-1], 2);     // Assumes ref_ratio = 2
        }

        for(int lev = 0; lev <= finest_level; ++lev)
        {
            // read in level 'lev' BoxArray from Header
            BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);

            // Set vector entries
            grids[lev] = ba;
            dmaps[lev] = DistributionMapping(ba, ParallelDescriptor::NProcs());
        }

        for(int lev = 0; lev <= finest_level; ++lev)
        {
           Print() << "  Level " << lev 
                   << " integrating " << grids[lev].numPts() 
                   << " cells on " << grids[lev].size() 
                   << " boxes for " << dt/std::pow(2,lev) << " seconds \n";
        }

    }

    // -----------------------------------------------------------------------------
    // Create MultiFabs with no ghost cells
    // -----------------------------------------------------------------------------
    int num_grow = 0;
    Vector<MultiFab> mf(finest_level+1);
    Vector<MultiFab> rY_source_ext(finest_level+1);
    Vector<MultiFab> mfE(finest_level+1);
    Vector<MultiFab> rY_source_energy_ext(finest_level+1);
    Vector<MultiFab> fctCount(finest_level+1);
    Vector<iMultiFab> dummyMask(finest_level+1);
    for(int lev = 0; lev <= finest_level; ++lev)
    {
       mf[lev].define(grids[lev],dmaps[lev],NUM_SPECIES+1,num_grow);
       rY_source_ext[lev].define(grids[lev],dmaps[lev],NUM_SPECIES,num_grow);
       mfE[lev].define(grids[lev],dmaps[lev],1,num_grow);
       rY_source_energy_ext[lev].define(grids[lev],dmaps[lev],1,num_grow);
       fctCount[lev].define(grids[lev],dmaps[lev],1,num_grow);
       dummyMask[lev].define(grids[lev],dmaps[lev],1,num_grow);
       dummyMask[lev].setVal(1);
    }

    FabArrayBase::mfiter_tile_size = IntVect(AMREX_D_DECL(1024,1024,1024));

    // -----------------------------------------------------------------------------
    // Initialize data
    // -----------------------------------------------------------------------------
    BL_PROFILE_VAR("main::initialize_data()", InitData);
    if ( initFromChk == 0 ) {
       for(int lev = 0; lev <= finest_level; ++lev)
       {
          const auto geomdata = geoms[lev].data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& box = mfi.tilebox();

            Array4<Real> const& rY_a    = mf[lev].array(mfi);
            Array4<Real> const& rYs_a   = rY_source_ext[lev].array(mfi);
            Array4<Real> const& E_a     = mfE[lev].array(mfi);
            Array4<Real> const& rE_a    = rY_source_energy_ext[lev].array(mfi);

            ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              initialize_data(i, j, k, fuel_idx, ode_iE,
                              rY_a, rYs_a, E_a, rE_a,
                              geomdata);
            });
          }
       }
    } else {
       // Load the data from chkfile
       for(int lev = 0; lev <= finest_level; ++lev)
       {
          // Assuming here a PeleLM state: vel + rho + species + rhoh + T + rhoRT
          MultiFab stateIn(grids[lev],dmaps[lev],NUM_SPECIES+AMREX_SPACEDIM+4,0);
          MultiFab forcingIn(grids[lev],dmaps[lev],NUM_SPECIES+1,0);
          VisMF::Read(stateIn,
                      amrex::MultiFabFileFullPrefix(lev, chkfile, level_prefix, "OldState"));
          VisMF::Read(forcingIn,
                      amrex::MultiFabFileFullPrefix(lev, chkfile, level_prefix, "ChemForcing"));
          // Copy into our local data holders
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
             const Box& box = mfi.tilebox();

             auto const& rhoYs   = mf[lev].array(mfi);
             auto const& temp    = mf[lev].array(mfi,NUM_SPECIES);
             auto const& rhoH    = mfE[lev].array(mfi);
             auto const& force_Y = rY_source_ext[lev].array(mfi);
             auto const& force_E = rY_source_energy_ext[lev].array(mfi);
             auto const& rhoYs_in = stateIn.const_array(mfi,AMREX_SPACEDIM+1); 
             auto const& temp_in  = stateIn.const_array(mfi,AMREX_SPACEDIM+NUM_SPECIES+2); 
             auto const& rhoh_in  = stateIn.const_array(mfi,AMREX_SPACEDIM+NUM_SPECIES+1); 
             auto const& force_a = forcingIn.const_array(mfi);
             ParallelFor(box,
             [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                for (int n = 0; n < NUM_SPECIES; n++) {
                  rhoYs(i,j,k,n) = rhoYs_in(i,j,k,n) * 1.0e-3;    // with MKS -> CGS conversion
                  force_Y(i,j,k,n) = force_a(i,j,k,n) * 1.0e-3;
                }
                temp(i,j,k) = temp_in(i,j,k);
                rhoH(i,j,k) = rhoh_in(i,j,k) * 10.0;              // with MKS -> CGS conversion
                force_E(i,j,k) = force_a(i,j,k,NUM_SPECIES) * 10.0;   
             });
          }
       }
    }
    BL_PROFILE_VAR_STOP(InitData);

    // Print initial data if needed
    BL_PROFILE_VAR_NS("PlotFile",PlotFile);
    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = Concatenate(pltfile,0);
      Vector<int> isteps(finest_level + 1, 0);
      Vector<IntVect> refRatios(finest_level, {AMREX_D_DECL(2,2,2)});
      Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
         plt_VarsName.push_back("SPEC"+std::to_string(k));
      }
      plt_VarsName.push_back("TEMP");

      amrex::WriteMultiLevelPlotfile(outfile, finest_level + 1, GetVecOfConstPtrs(mf),
                                     plt_VarsName, geoms, 0.0, isteps, refRatios);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    // -----------------------------------------------------------------------------
    // Set typical values
    // -----------------------------------------------------------------------------
#ifndef USE_ARKODE_PP

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
#ifndef USE_ARKODE_PP
      if (use_typ_vals) {
        Print() << "Using user-defined typical values for the absolute tolerances of the ode solver.\n";
        Vector<double> typ_vals(NUM_SPECIES+1);
        ppode.getarr("typ_vals", typ_vals,0,NUM_SPECIES+1);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          typ_vals[i] = std::max(typ_vals[i],1.e-10);
        }
        SetTypValsODE(typ_vals);
      }
#endif
    }

#endif

    Print() << " \n STARTING THE ADVANCE \n";

    /* REACT */
    BL_PROFILE_VAR_NS("React",ReactInLoop);
    BL_PROFILE_VAR_NS("Allocs",Allocs);
    BL_PROFILE_VAR_NS("Flatten",mainflatten);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
       BL_PROFILE_VAR("Advance"+std::to_string(lev),Advance);
#ifdef AMREX_USE_OMP
       const auto tiling = MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
       const bool tiling = TilingIfNotGPU();
#endif
       for ( MFIter mfi(mf[lev],tiling); mfi.isValid(); ++mfi) {

          const Box& box  = mfi.tilebox();
          int nc          = box.numPts();
          int extra_cells = 0;

          const auto len     = length(box);
          const auto lo      = lbound(box);

          auto const& rhoY    = mf[lev].array(mfi);
          auto const& T       = mf[lev].array(mfi, NUM_SPECIES);
          auto const& rhoE    = mfE[lev].array(mfi);
          auto const& frcExt  = rY_source_ext[lev].array(mfi);
          auto const& frcEExt = rY_source_energy_ext[lev].array(mfi);
          auto const& fc      = fctCount[lev].array(mfi);
          auto const& mask    = dummyMask[lev].array(mfi);

          // -------------------------------------------------------------
          // Integration with Array4 react function
          if (reactFunc == 1) {
            Real time = 0.0;
            Real dt_lev = dt/std::pow(2,lev);
            Real dt_incr = dt_lev/ndt;

            Print() << "  [" << lev << "]" << " integrating " << nc << " cells \n";
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
              dt_incr =  dt_lev/ndt;
            }
            BL_PROFILE_VAR_STOP(ReactInLoop);

          // -------------------------------------------------------------
          // Integration with 1dArray raw pointer react function
          }
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
                          rhoY, frcExt, T, rhoE, frcEExt,
                          tmp_vect, tmp_src_vect, tmp_vect_energy, tmp_src_vect_energy);
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
               Real dt_lev = dt/std::pow(2,lev);
               Real dt_incr = dt_lev/ndt;
               for (int ii = 0; ii < ndt; ++ii) {
                  tmp_fc[i] += react(&tmp_vect[i*(NUM_SPECIES+1)], &tmp_src_vect[i*NUM_SPECIES],
                                     &tmp_vect_energy[i], &tmp_src_vect_energy[i],
                                     dt_incr,time,ode_iE, ode_ncells
#ifdef AMREX_USE_GPU
                                     , amrex::Gpu::gpuStream()
#endif
                  );

                  dt_incr =  dt_lev/ndt;
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
                            tmp_vect, tmp_vect_energy, tmp_fc[icell], dt);
            });
            BL_PROFILE_VAR_STOP(mainflatten);

            delete[] tmp_vect;
            delete[] tmp_src_vect;
            delete[] tmp_vect_energy;
            delete[] tmp_src_vect_energy;
            delete[] tmp_fc;
            delete[] tmp_mask;
          }
       }
       BL_PROFILE_VAR_STOP(Advance);
    }

    // TODO multilevel max.
    {
      Vector<double> typ_vals(NUM_SPECIES+1);
      Print() << "ode.typ_vals= ";
      for (int i = 0; i < NUM_SPECIES+1; ++i) {
        Print() << std::max(1.e-10,mf[0].max(i)) << " ";
      }
      Print() << std::endl;
    }

    if (do_plt) {
        BL_PROFILE_VAR_START(PlotFile);
        std::string outfile = Concatenate(pltfile,1);
        // TODO: add fct count to this output
        Vector<int> isteps(finest_level + 1, 0);
        Vector<IntVect> refRatios(finest_level, {AMREX_D_DECL(2,2,2)});
        Vector<std::string> plt_VarsName;
        for (int k = 0; k < NUM_SPECIES; ++k) {
           plt_VarsName.push_back("SPEC"+std::to_string(k));
        }
        plt_VarsName.push_back("TEMP");

        amrex::WriteMultiLevelPlotfile(outfile, finest_level + 1, GetVecOfConstPtrs(mf),
                                       plt_VarsName, geoms, 0.0, isteps, refRatios);
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
