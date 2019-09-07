#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <AMReX_GpuDevice.H>
#include <kernel.H>

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {
      ParmParse pp;

      Vector<int> n_cells(BL_SPACEDIM,256);
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(n_cells[0]-1,n_cells[1]-1,n_cells[2]-1)));

      int max_size = 64;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_spec = NUM_SPECIES;
      int num_reac = NUM_REACTIONS;

      DistributionMapping dm(ba);

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,num_spec,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("INIT");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          For(bx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
            for (int n=0; n<num_spec; ++n) {
              Y(i,j,k,n) = 1./num_spec;
            }
            temp(i,j,k) = 450;
            rho(i,j,k) = 0.75;
          });
        }
      }

      MultiFab wdot(ba,dm,num_spec,num_grow);
      wdot.setVal(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("COMPUTE_W");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& box   = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          const auto& w    = wdot.array(mfi);
          int numPts = box.numPts();

          //printf("That many pts: %d \n", numPts);

	  amrex::ParallelFor(box,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
	    W_spec_d(i, j, k, rho, temp, Y, w);
          });
        }
      }

      VisMF::Write(wdot,"WDOT");

    }

    Finalize();

    return 0;
}

