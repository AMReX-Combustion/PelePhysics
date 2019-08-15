#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <AMReX_GpuDevice.H>
#include <kernel.H>

template <typename L>
void ForMarc (Box const& box, int nc, L f) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    Gpu::ExecutionConfig ec;
    ec.numBlocks.x = AMREX_GPU_MAX_THREADS;
    ec.numBlocks.y = 1;
    ec.numBlocks.z = 1;
    ec.numThreads.x = nc;
    ec.numThreads.y = 1;
    ec.numThreads.z = 1;
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
      for (int icell = blockIdx.x, stride = blockDim.x; icell < ncells; icell += stride)
      {
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        int n = threadIdx.x;
        f(i,j,k,n);
      }
    });
    AMREX_GPU_ERROR_CHECK();
}

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {
      ParmParse pp;

      Vector<int> n_cells(BL_SPACEDIM,64);
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(n_cells[0]-1,n_cells[1]-1,n_cells[2]-1)));

      int max_size = 32;
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

      kinit();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("MARC");
        FArrayBox Q;
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& box = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          const auto& w    = wdot.array(mfi);
          int numPts = box.numPts();

          ForMarc(box, num_reac,
          [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
          {
            Real wtmp; W_spec_d(rho(i,j,k),temp(i,j,k),&(Y(i,j,k,0)),numPts,&wtmp);
            if (n<num_spec) {
              w(i,j,k,n) = wtmp;
            }
          });
        }
      }

      VisMF::Write(wdot,"WDOT");

    }

    Finalize();

    return 0;
}

