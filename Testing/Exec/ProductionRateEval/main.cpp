#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <mechanism.h>
#include <EOS.H>
#include <AMReX_GpuDevice.H>
#include <GPU_misc.H>

#include <PlotFileFromMF.H>
#ifdef USE_SUNDIALS_PP
#include <reactor.h>
#else
#include <reactor.H> 
#endif

using namespace amrex;




std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {

  amrex::Print() << " Initialization of reactor... \n";
  int reactor_type = 2;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
  {
#ifdef USE_SUNDIALS_PP
    SetTolFactODE(relative_tol_chem,absolute_tol_chem);
#endif

    int ncells_chem = 1; // For CPU, this is how many cells integrated at once, for GPU this is ignored.
#ifdef USE_CUDA_SUNDIALS_PP
    reactor_info(&reactor_type,&ncells_chem);
#else
    reactor_init(&reactor_type,&ncells_chem);
#endif
  }

  amrex::Print() << " Initialization of EOS (CPP)... \n";
  EOS::init();

      
  std::vector<int> npts(3,32);

  Box domain(IntVect(D_DECL(0,0,0)),
             IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

  std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
  for (int i=0; i<BL_SPACEDIM; ++i) {
    phi[i] = domain.length(i);
    dx[i] = (phi[i] - plo[i])/domain.length(i);
  }
  ParmParse pp;
  int max_size = 32;
  pp.query("max_size",max_size);
  BoxArray ba(domain);
  ba.maxSize(max_size);

  int num_spec;
  num_spec = NUM_SPECIES;

  DistributionMapping dm{ba};

  int num_grow = 0;
  MultiFab mass_frac(ba,dm,num_spec,num_grow);
  MultiFab temperature(ba,dm,1,num_grow);
  MultiFab density(ba,dm,1,num_grow);

  IntVect tilesize(D_DECL(10240,8,32));
    
  int box_count =0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(mass_frac,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& gbox = mfi.tilebox();

    Array4<Real> const& Y_a    = mass_frac.array(mfi);
    Array4<Real> const& T_a    = temperature.array(mfi);
    Array4<Real> const& rho_a  = density.array(mfi);

    amrex::ParallelFor(gbox, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      initialize_data(i, j, k, Y_a, T_a, rho_a, dx, plo, phi);
    });
  }

  
  ParmParse ppa("amr");
  std::string pltfile("plt");  
  ppa.query("plot_file",pltfile);
  MultiFab wdots(ba,dm,num_spec,num_grow);
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& box = mfi.tilebox();

    const auto  mf     = mass_frac.array(mfi);
    const auto  temp   = temperature.array(mfi);
    const auto  rho    = density.array(mfi); 
    const auto  wdot   = wdots.array(mfi);

    ParallelFor(box,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real Yl[NUM_SPECIES];
      Real Wl[NUM_SPECIES];
      for (int n=0; n<NUM_SPECIES; ++n) Yl[n] = mf(i,j,k,n);
      EOS::RTY2WDOT(rho(i,j,k), temp(i,j,k), Yl, Wl);
      for (int n=0; n<NUM_SPECIES; ++n) wdot(i,j,k,n) = Wl[n];
    });

  }


  std::string outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
  PlotFileFromMF(wdots,outfile);

  }
  
  Finalize();

  return 0;
}
