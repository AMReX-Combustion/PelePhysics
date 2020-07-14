#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <mechanism.h>
#include <EOS.H>
#include <AMReX_GpuDevice.H>

#include <PlotFileFromMF.H>
#ifdef USE_SUNDIALS_PP
#include <reactor.h>
#else
#include <reactor.H> 
#endif

using namespace amrex;

AMREX_GPU_DEVICE
inline
void
initialize_data(int i, int j, int k,  
		Array4<Real> const& mf,
		Array4<Real> const& temp,
		Array4<Real> const& rho,
                const Real* dx, 
		const Real* plo, 
		const Real* phi ) noexcept
{
  Real dTemp = 5.0;
  Real dRho = 0.005;
  Real y = plo[1] + (j+0.5)*dx[1];
  Real x = plo[0] + (i+0.5)*dx[0];
  Real pi = 3.1415926535897932;
  Real L[3];
  Real P[3];
  Real Y_lo[NUM_SPECIES];
  Real Y_hi[NUM_SPECIES];

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_lo[n] = 0.0;
    Y_hi[n] = 1.0 / NUM_SPECIES ;
  }
  Y_lo[0] = 1.0;

  // T, Yk, rho
  Real Tavg = 500;
  Real Ravg = 0.01;
#if ( AMREX_SPACEDIM == 1 )
  temp(i,j,k) = Tavg;
  rho(i,j,k) = Ravg;
#else
  temp(i,j,k) = Tavg + dTemp * std::sin(2.0*pi*y/P[1]);
  rho(i,j,k) = Ravg + dRho * std::sin(2.0*pi*y/P[1]);
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf(i,j,k,n) = Y_lo[n] + (Y_hi[n]-Y_lo[n]) * x / L[0];
  }
  // corr Yk
  Real dummy = 0.0;
  for (int n = 0; n < NUM_SPECIES-1; n++) {
    dummy = dummy + mf(i,j,k,n);
  }
  mf(i,j,k,NUM_SPECIES-1) = 1.0 - dummy;
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    amrex::Print() << " Initialization of EOS (CPP)... \n";
    EOS::init();
      
    ParmParse pp;
    int nc = 512;
    pp.query("nc",nc);
    std::vector<int> npts(3,nc);

    Box domain(IntVect(D_DECL(0,0,0)),
               IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

    GpuArray<Real,3> plo, phi, dx;
    for (int i=0; i<BL_SPACEDIM; ++i) {
      phi[i] = domain.length(i);
      dx[i] = (phi[i] - plo[i])/domain.length(i);
    }
    int max_size = 128;
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
                           initialize_data(i, j, k, Y_a, T_a, rho_a, dx.data(), plo.data(), phi.data());
                         });
    }

  
    bool write_output = pp.countval("plot_file") > 0;
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

    if (write_output) {
      std::string pltfile;
      pp.get("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(wdots,outfile);
    }
  }
  
  Finalize();

  return 0;
}
