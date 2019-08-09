#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include <GPU_misc.H>
#include <AMReX_GpuDevice.H>
#include <chemistry_file.H>

#include <main_F.H>
#include <PlotFileFromMF.H>

std::string inputs_name = "";

using namespace amrex;

extern "C" {

  void CUDART_CB foo (cudaStream_t stream, cudaError_t error, void* p)
  {
    int idx = *( (int *) p);
    printf("Index %d\n",idx);
  }

  void CUDA_CB foo1( void*  userData )
  {
    int idx = *( (int *)(userData));
    printf("Index: %d",idx);
  }
}

template <typename L>
void For (Box const& box, int nc, L f) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    //ec.numBlocks.y = nc;
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
      for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
           icell < ncells; icell += stride)
      {
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        for (int n = blockDim.y*blockIdx.y+threadIdx.y, ystride = blockDim.y*gridDim.y; n < nc; n += ystride) {
          f(i,j,k,n);
        }
      }
    });
    AMREX_GPU_ERROR_CHECK();
}
template <unsigned int blockSize>
__device__ void warpReduce(volatile int *sdata, unsigned int tid) {
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
  if (blockSize >=  8) sdata[tid] += sdata[tid + 4];
  if (blockSize >=  4) sdata[tid] += sdata[tid + 2];
  if (blockSize >=  2) sdata[tid] += sdata[tid + 1];
}

template <unsigned int blockSize>
__global__ void reduce6(int *g_idata, int *g_odata, unsigned int n) {
  extern __shared__ int sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + tid;
  unsigned int gridSize = blockSize*2*gridDim.x;
  sdata[tid] = 0;

  while (i < n) { sdata[tid] += g_idata[i] + g_idata[i+blockSize]; i += gridSize; }
  __syncthreads();

  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }

  if (tid < 32) warpReduce<blockSize>(sdata, tid);
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}



int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {
      unsigned int len=1e6;
      int *test = (int*) malloc(len*sizeof(int));
      for (unsigned int i=0; i<len; ++i) {
        test[i] = i*i;
      }
      int *dev_test;
      cudaMalloc( (void**)&dev_test, len*sizeof(int) );
      cudaMemcpy( dev_test, test, len, cudaMemcpyHostToDevice );
      
      dim3 grid, block;
      grid.x = 512; grid.y = 1; block.x = 8; block.y = 16;
      
      int valTot;
      reduce6<512><<<grid,block>>>(dev_test,&valTot,len);
      Print() << "Total is " << valTot << std::endl;

      ParmParse pp;
    
      std::string probin_file = "probin";
      pp.query("probin_file",probin_file);
      int probin_file_length = probin_file.length();
      std::vector<int> probin_file_name(probin_file_length);

      for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

      //int fuel_idx = FUEL_ID;
      int fuel_idx = 10;
      int oxy_idx  = OXY_ID;
      int bath_idx = BATH_ID;
      Print() << "fuel_idx: " << fuel_idx << std::endl;
      Print() << "oxy_idx: " << oxy_idx << std::endl;
      Print() << "bath_idx: " << bath_idx << std::endl;

      extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx);
    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 128;
      }
      npts[1] = 256;
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
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

      for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
 	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
      }

      ParmParse ppa("amr");
      std::string pltfile("plt");  
      ppa.query("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,0);
      PlotFileFromMF(temperature,outfile);

      MultiFab wdots(ba,dm,num_spec,num_grow);
      wdots.setVal(0);
      MultiFab Qmfab(ba,dm,num_reac,num_grow);
      Qmfab.setVal(0);

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
          const auto& wdot = wdots.array(mfi);

          int numPts = box.numPts();

#if 0
          Q.resize(box,NUM_REACTIONS);
          Elixir iQ = Q.elixir();
          auto iQa = Q.array();
#else
          const auto& iQa = Qmfab.array(mfi);
#endif

#if 1
          int nR = NUM_REACTIONS;
          For(box, nR,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
              {
                //Real Yloc[NUM_SPECIES];
                //for (int L=0; L<NUM_SPECIES; ++L) {
                //  Yloc[L] = Y(i,j,k,L);
                //}
                iQa(i,j,k,n) = Q_reac1_d(rho(i,j,k),temp(i,j,k),&(Y(i,j,k,0)),numPts,n);
                //iQa(i,j,k,n) = Q_reac_d(rho(i,j,k),temp(i,j,k),Yloc,n);
              });
#endif
#if 1
          For(box, num_spec,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
              {
                //Real qloc[NUM_REACTIONS];
                //for (int L=0; L<NUM_REACTIONS; ++L) {
                //  qloc[L] = iQa(i,j,k,L);
                //}
                //wdot(i,j,k,n) = W_spec_d(qloc,n);
                wdot(i,j,k,n) = W_spec1_d(&(iQa(i,j,k,0)),numPts,n);
              });
#endif
          //void* p = nullptr;
          //cudaStreamAddCallback(Gpu::gpuStream(), foo, p, 0);
          //int idx = mfi.index();
          //cudaLaunchHostFunc(Gpu::gpuStream(), foo1, (void*)(&idx));
        }
      }

      outfile = amrex::Concatenate(pltfile,1);
      PlotFileFromMF(wdots,outfile);

      CKFINALIZE();
      extern_close();

    }

    Finalize();

    return 0;
}

