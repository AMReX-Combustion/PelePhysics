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

#if defined(USE_SUNDIALS_PP)
// ARKODE or CVODE
  #if defined(USE_CUDA_SUNDIALS_PP)
    #include <GPU_misc.H> 
  #endif
  #ifdef USE_ARKODE_PP
    #include <actual_CARKODE.h>
  #else
    #include <actual_Creactor.h>
  #endif
#else
  #if defined(USE_RK64_PP)
// Expl RK solver
    #include <actual_CRK64.h>
  #else
// DVODE
    #include <actual_reactor.H> 
  #endif
#endif

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

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
    std::string probin_file="probin";
    std::string fuel_name="none";
    std::string pltfile("plt");
    /* Mixture info */
    int fuel_idx   = -1;
    int oxy_idx    = -1;
    int bath_idx   = -1;
    /* ODE inputs */
    int ode_ncells = 1;
    int ode_iE     = -1;
    int third_dim  = 1024;
    int ndt        = 1; 
    Real dt        = 1.e-5;
#ifdef USE_SUNDIALS_PP
    /* ARKODE parameters for now but should be for all solvers */
    Real rtol=1e-9;
    Real atol=1e-9;
    int set_typ_vals = 0;
#endif

    {
      /* ParmParse from the inputs file */
      ParmParse pp;
      
      // probin file
      pp.query("probin_file",probin_file);

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

#if defined(USE_CUDA_SUNDIALS_PP)
      // nb of cells to integrate
      ppode.query("ode_ncells",ode_ncells);
#endif

#ifdef USE_SUNDIALS_PP
      /* Additional ARKODE queries */
      ppode.query("rtol",rtol);
      ppode.query("atol",atol);
      ppode.query("set_typ_vals",set_typ_vals);
#endif

    }

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
#else
    amrex::Print()<<"Using DVODE (implicit solver)";
#endif
#endif
    amrex::Print() << std::endl;

    amrex::Print() << "Type of reactor: ";
        amrex::Print() << ode_iE;
    amrex::Print() << std::endl;

    amrex::Print() << "Fuel: ";

    /* take care of probin init to initialize problem */
    int probin_file_length = probin_file.length();
    std::vector<int> probin_file_name(probin_file_length);
    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];
    if (fuel_name == "H2") {
        fuel_idx  = H2_ID;
        amrex::Print() << fuel_name << ", Oxy: O2";
        amrex::Print() << std::endl;
#ifdef CH4_ID
    } else if (fuel_name == "CH4") {
        fuel_idx  = CH4_ID;
        amrex::Print() << fuel_name << ", Oxy: O2";
        amrex::Print() << std::endl;
#endif
#ifdef NC12H26_ID
    } else if (fuel_name == "NC12H26") {
        fuel_idx  = NC12H26_ID;
        amrex::Print() << fuel_name << ", Oxy: O2";
        amrex::Print() << std::endl;
#endif
    }
    oxy_idx   = O2_ID;
    bath_idx  = N2_ID;
    extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx,&ode_iE);

    BL_PROFILE_VAR("reactor_info()", reactInfo);

    /* Initialize D/CVODE reactor */
#ifdef USE_SUNDIALS_PP
#ifdef USE_CUDA_SUNDIALS_PP
    reactor_info(&ode_iE, &ode_ncells);
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    // init species-specific abs tolerances
    if (set_typ_vals) {
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
        SetTypValsCVODE(typ_vals);
    }
    reactor_init(&ode_iE, &ode_ncells,rtol,atol);
#endif
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    reactor_init(&ode_iE, &ode_ncells);
#endif

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

    BL_PROFILE_VAR("initialize_data()", InitData);

    /* INITIALIZE DATA */
#ifndef USE_CUDA_SUNDIALS_PP
    IntVect tilesize(D_DECL(10240,8,32));
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi ){
#else
    int count_mf = 0;
    for (MFIter mfi(mf,false); mfi.isValid(); ++mfi ){
#endif
        const Box& box = mfi.tilebox();
        initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
               		BL_TO_FORTRAN_N_3D(mf[mfi],0),
               		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
		        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
		        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
#ifdef USE_CUDA_SUNDIALS_PP
	count_mf = count_mf + 1;
#endif
    }

    BL_PROFILE_VAR_STOP(InitData);

#ifdef USE_CUDA_SUNDIALS_PP
    amrex::Print() << "That many boxes: " << count_mf<< "\n";
#endif

    timer_initialize_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_initialize_stop,IOProc);

    BL_PROFILE_VAR("PlotFileFromMF()", PlotFile);

    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);
    std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);

    BL_PROFILE_VAR_STOP(PlotFile);

    /* EVALUATE */
    amrex::Print() << " \n STARTING THE ADVANCE \n";

#ifdef USE_CUDA_SUNDIALS_PP
    BL_PROFILE_VAR("Malloc()", Allocs);
    BL_PROFILE_VAR_STOP(Allocs);

    BL_PROFILE_VAR("React()", ReactInLoop);
    BL_PROFILE_VAR_STOP(ReactInLoop);

    BL_PROFILE_VAR("(un)flatten()", FlatStuff);
    BL_PROFILE_VAR_STOP(FlatStuff);
#endif

    BL_PROFILE_VAR("advance()", Advance);

    timer_adv = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_adv,IOProc);

#ifndef USE_CUDA_SUNDIALS_PP
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi )
#else
    for ( MFIter mfi(mf,false); mfi.isValid(); ++mfi )
#endif
    {
	/* Prints to follow the computation */
        /* ADVANCE */
        Real time = 0.0;
        Real dt_incr   = dt/ndt;
	Real fc_tmp;

        const Box& box = mfi.tilebox();
	int ncells     = box.numPts();

	const auto len     = amrex::length(box);
	const auto lo      = amrex::lbound(box);

#ifdef USE_CUDA_SUNDIALS_PP
        std::cout<<"stream "<<amrex::Gpu::gpuStream()<<std::endl;

        cudaError_t cuda_status = cudaSuccess;
        const auto ec = Gpu::ExecutionConfig(ncells);
#endif
        /* VERSION 2 */
        const auto rhoY    = mf.array(mfi);
        const auto rhoE    = mfE.array(mfi);
        const auto frcExt  = rY_source_ext.array(mfi);
        const auto frcEExt = rY_source_energy_ext.array(mfi);
        const auto fc      = fctCount.array(mfi);


#ifdef USE_CUDA_SUNDIALS_PP
	amrex::Print() << " Integrating " << ncells <<" cells \n";
        /* Pack the data */
	// rhoY,T
        amrex::Real *tmp_vect; 
	// rhoY_src_ext
        amrex::Real *tmp_src_vect;
	// rhoE/rhoH
        amrex::Real *tmp_vect_energy;
	amrex::Real *tmp_src_vect_energy;
        
	BL_PROFILE_VAR_START(Allocs);
        cudaMallocManaged(&tmp_vect, (Ncomp+1)*ncells*sizeof(amrex::Real));
        cudaMallocManaged(&tmp_src_vect, Ncomp*ncells*sizeof(amrex::Real));
        cudaMallocManaged(&tmp_vect_energy, ncells*sizeof(amrex::Real));
        cudaMallocManaged(&tmp_src_vect_energy, ncells*sizeof(amrex::Real));
	BL_PROFILE_VAR_STOP(Allocs);

        BL_PROFILE_VAR_START(FlatStuff);
        /* Packing of data */
	amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	[=] AMREX_GPU_DEVICE () noexcept {
	    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	        icell < ncells; icell += stride) {
	        int k =  icell /   (len.x*len.y);
		int j = (icell - k*(len.x*len.y)) /   len.x;
		int i = (icell - k*(len.x*len.y)) - j*len.x;
		i += lo.x;
		j += lo.y;
		k += lo.z;
                gpu_flatten(icell, i, j, k, rhoY, frcExt, rhoE, frcEExt, 
                                            tmp_vect, tmp_src_vect, tmp_vect_energy, tmp_src_vect_energy);
            }
        });
	BL_PROFILE_VAR_STOP(FlatStuff);
        
        cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());  

        /* Solve */
        BL_PROFILE_VAR_START(ReactInLoop);
	time = 0.0;
	for (int ii = 0; ii < ndt; ++ii) {
	    fc_tmp = react(tmp_vect, tmp_src_vect,
	                    tmp_vect_energy, tmp_src_vect_energy,
	                    &dt_incr, &time,
                            &ode_iE, &ncells, amrex::Gpu::gpuStream());
	    //printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp + (NUM_SPECIES + 1)]);
	    dt_incr =  dt/ndt;
        }
        BL_PROFILE_VAR_STOP(ReactInLoop);

        BL_PROFILE_VAR_START(FlatStuff);
        /* Unpacking of data */
	amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	[=] AMREX_GPU_DEVICE () noexcept {
	    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	        icell < ncells; icell += stride) {
	        int k =  icell /   (len.x*len.y);
		int j = (icell - k*(len.x*len.y)) /   len.x;
		int i = (icell - k*(len.x*len.y)) - j*len.x;
		i += lo.x;
		j += lo.y;
		k += lo.z;
                gpu_unflatten(icell, i, j, k, rhoY, rhoE, 
                                            tmp_vect, tmp_vect_energy);
            }
        });
        BL_PROFILE_VAR_STOP(FlatStuff);

        cudaFree(tmp_vect);
        cudaFree(tmp_src_vect);
        cudaFree(tmp_vect_energy);
        cudaFree(tmp_src_vect_energy);
       
        cuda_status = cudaStreamSynchronize(amrex::Gpu::gpuStream());  
#else

	amrex::Print() << " Integrating " << ncells << " cells with a "<<ode_ncells<< " ode cell buffer \n";
        /* Pack the data */
	// rhoY,T
	double tmp_vect[ode_ncells*(Ncomp+1)];
	// rhoY_src_ext
	double tmp_src_vect[ode_ncells*(Ncomp)];
	// rhoE/rhoH
	double tmp_vect_energy[ode_ncells];
	double tmp_src_vect_energy[ode_ncells];

	int indx_i[ode_ncells];
	int indx_j[ode_ncells];
	int indx_k[ode_ncells];
	
	int nc = 0;
	int num_cell_ode_int = 0;
	for         (int k = 0; k < len.z; ++k) {
	    for         (int j = 0; j < len.y; ++j) {
	        for         (int i = 0; i < len.x; ++i) {
		    /* Fill the vectors */
	            for (int sp=0;sp<Ncomp; sp++){
	                tmp_vect[nc*(Ncomp+1) + sp]   = rhoY(i+lo.x,j+lo.y,k+lo.z,sp);
		        tmp_src_vect[nc*Ncomp + sp]   = frcExt(i+lo.x,j+lo.y,k+lo.z,sp);
		    }
		    tmp_vect[nc*(Ncomp+1) + Ncomp]    = rhoY(i+lo.x,j+lo.y,k+lo.z,Ncomp);
		    tmp_vect_energy[nc]               = rhoE(i+lo.x,j+lo.y,k+lo.z,0);
		    tmp_src_vect_energy[nc]           = frcEExt(i+lo.x,j+lo.y,k+lo.z,0);
		    //
		    indx_i[nc] = i+lo.x;
		    indx_j[nc] = j+lo.y;
		    indx_k[nc] = k+lo.z;
		    //
		    nc = nc+1;
		    //
		    num_cell_ode_int = num_cell_ode_int + 1;
		    if (nc == ode_ncells) {
			time = 0.0;
			dt_incr =  dt/ndt;
			for (int ii = 0; ii < ndt; ++ii) {
#if defined(USE_SUNDIALS_PP) || defined(USE_RK64_PP)
	                    fc_tmp = react(tmp_vect, tmp_src_vect,
		                tmp_vect_energy, tmp_src_vect_energy,
		                &dt_incr, &time);
#else
                            double pressure = 1013250.0;
	                    fc_tmp = react(tmp_vect, tmp_src_vect,
		                tmp_vect_energy, tmp_src_vect_energy,
				&pressure,
		                &dt_incr, &time);
#endif
		            dt_incr =  dt/ndt;
			    printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp]);
			}
		        nc = 0;
		        for (int l = 0; l < ode_ncells ; ++l){
		            for (int sp=0;sp<Ncomp; sp++){
		                rhoY(indx_i[l],indx_j[l],indx_k[l],sp) = tmp_vect[l*(Ncomp+1) + sp];
		            }
		            rhoY(indx_i[l],indx_j[l],indx_k[l],Ncomp)  = tmp_vect[l*(Ncomp+1) + Ncomp];
		            rhoE(indx_i[l],indx_j[l],indx_k[l],0)      = tmp_vect_energy[l];
			    fc(indx_i[l],indx_j[l],indx_k[l],0)        = fc_tmp;
		        }
		    }
		}
	    }
	}
	if (nc != 0) {
		printf(" WARNING !! Not enough cells (%d) to fill %d \n", nc, ode_ncells);
	} else {
		printf(" Integrated %d cells \n",num_cell_ode_int);
	}
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
    
    extern_close();

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
