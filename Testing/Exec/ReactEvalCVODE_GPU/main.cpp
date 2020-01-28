#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include "mechanism.h"
#include <GPU_misc.H>

using namespace amrex;

#include <Transport_F.H>
#include <main_F.H>
#include <PlotFileFromMF.H>
#include <actual_Creactor_GPU.h>

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //TOTAL TIME
    Real timer_init = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_init,IOProc);
    Real timer_tot = 0.;
    //INITIALIZATIOn
    Real timer_initializ_stop = 0.;
    // ADVANCE
    Real timer_adv_init = 0.;
    Real timer_adv_stop = 0.;
    //Print
    Real timer_print_init = 0.;
    Real timer_print_stop = 0.;

    {

    int max_grid_size = 16;
    std::string probin_file="probin";
    std::string fuel_name="none";
    std::string pltfile("plt");
    /* CVODE inputs */
    int cvode_ncells = 1;
    int cvode_iE = 1;
    int fuel_idx = -1;
    int oxy_idx = -1;
    int bath_idx = -1;
    int third_dim = 1024;
    int ndt = 1; 
    Real dt = 1.e-5; 

    {
      /* ParmParse from the inputs file */
      ParmParse pp;
      
      // probin file
      pp.query("probin_file",probin_file);

      // domain size
      pp.query("max_grid_size",max_grid_size);

      // final time
      pp.query("dt",dt);

      // time stepping
      pp.query("ndt",ndt); 

      // third dim
      pp.query("third_dim",third_dim);

      pp.query("reactor_type",cvode_iE);
      // Select CVODE type of energy employed.
      //1 for UV, 2 for HP
      //   1 = Internal energy
      //   anything else = enthalpy (PeleLM restart)
         
      // nb of cells to integrate
      pp.query("cvode_ncells",cvode_ncells);

      // Get name of fuel 
      pp.get("fuel_name", fuel_name);

    }

    //if (fuel_name != FUEL_NAME) {
    //    amrex::Print() << fuel_name << "!=" <<FUEL_NAME << std::endl;
    //    amrex::Abort("fuel_name is inconsistent with chosen mechanism");
    //} else {
        amrex::Print() << "Fuel: ";
            amrex::Print() << fuel_name << ", Oxy: O2";
        amrex::Print() << std::endl;
    //}

    amrex::Print() << "Integration method: ";
        amrex::Print() << "BDF (stiff)";
    amrex::Print() << std::endl;

    amrex::Print() << "Integration iteration method: ";
        amrex::Print() << "Newton";
    amrex::Print() << std::endl;

    amrex::Print() << "Type of reactor: ";
        amrex::Print() << cvode_iE;
    amrex::Print() << std::endl;


    /* take care of probin init to initialize problem */
    int probin_file_length = probin_file.length();
    std::vector<int> probin_file_name(probin_file_length);
    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];
    if (fuel_name == "H2") {
        fuel_idx  = H2_ID;
        amrex::Print() << "FUEL IS H2 \n";
    } else if (fuel_name == "CH4") {
        fuel_idx  = CH4_ID;
        amrex::Print() << "FUEL IS CH4 \n";
#ifdef NC12H26_ID
    } else if (fuel_name == "NC12H26") {
        fuel_idx  = NC12H26_ID;
        amrex::Print() << "FUEL IS NC12H26 \n";
#endif
    }
    oxy_idx   = O2_ID;
    bath_idx  = N2_ID;
    extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx,&cvode_iE);

    BL_PROFILE_VAR("reactor_info()", reactInfo);

    /* Initialize D/CVODE reactor */
    reactor_info(&cvode_iE, &cvode_ncells);

    BL_PROFILE_VAR_STOP(reactInfo);

    /* make domain and BoxArray */
    std::vector<int> npts(3,1);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = third_dim;
    }

    amrex::Print() << "Integrating "<<npts[0]<< "x"<<npts[1]<< "x"<<npts[2]<< "  box for: ";
        amrex::Print() << dt << " seconds";
    amrex::Print() << std::endl;


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

    //IntVect tilesize(D_DECL(10240,8,32));

    BL_PROFILE_VAR("initialize_data()", InitData);

    int count_mf = 0;
    /* INITIALIZE DATA */
    for (MFIter mfi(mf,false); mfi.isValid(); ++mfi ){
        const Box& box = mfi.tilebox();
        initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
               		BL_TO_FORTRAN_N_3D(mf[mfi],0),
               		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
		        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
		        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
	count_mf = count_mf + 1;
    }

    BL_PROFILE_VAR_STOP(InitData);

    amrex::Print() << "That many boxes: " << count_mf<< "\n";

    timer_initializ_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_initializ_stop,IOProc);

    //timer_print_init = ParallelDescriptor::second();
    //ParallelDescriptor::ReduceRealMax(timer_print_init,IOProc);

    BL_PROFILE_VAR("PlotFileFromMF()", PlotFile);

    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);
    std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);

    BL_PROFILE_VAR_STOP(PlotFile);

    //timer_print_stop = ParallelDescriptor::second();
    //ParallelDescriptor::ReduceRealMax(timer_print_stop,IOProc);
     
    /* EVALUATE */
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    BL_PROFILE_VAR("Malloc()", Allocs);
    BL_PROFILE_VAR_STOP(Allocs);

    BL_PROFILE_VAR("React()", ReactInLoop);
    BL_PROFILE_VAR_STOP(ReactInLoop);

    BL_PROFILE_VAR("(un)flatten()", FlatStuff);
    BL_PROFILE_VAR_STOP(FlatStuff);

    BL_PROFILE_VAR("advance()", Advance);

    timer_adv_init = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_adv_init,IOProc);

    for ( MFIter mfi(mf,false); mfi.isValid(); ++mfi )
    {
        std::cout<<"stream "<<amrex::Gpu::gpuStream()<<std::endl;
	/* Prints to follow the computation */
        /* ADVANCE */
        Real time = 0.0;
        Real dt_incr   = dt/ndt;
        amrex::Real fc_pt;

        cudaError_t cuda_status = cudaSuccess;

        const Box& box = mfi.tilebox();
	int ncells     = box.numPts();

	const auto len     = amrex::length(box);
	const auto lo      = amrex::lbound(box);

        const auto ec = Gpu::ExecutionConfig(ncells);

        /* VERSION 2 */
        const auto rhoY    = mf.array(mfi);
        const auto rhoE    = mfE.array(mfi);
        const auto frcExt  = rY_source_ext.array(mfi);
        const auto frcEExt = rY_source_energy_ext.array(mfi);
        const auto fc      = fctCount.array(mfi);


        /* Pack the data NEED THOSE TO BE DEF ALWAYS */
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
        /* SECOND VERSION */
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
        //for (int ii = 0; ii < ncells; ++ii) {
	//    printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp + ii*(NUM_SPECIES + 1)]);
        //}
        //printf("\n \n ");
	for (int ii = 0; ii < ndt; ++ii) {
	    fc_pt = react(tmp_vect, tmp_src_vect,
	                    tmp_vect_energy, tmp_src_vect_energy,
	                    &dt_incr, &time,
                            &cvode_iE, &ncells, amrex::Gpu::gpuStream());
	    //printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp + (NUM_SPECIES + 1)]);
	    dt_incr =  dt/ndt;
        }
        //for (int ii = 0; ii < ncells; ++ii) {
        //    printf("%14.6e %14.6e \n", time, tmp_vect[Ncomp + (NUM_SPECIES + 1)*ii]);
        //}
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

    }

    BL_PROFILE_VAR_STOP(Advance);

    timer_adv_stop = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_adv_stop,IOProc);


    timer_print_init = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(timer_print_init,IOProc);

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


    amrex::Print() << "Run Time total        = " << timer_tot            - timer_init    << "\n"
                   << "Run Time init         = " << timer_initializ_stop - timer_init    << "\n"
                   << "Run Time advance      = " << timer_adv_stop       - timer_adv_init << "\n"
                   << "Run Time print plt    = " << timer_print_stop     - timer_print_init << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
