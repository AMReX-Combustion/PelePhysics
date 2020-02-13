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
  #ifdef USE_ARKODE_PP
    #include <actual_CARKODE.h>
  #else
    #include <actual_Creactor.h>
  #endif
#elif defined(USE_RK64_PP)
// Expl RK solver
  #include <actual_CRK64.h>
#else
// DVODE
  #include <actual_reactor.H> 
#endif

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    Real timer_tot = amrex::second();
    Real timer_init = 0.;
    Real timer_advance = 0.;
    Real timer_print = 0.;
    Real timer_print_tmp = 0.;


    {

    timer_init = amrex::second();

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
    int ndt        = 1; 
    Real dt        = 1.e-5;
#ifdef USE_ARKODE_PP 
    /* ARKODE parameters for now but should be for all solvers */
    Real rtol=1e-9;
    Real atol=1e-9;
#endif

    {
      /* ParmParse from the inputs file */
      ParmParse pp;
      
      // probin file
      pp.query("probin_file",probin_file);

      // domain size
      pp.query("max_grid_size",max_grid_size);

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

#ifdef USE_ARKODE_PP 
      /* Additional ARKODE queries */
      ppode.query("rtol",rtol);
      ppode.query("atol",atol);
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
        amrex::Print() << fuel_name << ", Oxy: O2";
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
    extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx,&ode_iE);

    /* Initialize D/CVODE reactor */
#ifdef _OPENMP
#pragma omp parallel
#endif
#ifdef USE_SUNDIALS_PP
#ifdef USE_ARKODE_PP
    reactor_init(&ode_iE, &ode_ncells,rtol,atol);
#else
    reactor_init(&ode_iE, &ode_ncells);
#endif
#else
    reactor_init(&ode_iE, &ode_ncells);
#endif

    /* make domain and BoxArray */
    std::vector<int> npts(3,1);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 1;
    }
    npts[1] = 1;

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

    IntVect tilesize(D_DECL(10240,8,32));

    /* INITIALIZE DATA */
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi ){
        const Box& box = mfi.tilebox();
        initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
               		BL_TO_FORTRAN_N_3D(mf[mfi],0),
               		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
		        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
		        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
    }

    timer_init = amrex::second() - timer_init; 

    timer_print = amrex::second();

    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);
    std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);

    timer_print = amrex::second() - timer_print;
     
    /* EVALUATE */
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    timer_advance = amrex::second();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi )
    {
	/* Prints to follow the computation */
        /* ADVANCE */
        Real time = 0.0;
        Real dt_incr   = dt/ndt;
	Real fc_tmp;

        const Box& box = mfi.tilebox();
	int ncells = box.numPts();
	amrex::Print() << " Integrating " << ncells << " cells with a "<<ode_ncells<< " ode cell buffer \n";

	const auto len     = amrex::length(box);
	const auto lo      = amrex::lbound(box);

	const auto rhoY    = mf.array(mfi);
	const auto rhoE    = mfE.array(mfi);
	const auto frcExt  = rY_source_ext.array(mfi); 
	const auto frcEExt = rY_source_energy_ext.array(mfi);
	const auto fc      = fctCount.array(mfi); 

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
    }

    timer_advance = amrex::second() - timer_advance;


    timer_print_tmp = amrex::second();

    outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);

    timer_print = amrex::second() - timer_print_tmp + timer_print;
    
    extern_close();

    }

    timer_tot = amrex::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance, timer_print},
                                     ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
                   << "Run Time init         = " << timer_init    << "\n"
                   << "Run Time advance      = " << timer_advance << "\n"
                   << "Run Time print plt    = " << timer_print << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
