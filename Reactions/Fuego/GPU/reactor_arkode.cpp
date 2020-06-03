#include <GPU/reactor_arkode.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <EOS.H>
#include <AMReX_Gpu.H>

using namespace amrex;

AMREX_GPU_DEVICE_MANAGED  int eint_rho = 1; // in/out = rhoE/rhoY
AMREX_GPU_DEVICE_MANAGED  int enth_rho = 2; // in/out = rhoH/rhoY 
AMREX_GPU_DEVICE_MANAGED int use_erkode=0;

/******************************************************************************************/
/* Initialization routine, called once at the begining of the problem */
int reactor_info(const int* reactor_type, const int* Ncells)
{
    amrex::ParmParse pp("ode");
    pp.query("use_erkode", use_erkode);
    if(use_erkode==1)
    {
        amrex::Print()<<"Using ERK ODE\n";
    }
    else
    {
        amrex::Print()<<"Using ARK ODE\n";
    }
    return(0);
}
/******************************************************************************************/
/* Main call routine */
int react(realtype *rY_in, realtype *rY_src_in, 
        realtype *rX_in, realtype *rX_src_in,
        realtype *dt_react, realtype *time,
        const int* reactor_type,const int* Ncells, cudaStream_t stream,
        double reltol,double abstol)
{

    int NCELLS, NEQ, neq_tot,flag;
    realtype time_init, time_out;
    void *arkode_mem    = NULL;
    N_Vector y         = NULL;

    NEQ = NUM_SPECIES;
    NCELLS         = *Ncells;
    neq_tot        = (NEQ + 1) * NCELLS;

    /* User data */
    UserData user_data;
    BL_PROFILE_VAR("AllocsInARKODE", AllocsARKODE);
    cudaMallocManaged(&user_data, sizeof(struct ARKODEUserData));
    BL_PROFILE_VAR_STOP(AllocsARKODE);
    user_data->ncells_d[0]             = NCELLS;
    user_data->neqs_per_cell[0]        = NEQ;
    user_data->ireactor_type           = *reactor_type; 
    user_data->iverbose                = 1;
    user_data->stream                  = stream;
    user_data->nbBlocks                = NCELLS/32;
    user_data->nbThreads               = 32;

    y = N_VNewManaged_Cuda(neq_tot);
    N_VSetCudaStream_Cuda(y, &stream);


    BL_PROFILE_VAR_START(AllocsARKODE);
    /* Define vectors to be used later in creact */
    cudaMalloc(&(user_data->rhoe_init), NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rhoesrc_ext), NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rYsrc), (NCELLS*NEQ)*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsARKODE);

    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
    BL_PROFILE_VAR("AsyncCpy", AsyncCpy);
    // rhoY,T
    cudaMemcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS), cudaMemcpyHostToDevice);
    // rhoY_src_ext
    cudaMemcpy(user_data->rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double), cudaMemcpyHostToDevice);
    // rhoE/rhoH
    cudaMemcpy(user_data->rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
    cudaMemcpy(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
    BL_PROFILE_VAR_STOP(AsyncCpy);

    /* Initial time and time to reach after integration */
    time_init = *time;
    time_out  = *time + (*dt_react);
    
    if(use_erkode == 0)
    {
       arkode_mem = ARKStepCreate(cF_RHS, NULL, *time, y);
        flag = ARKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
        flag = ARKStepSStolerances(arkode_mem, reltol, abstol); 
        flag = ARKStepResStolerance(arkode_mem, abstol);
        flag = ARKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);      /* call integrator */
    }
    else
    {
       arkode_mem = ERKStepCreate(cF_RHS, *time, y);
       flag = ERKStepSetUserData(arkode_mem, static_cast<void*>(user_data));
       flag = ERKStepSStolerances(arkode_mem, reltol, abstol); 
       flag = ERKStepEvolve(arkode_mem, time_out, y, &time_init, ARK_NORMAL);      /* call integrator */
    }

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    *dt_react = time_init - *time;
    *time  = time_init + (*dt_react);
#endif
    /* Pack data to return in main routine external */
    BL_PROFILE_VAR_START(AsyncCpy);
    cudaMemcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost);

    for  (int i = 0; i < NCELLS; i++) {
        rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
    }
    BL_PROFILE_VAR_STOP(AsyncCpy);


    /* Get estimate of how hard the integration process was */
    long int nfe,nfi;
    if(use_erkode==0)
    {
        flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
    }
    else
    {
        flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
    }

    //cleanup
    N_VDestroy(y);          /* Free the y vector */
    if(use_erkode==0)
    {
        ARKStepFree(&arkode_mem);
    }
    else
    {
        ERKStepFree(&arkode_mem);
    }

    cudaFree(user_data->rhoe_init);
    cudaFree(user_data->rhoesrc_ext);
    cudaFree(user_data->rYsrc);
    cudaFree(user_data);

    return nfe;
}
/******************************************************************************************/
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
        void *user_data)
{

    BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);

    cudaError_t cuda_status = cudaSuccess;

    /* Get Device pointers for Kernel call */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
    realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);

    // allocate working space 
    UserData udata = static_cast<ARKODEUserData*>(user_data);
    udata->dt_save = t;

    const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
    //amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
    amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
            [=] AMREX_GPU_DEVICE () noexcept {
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                    icell < udata->ncells_d[0]; icell += stride) {
            fKernelSpec(icell, user_data, yvec_d, ydot_d, udata->rhoe_init, 
                    udata->rhoesrc_ext, udata->rYsrc);    
            }
            }); 

    cuda_status = cudaStreamSynchronize(udata->stream);  
    assert(cuda_status == cudaSuccess);

    BL_PROFILE_VAR_STOP(fKernelSpec);

    return(0);
}
/******************************************************************************************/
/**********************************/
/*
 * CUDA kernels
 */
AMREX_GPU_DEVICE
inline void 
fKernelSpec(int icell, void *user_data, 
        realtype *yvec_d, realtype *ydot_d,  
        double *rhoe_init, double *rhoesrc_ext, double *rYs)
{
    UserData udata = static_cast<ARKODEUserData*>(user_data);

    amrex::GpuArray<amrex::Real,NUM_SPECIES> mw;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> massfrac;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> ei_pt;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> cdots_pt;
    amrex::Real Cv_pt, rho_pt, temp_pt, nrg_pt;

    int offset = icell * (NUM_SPECIES + 1); 

    /* MW CGS */
    get_mw(mw.arr);

    /* rho */ 
    rho_pt = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
        rho_pt = rho_pt + yvec_d[offset + n];
    }

    /* Yks, C CGS*/
    for (int i = 0; i < NUM_SPECIES; i++){
        massfrac[i] = yvec_d[offset + i] / rho_pt;
    }

    /* NRG CGS */
    nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell]*(udata->dt_save)) /rho_pt;

    /* temp */
    temp_pt = yvec_d[offset + NUM_SPECIES];

    /* Additional var needed */
    if (udata->ireactor_type == 1){
        /* UV REACTOR */
        EOS::EY2T(nrg_pt, massfrac.arr, temp_pt);
        EOS::T2Ei(temp_pt, ei_pt.arr);
        EOS::TY2Cv(temp_pt, massfrac.arr, Cv_pt);
    }else {
        /* HP REACTOR */
        EOS::HY2T(nrg_pt, massfrac.arr, temp_pt);
        EOS::TY2Cp(temp_pt, massfrac.arr, Cv_pt);
        EOS::T2Hi(temp_pt, ei_pt.arr);
    }

    EOS::RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

    /* Fill ydot vect */
    ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
    for (int i = 0; i < NUM_SPECIES; i++){
        ydot_d[offset + i]           = cdots_pt[i] + rYs[icell * NUM_SPECIES + i];
        ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * ei_pt[i];
    }
    ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho_pt * Cv_pt);
    
    /* Fill ydot vect */
    /*ydot_d[offset + NUM_SPECIES] = 1e6;
    for (int i = 0; i < NUM_SPECIES; i++){
        ydot_d[offset + i]           = 0.0;
    }*/
}
