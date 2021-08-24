#include "reactor.H"

using namespace amrex;

#define SUN_CUSP_CONTENT(S) ((SUNLinearSolverContent_Dense_custom)(S->content))
#define SUN_CUSP_SUBSYS_SIZE(S) (SUN_CUSP_CONTENT(S)->subsys_size)
#define SUN_CUSP_NUM_SUBSYS(S) (SUN_CUSP_CONTENT(S)->nsubsys)
#define SUN_CUSP_SUBSYS_NNZ(S) (SUN_CUSP_CONTENT(S)->subsys_nnz)

#define SUN_CUSP_LASTFLAG(S) (SUN_CUSP_CONTENT(S)->last_flag)
#define SUN_CUSP_STREAM(S) (SUN_CUSP_CONTENT(S)->stream)
#define SUN_CUSP_NBLOCK(S) (SUN_CUSP_CONTENT(S)->nbBlocks)
#define SUN_CUSP_NTHREAD(S) (SUN_CUSP_CONTENT(S)->nbThreads)

int reactor_init(int reactor_type, int Ncells)
{
  ParmParse pp("ode");
  int ianalytical_jacobian = 0;
  pp.query("analytical_jacobian", ianalytical_jacobian);
  int iverbose = 1;
  pp.query("verbose", iverbose);

  if (iverbose > 0) {
    Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
    Print() << "Number of cells in one solve is " << Ncells << "\n";
  }

  std::string solve_type_str = "none";
  ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  /*
  sparse_custom_solve          = 1;
  sparse_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {
    isolve_type = sparse_custom_solve;
  } else if (solve_type_str == "sparse") {
    isolve_type = sparse_solve;
  } else if (solve_type_str == "GMRES") {
    isolve_type = iterative_gmres_solve;
  } else {
    Abort("Wrong solve_type. Options are: sparse, sparse_custom, GMRES");
  }

  if (isolve_type == iterative_gmres_solve) {
    if (ianalytical_jacobian == 1) {
#if defined(AMREX_USE_CUDA)
      if (iverbose > 0) {
        Print() << "Using an Iterative GMRES Solver with sparse simplified "
                   "preconditioning \n";
      }
      int nJdata;
      int HP;
      if (reactor_type == eint_rho) {
        HP = 0;
      } else {
        HP = 1;
      }
      // Precond data
      SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata, &HP);
      if (iverbose > 0) {
        Print() << "--> SPARSE Preconditioner -- non zero entries: " << nJdata
                << ", which represents "
                << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
                << " % fill-in pattern\n";
      }
#elif defined(AMREX_USE_HIP)
      Abort("\n--> With HIP the only option is NP GMRES \n");
#else
      Abort("\n--> Not sure what do with analytic Jacobian in this case \n");
#endif
    } else {
      if (iverbose > 0) {
        Print()
          << "Using an Iterative GMRES Solver without preconditionning \n";
      }
    }

  } else if (isolve_type == sparse_custom_solve) {
#if defined(AMREX_USE_CUDA)
    if (ianalytical_jacobian == 1) {
      if (iverbose > 0) {
        Print() << "Using a custom sparse direct solver (with an analytical "
                   "Jacobian) \n";
      }
      int nJdata;
      int HP;
      if (reactor_type == eint_rho) {
        HP = 0;
      } else {
        HP = 1;
      }
      SPARSITY_INFO_SYST(&nJdata, &HP, Ncells);
      if (iverbose > 0) {
        Print() << "--> SPARSE Solver -- non zero entries: " << nJdata
                << ", which represents "
                << nJdata /
                     float(Ncells * (NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                     100.0
                << " % fill-in pattern\n";
      }
    } else {
      Abort(
        "\n--> Custom direct sparse solver requires an analytic Jacobian \n");
    }
#elif defined(AMREX_USE_HIP)
    Abort("\n--> With HIP the only option is NP GMRES \n");
#else
    Abort("\n--> No sparse solve for this case\n");
#endif

  } else if (isolve_type == sparse_solve) {
#if defined(AMREX_USE_CUDA)
    if (ianalytical_jacobian == 1) {
      Print() << "Using a Sparse Direct Solver based on cuSolver \n";
      int nJdata;
      int HP;
      if (reactor_type == eint_rho) {
        HP = 0;
      } else {
        HP = 1;
      }
      SPARSITY_INFO_SYST(&nJdata, &HP, Ncells);
      if (iverbose > 0) {
        Print() << "--> SPARSE Solver -- non zero entries: " << nJdata
                << ", which represents "
                << nJdata /
                     float(Ncells * (NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                     100.0
                << " % fill-in pattern\n";
      }
    } else {
      Abort("\n--> Sparse direct solvers requires an analytic Jacobian \n");
    }
#elif defined(AMREX_USE_HIP)
    Abort("\n--> No HIP equivalent to cuSolver implemented yet \n");
#else
    Abort("\n--> No batch solver for this case\n");
#endif

  } else {
    Abort("\n--> Bad linear solver choice for CVODE \n");
  }

  Print() << "\n--> DONE WITH INITIALIZATION (GPU)" << reactor_type << "\n";

  return (0);
}

int
react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rY_src_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time,
  const int& reactor_type
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  N_Vector y = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
  void* cvode_mem = NULL;
  int flag;
  int NCELLS, neq_tot;
  N_Vector atol;
  realtype* ratol;

  ParmParse pp("ode");
  int ianalytical_jacobian = 0;
  pp.query("analytical_jacobian", ianalytical_jacobian);
  int iverbose = 1;
  pp.query("verbose", iverbose);

  std::string solve_type_str = "none";
  ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  /*
  sparse_custom_solve          = 1;
  sparse_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {
    isolve_type = sparse_custom_solve;
  } else if (solve_type_str == "sparse") {
    isolve_type = sparse_solve;
  } else if (solve_type_str == "GMRES") {
    isolve_type = iterative_gmres_solve;
  }

  NCELLS = box.numPts();
  neq_tot = (NUM_SPECIES + 1) * NCELLS;

#ifdef AMREX_USE_OMP
  Gpu::streamSynchronize();
#endif

  CVODEUserData * user_data;

  BL_PROFILE_VAR("AllocsInCVODE", AllocsCVODE);
  user_data = (CVODEUserData*)The_Arena()->alloc(sizeof(struct CVODEUserData));
  BL_PROFILE_VAR_STOP(AllocsCVODE);

  user_data->ncells = NCELLS;
  user_data->neqs_per_cell = NUM_SPECIES;
  user_data->ireactor_type = reactor_type;
  user_data->ianalytical_jacobian = ianalytical_jacobian;
  user_data->isolve_type = isolve_type;
  user_data->iverbose = iverbose;
  user_data->nbBlocks = std::max(1, NCELLS / 32);
#ifdef AMREX_USE_GPU
  user_data->stream = stream;
#endif
  user_data->nbThreads = 32;

#ifdef AMREX_USE_CUDA
  if (user_data->ianalytical_jacobian == 1) {
    int HP;
    if (user_data->ireactor_type == 1) {
      HP = 0;
    } else {
      HP = 1;
    }

    BL_PROFILE_VAR("SparsityFuegoStuff", SparsityStuff);
    BL_PROFILE_VAR_STOP(SparsityStuff);
    if (user_data->isolve_type == iterative_gmres_solve) {
      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_INFO_SYST_SIMPLIFIED(&(user_data->NNZ), &HP);
      BL_PROFILE_VAR_STOP(SparsityStuff);

      BL_PROFILE_VAR_START(AllocsCVODE);
      user_data->csr_row_count_h =
        (int*)The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
      user_data->csr_col_index_h =
        (int*)The_Arena()->alloc(user_data->NNZ * sizeof(int));
      user_data->csr_jac_h =
        (double*)The_Arena()->alloc(user_data->NNZ * NCELLS * sizeof(double));
      user_data->csr_val_h =
        (double*)The_Arena()->alloc(user_data->NNZ * NCELLS * sizeof(double));

      user_data->csr_row_count_d =
        (int*)The_Device_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
      user_data->csr_col_index_d =
        (int*)The_Device_Arena()->alloc(user_data->NNZ * sizeof(int));
      user_data->csr_jac_d = (double*)The_Device_Arena()->alloc(
        user_data->NNZ * NCELLS * sizeof(double));
      user_data->csr_val_d = (double*)The_Device_Arena()->alloc(
        user_data->NNZ * NCELLS * sizeof(double));
      BL_PROFILE_VAR_STOP(AllocsCVODE);

      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
        user_data->csr_col_index_h, user_data->csr_row_count_h, &HP, 1);
      BL_PROFILE_VAR_STOP(SparsityStuff);

      amrex::Gpu::htod_memcpy(
        &user_data->csr_col_index_d, &user_data->csr_col_index_h,
        sizeof(user_data->NNZ * sizeof(int)));
      amrex::Gpu::htod_memcpy(
        &user_data->csr_row_count_d, &user_data->csr_row_count_h,
        sizeof((NUM_SPECIES + 2) * sizeof(int)));

      BL_PROFILE_VAR("CuSolverInit", CuSolverInit);
      size_t workspaceInBytes, internalDataInBytes;
      cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
      cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

      workspaceInBytes = 0;
      internalDataInBytes = 0;

      cusolver_status = cusolverSpCreate(&(user_data->cusolverHandle));
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusolver_status = cusolverSpSetStream(user_data->cusolverHandle, stream);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusparse_status = cusparseCreateMatDescr(&(user_data->descrA));
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      cusparse_status =
        cusparseSetMatType(user_data->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      cusparse_status =
        cusparseSetMatIndexBase(user_data->descrA, CUSPARSE_INDEX_BASE_ONE);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      cusolver_status = cusolverSpCreateCsrqrInfo(&(user_data->info));
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusolver_status = cusolverSpXcsrqrAnalysisBatched(
        user_data->cusolverHandle, NUM_SPECIES + 1, NUM_SPECIES + 1,
        user_data->NNZ, user_data->descrA, user_data->csr_row_count_h,
        user_data->csr_col_index_h, user_data->info);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      /*
      size_t free_mem = 0;
      size_t total_mem = 0;
      cudaStat1 = cudaMemGetInfo( &free_mem, &total_mem );
      assert( cudaSuccess == cudaStat1 );
      std::cout<<"(AFTER SA) Free: "<< free_mem<< " Tot:
      "<<total_mem<<std::endl;
      */

      cusolver_status = cusolverSpDcsrqrBufferInfoBatched(
        user_data->cusolverHandle, NUM_SPECIES + 1, NUM_SPECIES + 1,
        user_data->NNZ, user_data->descrA, user_data->csr_val_h,
        user_data->csr_row_count_h, user_data->csr_col_index_h, NCELLS,
        user_data->info, &internalDataInBytes, &workspaceInBytes);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cudaError_t cudaStat1 = cudaSuccess;
      cudaStat1 = cudaMalloc((void**)&(user_data->buffer_qr), workspaceInBytes);
      assert(cudaStat1 == cudaSuccess);
      BL_PROFILE_VAR_STOP(CuSolverInit);

    } else if (isolve_type == sparse_solve) {
      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_INFO_SYST(&(user_data->NNZ), &HP, 1);
      BL_PROFILE_VAR_STOP(SparsityStuff);

      BL_PROFILE_VAR_START(AllocsCVODE);
      user_data->csr_row_count_h =
        (int*)The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
      user_data->csr_col_index_h =
        (int*)The_Arena()->alloc(user_data->NNZ * sizeof(int));

      BL_PROFILE_VAR_STOP(AllocsCVODE);

      cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
      cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
      cudaError_t cuda_status = cudaSuccess;

      cusolver_status = cusolverSpCreate(&(user_data->cusolverHandle));
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusolver_status = cusolverSpSetStream(user_data->cusolverHandle, stream);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusparse_status = cusparseCreate(&(user_data->cuSPHandle));
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      cusparse_status = cusparseSetStream(user_data->cuSPHandle, stream);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      A = SUNMatrix_cuSparse_NewBlockCSR(
        NCELLS, (NUM_SPECIES + 1), (NUM_SPECIES + 1), user_data->NNZ,
        user_data->cuSPHandle);
      if (check_flag((void*)A, "SUNMatrix_cuSparse_NewBlockCSR", 0))
        return (1);

      int retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);
      if (check_flag(&retval, "SUNMatrix_cuSparse_SetFixedPattern", 1))
        return (1);

      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_PREPROC_SYST_CSR(
        user_data->csr_col_index_h, user_data->csr_row_count_h, &HP, 1, 0);
      amrex::Gpu::htod_memcpy(
        &user_data->csr_col_index_d, &user_data->csr_col_index_h,
        sizeof(user_data->csr_col_index_h));
      amrex::Gpu::htod_memcpy(
        &user_data->csr_row_count_d, &user_data->csr_row_count_h,
        sizeof(user_data->csr_row_count_h));
      SUNMatrix_cuSparse_CopyToDevice(
        A, NULL, user_data->csr_row_count_h, user_data->csr_col_index_h);
      BL_PROFILE_VAR_STOP(SparsityStuff);
    } else {
      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_INFO_SYST(&(user_data->NNZ), &HP, 1);
      BL_PROFILE_VAR_STOP(SparsityStuff);

      BL_PROFILE_VAR_START(AllocsCVODE);
      user_data->csr_row_count_h =
        (int*)The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
      user_data->csr_col_index_h =
        (int*)The_Arena()->alloc(user_data->NNZ * sizeof(int));

      user_data->csr_row_count_d =
        (int*)The_Device_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
      user_data->csr_col_index_d =
        (int*)The_Device_Arena()->alloc(user_data->NNZ * sizeof(int));
      BL_PROFILE_VAR_STOP(AllocsCVODE);

      cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

      cusparse_status = cusparseCreate(&(user_data->cuSPHandle));
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      cusparse_status = cusparseSetStream(user_data->cuSPHandle, stream);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

      A = SUNMatrix_cuSparse_NewBlockCSR(
        NCELLS, (NUM_SPECIES + 1), (NUM_SPECIES + 1), user_data->NNZ,
        user_data->cuSPHandle);
      if (check_flag((void*)A, "SUNMatrix_cuSparse_NewBlockCSR", 0))
        return (1);

      int retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);
      if (check_flag(&retval, "SUNMatrix_cuSparse_SetFixedPattern", 1))
        return (1);

      BL_PROFILE_VAR_START(SparsityStuff);
      SPARSITY_PREPROC_SYST_CSR(
        user_data->csr_col_index_h, user_data->csr_row_count_h, &HP, 1, 0);
      amrex::Gpu::htod_memcpy(
        &user_data->csr_col_index_d, &user_data->csr_col_index_h,
        sizeof(user_data->csr_col_index_h));
      amrex::Gpu::htod_memcpy(
        &user_data->csr_row_count_d, &user_data->csr_row_count_h,
        sizeof(user_data->csr_row_count_h));
      SUNMatrix_cuSparse_CopyToDevice(
        A, NULL, user_data->csr_row_count_h, user_data->csr_col_index_h);
      BL_PROFILE_VAR_STOP(SparsityStuff);
    }
  }
#endif

#if defined(AMREX_USE_CUDA)
  y = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
#elif defined(AMREX_USE_HIP)
  y = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/true,
    *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
#else
  y = N_VNew_Serial(neq_tot);
  if (check_flag((void*)y, "N_VNew_Serial", 0))
    return (1);
#endif

#if defined(AMREX_USE_CUDA)
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#endif

  cvode_mem = CVodeCreate(CV_BDF);
  if (check_flag((void*)cvode_mem, "CVodeCreate", 0))
    return (1);

  flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

  BL_PROFILE_VAR_START(AllocsCVODE);
  user_data->rhoe_init_d =
    (double*)The_Device_Arena()->alloc(NCELLS * sizeof(double));
  user_data->rhoesrc_ext_d =
    (double*)The_Device_Arena()->alloc(NCELLS * sizeof(double));
  user_data->rYsrc_d =
    (double*)The_Device_Arena()->alloc(NCELLS * NUM_SPECIES * sizeof(double));
  BL_PROFILE_VAR_STOP(AllocsCVODE);

#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
  const auto len = amrex::length(box);
  const auto lo = amrex::lbound(box);
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);

    box_flatten(
      icell, NCELLS, i, j, k, user_data->ireactor_type, rY_in, rY_src_in, T_in,
      rEner_in, rEner_src_in, yvec_d, user_data->rYsrc_d,
      user_data->rhoe_init_d, user_data->rhoesrc_ext_d);
  });
  BL_PROFILE_VAR_STOP(FlatStuff);

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  amrex::Real time_init = time;
  amrex::Real time_out = time + dt_react;

  // Call CVodeInit to initialize the integrator memory and specify the
  //  user's right hand side function, the inital time, and
  //  initial dependent variable vector y.
  flag = CVodeInit(cvode_mem, cF_RHS, time_init, y);
  if (check_flag(&flag, "CVodeInit", 1))
    return (1);

  atol = N_VClone(y);
#if defined(AMREX_USE_CUDA)
  ratol = N_VGetHostArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  ratol = N_VGetHostArrayPointer_Hip(atol);
#else
  ratol = N_VGetArrayPointer(atol);
#endif
  if (typVals[0] > 0) {
    printf(
      "Setting CVODE tolerances rtol = %14.8e atolfact = %14.8e in PelePhysics "
      "\n",
      relTol, absTol);
    for (int i = 0; i < NCELLS; i++) {
      int offset = i * (NUM_SPECIES + 1);
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        ratol[offset + k] = typVals[k] * absTol;
      }
    }
  } else {
    for (int i = 0; i < neq_tot; i++) {
      ratol[i] = absTol;
    }
  }
#if defined(AMREX_USE_CUDA)
  N_VCopyToDevice_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_VCopyToDevice_Hip(atol);
#endif
  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  flag = CVodeSVtolerances(cvode_mem, relTol, atol);
  if (check_flag(&flag, "CVodeSVtolerances", 1))
    return (1);

  if (user_data->isolve_type == iterative_gmres_solve) {
    if (user_data->ianalytical_jacobian == 0) {
      LS = SUNLinSol_SPGMR(y, PREC_NONE, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    } else {
      LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    }

    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);

    if (user_data->ianalytical_jacobian == 1) {
#if defined(AMREX_USE_CUDA)
      flag = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
      if (check_flag(&flag, "CVodeSetPreconditioner", 1))
        return (1);
#else
      Abort("No options for preconditioning on non-CUDA GPUs");
#endif
    }

#if defined(AMREX_USE_CUDA)
  } else if (user_data->isolve_type == sparse_solve) {
    LS = SUNLinSol_cuSolverSp_batchQR(y, A, user_data->cusolverHandle);
    if (check_flag((void*)LS, "SUNLinSol_cuSolverSp_batchQR", 0))
      return (1);

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacFn(cvode_mem, cJac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);

  } else {
    LS = SUNLinSol_dense_custom(y, A, stream);
    if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
      return (1);

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacFn(cvode_mem, cJac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
#endif
  }

  flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    return (1);

  flag = CVodeSetMaxOrd(cvode_mem, 2);
  if (check_flag(&flag, "CVodeSetMaxOrd", 1))
    return (1);

  BL_PROFILE_VAR("AroundCVODE", AroundCVODE);
  flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
  if (check_flag(&flag, "CVode", 1))
    return (1);
  BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
  dt_react = time_init - time;
  time = time_init;
#endif

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

  BL_PROFILE_VAR_START(FlatStuff);
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);

    box_unflatten(
      icell, NCELLS, i, j, k, user_data->ireactor_type, rY_in, T_in, rEner_in,
      rEner_src_in, FC_in, yvec_d, user_data->rhoe_init_d, nfe, dt_react);
  });
  BL_PROFILE_VAR_STOP(FlatStuff);

  if (user_data->iverbose > 1) {
    PrintFinalStats(cvode_mem);
  }

  N_VDestroy(y);
  CVodeFree(&cvode_mem);

  SUNLinSolFree(LS);

  if (user_data->isolve_type != iterative_gmres_solve) {
    SUNMatDestroy(A);
  }

  The_Device_Arena()->free(user_data->rhoe_init_d);
  The_Device_Arena()->free(user_data->rhoesrc_ext_d);
  The_Device_Arena()->free(user_data->rYsrc_d);

  if (user_data->ianalytical_jacobian == 1) {
#ifdef AMREX_USE_CUDA
    The_Arena()->free(user_data->csr_row_count_h);
    The_Arena()->free(user_data->csr_col_index_h);
    if (user_data->isolve_type == iterative_gmres_solve) {
      The_Arena()->free(user_data->csr_val_h);
      The_Arena()->free(user_data->csr_jac_h);
      The_Device_Arena()->free(user_data->csr_val_d);
      The_Device_Arena()->free(user_data->csr_jac_d);

      cusolverStatus_t cusolver_status =
        cusolverSpDestroy(user_data->cusolverHandle);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusolver_status = cusolverSpDestroyCsrqrInfo(user_data->info);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cudaFree(user_data->buffer_qr);
    } else if (user_data->isolve_type == sparse_solve) {
      cusolverStatus_t cusolver_status =
        cusolverSpDestroy(user_data->cusolverHandle);
      assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

      cusparseStatus_t cusparse_status = cusparseDestroy(user_data->cuSPHandle);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    } else {
      cusparseStatus_t cusparse_status = cusparseDestroy(user_data->cuSPHandle);
      assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    }
#else
    Abort("Shoudn't be there. Analytical Jacobian only available with CUDA !");
#endif
  }

  The_Arena()->free(user_data);

  N_VDestroy(atol);

  return nfe;
}

static int
cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);

#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot_in);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Hip(ydot_in);
#else
  realtype* yvec_d = N_VGetArrayPointer(y_in);
  realtype* ydot_d = N_VGetArrayPointer(ydot_in);
#endif

  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);
  udata->dt_save = t;

  auto ncells = udata->ncells;
  auto dt_save = udata->dt_save;
  auto reactor_type = udata->ireactor_type;
  auto rhoe_init = udata->rhoe_init_d;
  auto rhoesrc_ext = udata->rhoesrc_ext_d;
  auto rYsrc = udata->rYsrc_d;
#ifdef AMREX_USE_GPU
  const auto ec = Gpu::ExecutionConfig(udata->ncells);
  // launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem,
  // udata->stream>>>(
  launch_global<<<
    udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
    [=] AMREX_GPU_DEVICE() noexcept {
      for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
               stride = blockDim.x * gridDim.x;
           icell < udata->ncells; icell += stride) {
        fKernelSpec(
          icell, ncells, dt_save, reactor_type, yvec_d, ydot_d,
          rhoe_init, rhoesrc_ext, rYsrc);
      }
    });
  Gpu::Device::streamSynchronize();
#else
  for (int icell = 0; icell < udata->ncells; icell++) {
    fKernelSpec(
      icell, ncells, dt_save, reactor_type, yvec_d, ydot_d,
      rhoe_init, rhoesrc_ext, rYsrc);
  }
#endif

  BL_PROFILE_VAR_STOP(fKernelSpec);

  return (0);
}

// React for 1d arrays
int
react(
  realtype *rY_in,
  realtype *rY_src_in, 
  realtype *rX_in,
  realtype *rX_src_in,
  realtype &dt_react,
  realtype &time,
  int reactor_type,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  // -------------------------------------------------------------
  // CVODE structs
  N_Vector y = NULL;             // Solution vector
  SUNLinearSolver LS = NULL;     // Linear solver used in Newton solve
  SUNMatrix A = NULL;            // Chem. jacobian, if provided
  void* cvode_mem = NULL;        // CVODE internal memory pointer
  N_Vector atol;                 // Component-wise abs. tolerance
  realtype* ratol;               // Rel. tolerance

  // temporary to catch CVODE exceptions
  int flag;

  // -------------------------------------------------------------
  // ParmParse integrator options
  ParmParse pp("ode");
  int ianalytical_jacobian = 0;
  pp.query("analytical_jacobian", ianalytical_jacobian);
  int iverbose = 1;
  pp.query("verbose", iverbose);

  std::string solve_type_str = "none";
  ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  /*
  sparse_custom_solve          = 1;
  sparse_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {        // Direct solve using custom GaussJordan routine
    isolve_type = sparse_custom_solve;
  } else if (solve_type_str == "sparse") {        // Direct sparse solve using cuSolver (CUDA only !)
    isolve_type = sparse_solve;
  } else if (solve_type_str == "GMRES") {         // Iter. solve using JFNK
    isolve_type = iterative_gmres_solve;
  }

  // -------------------------------------------------------------
  // Build CVODE data

  // System sizes
  int NCELLS  = Ncells;
  int NEQ     = NUM_SPECIES+1;            // rhoYs + rhoH
  int NEQSTOT = NEQ * NCELLS;

#ifdef AMREX_USE_GPU
  Gpu::streamSynchronize();
#endif

  // Allocate CVODE data
  BL_PROFILE_VAR("AllocsInCVODE", AllocsCVODE);
  CVODEUserData * user_data;
  user_data = (CVODEUserData*)The_Arena()->alloc(sizeof(struct CVODEUserData));
  user_data->rhoe_init_d =
    (double*)The_Device_Arena()->alloc(NCELLS * sizeof(double));
  user_data->rhoesrc_ext_d =
    (double*)The_Device_Arena()->alloc(NCELLS * sizeof(double));
  user_data->rYsrc_d =
    (double*)The_Device_Arena()->alloc(NCELLS * NUM_SPECIES * sizeof(double));
  BL_PROFILE_VAR_STOP(AllocsCVODE);

  // Initialize CVODE options
  user_data->ncells = NCELLS;
  user_data->neqs_per_cell = NUM_SPECIES;
  user_data->ireactor_type = reactor_type;
  user_data->ianalytical_jacobian = ianalytical_jacobian;
  user_data->isolve_type = isolve_type;
  user_data->iverbose = iverbose;
  user_data->nbBlocks = std::max(1, NCELLS / 32);
#ifdef AMREX_USE_GPU
  user_data->stream = stream;
#endif
  user_data->nbThreads = 32;
 
#ifdef AMREX_USE_GPU
  // For now, Abort() is analytical_jac required on GPU.
  if (user_data->ianalytical_jacobian == 1) {
     amrex::Abort("ode.analytical_jacobian = 1 currently unavailable on GPU");
  }
#endif

  // create solution Nvector
#if defined(AMREX_USE_CUDA)
  y = N_VNewWithMemHelp_Cuda(
    NEQSTOT, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
#elif defined(AMREX_USE_HIP)
  y = N_VNewWithMemHelp_Hip(
    NEQSTOT, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
#else
  y = N_VNew_Serial(NEQSTOT);
  if (check_flag((void*)y, "N_VNew_Serial", 0))
    return (1);
#endif

  // On GPU, define execution policy
#if defined(AMREX_USE_CUDA)
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#endif

  // Create CVODE object, specifying a Backward-Difference Formula method
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_flag((void*)cvode_mem, "CVodeCreate", 0))
    return (1);
  flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

  // -------------------------------------------------------------
  // Initialize solution vector
#if defined(AMREX_USE_CUDA)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y);
#elif defined(AMREX_USE_HIP)
  realtype* yvec_d = N_VGetDeviceArrayPointer_Hip(y);
#else
  realtype* yvec_d = N_VGetArrayPointer(y);
#endif

  BL_PROFILE_VAR("AsyncCpy", AsyncCpy);
#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(
    yvec_d, rY_in, sizeof(realtype) * (NEQSTOT));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc_d, rY_src_in, (NUM_SPECIES * NCELLS) * sizeof(realtype));
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init_d, rX_in, sizeof(realtype) * NCELLS);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext_d, rX_src_in, sizeof(realtype) * NCELLS);
#else
  std::memcpy(yvec_d, rY_in, sizeof(realtype) * (NEQSTOT));
  std::memcpy(
    user_data->rYsrc_d, rY_src_in, (NUM_SPECIES * NCELLS) * sizeof(realtype));
  std::memcpy(user_data->rhoe_init_d, rX_in, sizeof(realtype) * NCELLS);
  std::memcpy(user_data->rhoesrc_ext_d, rX_src_in, sizeof(realtype) * NCELLS);
#endif
  BL_PROFILE_VAR_STOP(AsyncCpy);

  // -------------------------------------------------------------
  // Initialize integrator
  amrex::Real time_init = time;
  amrex::Real time_out  = time + dt_react;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function, the inital time, and
   * initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, cF_RHS, time_init, y);
  if (check_flag(&flag, "CVodeInit", 1)) {
     return(1);
  }

  // -------------------------------------------------------------
  // Define CVODE tolerances
  atol = N_VClone(y);
#if defined(AMREX_USE_CUDA)
  ratol = N_VGetHostArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  ratol = N_VGetHostArrayPointer_Hip(atol);
#else
  ratol = N_VGetArrayPointer(atol);
#endif
  if (typVals[0] > 0) {
    printf(
      "Setting CVODE tolerances rtol = %14.8e atolfact = %14.8e in PelePhysics "
      "\n",
      relTol, absTol);
    for (int i = 0; i < NCELLS; i++) {
      int offset = i * (NUM_SPECIES + 1);
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        ratol[offset + k] = typVals[k] * absTol;
      }
    }
  } else {
    for (int i = 0; i < NEQSTOT; i++) {
      ratol[i] = absTol;
    }
  }
#if defined(AMREX_USE_CUDA)
  N_VCopyToDevice_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_VCopyToDevice_Hip(atol);
#endif
  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  flag = CVodeSVtolerances(cvode_mem, relTol, atol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) {
    return (1);
  }

  // -------------------------------------------------------------
  // Define and set CVODE linear solver
  if (user_data->isolve_type == iterative_gmres_solve) {
    if (user_data->ianalytical_jacobian == 0) {
      LS = SUNLinSol_SPGMR(y, PREC_NONE, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    } else {
#ifdef AMREX_USE_GPU
     amrex::Abort("ode.analytical_jacobian = 1 currently unavailable on GPU"),
#endif
      LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    }

    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);

    if (user_data->ianalytical_jacobian == 1) {
#if defined(AMREX_USE_CUDA)
      flag = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
      if (check_flag(&flag, "CVodeSetPreconditioner", 1))
        return (1);
#else
      Abort("No options for preconditioning on non-CUDA GPUs");
#endif
    }

#if defined(AMREX_USE_CUDA)
  } else if (user_data->isolve_type == sparse_solve) {
    LS = SUNLinSol_cuSolverSp_batchQR(y, A, user_data->cusolverHandle);
    if (check_flag((void*)LS, "SUNLinSol_cuSolverSp_batchQR", 0))
      return (1);

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacFn(cvode_mem, cJac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);

  } else {
    LS = SUNLinSol_dense_custom(y, A, stream);
    if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
      return (1);

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacFn(cvode_mem, cJac);
    if (check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
#endif
  }

  // -------------------------------------------------------------
  // Other CVODE options
  flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    return (1);

  flag = CVodeSetMaxOrd(cvode_mem, 2);
  if (check_flag(&flag, "CVodeSetMaxOrd", 1))
    return (1);

  // -------------------------------------------------------------
  // Actual CVODE solve
  BL_PROFILE_VAR("AroundCVODE", AroundCVODE);
  flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
  if (check_flag(&flag, "CVode", 1))
    return (1);
  BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
  dt_react = time_init - time;
  time = time_init;
#endif

  // -------------------------------------------------------------
  // Get the result back
  BL_PROFILE_VAR_START(AsyncCpy);
#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy_async(
#else
  std::memcpy(
#endif
    rY_in, yvec_d, (NEQ * NCELLS) * sizeof(realtype));
  for (int i = 0; i < NCELLS; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }
  BL_PROFILE_VAR_STOP(AsyncCpy);

  // -------------------------------------------------------------
  // Get the number of RHS evaluations
  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if (user_data->iverbose > 1) {
    PrintFinalStats(cvode_mem);
  }

  // -------------------------------------------------------------
  // Clean up memory
  N_VDestroy(y);
  CVodeFree(&cvode_mem);

  SUNLinSolFree(LS);

  if (user_data->isolve_type != iterative_gmres_solve) {
    SUNMatDestroy(A);
  }

  The_Device_Arena()->free(user_data->rhoe_init_d);
  The_Device_Arena()->free(user_data->rhoesrc_ext_d);
  The_Device_Arena()->free(user_data->rYsrc_d);

  The_Arena()->free(user_data);

  N_VDestroy(atol);

  return nfe;
}

#ifdef AMREX_USE_CUDA
static int
Precond(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data)
{

  BL_PROFILE_VAR("Precond()", Precond);

  cudaError_t cuda_status = cudaSuccess;
  size_t workspaceInBytes, internalDataInBytes;
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

  workspaceInBytes = 0;
  internalDataInBytes = 0;

  realtype* u_d = N_VGetDeviceArrayPointer_Cuda(u);
  realtype* udot_d = N_VGetDeviceArrayPointer_Cuda(fu);

  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);
  udata->gamma = gamma;

  BL_PROFILE_VAR("fKernelComputeAJ()", fKernelComputeAJ);
  if (jok) {
    const auto ec = Gpu::ExecutionConfig(udata->ncells);
    launch_global<<<
      udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
      [=] AMREX_GPU_DEVICE() noexcept {
        for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
                 stride = blockDim.x * gridDim.x;
             icell < udata->ncells; icell += stride) {
          fKernelComputeAJsys(icell, user_data, u_d, udata->csr_val_d);
        }
      });
    *jcurPtr = SUNFALSE;
  } else {
    const auto ec = Gpu::ExecutionConfig(udata->ncells);
    launch_global<<<
      udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
      [=] AMREX_GPU_DEVICE() noexcept {
        for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
                 stride = blockDim.x * gridDim.x;
             icell < udata->ncells; icell += stride) {
          fKernelComputeallAJ(icell, user_data, u_d, udata->csr_val_d);
        }
      });
    *jcurPtr = SUNTRUE;
  }

  cuda_status = cudaStreamSynchronize(udata->stream);
  assert(cuda_status == cudaSuccess);
  BL_PROFILE_VAR_STOP(fKernelComputeAJ);

  BL_PROFILE_VAR("InfoBatched(inPrecond)", InfoBatched);
  cusolver_status = cusolverSpDcsrqrBufferInfoBatched(
    udata->cusolverHandle, udata->neqs_per_cell + 1, udata->neqs_per_cell + 1,
    (udata->NNZ), udata->descrA, udata->csr_val_d, udata->csr_row_count_d,
    udata->csr_col_index_d, udata->ncells, udata->info, &internalDataInBytes,
    &workspaceInBytes);

  assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

  cuda_status = cudaDeviceSynchronize();
  assert(cuda_status == cudaSuccess);

  BL_PROFILE_VAR_STOP(InfoBatched);

  BL_PROFILE_VAR_STOP(Precond);

  return (0);
}
#endif

#ifdef AMREX_USE_CUDA
static int
PSolve(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  realtype gamma,
  realtype delta,
  int lr,
  void* user_data)
{
  BL_PROFILE_VAR("Psolve()", cusolverPsolve);

  cudaError_t cuda_status = cudaSuccess;
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);

  realtype* z_d = N_VGetDeviceArrayPointer_Cuda(z);
  realtype* r_d = N_VGetDeviceArrayPointer_Cuda(r);

  cusolver_status = cusolverSpDcsrqrsvBatched(
    udata->cusolverHandle, udata->neqs_per_cell + 1, udata->neqs_per_cell + 1,
    (udata->NNZ), udata->descrA, udata->csr_val_d, udata->csr_row_count_d,
    udata->csr_col_index_d, r_d, z_d, udata->ncells, udata->info,
    udata->buffer_qr);

  assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

  cuda_status = cudaDeviceSynchronize();
  assert(cuda_status == cudaSuccess);

  N_VCopyFromDevice_Cuda(z);
  N_VCopyFromDevice_Cuda(r);

  BL_PROFILE_VAR_STOP(cusolverPsolve);

  // Checks
  // if (udata->iverbose > 4) {
  //    for(int batchId = 0 ; batchId < udata->ncells; batchId++){
  //        // measure |bj - Aj*xj|
  //        realtype *csrValAj = (udata->csr_val_d) + batchId * (udata->NNZ);
  //        double *xj       = N_VGetHostArrayPointer_Cuda(z) + batchId *
  //        (udata->neqs_per_cell+1); double *bj       =
  //        N_VGetHostArrayPointer_Cuda(r) + batchId * (udata->neqs_per_cell+1);
  //        // sup| bj - Aj*xj|
  //        double sup_res = 0;
  //        for(int row = 0 ; row < (udata->neqs_per_cell+1) ; row++){
  //            printf("\n     row %d: ", row);
  //            const int start = udata->csr_row_count_d[row] - 1;
  //            const int end = udata->csr_row_count_d[row +1] - 1;
  //            double Ax = 0.0; // Aj(row,:)*xj
  //            for(int colidx = start ; colidx < end ; colidx++){
  //                const int col = udata->csr_col_index_d[colidx] - 1;
  //                const double Areg = csrValAj[colidx];
  //                const double xreg = xj[col];
  //                printf("  (%d, %14.8e, %14.8e, %14.8e) ",
  //                col,Areg,xreg,bj[row] ); Ax = Ax + Areg * xreg;
  //            }
  //            double rresidi = bj[row] - Ax;
  //            sup_res = (sup_res > fabs(rresidi))? sup_res : fabs(rresidi);
  //        }
  //        printf("batchId %d: sup|bj - Aj*xj| = %E \n", batchId, sup_res);
  //    }
  //}

  return (0);
}
#endif

#ifdef AMREX_USE_CUDA
// Will not work for cuSolver_sparse_solve right now
static int
cJac(
  realtype t,
  N_Vector y_in,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3)
{
  cudaError_t cuda_status = cudaSuccess;

  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);

  realtype* yvec_d = N_VGetDeviceArrayPointer_Cuda(y_in);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(fy);

  realtype* Jdata;

  BL_PROFILE_VAR("cJac::SparsityStuff", cJacSparsityStuff);
  Jdata = SUNMatrix_cuSparse_Data(J);
  if (
    (SUNMatrix_cuSparse_Rows(J) !=
     (udata->neqs_per_cell + 1) * (udata->ncells)) ||
    (SUNMatrix_cuSparse_Columns(J) !=
     (udata->neqs_per_cell + 1) * (udata->ncells)) ||
    (SUNMatrix_cuSparse_NNZ(J) != udata->ncells * udata->NNZ)) {
    Print() << "Jac error: matrix is wrong size!\n";
    return 1;
  }
  BL_PROFILE_VAR_STOP(cJacSparsityStuff);

  BL_PROFILE_VAR("Jacobian()", fKernelJac);
  const auto ec = Gpu::ExecutionConfig(udata->ncells);
  launch_global<<<
    udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
    [=] AMREX_GPU_DEVICE() noexcept {
      for (int icell = blockDim.x * blockIdx.x + threadIdx.x,
               stride = blockDim.x * gridDim.x;
           icell < udata->ncells; icell += stride) {
        fKernelComputeAJchem(icell, user_data, yvec_d, Jdata);
      }
    });
  cuda_status = cudaStreamSynchronize(udata->stream);
  assert(cuda_status == cudaSuccess);
  BL_PROFILE_VAR_STOP(fKernelJac);

  return (0);
}
#endif

#ifdef AMREX_USE_CUDA
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeAJchem(int ncell, void* user_data, realtype* u_d, realtype* Jdata)
{
  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);

  int u_offset = ncell * (NUM_SPECIES + 1);
  int jac_offset = ncell * (udata->NNZ);

  realtype* u_curr = u_d + u_offset;
  realtype* csr_jac_cell = Jdata + jac_offset;

  Real mw[NUM_SPECIES] = {0.0};
  get_mw(mw);

  Real rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + u_curr[n];
  }

  GpuArray<Real, NUM_SPECIES> massfrac;
  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac[i] = u_curr[i] / rho_pt;
  }

  Real temp_pt = u_curr[NUM_SPECIES];

  int consP;
  if (udata->ireactor_type == 1) {
    consP = 0;
  } else {
    consP = 1;
  }
  GpuArray<Real, (NUM_SPECIES + 1) * (NUM_SPECIES + 1)> Jmat_pt;
  auto eos = pele::physics::PhysicsType::eos();
  eos.RTY2JAC(rho_pt, temp_pt, massfrac.arr, Jmat_pt.arr, consP);

  for (int i = 0; i < udata->neqs_per_cell; i++) {
    for (int k = 0; k < udata->neqs_per_cell; k++) {
      Jmat_pt[k * (udata->neqs_per_cell + 1) + i] =
        Jmat_pt[k * (udata->neqs_per_cell + 1) + i] * mw[i] / mw[k];
    }
    Jmat_pt[i * (udata->neqs_per_cell + 1) + udata->neqs_per_cell] =
      Jmat_pt[i * (udata->neqs_per_cell + 1) + udata->neqs_per_cell] / mw[i];
    Jmat_pt[udata->neqs_per_cell * (udata->neqs_per_cell + 1) + i] =
      Jmat_pt[udata->neqs_per_cell * (udata->neqs_per_cell + 1) + i] * mw[i];
  }
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell + 2; i++) {
    nbVals = udata->csr_row_count_d[i] - udata->csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx_cell = udata->csr_col_index_d[udata->csr_row_count_d[i - 1] + j];
      csr_jac_cell[udata->csr_row_count_d[i - 1] + j] =
        Jmat_pt[idx_cell * (udata->neqs_per_cell + 1) + i - 1];
    }
  }
}
#endif

#ifdef AMREX_USE_CUDA
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeallAJ(
  int ncell, void* user_data, realtype* u_d, double* csr_val_arg)
{
  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);

  Real mw[NUM_SPECIES];
  GpuArray<Real, NUM_SPECIES> massfrac, activity;
  GpuArray<Real, (NUM_SPECIES + 1) * (NUM_SPECIES + 1)> Jmat_pt;
  Real rho_pt, temp_pt;

  int u_offset = ncell * (NUM_SPECIES + 1);
  int jac_offset = ncell * (udata->NNZ);

  realtype* u_curr = u_d + u_offset;
  realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
  realtype* csr_val_cell = csr_val_arg + jac_offset;

  get_mw(mw);

  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + u_curr[n];
  }

  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac[i] = u_curr[i] / rho_pt;
  }

  temp_pt = u_curr[NUM_SPECIES];

  auto eos = pele::physics::PhysicsType::eos();
  eos.RTY2C(rho_pt, temp_pt, massfrac.arr, activity.arr);

  int consP;
  if (udata->ireactor_type == 1) {
    consP = 0;
  } else {
    consP = 1;
  }
  DWDOT_SIMPLIFIED(Jmat_pt.arr, activity.arr, &temp_pt, &consP);

  for (int i = 0; i < udata->neqs_per_cell; i++) {
    for (int k = 0; k < udata->neqs_per_cell; k++) {
      Jmat_pt[k * (udata->neqs_per_cell + 1) + i] =
        Jmat_pt[k * (udata->neqs_per_cell + 1) + i] * mw[i] / mw[k];
    }
    Jmat_pt[i * (udata->neqs_per_cell + 1) + udata->neqs_per_cell] =
      Jmat_pt[i * (udata->neqs_per_cell + 1) + udata->neqs_per_cell] / mw[i];
    Jmat_pt[udata->neqs_per_cell * (udata->neqs_per_cell + 1) + i] =
      Jmat_pt[udata->neqs_per_cell * (udata->neqs_per_cell + 1) + i] * mw[i];
  }
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell + 2; i++) {
    nbVals = udata->csr_row_count_d[i] - udata->csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx =
        udata->csr_col_index_d[udata->csr_row_count_d[i - 1] + j - 1] - 1;
      if (idx == (i - 1)) {
        csr_val_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          1.0 -
          (udata->gamma) * Jmat_pt[idx * (udata->neqs_per_cell + 1) + idx];
        csr_jac_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          Jmat_pt[idx * (udata->neqs_per_cell + 1) + idx];
      } else {
        csr_val_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          -(udata->gamma) * Jmat_pt[idx * (udata->neqs_per_cell + 1) + i - 1];
        csr_jac_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          Jmat_pt[idx * (udata->neqs_per_cell + 1) + i - 1];
      }
    }
  }
}
#endif

#ifdef AMREX_USE_CUDA
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeAJsys(
  int ncell, void* user_data, realtype* u_d, double* csr_val_arg)
{
  CVODEUserData * udata = static_cast<CVODEUserData*>(user_data);

  int jac_offset = ncell * (udata->NNZ);

  realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
  realtype* csr_val_cell = csr_val_arg + jac_offset;

  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell + 2; i++) {
    nbVals = udata->csr_row_count_d[i] - udata->csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx =
        udata->csr_col_index_d[udata->csr_row_count_d[i - 1] + j - 1] - 1;
      if (idx == (i - 1)) {
        csr_val_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          1.0 -
          (udata->gamma) * csr_jac_cell[udata->csr_row_count_d[i - 1] + j - 1];
      } else {
        csr_val_cell[udata->csr_row_count_d[i - 1] + j - 1] =
          -(udata->gamma) * csr_jac_cell[udata->csr_row_count_d[i - 1] + j - 1];
      }
    }
  }
}
#endif

#ifdef AMREX_USE_CUDA
__global__ void
fKernelDenseSolve(
  int ncell,
  realtype* x_d,
  realtype* b_d,
  int subsys_size,
  int subsys_nnz,
  realtype* csr_val)
{
  int stride = blockDim.x * gridDim.x;

  for (int icell = blockDim.x * blockIdx.x + threadIdx.x; icell < ncell;
       icell += stride) {
    int offset = icell * subsys_size;
    int offset_A = icell * subsys_nnz;

    realtype* csr_val_cell = csr_val + offset_A;
    realtype* x_cell = x_d + offset;
    realtype* b_cell = b_d + offset;

    sgjsolve(csr_val_cell, x_cell, b_cell);
  }
}
#endif

#ifdef AMREX_USE_CUDA
SUNLinearSolver
SUNLinSol_dense_custom(N_Vector y, SUNMatrix A, cudaStream_t stream)
{
  if (y == NULL || A == NULL)
    return (NULL);

  if (N_VGetVectorID(y) != SUNDIALS_NVEC_CUDA)
    return (NULL);
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE)
    return (NULL);

  if (N_VGetLength(y) != SUNMatrix_cuSparse_Columns(A))
    return (NULL);

  if (!N_VIsManagedMemory_Cuda(y))
    return (NULL);

  SUNLinearSolver S;
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) {
    return (NULL);
  }

  S->ops->gettype = SUNLinSolGetType_Dense_custom;
  S->ops->setup = SUNLinSolSetup_Dense_custom;
  S->ops->solve = SUNLinSolSolve_Dense_custom;
  S->ops->free = SUNLinSolFree_Dense_custom;

  SUNLinearSolverContent_Dense_custom content;
  content = NULL;
  content = (SUNLinearSolverContent_Dense_custom)malloc(sizeof *content);
  if (content == NULL) {
    SUNLinSolFree(S);
    return (NULL);
  }

  S->content = content;

  content->last_flag = 0;
  content->nsubsys = SUNMatrix_cuSparse_NumBlocks(A);
  content->subsys_size = SUNMatrix_cuSparse_BlockRows(A);
  content->subsys_nnz = SUNMatrix_cuSparse_BlockNNZ(A);
  content->nbBlocks = std::max(1, content->nsubsys / 32);
  content->nbThreads = 32;
  content->stream = stream;

  return (S);
}
#endif

#ifdef AMREX_USE_CUDA
SUNLinearSolver_Type
SUNLinSolGetType_Dense_custom(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}
#endif

#ifdef AMREX_USE_CUDA
int
SUNLinSolSetup_Dense_custom(SUNLinearSolver S, SUNMatrix A)
{
  return (SUNLS_SUCCESS);
}
#endif

#ifdef AMREX_USE_CUDA
int
SUNLinSolSolve_Dense_custom(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol)
{
  cudaError_t cuda_status = cudaSuccess;

  realtype* x_d = N_VGetDeviceArrayPointer_Cuda(x);
  realtype* b_d = N_VGetDeviceArrayPointer_Cuda(b);

  realtype* d_data = SUNMatrix_cuSparse_Data(A);

  BL_PROFILE_VAR("fKernelDenseSolve()", fKernelDenseSolve);
  const auto ec = Gpu::ExecutionConfig(SUN_CUSP_NUM_SUBSYS(S));
  // TODO: why is this AMREX version NOT working ?
  // launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem,
  // SUN_CUSP_STREAM(S)>>>(
  //    [=] AMREX_GPU_DEVICE () noexcept {
  //        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride =
  //        blockDim.x*gridDim.x;
  //            icell < SUN_CUSP_NUM_SUBSYS(S); icell += stride) {
  //            fKernelDenseSolve(icell, x_d, b_d,
  //                  SUN_CUSP_SUBSYS_SIZE(S), SUN_CUSP_SUBSYS_NNZ(S), data_d);
  //        }
  //    });
  fKernelDenseSolve<<<
    SUN_CUSP_NBLOCK(S), SUN_CUSP_NTHREAD(S), ec.sharedMem,
    SUN_CUSP_STREAM(S)>>>(
    SUN_CUSP_NUM_SUBSYS(S), x_d, b_d, SUN_CUSP_SUBSYS_SIZE(S),
    SUN_CUSP_SUBSYS_NNZ(S), d_data);

  cuda_status = cudaStreamSynchronize(SUN_CUSP_STREAM(S));
  assert(cuda_status == cudaSuccess);

  BL_PROFILE_VAR_STOP(fKernelDenseSolve);

  return (SUNLS_SUCCESS);
}
#endif

#ifdef AMREX_USE_CUDA
int
SUNLinSolFree_Dense_custom(SUNLinearSolver S)
{
  if (S == NULL)
    return (SUNLS_SUCCESS);

  if (S->content) {
    free(S->content);
    S->content = NULL;
  }

  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }

  free(S);
  S = NULL;

  return (SUNLS_SUCCESS);
}
#endif

static void
PrintFinalStats(void* cvodeMem)
{
  long lenrw, leniw;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = CVodeGetWorkSpace(cvodeMem, &lenrw, &leniw);
  check_flag(&flag, "CVodeGetWorkSpace", 1);
  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetLinWorkSpace(cvodeMem, &lenrwLS, &leniwLS);
  check_flag(&flag, "CVodeGetLinWorkSpace", 1);
  flag = CVodeGetNumLinIters(cvodeMem, &nli);
  check_flag(&flag, "CVodeGetNumLinIters", 1);
  // flag = CVodeGetNumJacEvals(cvodeMem, &nje);
  // check_flag(&flag, "CVodeGetNumJacEvals", 1);
  flag = CVodeGetNumLinRhsEvals(cvodeMem, &nfeLS);
  check_flag(&flag, "CVodeGetNumLinRhsEvals", 1);

  flag = CVodeGetNumPrecEvals(cvodeMem, &npe);
  check_flag(&flag, "CVodeGetNumPrecEvals", 1);
  flag = CVodeGetNumPrecSolves(cvodeMem, &nps);
  check_flag(&flag, "CVodeGetNumPrecSolves", 1);

  flag = CVodeGetNumLinConvFails(cvodeMem, &ncfl);
  check_flag(&flag, "CVodeGetNumLinConvFails", 1);

#ifdef _OPENMP
  Print() << "\nFinal Statistics: "
          << "(thread:" << omp_get_thread_num() << ", ";
  Print() << "cvodeMem:" << cvodeMem << ")\n";
#else
  Print() << "\nFinal Statistics:\n";
#endif
  Print() << "lenrw      = " << lenrw << "    leniw         = " << leniw
          << "\n";
  Print() << "lenrwLS    = " << lenrwLS << "    leniwLS       = " << leniwLS
          << "\n";
  Print() << "nSteps     = " << nst << "\n";
  Print() << "nRHSeval   = " << nfe << "    nLinRHSeval   = " << nfeLS << "\n";
  Print() << "nnLinIt    = " << nni << "    nLinIt        = " << nli << "\n";
  Print() << "nLinsetups = " << nsetups << "    nErrtf        = " << netf
          << "\n";
  Print() << "nPreceval  = " << npe << "    nPrecsolve    = " << nps << "\n";
  Print() << "nConvfail  = " << ncfn << "    nLinConvfail  = " << ncfl
          << "\n\n";
}
