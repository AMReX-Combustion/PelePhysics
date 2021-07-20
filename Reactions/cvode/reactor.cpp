#include "reactor.H"

#ifdef AMREX_USE_GPU

using namespace amrex;

#define SUN_CUSP_CONTENT(S) ((SUNLinearSolverContent_Dense_custom)(S->content))
#define SUN_CUSP_SUBSYS_SIZE(S) (SUN_CUSP_CONTENT(S)->subsys_size)
#define SUN_CUSP_NUM_SUBSYS(S) (SUN_CUSP_CONTENT(S)->nsubsys)
#define SUN_CUSP_SUBSYS_NNZ(S) (SUN_CUSP_CONTENT(S)->subsys_nnz)

#define SUN_CUSP_LASTFLAG(S) (SUN_CUSP_CONTENT(S)->last_flag)
#define SUN_CUSP_STREAM(S) (SUN_CUSP_CONTENT(S)->stream)
#define SUN_CUSP_NBLOCK(S) (SUN_CUSP_CONTENT(S)->nbBlocks)
#define SUN_CUSP_NTHREAD(S) (SUN_CUSP_CONTENT(S)->nbThreads)

int sparse_solve = 1;
int sparse_cusolver_solve = 5;
int iterative_gmres_solve = 99;
int eint_rho = 1;
int enth_rho = 2;

Array<Real, NUM_SPECIES + 1> typVals = {-1};
Real relTol = 1.0e-10;
Real absTol = 1.0e-10;

void
SetTypValsODE(const std::vector<double>& ExtTypVals)
{
  Vector<std::string> kname;
  pele::physics::eos::speciesNames(kname);

  Print() << "Set the typVals in PelePhysics: \n  ";
  int size_ETV = ExtTypVals.size();
  AMREX_ASSERT(size_ETV == typVals.size());
  for (int i = 0; i < size_ETV - 1; i++) {
    typVals[i] = ExtTypVals[i];
    Print() << kname[i] << ":" << typVals[i] << "  ";
  }
  typVals[size_ETV - 1] = ExtTypVals[size_ETV - 1];
  Print() << "Temp:" << typVals[size_ETV - 1] << " \n";
}

void
SetTolFactODE(double relative_tol, double absolute_tol)
{
  relTol = relative_tol;
  absTol = absolute_tol;
  Print() << "Set RTOL, ATOL = " << relTol << " " << absTol
          << " in PelePhysics\n";
}

int
reactor_init(int reactor_type, int Ncells)
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
  sparse_solve          = 1;
  sparse_cusolver_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {
    isolve_type = sparse_solve;
  } else if (solve_type_str == "sparse") {
    isolve_type = sparse_cusolver_solve;
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

  } else if (isolve_type == sparse_solve) {
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

  } else if (isolve_type == sparse_cusolver_solve) {
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
  sparse_solve          = 1;
  sparse_cusolver_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {
    isolve_type = sparse_solve;
  } else if (solve_type_str == "sparse") {
    isolve_type = sparse_cusolver_solve;
  } else if (solve_type_str == "GMRES") {
    isolve_type = iterative_gmres_solve;
  }

  NCELLS = box.numPts();
  neq_tot = (NUM_SPECIES + 1) * NCELLS;

  UserData user_data;

  BL_PROFILE_VAR("AllocsInCVODE", AllocsCVODE);
  user_data = (CVodeUserData*)The_Arena()->alloc(sizeof(struct CVodeUserData));
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

    } else if (isolve_type == sparse_cusolver_solve) {
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
      icell, i, j, k, user_data->ireactor_type, rY_in, rY_src_in, T_in,
      rEner_in, rEner_src_in, yvec_d, user_data->rYsrc_d,
      user_data->rhoe_init_d, user_data->rhoesrc_ext_d);
  });
  BL_PROFILE_VAR_STOP(FlatStuff);

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
  } else if (user_data->isolve_type == sparse_cusolver_solve) {
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

  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

  BL_PROFILE_VAR_START(FlatStuff);
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);

    box_unflatten(
      icell, i, j, k, user_data->ireactor_type, rY_in, T_in, rEner_in,
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
    } else if (user_data->isolve_type == sparse_cusolver_solve) {
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

  UserData udata = static_cast<CVodeUserData*>(user_data);
  udata->dt_save = t;

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
          icell, udata->dt_save, udata->ireactor_type, yvec_d, ydot_d,
          udata->rhoe_init_d, udata->rhoesrc_ext_d, udata->rYsrc_d);
      }
    });
  Gpu::Device::streamSynchronize();
#else
  for (int icell = 0; icell < udata->ncells; icell++) {
    fKernelSpec(
      icell, udata->dt_save, udata->ireactor_type, yvec_d, ydot_d,
      udata->rhoe_init_d, udata->rhoesrc_ext_d, udata->rYsrc_d);
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
  sparse_solve          = 1;
  sparse_cusolver_solve = 5;
  iterative_gmres_solve = 99;
  */
  int isolve_type = -1;
  if (solve_type_str == "sparse_custom") {        // Direct solve using custom GaussJordan routine
    isolve_type = sparse_solve;
  } else if (solve_type_str == "sparse") {        // Direct sparse solve using cuSolver (CUDA only !)
    isolve_type = sparse_cusolver_solve;
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
  UserData user_data;
  user_data = (CVodeUserData*)The_Arena()->alloc(sizeof(struct CVodeUserData));
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
  } else if (user_data->isolve_type == sparse_cusolver_solve) {
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

  UserData udata = static_cast<CVodeUserData*>(user_data);
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

  UserData udata = static_cast<CVodeUserData*>(user_data);

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

  UserData udata = static_cast<CVodeUserData*>(user_data);

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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelSpec(
  int icell,
  double dt_save,
  int reactor_type,
  realtype* yvec_d,
  realtype* ydot_d,
  double* rhoe_init,
  double* rhoesrc_ext,
  double* rYs)
{
  int offset = icell * (NUM_SPECIES + 1);

  Real mw[NUM_SPECIES] = {0.0};
  get_mw(mw);

  Real rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + yvec_d[offset + n];
  }

  GpuArray<Real, NUM_SPECIES> massfrac;
  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac[i] = yvec_d[offset + i] / rho_pt;
  }

  Real nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) / rho_pt;

  Real temp_pt = yvec_d[offset + NUM_SPECIES];

  GpuArray<Real, NUM_SPECIES> ei_pt;
  Real Cv_pt = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == 1) {
    eos.EY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.TY2Cv(temp_pt, massfrac.arr, Cv_pt);
    eos.T2Ei(temp_pt, ei_pt.arr);
  } else {
    eos.HY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.TY2Cp(temp_pt, massfrac.arr, Cv_pt);
    eos.T2Hi(temp_pt, ei_pt.arr);
  }

  GpuArray<Real, NUM_SPECIES> cdots_pt;
  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
  for (int i = 0; i < NUM_SPECIES; i++) {
    ydot_d[offset + i] = cdots_pt[i] + rYs[icell * NUM_SPECIES + i];
    ydot_d[offset + NUM_SPECIES] =
      ydot_d[offset + NUM_SPECIES] - ydot_d[offset + i] * ei_pt[i];
  }
  ydot_d[offset + NUM_SPECIES] =
    ydot_d[offset + NUM_SPECIES] / (rho_pt * Cv_pt);
}

#ifdef AMREX_USE_CUDA
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeAJchem(int ncell, void* user_data, realtype* u_d, realtype* Jdata)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

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
  UserData udata = static_cast<CVodeUserData*>(user_data);

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
  UserData udata = static_cast<CVodeUserData*>(user_data);

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

// Check function return value...
//     opt == 0 means SUNDIALS function allocates memory so check if
//              returned NULL pointer
//     opt == 1 means SUNDIALS function returns a flag so check if
//              flag >= 0
//     opt == 2 means function allocates memory so check if returned
//              NULL pointer
static int
check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == NULL) {
    if (ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      Abort("abort");
    }
    return (1);
  } else if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      if (ParallelDescriptor::IOProcessor()) {
        fprintf(
          stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
          *errflag);
        Abort("abort");
      }
      return (1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    if (ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      Abort("abort");
    }
    return (1);
  }

  return (0);
}

#else

#define SUN_CUSP_CONTENT(S) ((SUNLinearSolverContent_Sparse_custom)(S->content))
#define SUN_CUSP_REACTYPE(S) (SUN_CUSP_CONTENT(S)->reactor_type)
#define SUN_CUSP_NUM_SUBSYS(S) (SUN_CUSP_CONTENT(S)->nsubsys)
#define SUN_CUSP_SUBSYS_NNZ(S) (SUN_CUSP_CONTENT(S)->subsys_nnz)
#define SUN_CUSP_SUBSYS_SIZE(S) (SUN_CUSP_CONTENT(S)->subsys_size)

using namespace amrex;

N_Vector y = NULL;
SUNLinearSolver LS = NULL;
SUNMatrix A = NULL;
void* cvode_mem = NULL;
UserData data = NULL;
Real time_init = 0.0;
Array<double, NUM_SPECIES + 1> typVals = {-1};
double relTol = 1.0e-10;
double absTol = 1.0e-10;
int dense_solve = 1;
int sparse_solve = 5;
int iterative_gmres_solve = 99;
int sparse_solve_custom = 101;
int iterative_gmres_solve_custom = 199;
int hack_dump_sparsity_pattern = -5;
int eint_rho = 1; // in/out = rhoE/rhoY
int enth_rho = 2; // in/out = rhoH/rhoY

#ifdef _OPENMP
#pragma omp threadprivate(y, LS, A)
#pragma omp threadprivate(cvode_mem, data)
#pragma omp threadprivate(time_init)
#pragma omp threadprivate(typVals)
#pragma omp threadprivate(relTol, absTol)
#endif

void
SetTypValsODE(const std::vector<double>& ExtTypVals)
{
  int size_ETV = (NUM_SPECIES + 1);
  Vector<std::string> kname;
  pele::physics::eos::speciesNames(kname);
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV - 1; i++) {
    typVals[i] = ExtTypVals[i];
  }
  typVals[size_ETV - 1] = ExtTypVals[size_ETV - 1];
  if (omp_thread == 0) {
    Print() << "Set the typVals in PelePhysics: \n  ";
    for (int i = 0; i < size_ETV - 1; i++) {
      Print() << kname[i] << ":" << typVals[i] << "  ";
    }
    Print() << "Temp:" << typVals[size_ETV - 1] << " \n";
  }
}

void
SetTolFactODE(double relative_tol, double absolute_tol)
{
  relTol = relative_tol;
  absTol = absolute_tol;
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  if (omp_thread == 0) {
    Print() << "Set RTOL, ATOL = " << relTol << " " << absTol
            << " in PelePhysics\n";
  }
}

void
ReSetTolODE()
{
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    data != NULL, "Reactor object is not initialized !!");

  int neq_tot = (NUM_SPECIES + 1) * data->ncells;
  N_Vector atol = N_VNew_Serial(neq_tot);
  realtype* ratol = N_VGetArrayPointer(atol);

  if (typVals[0] > 0) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "Setting CVODE tolerances rtol = " << relTol
              << " atolfact = " << absTol << " in PelePhysics \n";
    }
    for (int i = 0; i < data->ncells; i++) {
      int offset = i * (NUM_SPECIES + 1);
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        ratol[offset + k] = typVals[k] * absTol;
      }
    }
  } else {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "Setting CVODE tolerances rtol = " << relTol
              << " atol = " << absTol << " in PelePhysics \n";
    }
    for (int i = 0; i < neq_tot; i++) {
      ratol[i] = absTol;
    }
  }
  // Call CVodeSVtolerances to specify the scalar relative tolerance and vector
  // absolute tolerances
  int flag = CVodeSVtolerances(cvode_mem, relTol, atol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) {
    Abort("Problem in ReSetTolODE");
  }

  N_VDestroy(atol);
}

void
cvodeErrHandler(
  int error_code,
  const char* module,
  const char* function,
  char* msg,
  void* eh_data)
{
  if (error_code != CV_WARNING) {
    std::cout << "From CVODE: " << msg << std::endl;
    Abort("Bad CVODE");
  }
}

// Initialization routine, called once at the begining of the problem
int
reactor_init(int reactor_type, int ode_ncells)
{

  BL_PROFILE_VAR("reactInit", reactInit);

  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif
  // Total number of eq to integrate
  int neq_tot = (NUM_SPECIES + 1) * ode_ncells;

  y = N_VNew_Serial(neq_tot);
  if (check_flag((void*)y, "N_VNew_Serial", 0))
    return (1);

  // Call CVodeCreate to create the solver memory and specify the Backward
  // Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_flag((void*)cvode_mem, "CVodeCreate", 0))
    return (1);

  // Does not work for more than 1 cell right now
  data = AllocUserData(reactor_type, ode_ncells);
  if (check_flag((void*)data, "AllocUserData", 2))
    return (1);

  // Number of species and cells in mechanism
  if ((data->iverbose > 0) && (omp_thread == 0)) {
    Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
    Print() << "Number of cells in one solve is " << data->ncells << "\n";
  }

  // Set the pointer to user-defined data
  int flag = CVodeSetUserData(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetUserData", 1))
    return (1);

  realtype time = 0.0e+0;
  // Call CVodeInit to initialize the integrator memory and specify the user's
  // right hand side function, the inital time, and initial dependent variable
  // vector y.
  flag = CVodeInit(cvode_mem, cF_RHS, time, y);
  if (check_flag(&flag, "CVodeInit", 1))
    return (1);

  // Definition of tolerances: one for each species
  N_Vector atol = N_VNew_Serial(neq_tot);
  realtype* ratol = N_VGetArrayPointer(atol);
  if (typVals[0] > 0) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "Setting CVODE tolerances rtol = " << relTol
              << " atolfact = " << absTol << " in PelePhysics \n";
    }
    for (int i = 0; i < data->ncells; i++) {
      int offset = i * (NUM_SPECIES + 1);
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        // ratol[offset + k] = std::max(typVals[k]*absTol,relTol);
        ratol[offset + k] = typVals[k] * absTol;
      }
    }
  } else {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "Setting CVODE tolerances rtol = " << relTol
              << " atol = " << absTol << " in PelePhysics \n";
    }
    for (int i = 0; i < neq_tot; i++) {
      ratol[i] = absTol;
    }
  }
  // Call CVodeSVtolerances to specify the scalar relative tolerance and vector
  // absolute tolerances
  flag = CVodeSVtolerances(cvode_mem, relTol, atol);
  if (check_flag(&flag, "CVodeSVtolerances", 1))
    return (1);

  // flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-1);
  // if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

  flag = CVodeSetMaxNonlinIters(cvode_mem, 50);
  if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1))
    return (1);

  flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
  if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1))
    return (1);

  flag = CVodeSetErrHandlerFn(cvode_mem, cvodeErrHandler, 0);
  if (check_flag(&flag, "CVodeSetErrHandlerFn", 1))
    return (1);

  if (data->isolve_type == dense_solve) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "\n--> Using a Direct Dense Solver\n";
    }
    // Create dense SUNMatrix for use in linear solves
    A = SUNDenseMatrix(neq_tot, neq_tot);
    if (check_flag((void*)A, "SUNDenseMatrix", 0))
      return (1);

    // Create dense SUNLinearSolver object for use by CVode
    LS = SUNDenseLinearSolver(y, A);
    if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
      return (1);

    // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVDlsSetLinearSolver", 1))
      return (1);

  } else if (data->isolve_type == sparse_solve_custom) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "\n--> Using a custom Direct Sparse Solver\n";
    }
    // Create dense SUNMatrix for use in linear solves
    A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ) * data->ncells, CSR_MAT);
    if (check_flag((void*)A, "SUNDenseMatrix", 0))
      return (1);

    // Create dense SUNLinearSolver object for use by CVode
    LS = SUNLinSol_sparse_custom(
      y, A, reactor_type, data->ncells, (NUM_SPECIES + 1), data->NNZ);
    if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
      return (1);

    // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVDlsSetLinearSolver", 1))
      return (1);

  } else if (data->isolve_type == sparse_solve) {
#ifdef USE_KLU_PP
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "\n--> Using a Direct Sparse Solver\n";
    }
    // Create sparse SUNMatrix for use in linear solves
    A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ) * data->ncells, CSC_MAT);
    if (check_flag((void*)A, "SUNSparseMatrix", 0))
      return (1);

    // Create KLU solver object for use by CVode
    LS = SUNLinSol_KLU(y, A);
    if (check_flag((void*)LS, "SUNLinSol_KLU", 0))
      return (1);

    // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
#else
    Abort("Sparse solver not valid without KLU solver.");
#endif

  } else if (
    (data->isolve_type == iterative_gmres_solve) ||
    (data->isolve_type == iterative_gmres_solve_custom)) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "\n--> Using an Iterative Solver (" << data->isolve_type
              << ")\n";
    }

    // Create the linear solver object
    if (data->ianalytical_jacobian == 0) {
      LS = SUNLinSol_SPGMR(y, PREC_NONE, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    } else {
      LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
      if (check_flag((void*)LS, "SUNDenseLinearSolver", 0))
        return (1);
    }

    // Set CVSpils linear solver to LS
    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    if (check_flag(&flag, "CVSpilsSetLinearSolver", 1))
      return (1);
  } else {
    Abort("Wrong choice of linear solver...");
  }

  if (data->ianalytical_jacobian == 0) {
    if ((data->iverbose > 0) && (omp_thread == 0)) {
      Print() << "    Without Analytical J/Preconditioner\n";
    }
#ifdef USE_KLU_PP
    if (data->isolve_type == sparse_solve) {
      Abort("Sparse Solver requires an Analytical J");
    }
#endif
    if (data->isolve_type == sparse_solve_custom) {
      Abort("Custom sparse solver requires an Analytical J");
    }
  } else {
    if (data->isolve_type == iterative_gmres_solve_custom) {
      // Set the JAcobian-times-vector function
      flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
      if (check_flag(&flag, "CVSpilsSetJacTimes", 1))
        return (1);

      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With a custom Sparse Preconditioner\n";
      }
      // Set the preconditioner solve and setup functions
      flag = CVSpilsSetPreconditioner(cvode_mem, Precond_custom, PSolve_custom);
      if (check_flag(&flag, "CVSpilsSetPreconditioner", 1))
        return (1);

    } else if (data->isolve_type == iterative_gmres_solve) {
      // Set the JAcobian-times-vector function
      flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
      if (check_flag(&flag, "CVSpilsSetJacTimes", 1))
        return (1);
#ifdef USE_KLU_PP
      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With a Sparse Preconditioner\n";
      }
      // Set the preconditioner solve and setup functions
      flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
      if (check_flag(&flag, "CVSpilsSetPreconditioner", 1))
        return (1);
#else
      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With a Preconditioner\n";
      }
      // Set the preconditioner solve and setup functions
      flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
      if (check_flag(&flag, "CVSpilsSetPreconditioner", 1))
        return (1);
#endif
#ifdef USE_KLU_PP
    } else if (data->isolve_type == sparse_solve) {
      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With a Sparse Analytical J\n";
      }
      // Set the user-supplied Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cJac_KLU);
      if (check_flag(&flag, "CVodeSetJacFn", 1))
        return (1);
#endif
    } else if (data->isolve_type == dense_solve) {
      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With Analytical J\n";
      }
      // Set the user-supplied Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cJac);
      if (check_flag(&flag, "CVodeSetJacFn", 1))
        return (1);

    } else if (data->isolve_type == sparse_solve_custom) {
      if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "    With a Sparse Analytical J\n";
      }
      // Set the user-supplied Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cJac_sps);
      if (check_flag(&flag, "CVodeSetJacFn", 1))
        return (1);
    }
  }

  // Set the max number of time steps
  flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    return (1);

  // Set the max order
  flag = CVodeSetMaxOrd(cvode_mem, 2);
  if (check_flag(&flag, "CVodeSetMaxOrd", 1))
    return (1);

  // Set the num of steps to wait inbetween Jac evals
  flag = CVodeSetJacEvalFrequency(cvode_mem, 100);
  if (check_flag(&flag, "CVodeSetJacEvalFrequency", 1))
    return (1);

  N_VDestroy(atol);

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "\n--> DONE WITH INITIALIZATION (CPU)" << data->ireactor_type
            << "\n";
  }

  data->reactor_cvode_initialized = true;

  BL_PROFILE_VAR_STOP(reactInit);

  return (0);
}

// Main routine for CVode integration: integrate a Box version
int
react(
  const Box& box,
  Array4<Real> const& rY_in,
  Array4<Real> const& rY_src_in,
  Array4<Real> const& T_in,
  Array4<Real> const& rEner_in,
  Array4<Real> const& rEner_src_in,
  Array4<Real> const& FC_in,
  Array4<int> const& mask,
  Real& dt_react,
  Real& time,
  const int& reactor_type)
{

  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "\n -------------------------------------\n";
  }

  // Initial time and time to reach after integration
  time_init = time;

  if ((data->iverbose > 3) && (omp_thread == 0)) {
    Print() << "BEG : time curr is " << time_init << " and dt_react is "
            << dt_react << " and final time should be " << time_init + dt_react
            << "\n";
  }

  if (data->ncells != 1) {
    Abort("CVODE react can only integrate one cell at a time");
  }
  int box_ncells = box.numPts();
  data->boxcell = 0;

  if ((data->iverbose > 2) && (omp_thread == 0)) {
    Print() << "Ncells in the box = " << box_ncells << "\n";
  }

  // Perform integration one cell at a time
  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (mask(i, j, k) != -1) {
      Real mass_frac[NUM_SPECIES];
      Real rho = 0.0;
      Real rho_inv;
      Real Enrg_loc;
      Real temp;

      realtype* yvec_d = N_VGetArrayPointer(y);

      BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
      for (int n = 0; n < NUM_SPECIES; n++) {
        yvec_d[n] = rY_in(i, j, k, n);
        data->rYsrc[n] = rY_src_in(i, j, k, n);
        rho += yvec_d[n];
      }
      rho_inv = 1.0 / rho;
      temp = T_in(i, j, k, 0);
      data->rhoX_init[0] = rEner_in(i, j, k, 0);
      data->rhoXsrc_ext[0] = rEner_src_in(i, j, k, 0);

      // T update with energy and Y
      for (int n = 0; n < NUM_SPECIES; n++) {
        mass_frac[n] = yvec_d[n] * rho_inv;
      }
      Enrg_loc = data->rhoX_init[0] / rho;
      auto eos = pele::physics::PhysicsType::eos();
      if (data->ireactor_type == 1) {
        eos.EY2T(Enrg_loc, mass_frac, temp);
      } else {
        eos.HY2T(Enrg_loc, mass_frac, temp);
      }
      yvec_d[NUM_SPECIES] = temp;
      BL_PROFILE_VAR_STOP(FlatStuff);

      // ReInit CVODE is faster
      CVodeReInit(cvode_mem, time_init, y);

      // Time to reach after integration
      Real time_out_lcl = time_init + dt_react;

      // Integration
      Real dummy_time;
      BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
      CVode(cvode_mem, time_out_lcl, y, &dummy_time, CV_NORMAL);
      BL_PROFILE_VAR_STOP(AroundCVODE);

      if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() << "Additional verbose info --\n";
        PrintFinalStats(cvode_mem, yvec_d[NUM_SPECIES]);
        Print() << "\n -------------------------------------\n";
      }

      // Get estimate of how hard the integration process was
      long int nfe, nfeLS;
      CVodeGetNumRhsEvals(cvode_mem, &nfe);
      CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
      FC_in(i, j, k, 0) = nfe + nfeLS;

      BL_PROFILE_VAR_START(FlatStuff);
      rho = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        rY_in(i, j, k, n) = yvec_d[n];
        rho += yvec_d[n];
      }
      rho_inv = 1.0 / rho;
      temp = yvec_d[NUM_SPECIES];

      // T update with energy and Y
      for (int n = 0; n < NUM_SPECIES; n++) {
        mass_frac[n] = yvec_d[n] * rho_inv;
      }
      rEner_in(i, j, k, 0) =
        data->rhoX_init[0] + (dummy_time - time_init) * data->rhoXsrc_ext[0];
      Enrg_loc = rEner_in(i, j, k, 0) * rho_inv;
      if (data->ireactor_type == 1) {
        eos.EY2T(Enrg_loc, mass_frac, temp);
      } else {
        eos.HY2T(Enrg_loc, mass_frac, temp);
      }
      T_in(i, j, k, 0) = temp;
      BL_PROFILE_VAR_STOP(FlatStuff);

      if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() << "END : time curr is " << dummy_time
                << " and actual dt_react is " << (dummy_time - time_init)
                << "\n";
      }
    } else {
      FC_in(i, j, k, 0) = 0.0;
    }
  });

  // Update dt_react with real time step taken ...
  // should be very similar to input dt_react
  // dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
  // If reactor mode is activated, update time to perform subcycling
  time = time_init + dt_react;
#endif

  // Get estimate of how hard the integration process was
  return 20;
}

// Main routine for CVode integration: integrate a Box version 2
int
react_2(
  const Box& box,
  Array4<Real> const& rY_in,
  Array4<Real> const& rY_src_in,
  Array4<Real> const& T_in,
  Array4<Real> const& rEner_in,
  Array4<Real> const& rEner_src_in,
  Array4<Real> const& FC_in,
  Array4<int> const& mask_in,
  Real& dt_react,
  Real& time)
{

  realtype dummy_time;
  int flag;
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "\n -------------------------------------\n";
  }

  // Initial time and time to reach after integration
  time_init = time;
  realtype time_out = time + dt_react;

  if ((data->iverbose > 3) && (omp_thread == 0)) {
    Print() << "BEG : time curr is " << time_init << " and dt_react is "
            << dt_react << " and final time should be " << time_out << "\n";
  }

  // Define full box_ncells length vectors to be integrated piece by piece by
  // CVode
  int box_ncells = box.numPts();
  if ((data->iverbose > 2) && (omp_thread == 0)) {
    Print() << "Ncells in the box = " << box_ncells << "\n";
  }
  BL_PROFILE_VAR("reactor::ExtForcingAlloc", ExtForcingAlloc);
  (data->Yvect_full) = new amrex::Real[box_ncells * (NUM_SPECIES + 1)];
  (data->rYsrc) = new amrex::Real[box_ncells * (NUM_SPECIES)];
  (data->rhoX_init) = new amrex::Real[box_ncells];
  (data->rhoXsrc_ext) = new amrex::Real[box_ncells];
  (data->FCunt) = new int[box_ncells];
  (data->mask) = new int[box_ncells];
  BL_PROFILE_VAR_STOP(ExtForcingAlloc);

  BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
  // Fill the full box_ncells length vectors from input Array4
  const auto len = length(box);
  const auto lo = lbound(box);
  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    box_flatten(
      icell, i, j, k, data->ireactor_type, rY_in, rY_src_in, T_in, rEner_in,
      rEner_src_in, data->Yvect_full, data->rYsrc, data->rhoX_init,
      data->rhoXsrc_ext);
  });
  BL_PROFILE_VAR_STOP(FlatStuff);

  BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
  BL_PROFILE_VAR_STOP(AroundCVODE);

  // We may need extra cells to fill the fixed data->ncells in this case since
  // we do not Init each time
  int extra_cells = box_ncells - box_ncells / (data->ncells) * (data->ncells);
  if ((data->iverbose > 2) && (omp_thread == 0)) {
    Print() << " Extra cells = " << extra_cells << "\n";
  }

  // Integrate data->ncells at a time with CVode The extra cell machinery is not
  // ope yet and most likely produce out of bound errors
  realtype* yvec_d = N_VGetArrayPointer(y);
  long int nfe, nfeLS;
  for (int i = 0; i < box_ncells + extra_cells; i += data->ncells) {
    if (data->mask[i] != -1) //   TODO: this doesn't really work, but we are not
                             //   using react_2, let alone react_2 + EB.
    {
      // Print() <<" dealing with cell " << i <<  "\n";
      int offset = i * (NUM_SPECIES + 1);
      data->boxcell = i;
      for (int k = 0; k < data->ncells * (NUM_SPECIES + 1); k++) {
        yvec_d[k] = data->Yvect_full[offset + k];
      }

      // ReInit CVODE is faster
      CVodeReInit(cvode_mem, time_init, y);

      BL_PROFILE_VAR_START(AroundCVODE);
      // Integration
      flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
      if (check_flag(&flag, "CVode", 1))
        return (1);
      BL_PROFILE_VAR_STOP(AroundCVODE);

      // Update full box length vector
      for (int k = 0; k < data->ncells * (NUM_SPECIES + 1); k++) {
        data->Yvect_full[offset + k] = yvec_d[k];
      }

      // Get estimate of how hard the integration process was
      flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
      flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
      for (int k = 0; k < data->ncells; k++) {
        data->FCunt[i + k] = nfe + nfeLS;
      }

      if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() << "END : time curr is " << dummy_time
                << " and actual dt_react is " << (dummy_time - time_init)
                << "\n";
      }
    } else {
      for (int k = 0; k < data->ncells; k++) {
        data->FCunt[i + k] = 0.0;
      }
    }
  }

#ifdef MOD_REACTOR
  // If reactor mode is activated, update time to perform subcycling
  time = time_init + dt_react;
#endif

  BL_PROFILE_VAR_START(FlatStuff);
  // Update the input/output Array4 rY_in and rEner_in
  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int icell = (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
    box_unflatten(
      icell, i, j, k, data->ireactor_type, rY_in, T_in, rEner_in, rEner_src_in,
      FC_in, data->Yvect_full, data->rhoX_init, nfe, dt_react);
  });
  BL_PROFILE_VAR_STOP(FlatStuff);

  BL_PROFILE_VAR_START(ExtForcingAlloc);
  delete[](data->Yvect_full);
  delete[](data->rYsrc);
  delete[](data->rhoX_init);
  delete[](data->rhoXsrc_ext);
  delete[](data->FCunt);
  delete[](data->mask);
  BL_PROFILE_VAR_STOP(ExtForcingAlloc);

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "Additional verbose info --\n";
    PrintFinalStats(cvode_mem, yvec_d[NUM_SPECIES]);
    Print() << "\n -------------------------------------\n";
  }

  // Get estimate of how hard the integration process was
  //long int nfe, nfeLS;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  return nfe + nfeLS;
}

// Main routine for CVode integration: classic version
int
react(
  realtype* rY_in,
  realtype* rY_src_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int reactor_type,
  int Ncells)
{

  realtype dummy_time;
  int flag;
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "\n -------------------------------------\n";
  }

  // Initial time and time to reach after integration
  time_init = time;
  realtype time_out = time + dt_react;

  if ((data->iverbose > 3) && (omp_thread == 0)) {
    Print() << "BEG : time curr is " << time_init << " and dt_react is "
            << dt_react << " and final time should be " << time_out << "\n";
  }

  // Define full box_ncells length vectors to be integrated piece by piece by
  // CVode
  if ((data->iverbose > 2) && (omp_thread == 0)) {
    Print() << "Ncells in the box = " << Ncells << "\n";
  }

  BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
  // Get Device MemCpy of in arrays
  // Get Device pointer of solution vector
  realtype* yvec_d = N_VGetArrayPointer(y);
  // rhoY,T
  std::memcpy(yvec_d, rY_in, sizeof(Real) * ((NUM_SPECIES + 1) * Ncells));
  // rhoY_src_ext
  std::memcpy(
    data->rYsrc, rY_src_in, sizeof(Real) * (NUM_SPECIES * Ncells));
  // rhoE/rhoH
  std::memcpy(data->rhoX_init, rX_in, sizeof(Real) * Ncells);
  std::memcpy(data->rhoXsrc_ext, rX_src_in, sizeof(Real) * Ncells);
  BL_PROFILE_VAR_STOP(FlatStuff);

  // Check if y is within physical bounds we may remove that eventually
  check_state(y);
  if (!(data->actual_ok_to_react)) {
#ifdef MOD_REACTOR
    // If reactor mode is activated, update time
    time = time_out;
#endif
    return 0;
  }

  BL_PROFILE_VAR_START(FlatStuff);
  // T update with energy and Y
  int offset;
  realtype rho, rho_inv, nrg_loc, temp;
  for (int i = 0; i < Ncells; i++) {
    offset = i * (NUM_SPECIES + 1);
    realtype* mass_frac = rY_in + offset;
    // get rho
    rho = 0;
    for (int kk = 0; kk < NUM_SPECIES; kk++) {
      rho += mass_frac[kk];
    }
    rho_inv = 1 / rho;
    // get Yks
    for (int kk = 0; kk < NUM_SPECIES; kk++) {
      mass_frac[kk] = mass_frac[kk] * rho_inv;
    }
    // get energy
    nrg_loc = rX_in[i] * rho_inv;
    // recompute T
    temp = rY_in[offset + NUM_SPECIES];
    auto eos = pele::physics::PhysicsType::eos();
    if (reactor_type == eint_rho) {
      eos.EY2T(nrg_loc, mass_frac, temp);
    } else {
      eos.HY2T(nrg_loc, mass_frac, temp);
    }
    // store T in y
    yvec_d[offset + NUM_SPECIES] = temp;
  }
  BL_PROFILE_VAR_STOP(FlatStuff);

  // ReInit CVODE is faster
  CVodeReInit(cvode_mem, time_init, y);

  // There should be no internal looping of CVOde
  data->boxcell = 0;

  BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
  flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
  // ONE STEP MODE FOR DEBUGGING
  // flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
  if (check_flag(&flag, "CVode", 1))
    return (1);
  BL_PROFILE_VAR_STOP(AroundCVODE);

  // Update dt_react with real time step taken ...  should be very similar to
  // input dt_react
  dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
  // If reactor mode is activated, update time
  time = time_init + dt_react;
#endif

  if ((data->iverbose > 3) && (omp_thread == 0)) {
    Print() << "END : time curr is " << dummy_time << " and actual dt_react is "
            << dt_react << "\n";
  }

  BL_PROFILE_VAR_START(FlatStuff);
  // Pack data to return in main routine external
  std::memcpy(
    rY_in, yvec_d, ((NUM_SPECIES + 1) * Ncells) * sizeof(realtype));
  for (int i = 0; i < Ncells; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }

  // T update with energy and Y
  for (int i = 0; i < Ncells; i++) {
    offset = i * (NUM_SPECIES + 1);
    realtype* mass_frac = yvec_d + offset;
    // get rho
    rho = 0;
    for (int kk = 0; kk < NUM_SPECIES; kk++) {
      rho += mass_frac[kk];
    }
    rho_inv = 1 / rho;
    // get Yks
    for (int kk = 0; kk < NUM_SPECIES; kk++) {
      mass_frac[kk] = mass_frac[kk] * rho_inv;
    }
    // get energy
    nrg_loc = rX_in[i] * rho_inv;
    // recompute T
    auto eos = pele::physics::PhysicsType::eos();
    if (reactor_type == eint_rho) {
      eos.EY2T(nrg_loc, mass_frac, temp);
    } else {
      eos.HY2T(nrg_loc, mass_frac, temp);
    }
    // store T in rY_in
    rY_in[offset + NUM_SPECIES] = temp;
  }
  BL_PROFILE_VAR_STOP(FlatStuff);

  if ((data->iverbose > 1) && (omp_thread == 0)) {
    Print() << "Additional verbose info --\n";
    PrintFinalStats(cvode_mem, rY_in[NUM_SPECIES]);
    Print() << "\n -------------------------------------\n";
  }

  // Get estimate of how hard the integration process was
  long int nfe, nfeLS;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  return nfe + nfeLS;
}

// RHS routine
int
cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{

  realtype* y_d = N_VGetArrayPointer(y_in);
  realtype* ydot_d = N_VGetArrayPointer(ydot_in);

  BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);
  fKernelSpec(&t, y_d, ydot_d, user_data);
  BL_PROFILE_VAR_STOP(fKernelSpec);

  return (0);
}

// RHS source terms evaluation
void
fKernelSpec(realtype* t, realtype* yvec_d, realtype* ydot_d, void* user_data)
{
  // Make local copies of pointers in user_data (cell M)
  UserData data_wk = (UserData)user_data;

  // Loop on packed cells
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    // Tmp vars
    realtype massfrac[NUM_SPECIES];
    realtype Xi[NUM_SPECIES];
    realtype cdot[NUM_SPECIES], molecular_weight[NUM_SPECIES];
    realtype cX;
    realtype temp, energy;
    realtype dt;

    // dt is curr time - time init
    dt = *t - time_init;

    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);

    // MW CGS
    CKWT(molecular_weight);

    // rho MKS
    realtype rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + yvec_d[offset + i];
    }

    temp = yvec_d[offset + NUM_SPECIES];

    // Yks
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = yvec_d[offset + i] / rho;
    }

    // NRG CGS
    energy = (data_wk->rhoX_init[data->boxcell + tid] +
              data_wk->rhoXsrc_ext[data_wk->boxcell + tid] * dt) /
             rho;

    auto eos = pele::physics::PhysicsType::eos();
    if (data_wk->ireactor_type == eint_rho) {
      // UV REACTOR
      eos.EY2T(energy, massfrac, temp);
      eos.TY2Cv(temp, massfrac, cX);
      eos.T2Ei(temp, Xi);
    } else if (data_wk->ireactor_type == enth_rho) {
      // HP REACTOR
      eos.HY2T(energy, massfrac, temp);
      eos.TY2Cp(temp, massfrac, cX);
      eos.T2Hi(temp, Xi);
    }
    eos.RTY2WDOT(rho, temp, massfrac, cdot);

    // Fill ydot vect
    ydot_d[offset + NUM_SPECIES] = data_wk->rhoXsrc_ext[data_wk->boxcell + tid];
    for (int i = 0; i < NUM_SPECIES; i++) {
      ydot_d[offset + i] =
        cdot[i] + data_wk->rYsrc[(data_wk->boxcell + tid) * (NUM_SPECIES) + i];
      ydot_d[offset + NUM_SPECIES] =
        ydot_d[offset + NUM_SPECIES] - ydot_d[offset + i] * Xi[i];
    }
    ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] / (rho * cX);
  }
}

// Analytical Jacobian evaluation
int
cJac(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{

  // Make local copies of pointers to input data (big M)
  realtype* ydata = N_VGetArrayPointer(u);

  // Make local copies of pointers in user_data (cell M)
  UserData data_wk = (UserData)user_data;

  BL_PROFILE_VAR("DenseJac", DenseJac);
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    // Tmp vars
    realtype* J_col_k;
    realtype massfrac[NUM_SPECIES], molecular_weight[NUM_SPECIES];
    realtype temp;
    realtype Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];

    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);

    // MW CGS
    CKWT(molecular_weight);

    // rho MKS
    realtype rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }

    temp = ydata[offset + NUM_SPECIES];

    // Yks
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] / rho;
    }

    // Jac
    int consP;
    if (data_wk->ireactor_type == eint_rho) {
      consP = 0;
    } else {
      consP = 1;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
    // fill the sunMat
    for (int k = 0; k < NUM_SPECIES; k++) {
      J_col_k = SM_COLUMN_D(J, offset + k);
      for (int i = 0; i < NUM_SPECIES; i++) {
        J_col_k[offset + i] = Jmat_tmp[k * (NUM_SPECIES + 1) + i] *
                              molecular_weight[i] / molecular_weight[k];
      }
      J_col_k[offset + NUM_SPECIES] =
        Jmat_tmp[k * (NUM_SPECIES + 1) + NUM_SPECIES] / molecular_weight[k];
    }
    J_col_k = SM_COLUMN_D(J, offset + NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++) {
      J_col_k[offset + i] =
        Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] * molecular_weight[i];
    }
    J_col_k = SM_COLUMN_D(J, offset);
  }
  BL_PROFILE_VAR_STOP(DenseJac);

  return (0);
}

// Analytical SPARSE Jacobian evaluation
int
cJac_sps(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{
  // Make local copies of pointers to input data (big M)
  realtype* ydata = N_VGetArrayPointer(u);
  sunindextype* rowPtrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype* colIndx_tmp = SUNSparseMatrix_IndexValues(J);
  realtype* Jdata = SUNSparseMatrix_Data(J);

  // Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData)user_data;

  // MW CGS
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  BL_PROFILE_VAR("FillSparseData", FillSpsData);
  // Fixed colVal
  for (int i = 0; i < data_wk->NNZ * data_wk->ncells; i++) {
    colIndx_tmp[i] = (sunindextype)data_wk->colVals_c[i];
  }
  rowPtrs_tmp[0] = (sunindextype)data_wk->rowPtrs_c[0];
  // Fixed rowPtrs
  for (int i = 0; i < data_wk->ncells * (NUM_SPECIES + 1); i++) {
    rowPtrs_tmp[i + 1] = (sunindextype)data_wk->rowPtrs_c[i + 1];
  }
  BL_PROFILE_VAR_STOP(FillSpsData);

  BL_PROFILE_VAR("SparseJac", SpsJac);
  // Temp vectors
  realtype temp_save_lcl, temp;
  realtype massfrac[NUM_SPECIES];
  realtype Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
  // Save Jac from cell to cell if more than one
  temp_save_lcl = 0.0;
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);
    int offset_J = tid * data_wk->NNZ;
    // rho MKS
    realtype rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }
    // Yks
    realtype rhoinv = 1.0 / rho;
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] * rhoinv;
    }
    temp = ydata[offset + NUM_SPECIES];
    // Do we recompute Jac ?
    if (fabs(temp - temp_save_lcl) > 1.0) {
      int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
      temp_save_lcl = temp;
      // rescale
      for (int i = 0; i < NUM_SPECIES; i++) {
        for (int k = 0; k < NUM_SPECIES; k++) {
          Jmat_tmp[k * (NUM_SPECIES + 1) + i] =
            Jmat_tmp[k * (NUM_SPECIES + 1) + i] * molecular_weight[i] /
            molecular_weight[k];
        }
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] =
          Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] / molecular_weight[i];
      }
      for (int i = 0; i < NUM_SPECIES; i++) {
        Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] =
          Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] * molecular_weight[i];
      }
    }
    // Go from Dense to Sparse
    for (int i = 1; i < NUM_SPECIES + 2; i++) {
      int nbVals = data_wk->rowPtrs_c[i] - data_wk->rowPtrs_c[i - 1];
      for (int j = 0; j < nbVals; j++) {
        int idx = data_wk->colVals_c[data_wk->rowPtrs_c[i - 1] + j];
        Jdata[offset_J + data_wk->rowPtrs_c[i - 1] + j] =
          Jmat_tmp[(i - 1) + (NUM_SPECIES + 1) * idx];
      }
    }
  }
  BL_PROFILE_VAR_STOP(SpsJac);

  return (0);
}

#ifdef USE_KLU_PP
// Analytical SPARSE Jacobian evaluation
int
cJac_KLU(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  SUNMatrix J,
  void* user_data,
  N_Vector /* tmp1 */,
  N_Vector /* tmp2 */,
  N_Vector /* tmp3 */)
{

  BL_PROFILE_VAR("SparseKLUJac", SpsKLUJac);
  // Make local copies of pointers to input data (big M)
  realtype* ydata = N_VGetArrayPointer(u);
  sunindextype* colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype* rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  realtype* Jdata = SUNSparseMatrix_Data(J);

  // Make local copies of pointers in user_data (cell M)
  UserData data_wk = (UserData)user_data;

  // MW CGS
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  // Fixed RowVals
  for (int i = 0; i < data_wk->NNZ; i++) {
    rowvals_tmp[i] = data_wk->rowVals[0][i];
  }
  // Fixed colPtrs
  colptrs_tmp[0] = data_wk->colPtrs[0][0];
  for (int i = 0; i < data_wk->ncells * (NUM_SPECIES + 1); i++) {
    colptrs_tmp[i + 1] = data_wk->colPtrs[0][i + 1];
  }

  // Temp vectors
  realtype temp_save_lcl, temp;
  realtype massfrac[NUM_SPECIES];
  realtype Jmat_tmp[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
  // Save Jac from cell to cell if more than one
  temp_save_lcl = 0.0;
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);
    // rho MKS
    realtype rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + ydata[offset + i];
    }
    // Yks
    realtype rhoinv = 1.0 / rho;
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = ydata[offset + i] * rhoinv;
    }
    temp = ydata[offset + NUM_SPECIES];
    // Do we recompute Jac ?
    if (fabs(temp - temp_save_lcl) > 1.0) {
      // NRG CGS
      int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
      temp_save_lcl = temp;
      // rescale
      for (int i = 0; i < NUM_SPECIES; i++) {
        for (int k = 0; k < NUM_SPECIES; k++) {
          Jmat_tmp[k * (NUM_SPECIES + 1) + i] =
            Jmat_tmp[k * (NUM_SPECIES + 1) + i] * molecular_weight[i] /
            molecular_weight[k];
        }
        Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] =
          Jmat_tmp[i * (NUM_SPECIES + 1) + NUM_SPECIES] / molecular_weight[i];
      }
      for (int i = 0; i < NUM_SPECIES; i++) {
        Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] =
          Jmat_tmp[NUM_SPECIES * (NUM_SPECIES + 1) + i] * molecular_weight[i];
      }
    }
    // Go from Dense to Sparse
    BL_PROFILE_VAR("DensetoSps", DtoS);
    for (int i = 1; i < NUM_SPECIES + 2; i++) {
      int nbVals = data_wk->colPtrs[0][i] - data_wk->colPtrs[0][i - 1];
      for (int j = 0; j < nbVals; j++) {
        int idx = data_wk->rowVals[0][data_wk->colPtrs[0][i - 1] + j];
        Jdata[data_wk->colPtrs[0][offset + i - 1] + j] =
          Jmat_tmp[(i - 1) * (NUM_SPECIES + 1) + idx];
      }
    }
    BL_PROFILE_VAR_STOP(DtoS);
  }

  BL_PROFILE_VAR_STOP(SpsKLUJac);

  return (0);
}
#endif

// Preconditioner setup routine for GMRES solver when custom sparse mode is
// activated Generate and preprocess P
int
Precond_custom(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data)
{
  // Make local copies of pointers to input data (big M)
  realtype* udata = N_VGetArrayPointer(u);
  // Make local copies of pointers in user_data (cell M)
  UserData data_wk = (UserData)user_data;

  // MW CGS
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  // Check if Jac is stale
  if (jok) {
    // jok = SUNTRUE: Copy Jbd to P
    *jcurPtr = SUNFALSE;
  } else {
    // Temp vectors
    realtype temp, temp_save_lcl;
    realtype activity[NUM_SPECIES], massfrac[NUM_SPECIES];
    // Save Jac from cell to cell if more than one
    temp_save_lcl = 0.0;
    for (int tid = 0; tid < data_wk->ncells; tid++) {
      // Offset in case several cells
      int offset = tid * (NUM_SPECIES + 1);
      // rho MKS
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++) {
        rho = rho + udata[offset + i];
      }
      // Yks
      for (int i = 0; i < NUM_SPECIES; i++) {
        massfrac[i] = udata[offset + i] / rho;
      }
      temp = udata[offset + NUM_SPECIES];
      // Activities
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2C(rho, temp, massfrac, activity);
      // Do we recompute Jac ?
      if (fabs(temp - temp_save_lcl) > 1.0) {
        // Formalism
        int consP;
        if (data_wk->ireactor_type == eint_rho) {
          consP = 0;
        } else {
          consP = 1;
        }
        DWDOT_SIMPLIFIED(data_wk->JSPSmat[tid], activity, &temp, &consP);

        for (int i = 0; i < NUM_SPECIES; i++) {
          for (int k = 0; k < NUM_SPECIES; k++) {
            (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] =
              (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] *
              molecular_weight[i] / molecular_weight[k];
          }
          (data_wk->JSPSmat[tid])[i * (NUM_SPECIES + 1) + NUM_SPECIES] =
            (data_wk->JSPSmat[tid])[i * (NUM_SPECIES + 1) + NUM_SPECIES] /
            molecular_weight[i];
        }
        for (int i = 0; i < NUM_SPECIES; i++) {
          (data_wk->JSPSmat[tid])[NUM_SPECIES * (NUM_SPECIES + 1) + i] =
            (data_wk->JSPSmat[tid])[NUM_SPECIES * (NUM_SPECIES + 1) + i] *
            molecular_weight[i];
        }
        temp_save_lcl = temp;
      } else {
        // if not: copy the one from prev cell
        for (int i = 0; i < NUM_SPECIES + 1; i++) {
          for (int k = 0; k < NUM_SPECIES + 1; k++) {
            (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] =
              (data_wk->JSPSmat[tid - 1])[k * (NUM_SPECIES + 1) + i];
          }
        }
      }
    }

    *jcurPtr = SUNTRUE;
  }

  int nbVals;
  for (int i = 1; i < NUM_SPECIES + 2; i++) {
    // nb non zeros elem should be the same for all cells
    nbVals = data_wk->rowPtrs[0][i] - data_wk->rowPtrs[0][i - 1];
    for (int j = 0; j < nbVals; j++) {
      // row of non zero elem should be the same for all cells
      int idx = data_wk->colVals[0][data_wk->rowPtrs[0][i - 1] + j];
      // Scale by -gamma
      // Add identity matrix
      for (int tid = 0; tid < data_wk->ncells; tid++) {
        if (idx == (i - 1)) {
          data_wk->Jdata[tid][data_wk->rowPtrs[tid][i - 1] + j] =
            1.0 -
            gamma * (data_wk->JSPSmat[tid])[idx * (NUM_SPECIES + 1) + idx];
        } else {
          data_wk->Jdata[tid][data_wk->rowPtrs[tid][i - 1] + j] =
            -gamma * (data_wk->JSPSmat[tid])[(i - 1) + (NUM_SPECIES + 1) * idx];
        }
      }
    }
  }

  return (0);
}

#ifdef USE_KLU_PP
// Preconditioner setup routine for GMRES solver when KLU sparse mode is
// activated Generate and preprocess P
int
Precond_sparse(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data)
{
  // Make local copies of pointers to input data (big M)
  realtype* udata = N_VGetArrayPointer(u);
  // Make local copies of pointers in user_data (cell M)
  UserData data_wk = (UserData)user_data;

  // MW CGS
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  // Check if Jac is stale
  if (jok) {
    // jok = SUNTRUE: Copy Jbd to P
    *jcurPtr = SUNFALSE;
  } else {
    // Temp vectors
    realtype activity[NUM_SPECIES], massfrac[NUM_SPECIES];
    // Save Jac from cell to cell if more than one
    realtype temp_save_lcl = 0.0;
    for (int tid = 0; tid < data_wk->ncells; tid++) {
      // Offset in case several cells
      int offset = tid * (NUM_SPECIES + 1);
      // rho MKS
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++) {
        rho = rho + udata[offset + i];
      }
      // Yks
      for (int i = 0; i < NUM_SPECIES; i++) {
        massfrac[i] = udata[offset + i] / rho;
      }
      realtype temp = udata[offset + NUM_SPECIES];
      // Activities
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2C(rho, temp, massfrac, activity);
      // Do we recompute Jac ?
      if (fabs(temp - temp_save_lcl) > 1.0) {
        // Formalism
        int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
        DWDOT_SIMPLIFIED(data_wk->JSPSmat[tid], activity, &temp, &consP);

        for (int i = 0; i < NUM_SPECIES; i++) {
          for (int k = 0; k < NUM_SPECIES; k++) {
            (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] =
              (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] *
              molecular_weight[i] / molecular_weight[k];
          }
          (data_wk->JSPSmat[tid])[i * (NUM_SPECIES + 1) + NUM_SPECIES] =
            (data_wk->JSPSmat[tid])[i * (NUM_SPECIES + 1) + NUM_SPECIES] /
            molecular_weight[i];
        }
        for (int i = 0; i < NUM_SPECIES; i++) {
          (data_wk->JSPSmat[tid])[NUM_SPECIES * (NUM_SPECIES + 1) + i] =
            (data_wk->JSPSmat[tid])[NUM_SPECIES * (NUM_SPECIES + 1) + i] *
            molecular_weight[i];
        }
        temp_save_lcl = temp;
      } else {
        // if not: copy the one from prev cell
        for (int i = 0; i < NUM_SPECIES + 1; i++) {
          for (int k = 0; k < NUM_SPECIES + 1; k++) {
            (data_wk->JSPSmat[tid])[k * (NUM_SPECIES + 1) + i] =
              (data_wk->JSPSmat[tid - 1])[k * (NUM_SPECIES + 1) + i];
          }
        }
      }
    }

    *jcurPtr = SUNTRUE;
  }

  for (int i = 1; i < NUM_SPECIES + 2; i++) {
    // nb non zeros elem should be the same for all cells
    int nbVals = data_wk->colPtrs[0][i] - data_wk->colPtrs[0][i - 1];
    for (int j = 0; j < nbVals; j++) {
      // row of non zero elem should be the same for all cells
      int idx = data_wk->rowVals[0][data_wk->colPtrs[0][i - 1] + j];
      // Scale by -gamma
      // Add identity matrix
      for (int tid = 0; tid < data_wk->ncells; tid++) {
        if (idx == (i - 1)) {
          data_wk->Jdata[tid][data_wk->colPtrs[tid][i - 1] + j] =
            1.0 -
            gamma * (data_wk->JSPSmat[tid])[idx * (NUM_SPECIES + 1) + idx];
        } else {
          data_wk->Jdata[tid][data_wk->colPtrs[tid][i - 1] + j] =
            -gamma * (data_wk->JSPSmat[tid])[(i - 1) * (NUM_SPECIES + 1) + idx];
        }
      }
    }
  }

  BL_PROFILE_VAR("KLU_factorization", KLU_factor);
  if (!(data_wk->FirstTimePrecond)) {
    for (int tid = 0; tid < data_wk->ncells; tid++) {
      klu_refactor(
        data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid],
        data_wk->Symbolic[tid], data_wk->Numeric[tid], &(data_wk->Common[tid]));
    }
  } else {
    for (int tid = 0; tid < data_wk->ncells; tid++) {
      data_wk->Numeric[tid] = klu_factor(
        data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid],
        data_wk->Symbolic[tid], &(data_wk->Common[tid]));
    }
    data_wk->FirstTimePrecond = false;
  }
  BL_PROFILE_VAR_STOP(KLU_factor);

  return (0);
}

#else
// Preconditioner setup routine for GMRES solver when no sparse mode is
// activated Generate and preprocess P
int
Precond(
  realtype /* tn */,
  N_Vector u,
  N_Vector /* fu */,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data)
{
  // Make local copies of pointers to input data (big M)
  realtype* udata = N_VGetArrayPointer(u);

  // Make local copies of pointers in user_data
  UserData data_wk = (UserData)user_data;
  realtype**(**P), **(**Jbd);
  sunindextype*(**pivot);
  P = (data_wk->P);
  Jbd = (data_wk->Jbd);
  pivot = (data_wk->pivot);

  // Tmp arrays
  realtype molecular_weight[NUM_SPECIES];
  realtype Jmat[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
  realtype massfrac[NUM_SPECIES], activity[NUM_SPECIES];
  realtype temp;
  sunindextype ierr;

  // MW CGS
  CKWT(molecular_weight);

  if (jok) {
    // jok = SUNTRUE: Copy Jbd to P
    denseCopy(Jbd[0][0], P[0][0], NUM_SPECIES + 1, NUM_SPECIES + 1);
    *jcurPtr = SUNFALSE;
  } else {
    // rho MKS
    realtype rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + udata[i];
    }
    // Yks
    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = udata[i] / rho;
    }
    temp = udata[NUM_SPECIES];
    // Activities
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2C(rho, temp, massfrac, activity);
    // jok = SUNFALSE: Generate Jbd from scratch and copy to P
    // Make local copies of problem variables, for efficiency.
    int consP;
    if (data_wk->ireactor_type == eint_rho) {
      consP = 0;
    } else {
      consP = 1;
    }
    DWDOT_SIMPLIFIED(Jmat, activity, &temp, &consP);

    // Compute Jacobian.  Load into P.
    denseScale(0.0, Jbd[0][0], NUM_SPECIES + 1, NUM_SPECIES + 1);
    for (int i = 0; i < NUM_SPECIES; i++) {
      for (int k = 0; k < NUM_SPECIES; k++) {
        (Jbd[0][0])[k][i] = Jmat[k * (NUM_SPECIES + 1) + i] *
                            molecular_weight[i] / molecular_weight[k];
      }
      (Jbd[0][0])[i][NUM_SPECIES] =
        Jmat[i * (NUM_SPECIES + 1) + NUM_SPECIES] / molecular_weight[i];
    }
    for (int i = 0; i < NUM_SPECIES; i++) {
      (Jbd[0][0])[NUM_SPECIES][i] =
        Jmat[NUM_SPECIES * (NUM_SPECIES + 1) + i] * molecular_weight[i];
    }
    (Jbd[0][0])[NUM_SPECIES][NUM_SPECIES] =
      Jmat[(NUM_SPECIES + 1) * (NUM_SPECIES + 1) - 1];

    denseCopy(Jbd[0][0], P[0][0], NUM_SPECIES + 1, NUM_SPECIES + 1);

    *jcurPtr = SUNTRUE;
  }

  // Scale by -gamma
  denseScale(-gamma, P[0][0], NUM_SPECIES + 1, NUM_SPECIES + 1);
  // denseScale(0.0, P[0][0], NUM_SPECIES+1, NUM_SPECIES+1);

  // Add identity matrix and do LU decompositions on blocks in place.
  denseAddIdentity(P[0][0], NUM_SPECIES + 1);
  ierr = denseGETRF(P[0][0], NUM_SPECIES + 1, NUM_SPECIES + 1, pivot[0][0]);
  if (ierr != 0)
    return (1);

  return (0);
}
#endif

// PSolve for GMRES solver when custom sparse mode is activated
int
PSolve_custom(
  realtype /* tn */,
  N_Vector /* u */,
  N_Vector /* fu */,
  N_Vector r,
  N_Vector z,
  realtype /* gamma */,
  realtype /* delta */,
  int /* lr */,
  void* user_data)
{
  // Make local copies of pointers in user_data
  UserData data_wk = (UserData)user_data;

  // Make local copies of pointers to input data (big M)
  realtype* zdata = N_VGetArrayPointer(z);
  realtype* rdata = N_VGetArrayPointer(r);

  N_VScale(1.0, r, z);

  // Solve the block-diagonal system Pz = r using LU factors stored
  // in P and pivot data in pivot, and return the solution in z.
  BL_PROFILE_VAR("GaussSolver", GaussSolver);
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    int offset = tid * (NUM_SPECIES + 1);
    double* z_d_offset = zdata + offset;
    double* r_d_offset = rdata + offset;
    sgjsolve_simplified(data_wk->Jdata[tid], z_d_offset, r_d_offset);
  }
  BL_PROFILE_VAR_STOP(GaussSolver);

  return (0);
}

#ifdef USE_KLU_PP
// PSolve for GMRES solver when KLU sparse mode is activated
int
PSolve_sparse(
  realtype /* tn */,
  N_Vector /* u */,
  N_Vector /* fu */,
  N_Vector r,
  N_Vector z,
  realtype /* gamma */,
  realtype /* delta */,
  int /* lr */,
  void* user_data)
{
  // Make local copies of pointers in user_data
  UserData data_wk = (UserData)user_data;

  // Make local copies of pointers to input data (big M)
  realtype* zdata = N_VGetArrayPointer(z);

  BL_PROFILE_VAR("KLU_inversion", PSolve_sparse);
  N_VScale(1.0, r, z);

  // Solve the block-diagonal system Pz = r using LU factors stored
  //   in P and pivot data in pivot, and return the solution in z.
  realtype zdata_cell[NUM_SPECIES + 1];
  for (int tid = 0; tid < data_wk->ncells; tid++) {
    int offset_beg = tid * (NUM_SPECIES + 1);
    std::memcpy(
      zdata_cell, zdata + offset_beg, (NUM_SPECIES + 1) * sizeof(realtype));
    klu_solve(
      data_wk->Symbolic[tid], data_wk->Numeric[tid], NUM_SPECIES + 1, 1,
      zdata_cell, &(data_wk->Common[tid]));
    std::memcpy(
      zdata + offset_beg, zdata_cell, (NUM_SPECIES + 1) * sizeof(realtype));
  }
  BL_PROFILE_VAR_STOP(PSolve_sparse);

  return (0);
}

#else
// PSolve for GMRES solver when no sparse mode is activated
int
PSolve(
  realtype /* tn */,
  N_Vector /* u */,
  N_Vector /* fu */,
  N_Vector r,
  N_Vector z,
  realtype /* gamma */,
  realtype /* delta */,
  int /* lr */,
  void* user_data)
{
  // Make local copies of pointers to input data (big M)
  realtype* zdata = N_VGetArrayPointer(z);

  // Extract the P and pivot arrays from user_data.
  UserData data_wk = (UserData)user_data;
  realtype**(**P);
  sunindextype*(**pivot);
  P = data_wk->P;
  pivot = data_wk->pivot;

  N_VScale(1.0, r, z);

  // Solve the block-diagonal system Pz = r using LU factors stored
  //   in P and pivot data in pivot, and return the solution in z.
  realtype* v = zdata;
  denseGETRS(P[0][0], NUM_SPECIES + 1, pivot[0][0], v);

  return (0);
}
#endif

SUNLinearSolver
SUNLinSol_sparse_custom(
  N_Vector a_y,
  SUNMatrix a_A,
  int reactor_type,
  int nsubsys,
  int subsys_size,
  int subsys_nnz)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Sparse_custom content;

  // Check that required arguments are not NULL
  if (a_y == NULL || a_A == NULL)
    return (NULL);
  if (SUNMatGetID(a_A) != SUNMATRIX_SPARSE)
    return (NULL);

  // Matrix should be square
  if (SUNSparseMatrix_Columns(a_A) != SUNSparseMatrix_Rows(a_A))
    return (NULL);

  // Check that it is a CSR matrix
  if (SUNSparseMatrix_SparseType(a_A) != CSR_MAT)
    return (NULL);

  // Matrix and vector dimensions must agree
  if (N_VGetLength(a_y) != SUNSparseMatrix_Columns(a_A))
    return (NULL);

  // All subsystems must be the same size
  if (SUNSparseMatrix_Columns(a_A) != (subsys_size * nsubsys))
    return (NULL);

  // Number of nonzeros per subsys must be the same
  if (SUNSparseMatrix_NNZ(a_A) != (subsys_nnz * nsubsys))
    return (NULL);

  // Create an empty linear solver
  S = SUNLinSolNewEmpty();
  if (S == NULL)
    return (NULL);

  // Attach operations
  S->ops->gettype = SUNLinSolGetType_Sparse_custom;
  S->ops->solve = SUNLinSolSolve_Sparse_custom;

  // Create content
  content = (SUNLinearSolverContent_Sparse_custom)malloc(sizeof *content);
  if (content == NULL) {
    SUNLinSolFree(S);
    return (NULL);
  }

  // Attach content
  S->content = content;

  // Fill content
  content->last_flag = 0;
  content->reactor_type = reactor_type;
  content->nsubsys = nsubsys;
  content->subsys_size = subsys_size;
  content->subsys_nnz = subsys_nnz;

  return (S);
}

SUNLinearSolver_Type SUNLinSolGetType_Sparse_custom(SUNLinearSolver /* S */)
{
  return (SUNLINEARSOLVER_DIRECT);
}

int
SUNLinSolSolve_Sparse_custom(
  SUNLinearSolver S, SUNMatrix a_A, N_Vector x, N_Vector b, realtype /* tol */)
{
  realtype* x_d = N_VGetArrayPointer(x);
  realtype* b_d = N_VGetArrayPointer(b);

  double* Data = (double*)SUNSparseMatrix_Data(a_A);

  BL_PROFILE_VAR("GaussSolver", GaussSolver);
  for (int tid = 0; tid < SUN_CUSP_NUM_SUBSYS(S); tid++) {
    int offset = tid * SUN_CUSP_SUBSYS_NNZ(S);
    int offset_RHS = tid * SUN_CUSP_SUBSYS_SIZE(S);
    double* Data_offset = Data + offset;
    double* x_d_offset = x_d + offset_RHS;
    double* b_d_offset = b_d + offset_RHS;
    sgjsolve(Data_offset, x_d_offset, b_d_offset);
  }
  BL_PROFILE_VAR_STOP(GaussSolver);

  return (SUNLS_SUCCESS);
}

void
check_state(N_Vector yvec)
{
  realtype* ydata = N_VGetArrayPointer(yvec);

  data->actual_ok_to_react = true;

  for (int tid = 0; tid < data->ncells; tid++) {
    // Offset in case several cells
    int offset = tid * (NUM_SPECIES + 1);
    // rho MKS
    realtype rho = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      rho = rho + ydata[offset + k];
    }
    realtype Temp = ydata[offset + NUM_SPECIES];
    if ((rho < 1.0e-10) || (rho > 1.e10)) {
      data->actual_ok_to_react = false;
      Print() << "rho " << rho << "\n";
    }
    if ((Temp < 200.0) || (Temp > 5000.0)) {
      data->actual_ok_to_react = false;
      Print() << "Temp " << Temp << "\n";
    }
  }
}

void
PrintFinalStats(void* cvodeMem, realtype Temp)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn;
  long int nli, npe, nps, ncfl, netfails;
  int flag;
  realtype hlast, hinused, hcur;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netfails);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetLastStep(cvodeMem, &hlast);
  check_flag(&flag, "CVodeGetLastStep", 1);
  flag = CVodeGetActualInitStep(cvodeMem, &hinused);
  check_flag(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetCurrentTime(cvodeMem, &hcur);
  check_flag(&flag, "CVodeGetCurrentTime", 1);

  if (data->isolve_type == dense_solve) {
    flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
    check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
    flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
    check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (data->isolve_type == iterative_gmres_solve) {
    flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
    check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
    flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
    // flag = CVSpilsGetNumJTSetupEvals(cvodeMem, &nje);
    check_flag(&flag, "CVSpilsGetNumJtimesEvals", 1);
    flag = CVSpilsGetNumPrecEvals(cvodeMem, &npe);
    check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
    flag = CVSpilsGetNumPrecSolves(cvodeMem, &nps);
    check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
    flag = CVSpilsGetNumLinIters(cvodeMem, &nli);
    check_flag(&flag, "CVSpilsGetNumLinIters", 1);
    flag = CVSpilsGetNumConvFails(cvodeMem, &ncfl);
    check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  }

  Print() << "-- Final Statistics --\n";
  Print() << "NonLinear (Newton) related --\n";
  Print() << Temp << " |DT(dt, dtcur) = " << nst << "(" << hlast << "," << hcur
          << "), RHS = " << nfe << ", Iterations = " << nni
          << ", ErrTestFails = " << netfails << ", LinSolvSetups = " << nsetups
          << "\n";
  if (data->isolve_type == dense_solve) {
    Print() << "Linear (Dense Direct Solve) related --\n";
    Print() << Temp << " |FD RHS = " << nfeLS << ", NumJacEvals = " << nje
            << " \n";
  } else if (data->isolve_type == iterative_gmres_solve) {
    // LinSolvSetups actually reflects the number of time the LinSolver has been
    // called. NonLinIterations can be taken without the need for LinItes
    Print() << "Linear (Krylov GMRES Solve) related --\n";
    Print() << Temp << " |RHSeval = " << nfeLS << ", jtvEval = " << nje
            << ", NumPrecEvals = " << npe << ", NumPrecSolves = " << nps
            << "\n";
    Print() << Temp << " |Iterations = " << nli << ", ConvFails = " << ncfl
            << "\n";
  }
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

int
check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    if (ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    }
    return (1);
  }

  // Check if flag < 0
  else if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      if (ParallelDescriptor::IOProcessor()) {
        fprintf(
          stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
          *errflag);
      }
      return (1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    if (ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    }
    return (1);
  }

  return (0);
}

UserData
AllocUserData(int reactor_type, int num_cells)
{
  // Make local copies of pointers in user_data
  UserData data_wk = (UserData)malloc(sizeof *data_wk);
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  // ParmParse from the inputs file: only done once
  ParmParse pp("ode");
  pp.query("analytical_jacobian", data_wk->ianalytical_jacobian);
  data_wk->iverbose = 1;
  pp.query("verbose", data_wk->iverbose);

  std::string solve_type_str = "none";
  ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  /* options are:
  dense_solve           = 1;
  sparse_solve          = 5;
  iterative_gmres_solve = 99;
  sparse_solve_custom   = 101;
  iterative_gmres_solve_custom = 199;
  hack_dump_sparsity_pattern = -5;
  */
  if (solve_type_str == "dense") {
    data_wk->isolve_type = dense_solve;
  } else if (solve_type_str == "sparse") {
    data_wk->isolve_type = sparse_solve;
  } else if (solve_type_str == "GMRES") {
    data_wk->isolve_type = iterative_gmres_solve;
  } else if (solve_type_str == "sparse_custom") {
    data_wk->isolve_type = sparse_solve_custom;
  } else if (solve_type_str == "GMRES_custom") {
    data_wk->isolve_type = iterative_gmres_solve_custom;
  } else if (solve_type_str == "diag") {
    data_wk->isolve_type = hack_dump_sparsity_pattern;
  } else {
    Abort("Wrong solve_type. Options are: dense, sparse, GMRES, sparse_custom, "
          "GMRES_custom");
  }

  (data_wk->ireactor_type) = reactor_type;

  (data_wk->ncells) = num_cells;

  // Not sure of the interest of doing that the following:
  // N_Vector Data = NULL;
  // Data = N_VNew_Serial(data_wk->ncells*(NUM_SPECIES+1));
  //(data_wk->Yvect_full)  = N_VGetArrayPointer_Serial(Data);
  // Data = N_VNew_Serial(data_wk->ncells*(NUM_SPECIES));
  //(data_wk->rYsrc)       = N_VGetArrayPointer_Serial(Data);
  // Data = N_VNew_Serial(data_wk->ncells);
  //(data_wk->rhoX_init)   = N_VGetArrayPointer_Serial(Data);
  // Data = N_VNew_Serial(data_wk->ncells);
  //(data_wk->rhoXsrc_ext) = N_VGetArrayPointer_Serial(Data);

  (data_wk->Yvect_full) = new amrex::Real[data_wk->ncells * (NUM_SPECIES + 1)];
  (data_wk->rYsrc) = new amrex::Real[data_wk->ncells * (NUM_SPECIES)];
  (data_wk->rhoX_init) = new amrex::Real[data_wk->ncells];
  (data_wk->rhoXsrc_ext) = new amrex::Real[data_wk->ncells];
  (data_wk->FCunt) = new int[data_wk->ncells];
  (data_wk->mask) = new int[data_wk->ncells];

  (data_wk->FirstTimePrecond) = true;
  (data_wk->reactor_cvode_initialized) = false;
  (data_wk->actual_ok_to_react) = true;

  int HP;
  if (data_wk->ireactor_type == eint_rho) {
    HP = 0;
  } else {
    HP = 1;
  }

  // Sparse Direct and Sparse (It) Precond data
  data_wk->colPtrs = new int*[data_wk->ncells];
  data_wk->rowVals = new int*[data_wk->ncells];
  data_wk->Jdata = new realtype*[data_wk->ncells];

#ifndef USE_KLU_PP
  if (data_wk->isolve_type == iterative_gmres_solve) {
    // Precond data
    (data_wk->P) = new realtype***[data_wk->ncells];
    (data_wk->Jbd) = new realtype***[data_wk->ncells];
    (data_wk->pivot) = new sunindextype**[data_wk->ncells];
    for (int i = 0; i < data_wk->ncells; ++i) {
      (data_wk->P)[i] = new realtype**[data_wk->ncells];
      (data_wk->Jbd)[i] = new realtype**[data_wk->ncells];
      (data_wk->pivot)[i] = new sunindextype*[data_wk->ncells];
    }

    for (int i = 0; i < data_wk->ncells; ++i) {
      (data_wk->P)[i][i] = newDenseMat(NUM_SPECIES + 1, NUM_SPECIES + 1);
      (data_wk->Jbd)[i][i] = newDenseMat(NUM_SPECIES + 1, NUM_SPECIES + 1);
      (data_wk->pivot)[i][i] = newIndexArray(NUM_SPECIES + 1);
    }
    //}

#else
  if (data_wk->isolve_type == sparse_solve) {
    // Sparse Matrix for Direct Sparse KLU solver
    (data_wk->PS) = new SUNMatrix[1];
    SPARSITY_INFO(&(data_wk->NNZ), &HP, data_wk->ncells);
    if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
      Print() << "--> SPARSE solver -- non zero entries: " << data_wk->NNZ
              << ", which represents "
              << data_wk->NNZ /
                   float(
                     (NUM_SPECIES + 1) * (NUM_SPECIES + 1) * (data_wk->ncells) *
                     (data_wk->ncells)) *
                   100.0
              << " % fill-in pattern\n";
    }
    (data_wk->PS)[0] = SUNSparseMatrix(
      (NUM_SPECIES + 1) * data_wk->ncells, (NUM_SPECIES + 1) * data_wk->ncells,
      data_wk->NNZ, CSC_MAT);
    data_wk->colPtrs[0] = (int*)SUNSparseMatrix_IndexPointers((data_wk->PS)[0]);
    data_wk->rowVals[0] = (int*)SUNSparseMatrix_IndexValues((data_wk->PS)[0]);
    data_wk->Jdata[0] = SUNSparseMatrix_Data((data_wk->PS)[0]);
    SPARSITY_PREPROC_CSC(
      data_wk->rowVals[0], data_wk->colPtrs[0], &HP, data_wk->ncells);

  } else if (data_wk->isolve_type == iterative_gmres_solve) {
    // KLU internal storage
    data_wk->Common = new klu_common[data_wk->ncells];
    data_wk->Symbolic = new klu_symbolic*[data_wk->ncells];
    data_wk->Numeric = new klu_numeric*[data_wk->ncells];
    // Sparse Matrices for It Sparse KLU block-solve
    data_wk->PS = new SUNMatrix[data_wk->ncells];
    // Number of non zero elements
    SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ), &HP);
    if (
      (data_wk->iverbose > 0) && (omp_thread == 0) &&
      (data_wk->ianalytical_jacobian != 0)) {
      Print() << "--> SPARSE Preconditioner -- non zero entries: "
              << data_wk->NNZ << ", which represents "
              << data_wk->NNZ / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                   100.0
              << " % fill-in pattern\n";
    }
    // Not used yet. TODO use to fetch sparse Mat
    data_wk->indx = new int[data_wk->NNZ];
    data_wk->JSPSmat = new realtype*[data_wk->ncells];
    for (int i = 0; i < data_wk->ncells; ++i) {
      (data_wk->PS)[i] = SUNSparseMatrix(
        NUM_SPECIES + 1, NUM_SPECIES + 1, data_wk->NNZ, CSC_MAT);
      data_wk->colPtrs[i] =
        (int*)SUNSparseMatrix_IndexPointers((data_wk->PS)[i]);
      data_wk->rowVals[i] = (int*)SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
      data_wk->Jdata[i] = SUNSparseMatrix_Data((data_wk->PS)[i]);
      // indx not used YET
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
        data_wk->rowVals[i], data_wk->colPtrs[i], data_wk->indx, &HP);
      data_wk->JSPSmat[i] = new realtype[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
      klu_defaults(&(data_wk->Common[i]));
      // data_wk->Common.btf = 0;
      //(data_wk->Common[i]).maxwork = 15;
      // data_wk->Common.ordering = 1;
      data_wk->Symbolic[i] = klu_analyze(
        NUM_SPECIES + 1, data_wk->colPtrs[i], data_wk->rowVals[i],
        &(data_wk->Common[i]));
    }
    //}
#endif

  } else if (data_wk->isolve_type == iterative_gmres_solve_custom) {
    // Sparse Direct and Sparse (It) Precond data
    data_wk->colVals = new int*[data_wk->ncells];
    data_wk->rowPtrs = new int*[data_wk->ncells];
    // Matrices for It Sparse custom block-solve
    data_wk->PS = new SUNMatrix[data_wk->ncells];
    data_wk->JSPSmat = new realtype*[data_wk->ncells];
    // Number of non zero elements
    SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ), &HP);
    for (int i = 0; i < data_wk->ncells; ++i) {
      (data_wk->PS)[i] = SUNSparseMatrix(
        NUM_SPECIES + 1, NUM_SPECIES + 1, data_wk->NNZ, CSR_MAT);
      data_wk->rowPtrs[i] =
        (int*)SUNSparseMatrix_IndexPointers((data_wk->PS)[i]);
      data_wk->colVals[i] = (int*)SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
      data_wk->Jdata[i] = SUNSparseMatrix_Data((data_wk->PS)[i]);
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
        data_wk->colVals[i], data_wk->rowPtrs[i], &HP, 0);
      data_wk->JSPSmat[i] = new realtype[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
    }
    if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
      Print() << "--> SPARSE Preconditioner -- non zero entries: "
              << data_wk->NNZ * data_wk->ncells << ", which represents "
              << data_wk->NNZ /
                   float(
                     (NUM_SPECIES + 1) * (NUM_SPECIES + 1) * data_wk->ncells) *
                   100.0
              << " % fill-in pattern\n";
    }
  } else if (data_wk->isolve_type == sparse_solve_custom) {
    // Number of non zero elements
    SPARSITY_INFO_SYST(&(data_wk->NNZ), &HP, 1);
    data_wk->PSc = SUNSparseMatrix(
      (NUM_SPECIES + 1) * data_wk->ncells, (NUM_SPECIES + 1) * data_wk->ncells,
      data_wk->NNZ * data_wk->ncells, CSR_MAT);
    data_wk->rowPtrs_c = (int*)SUNSparseMatrix_IndexPointers(data_wk->PSc);
    data_wk->colVals_c = (int*)SUNSparseMatrix_IndexValues(data_wk->PSc);
    SPARSITY_PREPROC_SYST_CSR(
      data_wk->colVals_c, data_wk->rowPtrs_c, &HP, data_wk->ncells, 0);
    if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
      Print() << "--> SPARSE solver -- non zero entries: "
              << data_wk->NNZ * data_wk->ncells << ", which represents "
              << data_wk->NNZ /
                   float(
                     (NUM_SPECIES + 1) * (NUM_SPECIES + 1) * data_wk->ncells) *
                   100.0
              << " % fill-in pattern\n";
    }
  } else if (data_wk->isolve_type == hack_dump_sparsity_pattern) {
    // Debug mode, makes no sense to call with OMP/MPI activated
    int counter;

    // CHEMISTRY JAC
    SPARSITY_INFO(&(data_wk->NNZ), &HP, 1);
    Print() << "--> Chem Jac -- non zero entries: " << data_wk->NNZ
            << ", which represents "
            << data_wk->NNZ / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                 100.0
            << " % fill-in pattern\n";
    SUNMatrix PS;
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), data_wk->NNZ, CSR_MAT);
    int *colIdx, *rowCount;
    rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_CSR(colIdx, rowCount, &HP, 1, 0);
    std::cout << " " << std::endl;
    std::cout << "*** Treating CHEM Jac (CSR symbolic analysis)***"
              << std::endl;
    std::cout << " " << std::endl;
    int nbVals;
    counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      nbVals = rowCount[i + 1] - rowCount[i];
      int* idx_arr = new int[nbVals];
      std::fill_n(idx_arr, nbVals, -1);
      std::memcpy(idx_arr, colIdx + rowCount[i], nbVals * sizeof(int));
      int idx = 0;
      for (int j = 0; j < NUM_SPECIES + 1; j++) {
        if ((j == idx_arr[idx]) && (nbVals > 0)) {
          std::cout << 1 << " ";
          idx = idx + 1;
          counter = counter + 1;
        } else {
          std::cout << 0 << " ";
        }
      }
      delete[] idx_arr;
      std::cout << std::endl;
    }
    std::cout << " There was " << counter << " non zero elems (compare to the "
              << data_wk->NNZ << " we need)" << std::endl;

    // SYST JAC
    SPARSITY_INFO_SYST(&(data_wk->NNZ), &HP, 1);
    Print() << "--> Syst Jac -- non zero entries: " << data_wk->NNZ
            << ", which represents "
            << data_wk->NNZ / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                 100.0
            << " % fill-in pattern\n";
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), data_wk->NNZ, CSR_MAT);
    rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_SYST_CSR(colIdx, rowCount, &HP, 1, 1);
    // CHEMISTRY JAC
    std::cout << " " << std::endl;
    std::cout << "*** Treating SYST Jac (CSR symbolic analysis)***"
              << std::endl;
    std::cout << " " << std::endl;
    counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      nbVals = rowCount[i + 1] - rowCount[i];
      int* idx_arr = new int[nbVals];
      std::fill_n(idx_arr, nbVals, -1);
      std::memcpy(idx_arr, colIdx + (rowCount[i] - 1), nbVals * sizeof(int));
      int idx = 0;
      for (int j = 0; j < NUM_SPECIES + 1; j++) {
        if ((j == idx_arr[idx] - 1) && ((nbVals - idx) > 0)) {
          std::cout << 1 << " ";
          idx = idx + 1;
          counter = counter + 1;
        } else {
          std::cout << 0 << " ";
        }
      }
      delete[] idx_arr;
      std::cout << std::endl;
    }
    std::cout << " There was " << counter << " non zero elems (compare to the "
              << data_wk->NNZ << " we need)" << std::endl;

    // SYST JAC SIMPLIFIED
    SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ), &HP);
    Print() << "--> Simplified Syst Jac (for Precond) -- non zero entries: "
            << data_wk->NNZ << ", which represents "
            << data_wk->NNZ / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                 100.0
            << " % fill-in pattern\n";
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), data_wk->NNZ, CSR_MAT);
    rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(colIdx, rowCount, &HP, 1);
    // CHEMISTRY JAC
    std::cout << " " << std::endl;
    std::cout << "*** Treating simplified SYST Jac (CSR symbolic analysis)***"
              << std::endl;
    std::cout << " " << std::endl;
    counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      nbVals = rowCount[i + 1] - rowCount[i];
      int* idx_arr = new int[nbVals];
      std::fill_n(idx_arr, nbVals, -1);
      std::memcpy(idx_arr, colIdx + (rowCount[i] - 1), nbVals * sizeof(int));
      int idx = 0;
      for (int j = 0; j < NUM_SPECIES + 1; j++) {
        if ((j == idx_arr[idx] - 1) && ((nbVals - idx) > 0)) {
          std::cout << 1 << " ";
          idx = idx + 1;
          counter = counter + 1;
        } else {
          std::cout << 0 << " ";
        }
      }
      delete[] idx_arr;
      std::cout << std::endl;
    }
    std::cout << " There was " << counter << " non zero elems (compare to the "
              << data_wk->NNZ << " we need)" << std::endl;

    Abort("Dump Sparsity Patern of different Jacobians in CSR format.");
  }

  return (data_wk);
}

void
reactor_close()
{

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  if (data->isolve_type == dense_solve) {
    SUNMatDestroy(A);
  }

  N_VDestroy(y);
  FreeUserData(data);
}

// Free data memory
// Probably not complete, how about the stuff allocated in KLU mode ?
void
FreeUserData(UserData data_wk)
{
  delete[](data_wk->Yvect_full);
  delete[](data_wk->rYsrc);
  delete[](data_wk->rhoX_init);
  delete[](data_wk->rhoXsrc_ext);
  delete[](data_wk->FCunt);
  delete[](data_wk->mask);

  delete[] data_wk->colPtrs;
  delete[] data_wk->rowVals;
  delete[] data_wk->Jdata;
#ifndef USE_KLU_PP
  if (data_wk->isolve_type == iterative_gmres_solve) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      destroyMat((data_wk->P)[i][i]);
      destroyMat((data_wk->Jbd)[i][i]);
      destroyArray((data_wk->pivot)[i][i]);
    }
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[](data_wk->P)[i];
      delete[](data_wk->Jbd)[i];
      delete[](data_wk->pivot)[i];
    }
    delete[](data_wk->P);
    delete[](data_wk->Jbd);
    delete[](data_wk->pivot);
    //}

#else
  if (data_wk->isolve_type == sparse_solve) {
    SUNMatDestroy(A);
    SUNMatDestroy((data_wk->PS)[0]);
    delete[](data_wk->PS);
  } else if (data_wk->isolve_type == iterative_gmres_solve) {
    delete[] data_wk->indx;
    for (int i = 0; i < data_wk->ncells; ++i) {
      klu_free_symbolic(&(data_wk->Symbolic[i]), &(data_wk->Common[i]));
      klu_free_numeric(&(data_wk->Numeric[i]), &(data_wk->Common[i]));
      delete[] data_wk->JSPSmat[i];
      SUNMatDestroy((data_wk->PS)[i]);
    }
    delete[] data_wk->JSPSmat;
    delete[] data_wk->Common;
    delete[] data_wk->Symbolic;
    delete[] data_wk->Numeric;
    delete[] data_wk->PS;
    //}
#endif

  } else if (data_wk->isolve_type == iterative_gmres_solve_custom) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[] data_wk->JSPSmat[i];
      SUNMatDestroy((data_wk->PS)[i]);
    }
    delete[] data_wk->colVals;
    delete[] data_wk->rowPtrs;
    delete[] data_wk->PS;
    delete[] data_wk->JSPSmat;
  } else if (data_wk->isolve_type == sparse_solve_custom) {
    SUNMatDestroy(A);
    SUNMatDestroy(data_wk->PSc);
  } else if (data_wk->isolve_type == hack_dump_sparsity_pattern) {
  }

  free(data_wk);
}

#endif
