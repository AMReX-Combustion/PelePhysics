#include "ReactorCvode.H"

#include <iostream>

namespace pele {
namespace physics {
namespace reactions {

int
ReactorCvode::init(int reactor_type, int ncells)
{
  BL_PROFILE("Pele::ReactorCvode::init()");
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
  amrex::ParmParse pp("ode");
  pp.query("rtol", relTol);
  pp.query("atol", absTol);
  checkCvodeOptions();

#ifndef AMREX_USE_GPU
  // ----------------------------------------------------------
  // On CPU, initialize cvode_mem/userData
  // ----------------------------------------------------------

  // Solution vector
  int neq_tot = (NUM_SPECIES + 1) * ncells;
  y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNew_Serial", 0)) {
    return (1);
  }

  // Call CVodeCreate to create the solver memory and specify the Backward
  // Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)cvode_mem, "CVodeCreate", 0)) {
    return (1);
  }

  udata_g = new CVODEUserData{};
  allocUserData(udata_g, ncells);
  if (utils::check_flag((void*)udata_g, "allocUserData", 2)) {
    return (1);
  }

  // Set the pointer to user-defined data
  int flag = CVodeSetUserData(cvode_mem, udata_g);
  if (utils::check_flag(&flag, "CVodeSetUserData", 1)) {
    return (1);
  }

  // Call CVodeInit to initialize the integrator memory and specify the user's
  // right hand side function, the inital time, and initial dependent variable
  // vector y.
  amrex::Real time = 0.0;
  flag = CVodeInit(cvode_mem, cF_RHS, time, y);
  if (utils::check_flag(&flag, "CVodeInit", 1)) {
    return (1);
  }

  // Setup tolerances
  set_sundials_solver_tols(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, udata_g->ncells,
    udata_g->verbose, relTol, absTol, "cvode");

  // Linear solver data
  if (
    udata_g->solve_type == cvode::denseFDDirect ||
    udata_g->solve_type == cvode::denseDirect) {
    // Create dense SUNMatrix for use in linear solves
    A = SUNDenseMatrix(
      neq_tot, neq_tot, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)A, "SUNDenseMatrix", 0)) {
      return (1);
    }

    // Create dense SUNLinearSolver object for use by CVode
    LS = SUNLinSol_Dense(y, A, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_Dense", 0)) {
      return (1);
    }

    // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1)) {
      return (1);
    }

  } else if (udata_g->solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
    // Create sparse SUNMatrix for use in linear solves
    A = SUNSparseMatrix(
      neq_tot, neq_tot, (udata_g->NNZ) * udata_g->ncells, CSC_MAT,
      *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)A, "SUNSparseMatrix", 0))
      return (1);

    // Create KLU solver object for use by CVode
    LS = SUNLinSol_KLU(y, A, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_KLU", 0))
      return (1);

    // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
#else
    amrex::Abort("sparseDirect solver_type not valid without KLU library.");
#endif

  } else if (udata_g->solve_type == cvode::customDirect) {
    // Create dense SUNMatrix for use in linear solves
    A = SUNSparseMatrix(
      neq_tot, neq_tot, (udata_g->NNZ) * udata_g->ncells, CSR_MAT,
      *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)A, "SUNDenseMatrix", 0)) {
      return (1);
    }

    // Create dense SUNLinearSolver object for use by CVode
    LS = cvode::SUNLinSol_sparse_custom(
      y, A, reactor_type, udata_g->ncells, (NUM_SPECIES + 1), udata_g->NNZ,
      *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_sparse_custom", 0)) {
      return (1);
    }

    // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1)) {
      return (1);
    }

  } else if (udata_g->solve_type == cvode::GMRES) {
    // Create the GMRES linear solver object
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_NONE, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) {
      return (1);
    }

    // Set CVode linear solver to LS
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1)) {
      return (1);
    }

  } else if (udata_g->solve_type == cvode::precGMRES) {
    // Create the GMRES linear solver object
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_LEFT, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) {
      return (1);
    }

    // Set CVode linear solver to LS
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1)) {
      return (1);
    }

  } else {
    amrex::Abort("Wrong choice of linear solver...");
  }

  // Analytical Jac. data for direct solver
  if (udata_g->analytical_jacobian == 1) {
    if (udata_g->solve_type == cvode::denseDirect) {
      // Set the user-supplied Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cvode::cJac);
      if (utils::check_flag(&flag, "CVodeSetJacFn", 1)) {
        return (1);
      }
    } else if (udata_g->solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
      // Set the user-supplied KLU Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cvode::cJac_KLU);
      if (utils::check_flag(&flag, "CVodeSetJacFn", 1))
        return (1);
#else
      amrex::Abort(
        "Shouldn't be there: sparseDirect solver_type not valid without "
        "KLU library.");
#endif
    } else if (udata_g->solve_type == cvode::customDirect) {
      // Set the user-supplied Jacobian routine Jac
      flag = CVodeSetJacFn(cvode_mem, cvode::cJac_sps);
      if (utils::check_flag(&flag, "CVodeSetJacFn", 1)) {
        return (1);
      }
    }
  }

  // Analytical Jac. data for iterative solver preconditioner
  if (udata_g->precond_type == cvode::denseSimpleAJac) {
    // Set the JAcobian-times-vector function
    flag = CVodeSetJacTimes(cvode_mem, nullptr, nullptr);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1)) {
      return (1);
    }
    // Set the preconditioner plain dense solve and setup functions
    flag = CVodeSetPreconditioner(cvode_mem, cvode::Precond, cvode::PSolve);
    if (utils::check_flag(&flag, "CVodeSetPreconditioner", 1)) {
      return (1);
    }
  } else if (udata_g->precond_type == cvode::sparseSimpleAJac) {
#ifdef PELE_USE_KLU
    // Set the JAcobian-times-vector function
    flag = CVodeSetJacTimes(cvode_mem, nullptr, nullptr);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);
    // Set the preconditioner KLU sparse solve and setup functions
    flag = CVodeSetPreconditioner(
      cvode_mem, cvode::Precond_sparse, cvode::PSolve_sparse);
    if (utils::check_flag(&flag, "CVodeSetPreconditioner", 1))
      return (1);
#else
    amrex::Abort(
      "sparseSimpleAJac precond_type not valid without KLU library.");
#endif
  } else if (udata_g->precond_type == cvode::customSimpleAJac) {
    // Set the JAcobian-times-vector function
    flag = CVodeSetJacTimes(cvode_mem, nullptr, nullptr);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1)) {
      return (1);
    }
    // Set the preconditioner to custom solve and setup functions
    flag = CVodeSetPreconditioner(
      cvode_mem, cvode::Precond_custom, cvode::PSolve_custom);
    if (utils::check_flag(&flag, "CVodeSetPreconditioner", 1)) {
      return (1);
    }
  }

  // CVODE runtime options
  flag = CVodeSetMaxNonlinIters(cvode_mem, 50); // Max newton iter.
  if (utils::check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) {
    return (1);
  }
  flag = CVodeSetMaxErrTestFails(cvode_mem, 100); // Max Err.test failure
  if (utils::check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) {
    return (1);
  }
  flag = CVodeSetErrHandlerFn(
    cvode_mem, cvode::cvodeErrHandler, nullptr); // Err. handler funct.
  if (utils::check_flag(&flag, "CVodeSetErrHandlerFn", 1)) {
    return (1);
  }
  flag = CVodeSetMaxNumSteps(cvode_mem, 10000); // Max substeps
  if (utils::check_flag(&flag, "CVodeSetMaxNumSteps", 1)) {
    return (1);
  }
  flag = CVodeSetMaxOrd(cvode_mem, udata_g->maxOrder); // Max order
  if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1)) {
    return (1);
  }
  flag = CVodeSetJacEvalFrequency(cvode_mem, 100); // Max Jac age
  if (utils::check_flag(&flag, "CVodeSetJacEvalFrequency", 1)) {
    return (1);
  }

  // End of CPU section
#endif

  return (0);
}

void
ReactorCvode::checkCvodeOptions() const
{
  // Query options
  amrex::ParmParse pp("ode");
  int verbose = 0;
  pp.query("verbose", verbose);

  if (verbose > 0) {
    amrex::Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
  }

  std::string solve_type_str = "none";
  amrex::ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  int solve_type = -1;
  int analytical_jacobian = 0;
  int precond_type = -1;

#ifdef AMREX_USE_GPU
  if (solve_type_str == "sparse_direct") {
    solve_type = cvode::sparseDirect;
    analytical_jacobian = 1;
#ifdef AMREX_USE_CUDA
    if (verbose > 0)
      amrex::Print()
        << " Using a cuSparse direct linear solve with analytical Jacobian\n";
#else
    amrex::Abort("solve_type 'sparse_direct' only available with CUDA");
#endif
  } else if (solve_type_str == "custom_direct") {
    solve_type = cvode::customDirect;
    analytical_jacobian = 1;
#ifdef AMREX_USE_CUDA
    if (verbose > 0)
      amrex::Print()
        << " Using a custom direct linear solve with analytical Jacobian\n";
#else
    amrex::Abort("solve_type 'custom_direct' only available with CUDA");
#endif
  } else if (solve_type_str == "magma_direct") {
    solve_type = cvode::magmaDirect;
    analytical_jacobian = 1;
#ifdef PELE_USE_MAGMA
    if (verbose > 0)
      amrex::Print() << " Using MAGMA direct linear solve\n";
#else
    amrex::Abort(
      "solve_type 'magma_direct' only available with if PELE_USE_MAGMA true");
#endif
  } else if (solve_type_str == "GMRES") {
    solve_type = cvode::GMRES;
    if (verbose > 0)
      amrex::Print() << " Using a JFNK GMRES linear solve\n";
  } else if (solve_type_str == "precGMRES") {
    solve_type = cvode::precGMRES;
    if (verbose > 0)
      amrex::Print() << " Using a JFNK GMRES linear solve";
    std::string prec_type_str = "cuSparse_simplified_AJacobian";
    ppcv.query("precond_type", prec_type_str);
    if (prec_type_str == "cuSparse_simplified_AJacobian") {
      precond_type = cvode::sparseSimpleAJac;
#ifdef AMREX_USE_CUDA
      amrex::Print() << " with a cuSparse simplified AJ-based preconditioner";
#else
      amrex::Abort(
        "precond_type 'cuSparse_simplified_AJacobian' only available with "
        "CUDA");
#endif
    } else {
      amrex::Abort(
        "Wrong precond_type. Only option is: 'cuSparse_simplified_AJacobian'");
    }
    amrex::Print() << "\n";
  } else {
    amrex::Abort(
      "Wrong solve_type. Options are: 'sparse_direct', 'custom_direct', "
      "'GMRES', 'precGMRES'");
  }

#else
  if (solve_type_str == "dense_direct") {
    solve_type = cvode::denseFDDirect;
    if (verbose > 0) {
      amrex::Print()
        << " Using a dense direct linear solve with Finite Difference "
           "Jacobian\n";
    }
  } else if (solve_type_str == "denseAJ_direct") {
    solve_type = cvode::denseDirect;
    analytical_jacobian = 1;
    if (verbose > 0) {
      amrex::Print()
        << " Using a dense direct linear solve with Analytical Jacobian\n";
    }
  } else if (solve_type_str == "sparse_direct") {
    solve_type = cvode::sparseDirect;
    analytical_jacobian = 1;
#ifndef PELE_USE_KLU
    amrex::Abort("solver_type sparse_direct requires the KLU library");
#endif
    if (verbose > 0) {
      amrex::Print()
        << " Using a sparse direct linear solve with KLU Analytical Jacobian\n";
    }
  } else if (solve_type_str == "custom_direct") {
    solve_type = cvode::customDirect;
    analytical_jacobian = 1;
    if (verbose > 0) {
      amrex::Print()
        << " Using a sparse custom direct linear solve with Analytical "
           "Jacobian\n";
    }
  } else if (solve_type_str == "GMRES") {
    solve_type = cvode::GMRES;
    if (verbose > 0) {
      amrex::Print() << " Using a JFNK GMRES linear solve\n";
    }
  } else if (solve_type_str == "precGMRES") {
    solve_type = cvode::precGMRES;
    if (verbose > 0) {
      amrex::Print() << " Using a JFNK GMRES linear solve";
    }
    std::string prec_type_str = "none";
    ppcv.query("precond_type", prec_type_str);
    if (prec_type_str == "dense_simplified_AJacobian") {
      precond_type = cvode::denseSimpleAJac;
      amrex::Print() << " with a dense simplified AJ-based preconditioner";
    } else if (prec_type_str == "sparse_simplified_AJacobian") {
      precond_type = cvode::sparseSimpleAJac;
#ifndef PELE_USE_KLU
      amrex::Abort(
        "precond_type sparse_simplified_AJacobian requires the KLU library");
#endif
      amrex::Print() << " with a sparse simplified AJ-based preconditioner";
    } else if (prec_type_str == "custom_simplified_AJacobian") {
      precond_type = cvode::customSimpleAJac;
      amrex::Print() << " with a custom simplified AJ-based preconditioner";
    } else {
      amrex::Abort(
        "Wrong precond_type. Options are: 'dense_simplified_AJacobian', "
        "'sparse_simplified_AJacobian', 'custom_simplified_AJacobian'");
    }
    amrex::Print() << "\n";
  } else if (solve_type_str == "diagnostic") {
    solve_type = cvode::hackDumpSparsePattern;
  } else {
    amrex::Abort(
      "Wrong solve_type. Options are: 'dense_direct', denseAJ_direct', "
      "'sparse_direct', 'custom_direct', 'GMRES', 'precGMRES'");
  }
#endif

  // Print additionnal information
  if (precond_type == cvode::sparseSimpleAJac) {
    int nJdata;
    const int HP = m_reactor_type == ReactorTypes::h_reactor_type;
    // Simplified AJ precond data
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA)
    SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata, &HP);
    if (verbose > 0) {
      amrex::Print()
        << "--> cuSparse AJ based matrix Preconditioner -- non zero entries: "
        << nJdata << ", which represents "
        << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
        << " % fill-in pattern\n";
    }
#elif defined(AMREX_USE_HIP)
    amrex::Abort(
      "\n--> precond_type sparse simplified_AJacobian not available with "
      "HIP \n");
#elif defined(AMREX_USE_DPCPP)
    amrex::Abort(
      "\n--> precond_type sparse simplified_AJacobian not available with "
      "DPCPP \n");
#endif

#else
    SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata, &HP);
    if (verbose > 0) {
      amrex::Print()
        << "--> KLU sparse AJ based matrix Preconditioner -- non zero entries: "
        << nJdata << ", which represents "
        << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
        << " % fill-in pattern\n";
    }
#endif
#ifndef AMREX_USE_GPU
  } else if (precond_type == cvode::customSimpleAJac) {
    int nJdata;
    const int HP = m_reactor_type == ReactorTypes::h_reactor_type;
    // Simplified AJ precond data
    SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata, &HP);
    if (verbose > 0) {
      amrex::Print()
        << "--> custom sparse AJ based matrix Preconditioner -- non zero "
           "entries: "
        << nJdata << ", which represents "
        << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
        << " % fill-in pattern\n";
    }
#endif
  }

  if (analytical_jacobian == 1) {
    int nJdata;
    const int HP = m_reactor_type == ReactorTypes::h_reactor_type;
    int ncells = 1; // Print the pattern of the diagonal block. ncells will
                    // actually vary on GPU.
#ifdef AMREX_USE_GPU
    if (solve_type == cvode::sparseDirect) {
#if defined(AMREX_USE_CUDA)
      SPARSITY_INFO_SYST(&nJdata, &HP, ncells);
      if (verbose > 0) {
        amrex::Print()
          << "--> cuSparse based matrix Solver -- non zero entries: " << nJdata
          << ", which represents "
          << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
          << " % fill-in pattern\n";
      }
#elif defined(AMREX_USE_HIP)
      amrex::Abort("\n--> Analytical Jacobian not available with HIP. Change "
                   "solve_type.\n");
#elif defined(AMREX_USE_DPCPP)
      amrex::Abort("\n--> Analytical Jacobian not available with DPCPP. Change "
                   "solve_type.\n");
#endif
    }

#else
    if (solve_type == cvode::customDirect) {
      SPARSITY_INFO_SYST(&nJdata, &HP, ncells);
      if (verbose > 0) {
        amrex::Print()
          << "--> sparse AJ-based matrix custom Solver -- non zero entries: "
          << nJdata << ", which represents "
          << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
          << " % fill-in pattern\n";
      }
    } else if (solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
      SPARSITY_INFO(&nJdata, &HP, ncells);
      if (verbose > 0) {
        amrex::Print()
          << "--> KLU sparse AJ-based matrix Solver -- non zero entries: "
          << nJdata << ", which represents "
          << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
          << " % fill-in pattern\n";
      }
#else
      amrex::Abort(
        "solver_type 'sparseDirect' uses a sparse KLU matrix and requires "
        "the KLU library");
#endif
    }
#endif
  }

#ifndef AMREX_USE_GPU
  if (solve_type == cvode::hackDumpSparsePattern) {
    // This is a diagnostic option -> dump sparsity pattern and abort.
    // Reactor type
    const int HP = m_reactor_type == ReactorTypes::h_reactor_type;

    // CHEMISTRY JAC
    int nJdata = 0;
    SPARSITY_INFO(&nJdata, &HP, 1);
    amrex::Print() << "--> Chem. Jac -- non zero entries: " << nJdata
                   << ", which represents "
                   << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                        100.0
                   << " % fill-in pattern\n";
    SUNMatrix PS;
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), nJdata, CSR_MAT,
      *amrex::sundials::The_Sundials_Context());
    int* rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    int* colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_CSR(colIdx, rowCount, &HP, 1, 0);
    amrex::Print()
      << "\n\n *** Treating CHEM Jac (CSR symbolic analysis)*** \n\n";
    int counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      int nbVals = rowCount[i + 1] - rowCount[i];
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
    amrex::Print() << " There was " << counter
                   << " non zero elems (compared to the " << nJdata
                   << " we need) \n";

    // SYST JAC
    SPARSITY_INFO_SYST(&nJdata, &HP, 1);
    amrex::Print() << "--> Syst. Jac -- non zero entries: " << nJdata
                   << ", which represents "
                   << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) *
                        100.0
                   << " % fill-in pattern\n";
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), nJdata, CSR_MAT,
      *amrex::sundials::The_Sundials_Context());
    rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_SYST_CSR(colIdx, rowCount, &HP, 1, 1);
    amrex::Print()
      << "\n\n *** Treating SYST Jac (CSR symbolic analysis)*** \n\n";
    counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      int nbVals = rowCount[i + 1] - rowCount[i];
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
    amrex::Print() << " There was " << counter
                   << " non zero elems (compared to the " << nJdata
                   << " we need) \n";

    // SYST JAC SIMPLIFIED
    SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata, &HP);
    amrex::Print()
      << "--> Simplified Syst Jac (for Precond) -- non zero entries: " << nJdata
      << ", which represents "
      << nJdata / float((NUM_SPECIES + 1) * (NUM_SPECIES + 1)) * 100.0
      << " % fill-in pattern\n";
    PS = SUNSparseMatrix(
      (NUM_SPECIES + 1), (NUM_SPECIES + 1), nJdata, CSR_MAT,
      *amrex::sundials::The_Sundials_Context());
    rowCount = (int*)SUNSparseMatrix_IndexPointers(PS);
    colIdx = (int*)SUNSparseMatrix_IndexValues(PS);
    SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(colIdx, rowCount, &HP, 1);
    amrex::Print() << "\n\n *** Treating simplified SYST Jac (CSR symbolic "
                      "analysis)*** \n\n";
    counter = 0;
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      int nbVals = rowCount[i + 1] - rowCount[i];
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
    amrex::Print() << " There was " << counter
                   << " non zero elems (compared to the " << nJdata
                   << " we need) \n";

    amrex::Abort("Chemistry sparsity pattern dumped -> exiting now !");
  }
#endif
}

void
ReactorCvode::allocUserData(
  CVODEUserData* udata,
  int a_ncells
#ifdef AMREX_USE_GPU
  ,
  SUNMatrix& a_A,
  amrex::gpuStream_t stream
#endif
) const
{
  // Query options
  amrex::ParmParse pp("ode");
  int verbose = 0;
  pp.query("verbose", verbose);

  std::string solve_type_str = "none";
  amrex::ParmParse ppcv("cvode");
  udata->maxOrder = 2;
  ppcv.query("max_order", udata->maxOrder);
  ppcv.query("solve_type", solve_type_str);
  // Defaults
  udata->solve_type = -1;
  udata->analytical_jacobian = 0;
  udata->precond_type = -1;

#ifdef AMREX_USE_GPU
  if (solve_type_str == "sparse_direct") {
    udata->solve_type = cvode::sparseDirect;
    udata->analytical_jacobian = 1;
  } else if (solve_type_str == "custom_direct") {
    udata->solve_type = cvode::customDirect;
    udata->analytical_jacobian = 1;
  } else if (solve_type_str == "magma_direct") {
    udata->solve_type = cvode::magmaDirect;
    udata->analytical_jacobian = 1;
  } else if (solve_type_str == "GMRES") {
    udata->solve_type = cvode::GMRES;
  } else if (solve_type_str == "precGMRES") {
    udata->solve_type = cvode::precGMRES;
    std::string prec_type_str = "cuSparse_simplified_AJacobian";
    ppcv.query("precond_type", prec_type_str);
    if (prec_type_str == "cuSparse_simplified_AJacobian") {
      udata->precond_type = cvode::sparseSimpleAJac;
    } else {
      amrex::Abort(
        "Wrong precond_type. Only option is: 'cuSparse_simplified_AJacobian'");
    }
    amrex::Print() << "\n";
  } else {
    amrex::Abort(
      "Wrong solve_type. Options are: 'sparse_direct', 'custom_direct', "
      "'GMRES', 'precGMRES'");
  }

#else
  if (solve_type_str == "dense_direct") {
    udata->solve_type = cvode::denseFDDirect;
  } else if (solve_type_str == "denseAJ_direct") {
    udata->solve_type = cvode::denseDirect;
    udata->analytical_jacobian = 1;
  } else if (solve_type_str == "sparse_direct") {
    udata->solve_type = cvode::sparseDirect;
    udata->analytical_jacobian = 1;
#ifndef PELE_USE_KLU
    amrex::Abort("solver_type sparse_direct requires the KLU library");
#endif
  } else if (solve_type_str == "custom_direct") {
    udata->solve_type = cvode::customDirect;
    udata->analytical_jacobian = 1;
  } else if (solve_type_str == "GMRES") {
    udata->solve_type = cvode::GMRES;
  } else if (solve_type_str == "precGMRES") {
    udata->solve_type = cvode::precGMRES;
    std::string prec_type_str = "sparse_simplified_AJacobian";
    ppcv.query("precond_type", prec_type_str);
    if (prec_type_str == "dense_simplified_AJacobian") {
      udata->precond_type = cvode::denseSimpleAJac;
    } else if (prec_type_str == "sparse_simplified_AJacobian") {
      udata->precond_type = cvode::sparseSimpleAJac;
#ifndef PELE_USE_KLU
      amrex::Abort(
        "precond_type sparse_simplified_AJacobian requires the KLU library");
#endif
    } else if (prec_type_str == "custom_simplified_AJacobian") {
      udata->precond_type = cvode::customSimpleAJac;
    } else {
      amrex::Abort(
        "Wrong precond_type. Options are: 'dense_simplified_AJacobian', "
        "'sparse_simplified_AJacobian', 'custom_simplified_AJacobian'");
    }
    amrex::Print() << "\n";
  } else {
    amrex::Abort(
      "Wrong solve_type. Options are: 'dense_direct', denseAJ_direct', "
      "'sparse_direct', 'custom_direct', 'GMRES', 'precGMRES'");
  }
#endif

  // Pass options to udata
  const int HP = m_reactor_type == ReactorTypes::h_reactor_type;
  int nspec_tot = (NUM_SPECIES)*a_ncells;
  udata->reactor_type = m_reactor_type;
  udata->ncells = a_ncells;
  udata->verbose = verbose;
#ifdef AMREX_USE_GPU
  udata->nbThreads = 32;
  udata->nbBlocks = std::max(1, a_ncells / udata->nbThreads);
  udata->stream = stream;
#endif

  // Alloc internal udata solution/forcing containers
  udata->rYsrc_ext =
    (amrex::Real*)amrex::The_Arena()->alloc(nspec_tot * sizeof(amrex::Real));
  udata->rhoe_init =
    (amrex::Real*)amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real));
  udata->rhoesrc_ext =
    (amrex::Real*)amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real));
  udata->mask = (int*)amrex::The_Arena()->alloc(a_ncells * sizeof(int));

#ifndef AMREX_USE_GPU
  udata->FCunt = (int*)amrex::The_Arena()->alloc(a_ncells * sizeof(int));
  udata->FirstTimePrecond = true;
#endif

  // Alloc internal udata Analytical Jacobian containers
#ifdef AMREX_USE_GPU
  if (udata->solve_type == cvode::sparseDirect) {
#ifdef AMREX_USE_CUDA
    SPARSITY_INFO_SYST(&(udata->NNZ), &HP, 1);
    udata->csr_row_count_h =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_h =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));
    udata->csr_row_count_d =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_d =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));

    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cusolver_status = cusolverSpCreate(&(udata->cusolverHandle));
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusolver_status = cusolverSpSetStream(udata->cusolverHandle, stream);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
    cusparse_status = cusparseCreate(&(udata->cuSPHandle));
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparse_status = cusparseSetStream(udata->cuSPHandle, stream);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);

    a_A = SUNMatrix_cuSparse_NewBlockCSR(
      a_ncells, (NUM_SPECIES + 1), (NUM_SPECIES + 1), udata->NNZ,
      udata->cuSPHandle, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)a_A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) {
      amrex::Abort("Something went wrong while creating cuSparse_NewBlockCSR");
    }

    int retval = SUNMatrix_cuSparse_SetFixedPattern(a_A, 1);
    // if (utils::check_flag(&retval, "SUNMatrix_cuSparse_SetFixedPattern", 1))
    // return(1);

    SPARSITY_PREPROC_SYST_CSR(
      udata->csr_col_index_h, udata->csr_row_count_h, &HP, 1, 0);
    amrex::Gpu::htod_memcpy(
      &udata->csr_col_index_d, &udata->csr_col_index_h,
      sizeof(udata->csr_col_index_h));
    amrex::Gpu::htod_memcpy(
      &udata->csr_row_count_d, &udata->csr_row_count_h,
      sizeof(udata->csr_row_count_h));
    SUNMatrix_cuSparse_CopyToDevice(
      a_A, NULL, udata->csr_row_count_h, udata->csr_col_index_h);
#else
    amrex::Abort(
      "Solver_type sparse_direct is only available with CUDA on GPU");
#endif
  } else if (udata->solve_type == cvode::customDirect) {
#ifdef AMREX_USE_CUDA
    SPARSITY_INFO_SYST(&(udata->NNZ), &HP, 1);
    udata->csr_row_count_h =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_h =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));
    udata->csr_row_count_d =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_d =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));

    cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
    cusparse_status = cusparseCreate(&(udata->cuSPHandle));
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparse_status = cusparseSetStream(udata->cuSPHandle, stream);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);

    a_A = SUNMatrix_cuSparse_NewBlockCSR(
      a_ncells, (NUM_SPECIES + 1), (NUM_SPECIES + 1), udata->NNZ,
      udata->cuSPHandle, *amrex::sundials::The_Sundials_Context());
    // if (utils::check_flag((void *)a_A, "SUNMatrix_cuSparse_NewBlockCSR", 0))
    // return(1);

    int retval = SUNMatrix_cuSparse_SetFixedPattern(a_A, 1);
    // if(utils::check_flag(&retval, "SUNMatrix_cuSparse_SetFixedPattern", 1))
    // return(1);

    SPARSITY_PREPROC_SYST_CSR(
      udata->csr_col_index_h, udata->csr_row_count_h, &HP, 1, 0);
    amrex::Gpu::htod_memcpy(
      &udata->csr_col_index_d, &udata->csr_col_index_h,
      sizeof(udata->csr_col_index_h));
    amrex::Gpu::htod_memcpy(
      &udata->csr_row_count_d, &udata->csr_row_count_h,
      sizeof(udata->csr_row_count_h));
    SUNMatrix_cuSparse_CopyToDevice(
      a_A, NULL, udata->csr_row_count_h, udata->csr_col_index_h);
#else
    amrex::Abort(
      "Solver_type custom_direct is only available with CUDA on GPU");
#endif
  } else if (udata->solve_type == cvode::magmaDirect) {
#ifdef PELE_USE_MAGMA
    a_A = SUNMatrix_MagmaDenseBlock(
      a_ncells, (NUM_SPECIES + 1), (NUM_SPECIES + 1), SUNMEMTYPE_DEVICE,
      *amrex::sundials::The_SUNMemory_Helper(), NULL);
#else
    amrex::Abort("Solver_type magma_direct reauires PELE_USE_MAGMA = TRUE");
#endif
  }

#else
  if (udata->solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
    // CSC matrices data -> one big matrix used for the direct solve
    udata->colPtrs = new int*[1];
    udata->rowVals = new int*[1];
    udata->Jdata = new amrex::Real*[1];

    // Number of non zero elements in ODE system
    SPARSITY_INFO(&(udata->NNZ), &HP, udata->ncells);
    // Build Sparse Matrix for direct sparse KLU solver
    (udata->PS) = new SUNMatrix[1];
    (udata->PS)[0] = SUNSparseMatrix(
      (NUM_SPECIES + 1) * udata->ncells, (NUM_SPECIES + 1) * udata->ncells,
      udata->NNZ * udata->ncells, CSC_MAT,
      *amrex::sundials::The_Sundials_Context());
    udata->colPtrs[0] = (int*)SUNSparseMatrix_IndexPointers((udata->PS)[0]);
    udata->rowVals[0] = (int*)SUNSparseMatrix_IndexValues((udata->PS)[0]);
    udata->Jdata[0] = SUNSparseMatrix_Data((udata->PS)[0]);
    SPARSITY_PREPROC_CSC(
      udata->rowVals[0], udata->colPtrs[0], &HP, udata->ncells);
#endif
  } else if (udata->solve_type == cvode::customDirect) {
    // Number of non zero elements in ODE system
    SPARSITY_INFO_SYST(&(udata->NNZ), &HP, udata->ncells);
    // Build the SUNmatrix as CSR sparse and fill ptrs to row/Vals
    udata->PSc = SUNSparseMatrix(
      (NUM_SPECIES + 1) * udata->ncells, (NUM_SPECIES + 1) * udata->ncells,
      udata->NNZ * udata->ncells, CSR_MAT,
      *amrex::sundials::The_Sundials_Context());
    udata->rowPtrs_c = (int*)SUNSparseMatrix_IndexPointers(udata->PSc);
    udata->colVals_c = (int*)SUNSparseMatrix_IndexValues(udata->PSc);
    SPARSITY_PREPROC_SYST_CSR(
      udata->colVals_c, udata->rowPtrs_c, &HP, udata->ncells, 0);
  }
#endif

  // Alloc internal udata Preconditioner containers
#ifdef AMREX_USE_GPU
  if (udata->precond_type == cvode::sparseSimpleAJac) {
#ifdef AMREX_USE_CUDA
    SPARSITY_INFO_SYST_SIMPLIFIED(&(udata->NNZ), &HP);
    udata->csr_row_count_h =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_h =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));

    udata->csr_row_count_d =
      (int*)amrex::The_Arena()->alloc((NUM_SPECIES + 2) * sizeof(int));
    udata->csr_col_index_d =
      (int*)amrex::The_Arena()->alloc(udata->NNZ * sizeof(int));
    udata->csr_jac_d = (amrex::Real*)amrex::The_Arena()->alloc(
      udata->NNZ * a_ncells * sizeof(amrex::Real));
    udata->csr_val_d = (amrex::Real*)amrex::The_Arena()->alloc(
      udata->NNZ * a_ncells * sizeof(amrex::Real));

    SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
      udata->csr_col_index_h, udata->csr_row_count_h, &HP, 1);

    amrex::Gpu::htod_memcpy(
      &udata->csr_col_index_d, &udata->csr_col_index_h,
      sizeof(udata->NNZ * sizeof(int)));
    amrex::Gpu::htod_memcpy(
      &udata->csr_row_count_d, &udata->csr_row_count_h,
      sizeof((NUM_SPECIES + 2) * sizeof(int)));

    size_t workspaceInBytes = 0;
    size_t internalDataInBytes = 0;

    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cusolver_status = cusolverSpCreate(&(udata->cusolverHandle));
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusolver_status = cusolverSpSetStream(udata->cusolverHandle, stream);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
    cusparse_status = cusparseCreateMatDescr(&(udata->descrA));
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparse_status =
      cusparseSetMatType(udata->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusparse_status =
      cusparseSetMatIndexBase(udata->descrA, CUSPARSE_INDEX_BASE_ONE);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
    cusolver_status = cusolverSpCreateCsrqrInfo(&(udata->info));
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    // symbolic analysis
    cusolver_status = cusolverSpXcsrqrAnalysisBatched(
      udata->cusolverHandle,
      NUM_SPECIES + 1, // size per subsystem
      NUM_SPECIES + 1, // size per subsystem
      udata->NNZ, udata->descrA, udata->csr_row_count_h, udata->csr_col_index_h,
      udata->info);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    /*
    size_t free_mem = 0;
    size_t total_mem = 0;
    cudaStat1 = cudaMemGetInfo( &free_mem, &total_mem );
    AMREX_ASSERT( cudaSuccess == cudaStat1 );
    std::cout<<"(AFTER SA) Free: "<< free_mem<< " Tot: "<<total_mem<<std::endl;
    */

    // allocate working space
    cusolver_status = cusolverSpDcsrqrBufferInfoBatched(
      udata->cusolverHandle,
      NUM_SPECIES + 1, // size per subsystem
      NUM_SPECIES + 1, // size per subsystem
      udata->NNZ, udata->descrA, udata->csr_val_d, udata->csr_row_count_h,
      udata->csr_col_index_h, a_ncells, udata->info, &internalDataInBytes,
      &workspaceInBytes);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    // amrex::Print() << " BufferInfo workspaceInBytes " << workspaceInBytes <<
    // "\n";

    cudaError_t cudaStat1 = cudaSuccess;
    cudaStat1 = cudaMalloc((void**)&(udata->buffer_qr), workspaceInBytes);
    AMREX_ASSERT(cudaStat1 == cudaSuccess);
#else
    amrex::Abort(
      "cuSparse_simplified_AJacobian is only available with CUDA on GPU");
#endif
  }

#else
  if (udata->precond_type == cvode::denseSimpleAJac) {
    // Matrix data : big bunch of dimensions, not sure why. Generally ncells ==
    // 1 so not too bad Simply create the space.
    (udata->P) = new amrex::Real***[udata->ncells];
    (udata->Jbd) = new amrex::Real***[udata->ncells];
    (udata->pivot) = new sunindextype**[udata->ncells];
    for (int i = 0; i < udata->ncells; ++i) {
      (udata->P)[i] = new amrex::Real**[udata->ncells];
      (udata->Jbd)[i] = new amrex::Real**[udata->ncells];
      (udata->pivot)[i] = new sunindextype*[udata->ncells];
    }
    for (int i = 0; i < udata->ncells; ++i) {
      (udata->P)[i][i] =
        SUNDlsMat_newDenseMat(NUM_SPECIES + 1, NUM_SPECIES + 1);
      (udata->Jbd)[i][i] =
        SUNDlsMat_newDenseMat(NUM_SPECIES + 1, NUM_SPECIES + 1);
      (udata->pivot)[i][i] = SUNDlsMat_newIndexArray(NUM_SPECIES + 1);
    }
  } else if (udata->precond_type == cvode::sparseSimpleAJac) {
#ifdef PELE_USE_KLU
    // CSC matrices data for each submatrix (cells)
    udata->colPtrs = new int*[udata->ncells];
    udata->rowVals = new int*[udata->ncells];
    udata->Jdata = new amrex::Real*[udata->ncells];

    // KLU internal storage
    udata->Common = new klu_common[udata->ncells];
    udata->Symbolic = new klu_symbolic*[udata->ncells];
    udata->Numeric = new klu_numeric*[udata->ncells];
    // Sparse Matrices for It Sparse KLU block-solve
    udata->PS = new SUNMatrix[udata->ncells];
    // Number of non zero elements
    SPARSITY_INFO_SYST_SIMPLIFIED(&(udata->NNZ), &HP);
    // Not used yet. TODO use to fetch sparse Mat
    udata->indx = new int[udata->NNZ];
    udata->JSPSmat = new amrex::Real*[udata->ncells];
    for (int i = 0; i < udata->ncells; ++i) {
      (udata->PS)[i] = SUNSparseMatrix(
        NUM_SPECIES + 1, NUM_SPECIES + 1, udata->NNZ, CSC_MAT,
        *amrex::sundials::The_Sundials_Context());
      udata->colPtrs[i] = (int*)SUNSparseMatrix_IndexPointers((udata->PS)[i]);
      udata->rowVals[i] = (int*)SUNSparseMatrix_IndexValues((udata->PS)[i]);
      udata->Jdata[i] = SUNSparseMatrix_Data((udata->PS)[i]);
      // indx not used YET
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
        udata->rowVals[i], udata->colPtrs[i], udata->indx, &HP);
      udata->JSPSmat[i] =
        new amrex::Real[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
      klu_defaults(&(udata->Common[i]));
      // udata->Common.btf = 0;
      //(udata->Common[i]).maxwork = 15;
      // udata->Common.ordering = 1;
      udata->Symbolic[i] = klu_analyze(
        NUM_SPECIES + 1, udata->colPtrs[i], udata->rowVals[i],
        &(udata->Common[i]));
    }
#endif
  } else if (udata->precond_type == cvode::customSimpleAJac) {
    // CSR matrices data for each submatrix (cells)
    udata->colVals = new int*[udata->ncells];
    udata->rowPtrs = new int*[udata->ncells];
    // Matrices for each sparse custom block-solve
    udata->PS = new SUNMatrix[udata->ncells];
    udata->JSPSmat = new amrex::Real*[udata->ncells];
    // Number of non zero elements
    SPARSITY_INFO_SYST_SIMPLIFIED(&(udata->NNZ), &HP);
    for (int i = 0; i < udata->ncells; ++i) {
      (udata->PS)[i] = SUNSparseMatrix(
        NUM_SPECIES + 1, NUM_SPECIES + 1, udata->NNZ, CSR_MAT,
        *amrex::sundials::The_Sundials_Context());
      udata->rowPtrs[i] = (int*)SUNSparseMatrix_IndexPointers((udata->PS)[i]);
      udata->colVals[i] = (int*)SUNSparseMatrix_IndexValues((udata->PS)[i]);
      udata->Jdata[i] = SUNSparseMatrix_Data((udata->PS)[i]);
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
        udata->colVals[i], udata->rowPtrs[i], &HP, 0);
      udata->JSPSmat[i] =
        new amrex::Real[(NUM_SPECIES + 1) * (NUM_SPECIES + 1)];
    }
  }
#endif
}

int
ReactorCvode::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rYsrc_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorCvode::react()");
  // CPU and GPU version are very different such that the entire function
  // is split between a GPU region and a CPU region

  amrex::Real time_start = time;
  amrex::Real time_final = time + dt_react;
  amrex::Real CvodeActual_time_final = 0.0;

  //----------------------------------------------------------
  // GPU Region
  //----------------------------------------------------------
#ifdef AMREX_USE_GPU
  const int ncells = box.numPts();
  const int neq_tot = (NUM_SPECIES + 1) * ncells;

  // Solution vector and execution policy
#if defined(AMREX_USE_CUDA)
  N_Vector y = N_VNewWithMemHelp_Cuda(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  N_Vector y = N_VNewWithMemHelp_Hip(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  // SUNHipExecPolicy* reduce_exec_policy =
  // new SUNHipBlockReduceAtomicExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_DPCPP)
  N_Vector y = N_VNewWithMemHelp_Sycl(
    neq_tot, false, *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Sycl", 0))
    return (1);
  SUNSyclExecPolicy* stream_exec_policy =
    new SUNSyclThreadDirectExecPolicy(256);
  SUNSyclExecPolicy* reduce_exec_policy =
    new SUNSyclBlockReduceExecPolicy(256, 0);
  N_VSetKernelExecPolicy_Sycl(y, stream_exec_policy, reduce_exec_policy);
#endif

  // amrex::Print() << "using atomics? "  << reduce_exec_policy->atomic() <<
  // "\n";

  // Solution data array
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y);

  amrex::Gpu::streamSynchronize();
  SUNMatrix A = NULL;
  CVODEUserData* user_data = new CVODEUserData{};
  allocUserData(user_data, ncells, A, stream);

  // Fill data
  flatten(
    box, ncells, rY_in, rYsrc_in, T_in, rEner_in, rEner_src_in, yvec_d,
    user_data->rYsrc_ext, user_data->rhoe_init, user_data->rhoesrc_ext);

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  // Setup Cvode object
  void* cvode_mem =
    CVodeCreate(CV_BDF, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)cvode_mem, "CVodeCreate", 0))
    return (1);
  int flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

  // Call CVodeInit to initialize the integrator memory and specify the user's
  // right hand side function, the inital time, and initial dependent variable
  // vector y.
  flag = CVodeInit(cvode_mem, cF_RHS, time_start, y);
  if (utils::check_flag(&flag, "CVodeInit", 1))
    return (1);

  // Setup tolerances with typical values
  set_sundials_solver_tols(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, user_data->ncells,
    user_data->verbose, relTol, absTol, "cvode");

  // Linear solver data
  SUNLinearSolver LS = NULL;
  if (user_data->solve_type == cvode::sparseDirect) {
#if defined(AMREX_USE_CUDA)
    LS = SUNLinSol_cuSolverSp_batchQR(
      y, A, user_data->cusolverHandle,
      *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_cuSolverSp_batchQR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
#else
    amrex::Abort(
      "Shoudn't be there. solve_type sparse_direct only available with CUDA");
#endif
  } else if (user_data->solve_type == cvode::customDirect) {
#if defined(AMREX_USE_CUDA)
    LS = cvode::SUNLinSol_dense_custom(
      y, A, stream, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_dense_custom", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
    flag = CVodeSetJacFn(cvode_mem, cvode::cJac);
    if (utils::check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
#else
    amrex::Abort(
      "Shoudn't be there. solve_type custom_direct only available with CUDA");
#endif
  } else if (user_data->solve_type == cvode::magmaDirect) {
#ifdef PELE_USE_MAGMA
    LS = SUNLinSol_MagmaDense(y, A, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_MagmaDense", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
#else
    amrex::Abort(
      "Shoudn't be there. solve_type magma_direct only available with "
      "PELE_USE_MAGMA = TRUE");
#endif
  } else if (user_data->solve_type == cvode::GMRES) {
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_NONE, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);
  } else if (user_data->solve_type == cvode::precGMRES) {
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_LEFT, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);
  }

  // Analytical Jac. data for direct solver
  // Sparse/custom/magma direct uses the same Jacobian functions
  if (user_data->analytical_jacobian == 1) {
    flag = CVodeSetJacFn(cvode_mem, cvode::cJac);
    if (utils::check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
  }

  // Analytical Jac. data for iterative solver preconditioner
  if (user_data->precond_type == cvode::sparseSimpleAJac) {
    flag = CVodeSetPreconditioner(cvode_mem, cvode::Precond, cvode::PSolve);
    if (utils::check_flag(&flag, "CVodeSetPreconditioner", 1))
      return (1);
  }

  // CVODE runtime options
  flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (utils::check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    return (1);
  flag = CVodeSetMaxOrd(cvode_mem, user_data->maxOrder);
  if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1))
    return (1);

  // Actual CVODE solve
  BL_PROFILE_VAR("Pele::ReactorCvode::react():CVode", AroundCVODE);
  flag = CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
  if (utils::check_flag(&flag, "CVode", 1))
    return (1);
  BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
  dt_react =
    time_start - CvodeActual_time_final; // Actual dt_react performed by Cvode
  time += dt_react;                      // Increment time in reactor mode
#endif

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  // Get workload estimate
  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

  amrex::Gpu::DeviceVector<long int> v_nfe(ncells, nfe);
  long int* d_nfe = v_nfe.data();
  unflatten(
    box, ncells, rY_in, T_in, rEner_in, rEner_src_in, FC_in, yvec_d,
    user_data->rhoe_init, d_nfe, dt_react);

  if (user_data->verbose > 1) {
    print_final_stats(cvode_mem);
  }

  // Clean up
  N_VDestroy(y);
  CVodeFree(&cvode_mem);

  SUNLinSolFree(LS);
  if (A != nullptr) {
    SUNMatDestroy(A);
  }
  freeUserData(user_data);

#else
  //----------------------------------------------------------
  // CPU Region
  //----------------------------------------------------------

  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  // Update TypicalValues
  set_sundials_solver_tols(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, udata_g->ncells,
    udata_g->verbose, relTol, absTol, "cvode");

  // Perform integration one cell at a time
  const int icell = 0;
  const int ncells = 1;
  const auto captured_reactor_type = m_reactor_type;
  ParallelFor(
    box, [=, &CvodeActual_time_final] AMREX_GPU_DEVICE(
           int i, int j, int k) noexcept {
      if (mask(i, j, k) != -1) {

        amrex::Real* yvec_d = N_VGetArrayPointer(y);
        utils::box_flatten<Ordering>(
          icell, i, j, k, ncells, captured_reactor_type, rY_in, rYsrc_in, T_in,
          rEner_in, rEner_src_in, yvec_d, udata_g->rYsrc_ext,
          udata_g->rhoe_init, udata_g->rhoesrc_ext);

        // ReInit CVODE is faster
        CVodeReInit(cvode_mem, time_start, y);

        BL_PROFILE_VAR("Pele::ReactorCvode::react():CVode", AroundCVODE);
        CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
        BL_PROFILE_VAR_STOP(AroundCVODE);

        if ((udata_g->verbose > 1) && (omp_thread == 0)) {
          amrex::Print() << "Additional verbose info --\n";
          print_final_stats(cvode_mem);
          amrex::Print() << "\n -------------------------------------\n";
        }

        amrex::Real actual_dt = CvodeActual_time_final - time_start;

        // Get estimate of how hard the integration process was
        long int nfe = 0;
        long int nfeLS = 0;
        CVodeGetNumRhsEvals(cvode_mem, &nfe);
        CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
        const long int nfe_tot = nfe + nfeLS;

        utils::box_unflatten<Ordering>(
          icell, i, j, k, ncells, captured_reactor_type, rY_in, T_in, rEner_in,
          rEner_src_in, FC_in, yvec_d, udata_g->rhoe_init, nfe_tot, dt_react);

        if ((udata_g->verbose > 3) && (omp_thread == 0)) {
          amrex::Print() << "END : time curr is " << CvodeActual_time_final
                         << " and actual dt_react is " << actual_dt << "\n";
        }
      } else {
        FC_in(i, j, k, 0) = 0.0;
      }
    });

#ifdef MOD_REACTOR
  dt_react =
    time_start -
    time_final; // In this case, assumes all individual CVODE calls nailed it.
  time += dt_react; // Increment time in reactor mode
#endif

  long int nfe =
    20; // Dummy, the return value is no longer used for this function.
#endif

  return nfe;
}

int
ReactorCvode::react(
  realtype* rY_in,
  realtype* rYsrc_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorCvode::react()");

  std::cout << "Reacting (flattened)\n";

  // CPU and GPU version are very different such that the entire file
  // is split between a GPU region and a CPU region

  amrex::Real time_start = time;
  amrex::Real time_final = time + dt_react;
  amrex::Real CvodeActual_time_final = 0.0;

  //----------------------------------------------------------
  // GPU Region
  //----------------------------------------------------------
#ifdef AMREX_USE_GPU
  int neq_tot = (NUM_SPECIES + 1) * ncells;
  N_Vector y = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
  void* cvode_mem = NULL;

  // Fill user_data
  amrex::Gpu::streamSynchronize();
  CVODEUserData* user_data = new CVODEUserData{};
  allocUserData(user_data, ncells, A, stream);

  // Solution vector and execution policy
#if defined(AMREX_USE_CUDA)
  y = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0))
    return (1);
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy =
    new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  y = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0))
    return (1);
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy =
    new SUNHipBlockReduceExecPolicy(256, 0, stream);
  // SUNHipExecPolicy* reduce_exec_policy =
  // new SUNHipBlockReduceAtomicExecPolicy(256, 0, stream);
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_DPCPP)
  y = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(),
    *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)y, "N_VNewWithMemHelp_Sycl", 0))
    return (1);
  SUNSyclExecPolicy* stream_exec_policy =
    new SUNSyclThreadDirectExecPolicy(256);
  SUNSyclExecPolicy* reduce_exec_policy =
    new SUNSyclBlockReduceExecPolicy(256, 0);
  N_VSetKernelExecPolicy_Sycl(y, stream_exec_policy, reduce_exec_policy);
#endif

  // amrex::Print() << "using atomics? "  << reduce_exec_policy->atomic() <<
  // "\n";

  // Solution data array
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y);

  // Fill data
  BL_PROFILE_VAR("Pele::ReactorCvode::react():ASyncCopy", AsyncCopy);
  amrex::Gpu::htod_memcpy_async(yvec_d, rY_in, sizeof(amrex::Real) * (neq_tot));
  amrex::Gpu::htod_memcpy_async(
    user_data->rYsrc_ext, rYsrc_in, sizeof(amrex::Real) * NUM_SPECIES * ncells);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoe_init, rX_in, sizeof(amrex::Real) * ncells);
  amrex::Gpu::htod_memcpy_async(
    user_data->rhoesrc_ext, rX_src_in, sizeof(amrex::Real) * ncells);
  BL_PROFILE_VAR_STOP(AsyncCopy);

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  // Initialize integrator
  cvode_mem = CVodeCreate(CV_BDF, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag((void*)cvode_mem, "CVodeCreate", 0))
    return (1);
  int flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

  // Call CVodeInit to initialize the integrator memory and specify the
  //  user's right hand side function, the inital time, and
  //  initial dependent variable vector y.
  flag = CVodeInit(cvode_mem, cF_RHS, time_start, y);
  if (utils::check_flag(&flag, "CVodeInit", 1))
    return (1);

  // Setup tolerances with typical values
  set_sundials_solver_tols(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, user_data->ncells,
    user_data->verbose, relTol, absTol, "cvode");

  // Linear solver data
  if (user_data->solve_type == cvode::sparseDirect) {
#if defined(AMREX_USE_CUDA)
    LS = SUNLinSol_cuSolverSp_batchQR(
      y, A, user_data->cusolverHandle,
      *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_cuSolverSp_batchQR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

#else
    amrex::Abort(
      "Shoudn't be there. solve_type sparse_direct only available with CUDA");
#endif
  } else if (user_data->solve_type == cvode::customDirect) {
#if defined(AMREX_USE_CUDA)
    LS = cvode::SUNLinSol_dense_custom(
      y, A, stream, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_dense_custom", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);

    flag = CVodeSetJacFn(cvode_mem, cvode::cJac);
    if (utils::check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
#else
    amrex::Abort(
      "Shoudn't be there. solve_type custom_direct only available with CUDA");
#endif
  } else if (user_data->solve_type == cvode::magmaDirect) {
#ifdef PELE_USE_MAGMA
    LS = SUNLinSol_MagmaDense(y, A, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_MagmaDense", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
#else
    amrex::Abort(
      "Shoudn't be there. solve_type magma_direct only available with "
      "PELE_USE_MAGMA = TRUE");
#endif
  } else if (user_data->solve_type == cvode::GMRES) {
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_NONE, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);
  } else if (user_data->solve_type == cvode::precGMRES) {
    LS = SUNLinSol_SPGMR(
      y, SUN_PREC_LEFT, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag((void*)LS, "SUNLinSol_SPGMR", 0))
      return (1);
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1))
      return (1);
    flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1))
      return (1);
  }

  // Analytical Jac. data for direct solver
  // Both sparse/custom direct uses the same Jacobian functions
  if (user_data->analytical_jacobian == 1) {
    flag = CVodeSetJacFn(cvode_mem, cvode::cJac);
    if (utils::check_flag(&flag, "CVodeSetJacFn", 1))
      return (1);
  }

  // Analytical Jac. data for iterative solver preconditioner
  if (user_data->precond_type == cvode::sparseSimpleAJac) {
    flag = CVodeSetPreconditioner(cvode_mem, cvode::Precond, cvode::PSolve);
    if (utils::check_flag(&flag, "CVodeSetPreconditioner", 1))
      return (1);
  }

  // CVODE runtime options
  flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (utils::check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    return (1);
  flag = CVodeSetMaxOrd(cvode_mem, user_data->maxOrder);
  if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1))
    return (1);

  // Actual CVODE solve
  BL_PROFILE_VAR("Pele::ReactorCvode::react():CVode", AroundCVODE);
  flag = CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
  if (utils::check_flag(&flag, "CVode", 1))
    return (1);
  BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
  dt_react =
    time_start - CvodeActual_time_final; // Actual dt_react performed by Cvode
  time += dt_react;                      // Increment time in reactor mode
#endif

#ifdef AMREX_USE_OMP
  Gpu::Device::streamSynchronize();
#endif

  // Get the result back
  BL_PROFILE_VAR_START(AsyncCopy);
  amrex::Gpu::dtoh_memcpy_async(rY_in, yvec_d, sizeof(amrex::Real) * neq_tot);
  for (int i = 0; i < ncells; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }
  BL_PROFILE_VAR_STOP(AsyncCopy);

  // Get the number of RHS evaluations
  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if (user_data->verbose > 1) {
    print_final_stats(cvode_mem);
  }

  // Clean up
  N_VDestroy(y);
  CVodeFree(&cvode_mem);

  SUNLinSolFree(LS);
  if (A != nullptr) {
    SUNMatDestroy(A);
  }
  freeUserData(user_data);

  //----------------------------------------------------------
  // CPU Region
  //----------------------------------------------------------
#else
  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  // Pointer of solution vector
  amrex::Real* yvec_d = N_VGetArrayPointer(y);
  std::memcpy(
    yvec_d, rY_in, sizeof(amrex::Real) * ((NUM_SPECIES + 1) * ncells));
  std::memcpy(
    udata_g->rYsrc_ext, rYsrc_in, sizeof(amrex::Real) * (NUM_SPECIES * ncells));
  std::memcpy(udata_g->rhoe_init, rX_in, sizeof(amrex::Real) * ncells);
  std::memcpy(udata_g->rhoesrc_ext, rX_src_in, sizeof(amrex::Real) * ncells);

  // ReInit CVODE is faster
  CVodeReInit(cvode_mem, time_start, y);

  BL_PROFILE_VAR("Pele::react():CVode", AroundCVODE);
  int flag =
    CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
  // ONE STEP MODE FOR DEBUGGING
  // flag = CVode(cvode_mem, time_final, y, &CvodeActual_time_final,
  // CV_ONE_STEP);
  if (utils::check_flag(&flag, "CVode", 1)) {
    return (1);
  }
  BL_PROFILE_VAR_STOP(AroundCVODE);

  // Update TypicalValues
  set_sundials_solver_tols(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, udata_g->ncells,
    udata_g->verbose, relTol, absTol, "cvode");

#ifdef MOD_REACTOR
  dt_react =
    time_start - CvodeActual_time_final; // Actual dt_react performed by Cvode
  time += dt_react;                      // Increment time in reactor mode
#endif

  // Pack data to return in main routine external
  std::memcpy(
    rY_in, yvec_d, sizeof(amrex::Real) * ((NUM_SPECIES + 1) * ncells));
  for (int i = 0; i < ncells; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }

  if ((udata_g->verbose > 1) && (omp_thread == 0)) {
    amrex::Print() << "Additional verbose info --\n";
    print_final_stats(cvode_mem);
    amrex::Print() << "\n -------------------------------------\n";
  }

  // Get estimate of how hard the integration process was
  long int nfe, nfeLS;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  nfe += nfeLS;
#endif

  return nfe;
}

int
ReactorCvode::cF_RHS(
  realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::ReactorCvode::cF_RHS()");
#if defined(AMREX_USE_CUDA)
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer_Cuda(y_in);
  amrex::Real* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot_in);
#elif defined(AMREX_USE_HIP)
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer_Hip(y_in);
  amrex::Real* ydot_d = N_VGetDeviceArrayPointer_Hip(ydot_in);
#elif defined(AMREX_USE_DPCPP)
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer_Sycl(y_in);
  amrex::Real* ydot_d = N_VGetDeviceArrayPointer_Sycl(ydot_in);
#else
  amrex::Real* yvec_d = N_VGetArrayPointer(y_in);
  amrex::Real* ydot_d = N_VGetArrayPointer(ydot_in);
#endif

  auto* udata = static_cast<CVODEUserData*>(user_data);
  udata->dt_save = t;

  const auto ncells = udata->ncells;
  const auto dt_save = udata->dt_save;
  const auto reactor_type = udata->reactor_type;
  auto* rhoe_init = udata->rhoe_init;
  auto* rhoesrc_ext = udata->rhoesrc_ext;
  auto* rYsrc_ext = udata->rYsrc_ext;
  amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    utils::fKernelSpec<Ordering>(
      icell, ncells, dt_save, reactor_type, yvec_d, ydot_d, rhoe_init,
      rhoesrc_ext, rYsrc_ext);
  });
  amrex::Gpu::Device::streamSynchronize();
  return 0;
}

void
ReactorCvode::freeUserData(CVODEUserData* data_wk)
{
  amrex::The_Arena()->free(data_wk->rYsrc_ext);
  amrex::The_Arena()->free(data_wk->rhoe_init);
  amrex::The_Arena()->free(data_wk->rhoesrc_ext);
  amrex::The_Arena()->free(data_wk->mask);

#ifdef AMREX_USE_GPU

  if (data_wk->solve_type == cvode::sparseDirect) {
#ifdef AMREX_USE_CUDA
    amrex::The_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Arena()->free(data_wk->csr_col_index_h);
    cusolverStatus_t cusolver_status =
      cusolverSpDestroy(data_wk->cusolverHandle);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusparseStatus_t cusparse_status = cusparseDestroy(data_wk->cuSPHandle);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
#endif
  } else if (data_wk->solve_type == cvode::customDirect) {
#ifdef AMREX_USE_CUDA
    amrex::The_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Arena()->free(data_wk->csr_col_index_h);
    cusparseStatus_t cusparse_status = cusparseDestroy(data_wk->cuSPHandle);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
#endif
  }
  // Preconditioner analytical Jacobian data
  if (data_wk->precond_type == cvode::sparseSimpleAJac) {
#ifdef AMREX_USE_CUDA
    amrex::The_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Arena()->free(data_wk->csr_col_index_h);
    amrex::The_Arena()->free(data_wk->csr_val_d);
    amrex::The_Arena()->free(data_wk->csr_jac_d);
    cusolverStatus_t cusolver_status =
      cusolverSpDestroy(data_wk->cusolverHandle);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusolver_status = cusolverSpDestroyCsrqrInfo(data_wk->info);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cudaFree(data_wk->buffer_qr);
#endif
  }
  delete data_wk;

#else
  amrex::The_Arena()->free(data_wk->FCunt);

  // Direct solver Jac. data
  if (data_wk->solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
    delete[] data_wk->colPtrs;
    delete[] data_wk->rowVals;
    delete[] data_wk->Jdata;
    SUNMatDestroy(A);
    SUNMatDestroy((data_wk->PS)[0]);
    delete[](data_wk->PS);
#endif
  } else if (data_wk->solve_type == cvode::customDirect) {
    SUNMatDestroy(A);
    SUNMatDestroy(data_wk->PSc);
  }

  // Preconditionner Jac. data
  if (data_wk->precond_type == cvode::denseSimpleAJac) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      SUNDlsMat_destroyMat((data_wk->P)[i][i]);
      SUNDlsMat_destroyMat((data_wk->Jbd)[i][i]);
      SUNDlsMat_destroyArray((data_wk->pivot)[i][i]);
    }
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[](data_wk->P)[i];
      delete[](data_wk->Jbd)[i];
      delete[](data_wk->pivot)[i];
    }
    delete[](data_wk->P);
    delete[](data_wk->Jbd);
    delete[](data_wk->pivot);
  } else if (data_wk->precond_type == cvode::sparseSimpleAJac) {
#ifdef PELE_USE_KLU
    delete[] data_wk->colPtrs;
    delete[] data_wk->rowVals;
    delete[] data_wk->Jdata;
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
#endif
  } else if (data_wk->precond_type == cvode::customSimpleAJac) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[] data_wk->JSPSmat[i];
      SUNMatDestroy((data_wk->PS)[i]);
    }
    delete[] data_wk->colVals;
    delete[] data_wk->rowPtrs;
    delete[] data_wk->PS;
    delete[] data_wk->JSPSmat;
  }

  delete data_wk;
#endif
}

void
ReactorCvode::close()
{
#ifndef AMREX_USE_GPU
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  if (udata_g->solve_type == cvode::denseDirect) {
    SUNMatDestroy(A);
  }

  N_VDestroy(y);
  freeUserData(udata_g);
#endif
}

void
ReactorCvode::print_final_stats(void* cvodemem)
{
  long lenrw, leniw;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = CVodeGetWorkSpace(cvodemem, &lenrw, &leniw);
  utils::check_flag(&flag, "CVodeGetWorkSpace", 1);
  flag = CVodeGetNumSteps(cvodemem, &nst);
  utils::check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodemem, &nfe);
  utils::check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodemem, &nsetups);
  utils::check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodemem, &netf);
  utils::check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodemem, &nni);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodemem, &ncfn);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetLinWorkSpace(cvodemem, &lenrwLS, &leniwLS);
  utils::check_flag(&flag, "CVodeGetLinWorkSpace", 1);
  flag = CVodeGetNumLinIters(cvodemem, &nli);
  utils::check_flag(&flag, "CVodeGetNumLinIters", 1);
  // flag = CVodeGetNumJacEvals(cvodemem, &nje);
  // utils::check_flag(&flag, "CVodeGetNumJacEvals", 1);
  flag = CVodeGetNumLinRhsEvals(cvodemem, &nfeLS);
  utils::check_flag(&flag, "CVodeGetNumLinRhsEvals", 1);

  flag = CVodeGetNumPrecEvals(cvodemem, &npe);
  utils::check_flag(&flag, "CVodeGetNumPrecEvals", 1);
  flag = CVodeGetNumPrecSolves(cvodemem, &nps);
  utils::check_flag(&flag, "CVodeGetNumPrecSolves", 1);

  flag = CVodeGetNumLinConvFails(cvodemem, &ncfl);
  utils::check_flag(&flag, "CVodeGetNumLinConvFails", 1);

#ifdef AMREX_USE_OMP
  amrex::Print() << "\nFinal Statistics: "
                 << "(thread:" << omp_get_thread_num() << ", ";
  amrex::Print() << "cvode_mem:" << cvodemem << ")\n";
#else
  amrex::Print() << "\nFinal Statistics:\n";
#endif
  amrex::Print() << "lenrw      = " << lenrw << "    leniw         = " << leniw
                 << "\n";
  amrex::Print() << "lenrwLS    = " << lenrwLS
                 << "    leniwLS       = " << leniwLS << "\n";
  amrex::Print() << "nSteps     = " << nst << "\n";
  amrex::Print() << "nRHSeval   = " << nfe << "    nLinRHSeval   = " << nfeLS
                 << "\n";
  amrex::Print() << "nnLinIt    = " << nni << "    nLinIt        = " << nli
                 << "\n";
  amrex::Print() << "nLinsetups = " << nsetups << "    nErrtf        = " << netf
                 << "\n";
  amrex::Print() << "nPreceval  = " << npe << "    nPrecsolve    = " << nps
                 << "\n";
  amrex::Print() << "nConvfail  = " << ncfn << "    nLinConvfail  = " << ncfl
                 << "\n\n";
}

} // namespace reactions
} // namespace physics
} // namespace pele
