#include "reactor.H"

amrex::Array<double, NUM_SPECIES + 1> typVals={-1};
int eint_rho=1; // in/out = rhoE/rhoY
int enth_rho=2; // in/out = rhoH/rhoY
amrex::Real relTol=1e-6;
amrex::Real absTol=1e-10;
amrex::Real time_init=0.0;
int dense_solve = 1;
int sparse_solve = 5;
int iterative_gmres_solve = 99;
int sparse_solve_custom = 101;
int iterative_gmres_solve_custom = 199;
int hack_dump_sparsity_pattern = -5;
int sparse_custom_solve = 2;
int use_erkstep = 0;
int rk_method = 40;
int rk_controller = 0;

//===================================================================
void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals)
{
  int size_ETV = ExtTypVals.size();
  amrex::Vector<std::string> kname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV - 1; i++) 
  {
    typVals[i] = ExtTypVals[i];
  }
  typVals[size_ETV - 1] = ExtTypVals[size_ETV - 1];

  if (omp_thread == 0) 
  {
    amrex::Print() << "Set the typVals in PelePhysics: \n  ";
    for (int i = 0; i < size_ETV - 1; i++) 
    {
      amrex::Print() << kname[i] << ":" << typVals[i] << "  ";
    }
    amrex::Print() << "Temp:" << typVals[size_ETV - 1] << " \n";
  }

}
//===================================================================
void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol)
{
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif
  relTol = relative_tol;
  absTol = absolute_tol;
  if(omp_thread == 0){
  amrex::Print() << "Set RTOL, ATOL = " << relTol << " " << absTol
          << " in PelePhysics\n";
  }

}
//===================================================================
// Check function return value...
//     opt == 0 means SUNDIALS function allocates memory so check if
//              returned NULL pointer
//     opt == 1 means SUNDIALS function returns a flag so check if
//              flag >= 0
//     opt == 2 means function allocates memory so check if returned
//              NULL pointer
int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == NULL) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  } else if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        fprintf(
          stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
          *errflag);
        amrex::Abort("abort");
      }
      return (1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  }

  return (0);
}
//===================================================================
