#include "ReactorCvodeUtils.H"

namespace pele {
namespace physics {
namespace reactions {
namespace cvode {

// Error function for CVODE
void
cvodeErrHandler(
  int error_code,
  const char* /*module*/,
  const char* /*function*/,
  char* msg,
  void* /*eh_data*/)
{
  if (error_code != CV_WARNING) {
    std::cout << "From CVODE: " << msg << std::endl;
    amrex::Abort("Aborting from CVODE");
  }
}

void
printFinalStats(void* cvodeMem)
{
  long lenrw, leniw;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = CVodeGetWorkSpace(cvodeMem, &lenrw, &leniw);
  utils::check_flag(&flag, "CVodeGetWorkSpace", 1);
  flag = CVodeGetNumSteps(cvodeMem, &nst);
  utils::check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  utils::check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  utils::check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netf);
  utils::check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetLinWorkSpace(cvodeMem, &lenrwLS, &leniwLS);
  utils::check_flag(&flag, "CVodeGetLinWorkSpace", 1);
  flag = CVodeGetNumLinIters(cvodeMem, &nli);
  utils::check_flag(&flag, "CVodeGetNumLinIters", 1);
  // flag = CVodeGetNumJacEvals(cvodeMem, &nje);
  // utils::check_flag(&flag, "CVodeGetNumJacEvals", 1);
  flag = CVodeGetNumLinRhsEvals(cvodeMem, &nfeLS);
  utils::check_flag(&flag, "CVodeGetNumLinRhsEvals", 1);

  flag = CVodeGetNumPrecEvals(cvodeMem, &npe);
  utils::check_flag(&flag, "CVodeGetNumPrecEvals", 1);
  flag = CVodeGetNumPrecSolves(cvodeMem, &nps);
  utils::check_flag(&flag, "CVodeGetNumPrecSolves", 1);

  flag = CVodeGetNumLinConvFails(cvodeMem, &ncfl);
  utils::check_flag(&flag, "CVodeGetNumLinConvFails", 1);

#ifdef AMREX_USE_OMP
  amrex::Print() << "\nFinal Statistics: "
                 << "(thread:" << omp_get_thread_num() << ", ";
  amrex::Print() << "cvodeMem:" << cvodeMem << ")\n";
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
} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
