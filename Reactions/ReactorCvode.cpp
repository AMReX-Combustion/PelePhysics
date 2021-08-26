#include "ReactorCvode.H"

namespace pele {
namespace physics {
namespace reactions {
int
ReactorCvode::init(int reactor_type, int Ncells)
{
  BL_PROFILE("Pele::ReactorCvode::init()");
  return 0;
}

int
ReactorCvode::react(
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
  BL_PROFILE("Pele::ReactorCvode::react()");
  return 0;
}

int
ReactorCvode::react(
  realtype* rY_in,
  realtype* rY_src_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int reactor_type,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorCvode::react()");
  return 0;
}

int
ReactorCvode::cF_RHS(
  realtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::ReactorCvode::cF_RHS()");
  return 0;
}

void
ReactorCvode::SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals)
{
  int size_ETV = ExtTypVals.size();
  amrex::Vector<std::string> kname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV - 1; i++) {
    typVals[i] = ExtTypVals[i];
  }
  typVals[size_ETV - 1] = ExtTypVals[size_ETV - 1];

  if (omp_thread == 0) {
    amrex::Print() << "Set the typVals in PelePhysics: \n  ";
    for (int i = 0; i < size_ETV - 1; i++) {
      amrex::Print() << kname[i] << ":" << typVals[i] << "  ";
    }
    amrex::Print() << "Temp:" << typVals[size_ETV - 1] << " \n";
  }
}

void
ReactorCvode::SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol)
{
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif
  relTol = relative_tol;
  absTol = absolute_tol;
  if (omp_thread == 0) {
    amrex::Print() << "Set RTOL, ATOL = " << relTol << " " << absTol
                   << " in PelePhysics\n";
  }
}

void
ReactorCvode::close()
{
}

} // namespace reactions
} // namespace physics
} // namespace pele
