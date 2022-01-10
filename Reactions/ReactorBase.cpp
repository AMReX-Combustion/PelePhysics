#include "ReactorBase.H"

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

#include <cvode/cvode.h>

#ifdef AMREX_USE_GPU
#include "AMReX_SUNMemory.H"
#endif

#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele {
namespace physics {
namespace reactions {
void
ReactorBase::set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals)
{
  const int size_ETV = ExtTypVals.size();
  amrex::Vector<std::string> kname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV; i++) {
    m_typ_vals[i] = ExtTypVals[i];
  }

  if (omp_thread == 0) {
    amrex::Print() << "Set the typVals in PelePhysics: \n  ";
    for (int i = 0; i < size_ETV - 1; i++) {
      amrex::Print() << kname[i] << ":" << m_typ_vals[i] << "  ";
    }
    amrex::Print() << "Temp:" << m_typ_vals[size_ETV - 1] << std::endl;
  }
}
void
ReactorBase::set_sundials_solver_tols(
  sundials::Context& sunctx,
  void* sundials_mem,
  const int ncells,
  const int verbose,
  const amrex::Real relTol,
  const amrex::Real absTol,
  const std::string solvername)
{
  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  const int neq_tot = (NUM_SPECIES + 1) * ncells;

#if defined(AMREX_USE_CUDA)
  N_Vector atol = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetHostArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_Vector atol = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetHostArrayPointer_Hip(atol);
#elif defined(AMREX_USE_DPCPP)
  N_Vector atol = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(), sunctx);
  amrex::Real* ratol = N_VGetHostArrayPointer_Sycl(atol);
#else
  N_Vector atol = N_VNew_Serial(neq_tot, sunctx);
  amrex::Real* ratol = N_VGetArrayPointer(atol);
#endif

  if (m_typ_vals[0] > 0.0) {
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances with TypVals rtol = " << relTol
                     << " atolfact = " << absTol << " in PelePhysics \n";
    }
    for (int i = 0; i < ncells; i++) {
      const int offset = i * (NUM_SPECIES + 1);
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        ratol[offset + k] = m_typ_vals[k] * absTol;
      }
    }
  } else {
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances rtol = " << relTol << " atol = " << absTol
                     << " in PelePhysics \n";
    }
    for (int i = 0; i < neq_tot; i++) {
      ratol[i] = absTol;
    }
  }

#if defined(AMREX_USE_CUDA)
  N_VCopyToDevice_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_VCopyToDevice_Hip(atol);
#elif defined(AMREX_USE_DPCPP)
  N_VCopyToDevice_Sycl(atol);
#endif

  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  int flag;
  if (solvername == "cvode") {
    flag = CVodeSVtolerances(sundials_mem, relTol, atol);
  } else if (solvername == "arkstep") {
    flag = ARKStepSVtolerances(sundials_mem, relTol, atol);
  } else if (solvername == "erkstep") {
    flag = ERKStepSVtolerances(sundials_mem, relTol, atol);
  } else {
    amrex::Abort("setSundialsSolverTols not implemented for this solver type");
  }
  if (utils::check_flag(&flag, "SVtolerances", 1)) {
    amrex::Abort("Problem in setSundialsSolverTols");
  }

  N_VDestroy(atol);
}

} // namespace reactions
} // namespace physics
} // namespace pele
