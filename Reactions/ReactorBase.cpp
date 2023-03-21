#include "ReactorBase.H"

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

#include <cvode/cvode.h>

#ifdef AMREX_USE_GPU
#include "AMReX_SUNMemory.H"
#endif

#ifdef AMREX_USE_SYCL
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele::physics::reactions {

void
ReactorBase::set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals)
{
  const int size_ETV = static_cast<int>(ExtTypVals.size());
  amrex::Vector<std::string> kname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(kname);
  int omp_thread = 0;

#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  for (int i = 0; i < size_ETV; i++) {
    m_typ_vals[i] = ExtTypVals[i];
  }

  // cppcheck-suppress knownConditionTrueFalse
  if ((omp_thread == 0) && (verbose > 0)) {
    amrex::Print() << "Set the typical values for PelePhysics ODE integration:"
                   << std::endl;
    for (int i = 0; i < size_ETV - 1; i++) {
      amrex::Print() << kname[i] << " : " << m_typ_vals[i] << std::endl;
    }
    amrex::Print() << "Temp : " << m_typ_vals[size_ETV - 1] << std::endl;
  }
#if defined(AMREX_USE_GPU)
  m_typ_vals_gpu =
    (amrex::Real*)amrex::The_Arena()->alloc(size_ETV * sizeof(amrex::Real));
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, &(m_typ_vals[0]), &(m_typ_vals[0]) + size_ETV,
    m_typ_vals_gpu);
#endif
}

void
ReactorBase::set_sundials_solver_tols(
  // cppcheck-suppress constParameter
  sundials::Context& sunctx,
  void* sundials_mem,
  const int ncells,
  const amrex::Real relTol,
  const amrex::Real absTol,
  const std::string& solvername)
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
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_Vector atol = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetHostArrayPointer_Hip(atol);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Hip(atol);
#elif defined(AMREX_USE_SYCL)
  N_Vector atol = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(), sunctx);
  amrex::Real* ratol = N_VGetHostArrayPointer_Sycl(atol);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Sycl(atol);
#else
  N_Vector atol = N_VNew_Serial(neq_tot, sunctx);
  amrex::Real* ratol = N_VGetArrayPointer(atol);
#endif

  if (m_typ_vals[0] > 0.0) {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances with TypVals rtol = " << relTol
                     << " atolfact = " << absTol << " in PelePhysics \n";
    }
#if defined(AMREX_USE_GPU)
    const int nbThreads = 256;
    const int nbBlocks = std::max(1, neq_tot / nbThreads);
    auto arr = m_typ_vals_gpu;
    amrex::launch_global<256>
      <<<nbBlocks, 256>>>([=] AMREX_GPU_DEVICE() noexcept {
        const int icell = blockDim.x * blockIdx.x + threadIdx.x;
        if (icell < neq_tot) {
          ratol[icell] = arr[icell / ncells] * absTol;
        }
      });
#else
    for (int k = 0; k < NUM_SPECIES + 1; k++) {
      amrex::Real T = m_typ_vals[k] * absTol;
      for (int i = 0; i < ncells; i++) {
        const int offset = k * ncells + i;
        ratol[offset] = T;
      }
    }
#endif
  } else {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances rtol = " << relTol << " atol = " << absTol
                     << " in PelePhysics \n";
    }
#if defined(AMREX_USE_GPU)
    const int nbThreads = 256;
    const int nbBlocks = std::max(1, neq_tot / nbThreads);
    amrex::launch_global<256>
      <<<nbBlocks, 256>>>([=] AMREX_GPU_DEVICE() noexcept {
        int icell = blockDim.x * blockIdx.x + threadIdx.x;
        if (icell < neq_tot) {
          ratol[icell] = absTol;
        }
      });
#else
    for (int i = 0; i < neq_tot; i++) {
      ratol[i] = absTol;
    }
#endif
  }

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
  if (utils::check_flag(&flag, "SVtolerances", 1) != 0) {
    amrex::Abort("Problem in setSundialsSolverTols");
  }

  N_VDestroy(atol);
}

} // namespace pele::physics::reactions
