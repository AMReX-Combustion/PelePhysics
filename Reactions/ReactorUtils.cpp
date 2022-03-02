#include "ReactorUtils.H"

namespace pele {
namespace physics {
namespace reactions {
namespace utils {

// Check function return value...
//     opt == 0 means SUNDIALS function allocates memory so check if
//              returned NULL pointer
//     opt == 1 means SUNDIALS function returns a flag so check if
//              flag >= 0
//     opt == 2 means function allocates memory so check if returned
//              NULL pointer
int
check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == nullptr) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  }
  if (opt == 1) {
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
  } else if (opt == 2 && flagvalue == nullptr) {
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

#ifdef AMREX_USE_GPU
N_Vector setNVectorGPU(int nvsize, int atomic_reductions, amrex::gpuStream_t stream) 
{
#if defined(AMREX_USE_CUDA)
  N_Vector y = N_VNewWithMemHelp_Cuda(
    nvsize, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Cuda", 0)) {
    amrex::Abort("Unable to create NVector Cuda");
  }
  SUNCudaExecPolicy* stream_exec_policy =
    new SUNCudaThreadDirectExecPolicy(256, stream);
  SUNCudaExecPolicy* reduce_exec_policy;
  if (atomic_reductions) {
    reduce_exec_policy = new SUNCudaBlockReduceAtomicExecPolicy(256, 0, stream);
  } else {
    reduce_exec_policy = new SUNCudaBlockReduceExecPolicy(256, 0, stream);
  }
  N_VSetKernelExecPolicy_Cuda(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_HIP)
  N_Vector y = N_VNewWithMemHelp_Hip(
    nvsize, false, *amrex::sundials::The_SUNMemory_Helper(),
    *amrex::sundials::The_Sundials_Context());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Hip", 0)) {
    amrex::Abort("Unable to create NVector Hip");
  }
  SUNHipExecPolicy* stream_exec_policy =
    new SUNHipThreadDirectExecPolicy(256, stream);
  SUNHipExecPolicy* reduce_exec_policy;
  if (atomic_reductions) {
    reduce_exec_policy = new SUNHipBlockReduceAtomicExecPolicy(256, 0, stream);
  } else {
    reduce_exec_policy = new SUNHipBlockReduceExecPolicy(256, 0, stream);
  }
  N_VSetKernelExecPolicy_Hip(y, stream_exec_policy, reduce_exec_policy);
#elif defined(AMREX_USE_DPCPP)
  N_Vector y = N_VNewWithMemHelp_Sycl(
    nvsize, false, *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(),
    *amrex::sundials::The_Sundials_Context());
  if (check_flag((void*)y, "N_VNewWithMemHelp_Sycl", 0)) {
    amrex::Abort("Unable to create NVector Sycl");
  }
  SUNSyclExecPolicy* stream_exec_policy =
    new SUNSyclThreadDirectExecPolicy(256);
  SUNSyclExecPolicy* reduce_exec_policy =
    new SUNSyclBlockReduceExecPolicy(256, 0);
  N_VSetKernelExecPolicy_Sycl(y, stream_exec_policy, reduce_exec_policy);
#endif
  return y;

  delete stream_exec_policy;
  delete reduce_exec_policy;
}
#endif
} // namespace utils
} // namespace reactions
} // namespace physics
} // namespace pele
