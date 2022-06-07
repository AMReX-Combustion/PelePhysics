#ifndef _NEURAL_NET_H_
#define _NEURAL_NET_H_

#include <torch/script.h>
#include "ManifoldFunc.H"

namespace pele {
namespace physics {

class NNFunc: public ManifoldFunc
{
public:

  NNFunc() {}

  NNFunc()
  {
    // Get filename and load
    amrex::ParmParse pp("nn");

    std::string nnfile;
    pp.get("filename", nnfile);

    try
    {
      // Deserialize the ScriptModule from a file using torch::jit::load().
      nnmodel = torch::jit::load(nnfile);
      // Convert to proper data type
      nnmodel.to(dtype);
      // Move to GPU if needed
#ifdef AMREX_USE_CUDA
      nnmodel.to(device);
#endif
    }
    catch(const c10::Error& e)
    {
      amrex::Error("Unable to load neural network model for manifold EOS.");
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_value(const int ivar, amrex::Real indata[], amrex::Real& out)
  {
    torch::Tensor intensor = torch::from_blob(indata, {Ndim}, tensoropt);
    get_value(ivar, intensor, out);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_derivs(const int ivar, amrex::Real indata[], amrex::Real derivs[])
  {
    torch::Tensor intensor = torch::from_blob(indata, {Ndim}, tensoropt);
    get_derivs(ivar, intensor, derivs);
  }

private:

  // Variables
  torch::jit::script::Module nnmodel;
  int Ndim;

  // Whether we are on an NVIDIA GPU or not
  // TODO: determine if we can support other GPUs -- need to edit constructor as well in that case
#ifdef AMREX_USE_CUDA
  static const torch::Device device = torch::kCUDA;
#else
  static const torch::Device device = torch::kCPU;
#endif

  // Set data type
#ifdef BL_USE_FLOAT
  static const torch::ScalarType dtype = torch::kFloat32;
#else
  static const torch::ScalarType dtype = torch::kFloat64;
#endif

  // Tensor creation options
  static const torch::TensorOpt tensoropt = torch::TensorOptions().dtype(dtype).device(device);

  // Private helper functions
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  torch::Tensor eval_model(const torch::Tensor& indata)
  {
    return nnmodel.forward({indata}).toTensor();
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  torch::Tensor eval_model(const int ivar, const torch::Tensor& indata)
  {
    return nnmodel.forward({indata}).toTensor()[ivar];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_value(const int ivar, const torch::Tensor& indata, amrex::Real& out)
  {
    out = eval_model(ivar, indata)[0];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_derivs(const int ivar, torch::Tensor& indata, amrex::Real derivs[])
  {
    torch::Tensor out_tensor = eval_model(ivar, indata);
    out_tensor.backward();

    auto grad = indata.grad();
    auto acc = grad.accessor<Real,1>();

    for(int i = 0; i < acc.size(0); i++)
    {
      derivs[i] = acc[i];
    }
  }

}; // class NNFunc
}  // namespace physics
}  // namespace pele

#endif
