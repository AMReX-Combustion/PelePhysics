#ifndef BLACK_BOX_FUNC_FACTORY_H
#define BLACK_BOX_FUNC_FACTORY_H

#include "Table.H"
#include "NeuralNetHomerolled.H"

namespace pele::physics {

// Factory struct to allow either compile time or runtime
// selection of BlackBoxFunction type
template <typename FuncType>
struct BlackBoxFunctionFactory
{
};

template <unsigned int Dimension>
struct BlackBoxFunctionFactory<TabulatedFunction<Dimension>>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory() = default;

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory(const BlackBoxFunctionData* mf_data)
    : func(static_cast<const TabulatedFunctionData*>(mf_data))
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      mf_data->bbmodel == BlackBoxModel::TABLE,
      "Runtime Table/Network must match what you compiled with");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabulatedFunction<Dimension>* get_func() { return &func; }

private:
  TabulatedFunction<Dimension> func;
};

template <>
struct BlackBoxFunctionFactory<NeuralNetFunction>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory() = default;

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory(const BlackBoxFunctionData* mf_data)
    : func(static_cast<const NeuralNetFunctionData*>(mf_data))
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      mf_data->bbmodel == BlackBoxModel::NEURAL_NET,
      "Runtime Table/Network must match what you compiled with");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NeuralNetFunction* get_func() { return &func; }

private:
  NeuralNetFunction func;
};

template <>
struct BlackBoxFunctionFactory<BlackBoxFunction>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory() = default;

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunctionFactory(const BlackBoxFunctionData* mf_data)
  {
    if (mf_data->bbmodel == BlackBoxModel::TABLE) {
      const auto* tf_data = static_cast<const TabulatedFunctionData*>(mf_data);
      // Note: Table dimension not specified a priori, slight performance loss
      func = new pele::physics::TabulatedFunction<>(tf_data);
    } else if (mf_data->bbmodel == BlackBoxModel::NEURAL_NET) {
      const auto* nnf_data = static_cast<const NeuralNetFunctionData*>(mf_data);
      func = new pele::physics::NeuralNetFunction(nnf_data);
    } else {
      amrex::Abort("invalid black box function type requested");
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  ~BlackBoxFunctionFactory() { delete func; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFunction* get_func() { return func; }

private:
  BlackBoxFunction* func;
};
} // namespace pele::physics
#endif
