#ifndef BLACK_BOX_FUNC_FACTORY_H
#define BLACK_BOX_FUNC_FACTORY_H

#include "Table.H"
#include "NeuralNetHomerolled.H"

namespace pele {
namespace physics {

// Factory struct to allow either compile time or runtime
// selection of BlackBoxFunc type
template <typename FuncType>
struct BlackBoxFuncFactory
{
};

template <unsigned int Dimension>
struct BlackBoxFuncFactory<TabFunc<Dimension>>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory(const ManFuncData* mf_data)
    : func(static_cast<const TabFuncData*>(mf_data))
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      mf_data->manmodel == ManifoldModel::TABLE,
      "Runtime Table/Network must match what you compiled with");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabFunc<Dimension>* get_func() { return &func; }

private:
  TabFunc<Dimension> func;
};

template <>
struct BlackBoxFuncFactory<NNFunc>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory(const ManFuncData* mf_data)
    : func(static_cast<const NNFuncData*>(mf_data))
  {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      mf_data->manmodel == ManifoldModel::NEURAL_NET,
      "Runtime Table/Network must match what you compiled with");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc* get_func() { return &func; }

private:
  NNFunc func;
};

template <>
struct BlackBoxFuncFactory<ManifoldFunc>
{

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  BlackBoxFuncFactory(const ManFuncData* mf_data)
  {
    if (mf_data->manmodel == ManifoldModel::TABLE) {
      const TabFuncData* tf_data = static_cast<const TabFuncData*>(mf_data);
      // Note: Table dimension not specified a priori, slight performance loss
      func = new pele::physics::TabFunc<>(tf_data);
    } else if (mf_data->manmodel == ManifoldModel::NEURAL_NET) {
      const NNFuncData* nnf_data = static_cast<const NNFuncData*>(mf_data);
      func = new pele::physics::NNFunc(nnf_data);
    } else {
      amrex::Abort("invalid black box function type requested");
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  ~BlackBoxFuncFactory() { delete func; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  ManifoldFunc* get_func() { return func; }

private:
  ManifoldFunc* func;
};
} // namespace physics
} // namespace pele
#endif
