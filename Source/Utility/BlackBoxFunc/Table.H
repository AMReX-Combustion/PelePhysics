#ifndef _TABLE_H_
#define _TABLE_H_

#define MAXD_TABLE 5
#include <iostream>
#include <fstream>
#include "ManifoldFunc.H"

namespace pele {
namespace physics {

struct TabFuncData : ManFuncData
{
  int varSpacing;
  int* dimLengths;
  int* dimDataSpacing;
  amrex::Real* grids;
  amrex::Real* values;
};

class TabFuncParams
  : public PeleParamsGeneric<ManFuncData>::
      Register<TabFuncParams, PeleParams<TabFuncData, ManFuncData>>
{
public:
  static std::string identifier() { return "Table"; }
};

template <>
struct InitParm<TabFuncData, ManFuncData>
{
  static void host_initialize(PeleParams<TabFuncData, ManFuncData>* parm_in)
  {
    std::string pp_pref =
      parm_in->m_host_only_parm.parm_parse_prefix.empty()
        ? "table"
        : parm_in->m_host_only_parm.parm_parse_prefix + ".table";
    amrex::ParmParse pp(pp_pref);
    std::string tablefile;
    int verbose = 2;

    // TODO: allow user-specific parmparse prefix to enable loading multiple
    // networks
    pp.get("filename", tablefile);
    pp.query("v", verbose);
    if (verbose > 0) {
      amrex::Print() << "Loading tabulated data from file: " << tablefile
                     << std::endl;
    }
    parm_in->m_h_parm.manmodel = ManifoldModel::TABLE;
    read_table(tablefile, parm_in->m_h_parm);
    if (verbose > 1) {
      print(parm_in->m_h_parm);
    }
  }

  static void host_deallocate(PeleParams<TabFuncData, ManFuncData>* parm_in)
  {
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.dimLengths);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.dimDataSpacing);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.grids);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.values);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.varnames);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.dimnames);
  }

  static void read_table(const std::string& tablefile, TabFuncData& m_h_tf_data)
  {
    std::ifstream fi(tablefile, std::ios::binary | std::ios::in);
    if (!fi.is_open()) {
      // AMREX ABORT
      std::cout << ("Unable to open input file " + tablefile) << std::endl;
    }
    fi.seekg(0, std::ios::beg);

    // Number of dimensions
    fi.read(reinterpret_cast<char*>(&m_h_tf_data.Ndim), sizeof(int));

    // Dimension names
    m_h_tf_data.dimnames = static_cast<char*>(amrex::The_Pinned_Arena()->alloc(
      m_h_tf_data.Ndim * m_h_tf_data.len_str * sizeof(char)));
    fi.read(
      m_h_tf_data.dimnames,
      m_h_tf_data.Ndim * m_h_tf_data.len_str * sizeof(char));

    // Dimension lengths
    m_h_tf_data.dimLengths = static_cast<int*>(
      amrex::The_Pinned_Arena()->alloc(m_h_tf_data.Ndim * sizeof(int)));
    fi.read(
      reinterpret_cast<char*>(m_h_tf_data.dimLengths),
      m_h_tf_data.Ndim * sizeof(int));
    // Number of variables
    fi.read(reinterpret_cast<char*>(&m_h_tf_data.Nvar), sizeof(int));
    // No separate passed variables for a table
    m_h_tf_data.Nmanpar = m_h_tf_data.Ndim;

    // Grid locations for all dimensions
    int Ngridtot = 0;
    for (int ii = 0; ii < m_h_tf_data.Ndim; ii++) {
      Ngridtot += m_h_tf_data.dimLengths[ii];
    }
    m_h_tf_data.grids = static_cast<amrex::Real*>(
      amrex::The_Pinned_Arena()->alloc(Ngridtot * sizeof(amrex::Real)));
    fi.read(
      reinterpret_cast<char*>(m_h_tf_data.grids),
      sizeof(amrex::Real) * Ngridtot);

    // Model Name
    fi.read(m_h_tf_data.model_name, sizeof(char) * m_h_tf_data.len_str);
    m_h_tf_data.model_name[m_h_tf_data.len_str] = '\0';

    // Variable Names
    m_h_tf_data.varnames = static_cast<char*>(amrex::The_Pinned_Arena()->alloc(
      m_h_tf_data.Nvar * m_h_tf_data.len_str * sizeof(char)));
    fi.read(
      m_h_tf_data.varnames,
      m_h_tf_data.Nvar * m_h_tf_data.len_str * sizeof(char));

    // Data
    m_h_tf_data.dimDataSpacing = static_cast<int*>(
      amrex::The_Pinned_Arena()->alloc(m_h_tf_data.Ndim * sizeof(int)));
    m_h_tf_data.dimDataSpacing[0] = 1;
    for (int ii = 1; ii < m_h_tf_data.Ndim; ii++) {
      m_h_tf_data.dimDataSpacing[ii] =
        m_h_tf_data.dimDataSpacing[ii - 1] * m_h_tf_data.dimLengths[ii - 1];
    }
    m_h_tf_data.varSpacing = m_h_tf_data.dimDataSpacing[m_h_tf_data.Ndim - 1] *
                             m_h_tf_data.dimLengths[m_h_tf_data.Ndim - 1];
    int Ndata = m_h_tf_data.varSpacing * m_h_tf_data.Nvar;
    m_h_tf_data.values = static_cast<amrex::Real*>(
      amrex::The_Pinned_Arena()->alloc(Ndata * sizeof(amrex::Real)));
    fi.read(
      reinterpret_cast<char*>(m_h_tf_data.values), sizeof(amrex::Real) * Ndata);
  }

  static void print(const TabFuncData& m_h_tf_data)
  {
    amrex::Print() << "TABULATED FUNCTION" << std::endl;
    amrex::Print() << "Ndim: " << m_h_tf_data.Ndim << std::endl;
    amrex::Print() << "Nxi : ";
    for (int ii = 0; ii < m_h_tf_data.Ndim; ii++) {
      std::string dimname(
        &m_h_tf_data.dimnames[ii * m_h_tf_data.len_str], m_h_tf_data.len_str);
      amrex::Print() << amrex::trim(dimname) << " "
                     << m_h_tf_data.dimLengths[ii] << " | ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "Nv: " << m_h_tf_data.Nvar << std::endl;
    amrex::Print() << "Grids: " << std::endl;
    int start = 0;
    for (int ii = 0; ii < m_h_tf_data.Ndim; ii++) {
      for (int jj = start; jj < start + m_h_tf_data.dimLengths[ii]; jj++) {
        amrex::Print() << m_h_tf_data.grids[jj] << " ";
      }
      start += m_h_tf_data.dimLengths[ii];
      amrex::Print() << std::endl << std::endl;
    }
    amrex::Print() << "Model name: " << m_h_tf_data.model_name << std::endl;

    amrex::Print() << "Index | Variable | Min | Max" << std::endl;
    for (int ii = 0; ii < m_h_tf_data.Nvar; ii++) {
      std::string varname(
        &m_h_tf_data.varnames[ii * m_h_tf_data.len_str], m_h_tf_data.len_str);
      amrex::Real* dat_beg = &m_h_tf_data.values[ii * m_h_tf_data.varSpacing];
      amrex::Real* dat_end =
        &m_h_tf_data.values[(ii + 1) * m_h_tf_data.varSpacing];
      amrex::Real* min_entry = std::min_element(dat_beg, dat_end);
      amrex::Real* max_entry = std::max_element(dat_beg, dat_end);
      amrex::Print() << ii << " | " << amrex::trim(varname) << " | "
                     << *min_entry << " | " << *max_entry << std::endl;
    }
  }
};

// Wrapper for interpolation process, chooses appropriate approach based on
// dimensionality. Table dimensionlity specified at compile time through
// template is done recursively and entirely inlined, with dimensionality
// automatically detected. Runtime table dimensionality must use select case
// to find the dimensionality, allows up to 5 dimensions
template <unsigned int InterpDim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
interpolate(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real data[],
  const int data_spacing[],
  const int numdim)
{
  AMREX_ASSERT(numdim == InterpDim);
  const int idx = indices[InterpDim - 1];
  const int dimDataSpac = data_spacing[InterpDim - 1];
  const amrex::Real out_value_left = interpolate<InterpDim - 1>(
    indices, alphas, &data[idx * dimDataSpac], data_spacing, numdim - 1);
  const amrex::Real out_value_right = interpolate<InterpDim - 1>(
    indices, alphas, &data[(idx + 1) * dimDataSpac], data_spacing, numdim - 1);
  const amrex::Real alpha = alphas[InterpDim - 1];
  return alpha * out_value_left + (1.0 - alpha) * out_value_right;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
interpolate<1>(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real data[],
  const int* /*data_spacing[]*/,
  const int numdim)
{
  AMREX_ASSERT(numdim == 1);
  const int idx = indices[0];
  const amrex::Real out_value_left = data[idx];
  const amrex::Real out_value_right = data[idx + 1];
  const amrex::Real alpha = alphas[0];
  return alpha * out_value_left + (1.0 - alpha) * out_value_right;
}

// Default: runtime selection of interpolation dimension
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
interpolate<0>(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real data[],
  const int data_spacing[],
  const int numdim)
{
  amrex::Real out_value;
  AMREX_ASSERT(numdim <= MAXD_TABLE);
  switch (numdim) {
  case 1:
    out_value = interpolate<1>(indices, alphas, data, data_spacing, numdim);
    break;
  case 2:
    out_value = interpolate<2>(indices, alphas, data, data_spacing, numdim);
    break;
  case 3:
    out_value = interpolate<3>(indices, alphas, data, data_spacing, numdim);
    break;
  case 4:
    out_value = interpolate<4>(indices, alphas, data, data_spacing, numdim);
    break;
  case 5:
    out_value = interpolate<5>(indices, alphas, data, data_spacing, numdim);
    break;
  default:
    amrex::Abort("Tabulated function: Failure: only 1D-5D tables are allowed "
                 "at runtime \n");
  }
  return out_value;
}

// Wrapper for differentiation process, chooses appropriate approach based on
// dimensionality same caveats apply as for the interpolation process
template <unsigned int DiffDim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
differentiate(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real dxinv[],
  const amrex::Real data[],
  const int data_spacing[],
  const int numdim,
  amrex::Real derivs[])
{
  constexpr int dim = DiffDim - 1;
  AMREX_ASSERT(DiffDim == numdim);
  const amrex::Real alpha = alphas[dim];
  const amrex::Real oneMinusAlpha = 1.0 - alpha;
  amrex::Real derivs_left[dim];
  amrex::Real derivs_right[dim];
  const amrex::Real value_left = differentiate<dim>(
    indices, alphas, dxinv, &data[indices[dim] * data_spacing[dim]],
    data_spacing, dim, derivs_left);
  const amrex::Real value_right = differentiate<dim>(
    indices, alphas, dxinv, &data[(indices[dim] + 1) * data_spacing[dim]],
    data_spacing, dim, derivs_right);
  // interpolate existing derivatives
  for (int ddim = 0; ddim < dim; ddim++) {
    derivs[ddim] =
      alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
  }
  derivs[dim] = (value_right - value_left) * dxinv[dim]; // take the derivative
  return alpha * value_left + oneMinusAlpha * value_right; // interpolate
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
differentiate<1>(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real dxinv[],
  const amrex::Real data[],
  const int* /*data_spacing[]*/,
  const int numdim,
  amrex::Real derivs[])
{
  AMREX_ASSERT(numdim == 1);
  const amrex::Real alpha = alphas[0];
  const amrex::Real oneMinusAlpha = 1.0 - alpha;
  const int idx = indices[0];
  const amrex::Real value_left = data[idx];
  const amrex::Real value_right = data[idx + 1];
  derivs[0] = (value_right - value_left) * dxinv[0]; // take the derivative
  return alpha * value_left + oneMinusAlpha * value_right; // interpolate
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
differentiate<0>(
  const int indices[],
  const amrex::Real alphas[],
  const amrex::Real dxinv[],
  const amrex::Real data[],
  const int data_spacing[],
  const int numdim,
  amrex::Real derivs[])
{
  AMREX_ASSERT(numdim <= MAXD_TABLE);
  amrex::Real out_value;
  switch (numdim) {
  case 1:
    out_value = differentiate<1>(
      indices, alphas, dxinv, data, data_spacing, numdim, derivs);
    break;
  case 2:
    out_value = differentiate<2>(
      indices, alphas, dxinv, data, data_spacing, numdim, derivs);
    break;
  case 3:
    out_value = differentiate<3>(
      indices, alphas, dxinv, data, data_spacing, numdim, derivs);
    break;
  case 4:
    out_value = differentiate<4>(
      indices, alphas, dxinv, data, data_spacing, numdim, derivs);
    break;
  case 5:
    out_value = differentiate<5>(
      indices, alphas, dxinv, data, data_spacing, numdim, derivs);
    break;
  default:
    amrex::Abort("Tabulated function: Failure: only 1D-5D tables are allowed "
                 "at runtime \n");
  }
  return out_value;
}

// Struct to contain tabulated data
// slight performance benefit possible use template to specify dimensionality
// default (0) allows arbitrary dimensionality
template <unsigned int TableDim = 0>
class TabFunc : public ManifoldFunc
{
public:
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabFunc() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabFunc(const TabFuncData* tf_data_in)
  {
    tf_data = tf_data_in;
    AMREX_ALWAYS_ASSERT(
      TableDim == tf_data->Ndim ||
      (TableDim == 0 && tf_data->Ndim <= MAXD_TABLE));
  }

  // Lookup directly if you already know index
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void
  get_value(const int ivar, const amrex::Real interpdata[], amrex::Real& out)
  {
    // Get alphas and indices
    int indices[maxdim];
    amrex::Real alphas[maxdim];
    amrex::Real dxinv[maxdim];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    out = interpolate<TableDim>(
      indices, alphas, &tf_data->values[ivar * tf_data->varSpacing],
      tf_data->dimDataSpacing, tf_data->Ndim);
  }

  // Lookup a few variables directly if you already know indices
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_values(
    const int nvar,
    const int ivar[],
    const amrex::Real interpdata[],
    amrex::Real out[])
  {
    // Get alphas and indices
    int indices[maxdim];
    amrex::Real alphas[maxdim];
    amrex::Real dxinv[maxdim];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    for (int i = 0; i < nvar; ++i) {
      out[i] =
        (ivar[i] >= 0)
          ? interpolate<TableDim>(
              indices, alphas, &tf_data->values[ivar[i] * tf_data->varSpacing],
              tf_data->dimDataSpacing, tf_data->Ndim)
          : 0.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_all_values(const amrex::Real interpdata[], amrex::Real out[])
  {
    // Get alphas and indices
    int indices[maxdim];
    amrex::Real alphas[maxdim];
    amrex::Real dxinv[maxdim];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    for (int i = 0; i < tf_data->Nvar; ++i) {
      out[i] = interpolate<TableDim>(
        indices, alphas, &tf_data->values[i * tf_data->varSpacing],
        tf_data->dimDataSpacing, tf_data->Ndim);
    }
  }

  // Get derivatives of indexed variable with respect to table dimensions using
  // finite difference
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void
  get_derivs(int ivar, const amrex::Real derivloc[], amrex::Real derivs[])
  {
    // Get alphas and indices
    int indices[maxdim];
    amrex::Real alphas[maxdim];
    amrex::Real dxinv[maxdim];
    get_indices_alphas_dxinv(derivloc, indices, alphas, dxinv);

    // finite differences from table
    amrex::Real out = differentiate<TableDim>(
      indices, alphas, dxinv, &tf_data->values[ivar * tf_data->varSpacing],
      tf_data->dimDataSpacing, tf_data->Ndim, derivs);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual ManifoldModel model() { return ManifoldModel::TABLE; }

private:
  const TabFuncData* tf_data;
  static constexpr int maxdim = TableDim != 0 ? TableDim : MAXD_TABLE;

  // -----------------------------------------------------------
  // Search for the closest index in an array to a given value
  // using the bisection technique.
  // INPUTS/OUTPUTS:
  // xtable(0:n-1) => array to search in (ascending order)
  // n             => array size
  // x             => x location
  // idxlo (return)=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
  // -----------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  int locate(const amrex::Real* xtable, const int n, const amrex::Real x)
  {
    int idxlo = 0;
    // If x is out of bounds, return boundary index
    if (x >= xtable[n - 2]) {
      idxlo = n - 2;
      return idxlo;
    }
    if (x <= xtable[0]) {
      idxlo = 0;
      return idxlo;
    }

    // Do the bisection
    int idxhi = n - 2;
    bool notdone = true;
    while (notdone) {
      if (idxhi - idxlo <= 1) {
        notdone = false;
      } else {
        const int idxmid = (idxhi + idxlo) / 2;
        if (x >= xtable[idxmid]) {
          idxlo = idxmid;
        } else {
          idxhi = idxmid;
        }
      }
    }
    return idxlo;
  }

  // Get quantities for use in interpolation
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_indices_alphas_dxinv(
    const amrex::Real interpdata[],
    int indices[],
    amrex::Real alphas[],
    amrex::Real dxinv[])
  {
    int start = 0;
    for (int dim = 0; dim < tf_data->Ndim; dim++) {
      int dimLen = tf_data->dimLengths[dim];
      // Find index to left of data
      // TODO: test alternatives: lower_bound, find, home rolled while loop,
      // home rolled fixed grid direct SOFAR: lower_bound is much faster than
      // home rolled while loop. Even then, finding indexes is ~50% of runtime

#ifndef AMREX_USE_GPU
      auto beg = &tf_data->grids[start];
      int idx =
        std::lower_bound(beg + 1, beg + dimLen - 1, interpdata[dim]) - beg - 1;
      idx = std::min(
        std::max(idx, 0),
        dimLen - 2); // don't let idx be negative or out of bounds
#else
      // home rolled
      // int idx = 1;
      // while (tf_data->grids[start+idx] < interpdata[dim] && idx < dimLen)
      //    {idx+=1;};
      // idx-=1;
      int idx = locate(&tf_data->grids[start], dimLen, interpdata[dim]);
#endif

      indices[dim] = idx;
      idx += start;
      amrex::Real dx_inv =
        1.0 / (tf_data->grids[idx + 1] - tf_data->grids[idx]);
      dxinv[dim] = dx_inv;
      amrex::Real alpha = (tf_data->grids[idx + 1] - interpdata[dim]) * dx_inv;
      alphas[dim] = std::min(std::max(alpha, 0.0), 1.0);
      start += dimLen;
    }
  }
};
} // namespace physics
} // namespace pele
#endif
