#ifndef _TABLE_H_
#define _TABLE_H_

#define MAXD_TABLE 5
#include <iostream>
#include <fstream>
#include "ManifoldFunc.H"

namespace pele {
namespace physics {

// Struct to contain tabulated data
class TabFuncParams : public ManFuncParams
{
public:
  TabFuncParams() {}
  virtual ~TabFuncParams() {}

  virtual void initialize()
  {
    m_h_tf_data.manmodel = ManifoldModel::TABLE;

    amrex::ParmParse pp("manifold");
    std::string tablefile;
    pp.get("filename", tablefile);
    pp.query("v", m_verbose);
    read_table(tablefile);
    if (m_verbose >= 1) {
      print();
    }
    allocate();
  }

  void read_table(std::string& tablefile)
  {
    std::ifstream fi(tablefile, std::ios::binary | std::ios::in);
    if (!fi.is_open()) {
      // AMREX ABORT
      std::cout << ("Unable to open input file " + tablefile) << std::endl;
    }
    fi.seekg(0, std::ios::beg);

    // Number of dimensions
    fi.read(reinterpret_cast<char*>(&m_h_tf_data.Ndim), sizeof(int));
    if (m_h_tf_data.Ndim > MAXD_TABLE) {
      amrex::Abort(
        "Requesting too many table dimensions, need to increase MAXD_TABLE");
    }

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

  void print()
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

  struct TabFuncData : ManFuncParams::ManFuncData
  {
    int varSpacing;
    int* dimLengths;
    int* dimDataSpacing;
    amrex::Real* grids;
    amrex::Real* values;
  };

  virtual void allocate()
  {
    if (!m_device_allocated) {
      m_d_tf_data = static_cast<TabFuncData*>(
        amrex::The_Device_Arena()->alloc(sizeof(m_h_tf_data)));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  virtual void deallocate()
  {
    amrex::The_Pinned_Arena()->free(m_h_tf_data.dimLengths);
    amrex::The_Pinned_Arena()->free(m_h_tf_data.dimDataSpacing);
    amrex::The_Pinned_Arena()->free(m_h_tf_data.grids);
    amrex::The_Pinned_Arena()->free(m_h_tf_data.values);
    amrex::The_Pinned_Arena()->free(m_h_tf_data.varnames);
    amrex::The_Pinned_Arena()->free(m_h_tf_data.dimnames);
    if (m_device_allocated)
      amrex::The_Device_Arena()->free(m_d_tf_data);
  }

  virtual void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_tf_data, &m_h_tf_data + 1, m_d_tf_data);
    }
  }

  virtual ManFuncData& host_manfunc_data() { return m_h_tf_data; }

  virtual ManFuncData* device_manfunc_data() { return m_d_tf_data; }

  TabFuncData& host_tabfunc_data() { return m_h_tf_data; }

  TabFuncData* device_tabfunc_data() { return m_d_tf_data; }

private:
  int m_verbose = 0;
  TabFuncData m_h_tf_data;
  TabFuncData* m_d_tf_data;
  bool m_device_allocated{false};
};

// Struct to contain tabulated data
class TabFunc : public ManifoldFunc
{
public:
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabFunc() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  TabFunc(TabFuncParams::TabFuncData* tf_data_in) { tf_data = tf_data_in; }

  // Lookup directly if you already know index
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void
  get_value(const int ivar, const amrex::Real interpdata[], amrex::Real& out)
  {
    // Get alphas and indices
    int indices[MAXD_TABLE];
    amrex::Real alphas[MAXD_TABLE];
    amrex::Real dxinv[MAXD_TABLE];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    out = interpolate(
      indices, alphas, &tf_data->values[ivar * tf_data->varSpacing]);
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
    int indices[MAXD_TABLE];
    amrex::Real alphas[MAXD_TABLE];
    amrex::Real dxinv[MAXD_TABLE];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    for (int i = 0; i < nvar; ++i) {
      out[i] =
        (ivar[i] >= 0)
          ? interpolate(
              indices, alphas, &tf_data->values[ivar[i] * tf_data->varSpacing])
          : 0.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_all_values(const amrex::Real interpdata[], amrex::Real out[])
  {
    // Get alphas and indices
    int indices[MAXD_TABLE];
    amrex::Real alphas[MAXD_TABLE];
    amrex::Real dxinv[MAXD_TABLE];
    get_indices_alphas_dxinv(interpdata, indices, alphas, dxinv);

    // interpolate down
    for (int i = 0; i < tf_data->Nvar; ++i) {
      out[i] =
        interpolate(indices, alphas, &tf_data->values[i * tf_data->varSpacing]);
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
    int indices[MAXD_TABLE];
    amrex::Real alphas[MAXD_TABLE];
    amrex::Real dxinv[MAXD_TABLE];
    get_indices_alphas_dxinv(derivloc, indices, alphas, dxinv);

    // finite differences from table
    amrex::Real out = differentiate(
      indices, alphas, dxinv, &tf_data->values[ivar * tf_data->varSpacing],
      derivs);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void calculate_Wdot(
    const int paramidx, const amrex::Real indata[], amrex::Real& out)
  {
    out = 0.0;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual ManifoldModel model() { return tf_data->manmodel; }

  /*
  int idx_density;
  int idx_T;
  int idx_Wdot;
  int idx_rhoD;
  int idx_mu;
  amrex::Real Pnom;
  */

private:
  TabFuncParams::TabFuncData* tf_data;

  // -----------------------------------------------------------
  // Search for the closest index in an array to a given value
  // using the bisection technique.
  // INPUTS/OUTPUTS:
  // xtable(0:n-1) => array to search in (ascending order)
  // n             => array size
  // x             => x location
  // idxlo (return)=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
  // -----------------------------------------------------------
  AMREX_GPU_DEVICE
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
      // {idx+=1;}; idx-=1;
      int idx = locate(&tf_data->grids[start], dimLen, interpdata[dim]);
#endif
      // test to see cost if no lookup
      // int idx=42;

      indices[dim] = idx;
      idx += start;
      amrex::Real dx_inv =
        1.0 / (tf_data->grids[idx + 1] - tf_data->grids[idx]);
      dxinv[dim] = dx_inv;
      amrex::Real alpha = (tf_data->grids[idx + 1] - interpdata[dim]) * dx_inv;
      alphas[dim] = std::min(
        std::max(alpha, 0.0), 1.0); // Clip to avoid extrapolating when out of
                                    // bounds // REPLACE with AMReX function
      // std::cout << "alpha " << dim << ": " << alphas[dim] << std::endl;
      start += dimLen;
    }
  }

  // Wrapper for interpolation process, chooses appropriate approach based on
  // dimensionality some notes: the process is entirely recursive, but for  Ndim
  // <= 5 is unrolled to allow inlining Inlining not supported for Ndim > 5, but
  // that is unlikely to be used anyway minor performance gain possible by
  // modifying code so selecting dimensionality isn't neeeded for each lookup
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    amrex::Real out_value;
    switch (tf_data->Ndim) {
    case 1:
      out_value = interpolate_1D(indices, alphas, data);
      break;
    case 2:
      out_value = interpolate_2D(indices, alphas, data);
      break;
    case 3:
      out_value = interpolate_3D(indices, alphas, data);
      break;
    case 4:
      out_value = interpolate_4D(indices, alphas, data);
      break;
    case 5:
      out_value = interpolate_5D(indices, alphas, data);
      break;
    default:
#if MAXD_TABLE > 5
      out_value = interpolate_ND(indices, alphas, data);
#else
      amrex::Abort(
        "Tabulated function: Failure: table must have 1-5 dimensions \n)");
#endif
    }
    return out_value;
  }

  // recursive function for interpolating each dimension - cannot be inlined -
  // use for Ndim >5 until Ndim=5
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_ND(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    return interpolate_ND(indices, alphas, data, tf_data->Ndim - 1);
  }

  AMREX_GPU_HOST_DEVICE
  amrex::Real interpolate_ND(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real data[],
    const int dim)
  {
    if (dim == 5) {
      return interpolate_5D(indices, alphas, data);
    } else if (dim > 0) {
      amrex::Real out_value_left = interpolate_ND(
        indices, alphas, &data[indices[dim] * tf_data->dimDataSpacing[dim]],
        dim - 1);
      amrex::Real out_value_right = interpolate_ND(
        indices, alphas,
        &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]], dim - 1);
      return alphas[dim] * out_value_left +
             (1.0 - alphas[dim]) * out_value_right;
    } else {
      amrex::Real out_value_left = data[indices[dim]];
      amrex::Real out_value_right = data[indices[dim] + 1];
      return alphas[dim] * out_value_left +
             (1.0 - alphas[dim]) * out_value_right;
    }
  }

  // unrolled versions of above (for inlining/performance)
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_5D(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    int idx = indices[4];
    int dimDataSpac = tf_data->dimDataSpacing[4];
    amrex::Real out_value_left =
      interpolate_4D(indices, alphas, &data[idx * dimDataSpac]);
    amrex::Real out_value_right =
      interpolate_4D(indices, alphas, &data[(idx + 1) * dimDataSpac]);
    amrex::Real alpha = alphas[4];
    return alpha * out_value_left + (1.0 - alpha) * out_value_right;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_4D(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    int idx = indices[3];
    int dimDataSpac = tf_data->dimDataSpacing[3];
    amrex::Real out_value_left =
      interpolate_3D(indices, alphas, &data[idx * dimDataSpac]);
    amrex::Real out_value_right =
      interpolate_3D(indices, alphas, &data[(idx + 1) * dimDataSpac]);
    amrex::Real alpha = alphas[3];
    return alpha * out_value_left + (1.0 - alpha) * out_value_right;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_3D(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    int idx = indices[2];
    int dimDataSpac = tf_data->dimDataSpacing[2];
    amrex::Real out_value_left =
      interpolate_2D(indices, alphas, &data[idx * dimDataSpac]);
    amrex::Real out_value_right =
      interpolate_2D(indices, alphas, &data[(idx + 1) * dimDataSpac]);
    amrex::Real alpha = alphas[2];
    return alpha * out_value_left + (1.0 - alpha) * out_value_right;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_2D(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    int idx = indices[1];
    int dimDataSpac = tf_data->dimDataSpacing[1];
    amrex::Real out_value_left =
      interpolate_1D(indices, alphas, &data[idx * dimDataSpac]);
    amrex::Real out_value_right =
      interpolate_1D(indices, alphas, &data[(idx + 1) * dimDataSpac]);
    amrex::Real alpha = alphas[1];
    return alpha * out_value_left + (1.0 - alpha) * out_value_right;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real interpolate_1D(
    const int indices[], const amrex::Real alphas[], const amrex::Real data[])
  {
    int idx = indices[0];
    amrex::Real out_value_left = data[idx];
    amrex::Real out_value_right = data[idx + 1];
    amrex::Real alpha = alphas[0];
    return alpha * out_value_left + (1.0 - alpha) * out_value_right;
  }

  // Wrapper for interpolation process, chooses appropriate approach based on
  // dimensionality same caveats apply as for the interpolation process
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    amrex::Real out_value;
    switch (tf_data->Ndim) {
    case 1:
      out_value = differentiate_1D(indices, alphas, dxinv, data, derivs);
      break;
    case 2:
      out_value = differentiate_2D(indices, alphas, dxinv, data, derivs);
      break;
    case 3:
      out_value = differentiate_3D(indices, alphas, dxinv, data, derivs);
      break;
    case 4:
      out_value = differentiate_4D(indices, alphas, dxinv, data, derivs);
      break;
    case 5:
      out_value = differentiate_5D(indices, alphas, dxinv, data, derivs);
      break;
    default:
#if MAXD_TABLE > 5
      out_value = differentiate_ND(indices, alphas, dxinv, data, derivs);
#else
      amrex::Abort(
        "Tabulated function: Failure: table must have 1-5 dimensions \n)");
#endif
    }
    return out_value;
  }

  // Recursive functions
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_ND(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    amrex::Real outval =
      differentiate_ND(indices, alphas, dxinv, data, derivs, tf_data->Ndim - 1);
    return outval;
  }

  AMREX_GPU_HOST_DEVICE
  amrex::Real differentiate_ND(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[],
    const int dim)
  {
    amrex::Real value_left, value_right;
    amrex::Real alpha = alphas[dim];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    if (dim == 0) {
      int idx = indices[0];
      value_left = data[idx];
      value_right = data[idx + 1];
    } else {
      amrex::Real derivs_left[MAXD_TABLE];
      amrex::Real derivs_right[MAXD_TABLE];
      // Exit recursion for performance when possible
      if (dim == 5) {
        value_left = differentiate_5D(
          indices, alphas, dxinv,
          &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left);
        value_right = differentiate_5D(
          indices, alphas, dxinv,
          &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]],
          derivs_right);
      } else {
        value_left = differentiate_ND(
          indices, alphas, dxinv,
          &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left,
          dim - 1);
        value_right = differentiate_ND(
          indices, alphas, dxinv,
          &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]],
          derivs_right, dim - 1);
      }
      // interpolate existing derivatives
      for (int ddim = 0; ddim < dim; ddim++) {
        derivs[ddim] =
          alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
      }
    }
    derivs[dim] =
      (value_right - value_left) * dxinv[dim]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }

  // inlinable functions
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_5D(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    const int dim = 4;
    amrex::Real alpha = alphas[dim];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    amrex::Real derivs_left[dim];
    amrex::Real derivs_right[dim];
    amrex::Real value_left = differentiate_4D(
      indices, alphas, dxinv,
      &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left);
    amrex::Real value_right = differentiate_4D(
      indices, alphas, dxinv,
      &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]], derivs_right);
    // interpolate existing derivatives
    for (int ddim = 0; ddim < dim; ddim++) {
      derivs[ddim] =
        alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
    }
    derivs[dim] =
      (value_right - value_left) * dxinv[dim]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_4D(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    const int dim = 3;
    amrex::Real alpha = alphas[dim];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    amrex::Real derivs_left[dim];
    amrex::Real derivs_right[dim];
    amrex::Real value_left = differentiate_3D(
      indices, alphas, dxinv,
      &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left);
    amrex::Real value_right = differentiate_3D(
      indices, alphas, dxinv,
      &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]], derivs_right);
    // interpolate existing derivatives
    for (int ddim = 0; ddim < dim; ddim++) {
      derivs[ddim] =
        alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
    }
    derivs[dim] =
      (value_right - value_left) * dxinv[dim]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_3D(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    const int dim = 2;
    amrex::Real alpha = alphas[dim];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    amrex::Real derivs_left[dim];
    amrex::Real derivs_right[dim];
    amrex::Real value_left = differentiate_2D(
      indices, alphas, dxinv,
      &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left);
    amrex::Real value_right = differentiate_2D(
      indices, alphas, dxinv,
      &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]], derivs_right);
    // interpolate existing derivatives
    for (int ddim = 0; ddim < dim; ddim++) {
      derivs[ddim] =
        alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
    }
    derivs[dim] =
      (value_right - value_left) * dxinv[dim]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_2D(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    const int dim = 1;
    amrex::Real alpha = alphas[dim];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    amrex::Real derivs_left[dim];
    amrex::Real derivs_right[dim];
    amrex::Real value_left = differentiate_1D(
      indices, alphas, dxinv,
      &data[indices[dim] * tf_data->dimDataSpacing[dim]], derivs_left);
    amrex::Real value_right = differentiate_1D(
      indices, alphas, dxinv,
      &data[(indices[dim] + 1) * tf_data->dimDataSpacing[dim]], derivs_right);
    // interpolate existing derivatives
    for (int ddim = 0; ddim < dim; ddim++) {
      derivs[ddim] =
        alpha * derivs_left[ddim] + oneMinusAlpha * derivs_right[ddim];
    }
    derivs[dim] =
      (value_right - value_left) * dxinv[dim]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real differentiate_1D(
    const int indices[],
    const amrex::Real alphas[],
    const amrex::Real dxinv[],
    const amrex::Real data[],
    amrex::Real derivs[])
  {
    amrex::Real alpha = alphas[0];
    amrex::Real oneMinusAlpha = 1.0 - alpha;
    int idx = indices[0];
    amrex::Real value_left = data[idx];
    amrex::Real value_right = data[idx + 1];
    derivs[0] = (value_right - value_left) * dxinv[0]; // take the derivative
    return alpha * value_left + oneMinusAlpha * value_right; // interpolate
  }
};
} // namespace physics
} // namespace pele
#endif
