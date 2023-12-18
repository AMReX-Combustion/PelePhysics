#ifndef AMR_PARAM_H
#define AMR_PARAM_H

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <string>

namespace PeleRad {

struct AMRParam
{
public:
  amrex::ParmParse pp_;

  int max_level_;
  int ref_ratio_;
  int n_cell_;
  int max_grid_size_;
  int prob_type_;
  std::string plot_file_name_;

  AMREX_GPU_HOST
  AMRParam(amrex::ParmParse const& pp) : pp_(pp)
  {
    pp_.query("max_level", max_level_);
    pp_.query("ref_ratio", ref_ratio_);
    pp_.query("n_cell", n_cell_);
    pp_.query("max_grid_size", max_grid_size_);
    pp_.query("prob_type", prob_type_);
    pp.query("plot_file_name", plot_file_name_);
  }
};

} // namespace PeleRad
#endif
