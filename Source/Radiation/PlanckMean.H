#ifndef PLANCK_MEAN_H
#define PLANCK_MEAN_H

#include <fstream>

#include <AMReX.H>
#include <AMReX_Gpu.H>

namespace PeleRad {

class PlanckMean
{
private:
  amrex::GpuArray<amrex::Real, 126UL> kpco2_;
  amrex::GpuArray<amrex::Real, 126UL> kph2o_;
  amrex::GpuArray<amrex::Real, 126UL> kpco_;
  amrex::GpuArray<amrex::Real, 126UL> kpch4_;
  amrex::GpuArray<amrex::Real, 126UL> kpc2h4_;
  amrex::GpuArray<amrex::Real, 126UL> kpsoot_;

  AMREX_GPU_HOST
  void read_kp_file(
    std::string file_path, amrex::GpuArray<amrex::Real, 126UL>& file_data)
  {
    std::ifstream datastream(file_path);
    if (!datastream.is_open()) {
      amrex::Abort("PeleRad: Failed to open data file: " + file_path);
    }
    int idx = 0;
    for (amrex::Real T_temp, kp_temp; datastream >> T_temp >> kp_temp;) {
      if (!amrex::almostEqual(
            T_temp, 300.0 + 20.0 * static_cast<amrex::Real>(idx))) {
        amrex::Abort(
          "PeleRad: Invalid temperatures in data file: " + file_path);
      }
      if (idx > 125) {
        amrex::Abort(
          "PeleRad: Invalid number of entries in data file, must have 126: " +
          file_path);
      }
      file_data[idx] = kp_temp;
      idx++;
    }
    if (idx < 126) {
      amrex::Abort(
        "PeleRad: Invalid number of entries in data file, must have 126: " +
        file_path);
    }
    datastream.close();
  }

public:
  AMREX_GPU_HOST
  PlanckMean() = default;

  AMREX_GPU_HOST
  PlanckMean(std::string data_path) { load(data_path); }

  AMREX_GPU_HOST
  PlanckMean(PlanckMean const&) = delete;

  AMREX_GPU_HOST
  PlanckMean& operator=(PlanckMean const&) = delete;

  AMREX_GPU_HOST
  void load(std::string data_path)
  {
    read_kp_file(data_path + "/kpl_co2.dat", kpco2_);
    read_kp_file(data_path + "/kpl_h2o.dat", kph2o_);
    read_kp_file(data_path + "/kpl_co.dat", kpco_);
    read_kp_file(data_path + "/kpl_ch4.dat", kpch4_);
    read_kp_file(data_path + "/kpl_c2h4.dat", kpc2h4_);
    read_kp_file(data_path + "/kpl_soot.dat", kpsoot_);
  }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kpco2() const { return kpco2_; }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kph2o() const { return kph2o_; }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kpco() const { return kpco_; }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kpch4() const { return kpch4_; }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kpc2h4() const { return kpc2h4_; }

  AMREX_GPU_HOST_DEVICE
  const amrex::GpuArray<amrex::Real, 126UL>& kpsoot() const { return kpsoot_; }
};

} // namespace PeleRad

#endif
