#ifndef PMF_DATA_H
#define PMF_DATA_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <mechanism.H>

// Define a class to hold/manage the PMF data,
namespace pele::physics::PMF {

class PmfData
{
public:
  PmfData() = default;

  ~PmfData() = default;

  void initialize()
  {
    amrex::ParmParse pp("pmf");
    std::string datafile;
    pp.query("v", m_verbose);
    int do_average = 1;
    pp.query("do_cellAverage", do_average);
    if (!pp.contains("datafile")) {
      amrex::Abort("pmf.datafile is required when using pmf");
    }
    pp.get("datafile", datafile);
    read_pmf(datafile, do_average, m_verbose);
    allocate();
  }

  void read_pmf(const std::string& fname, int a_doAverage, int a_verbose);

  struct DataContainer
  {
    int m_nPoint;
    int m_nVar;
    int m_doAverage = 0;
    amrex::Real* pmf_X;
    amrex::Real* pmf_Y;
  };

  void allocate()
  {
    if (!m_device_allocated) {
      m_data_d = (DataContainer*)amrex::The_Arena()->alloc(sizeof(m_data_h));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void deallocate()
  {
    if (m_host_allocated) {
      amrex::The_Pinned_Arena()->free(m_data_h.pmf_X);
      amrex::The_Pinned_Arena()->free(m_data_h.pmf_Y);
    }
    if (m_device_allocated) {
      amrex::The_Arena()->free(m_data_d);
    }
  }

  void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_data_h, &m_data_h + 1, m_data_d);
    }
  }

  // Variables names
  amrex::Vector<std::string> pmf_names;

  // Accessors
  DataContainer getHostData() { return m_data_h; }

  DataContainer* getDeviceData() { return m_data_d; }

private:
  int m_verbose = 0;
  DataContainer m_data_h;
  DataContainer* m_data_d = nullptr;
  bool m_device_allocated{false};
  bool m_host_allocated{false};
};
} // namespace pele::physics::PMF
#endif
