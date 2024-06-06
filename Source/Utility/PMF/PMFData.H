#ifndef PMF_DATA_H
#define PMF_DATA_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <mechanism.H>
#include "PeleParamsGeneric.H"

// Define a class to hold/manage the PMF data,
namespace pele::physics {
namespace PMF::PmfData {

struct DataContainer
{
  int m_nPoint;
  int m_nVar;
  int m_doAverage = 0;
  amrex::Real* pmf_X;
  amrex::Real* pmf_Y;
};
} // namespace PMF::PmfData

template <>
struct InitParm<PMF::PmfData::DataContainer>
{

  static void read_pmf(
    const std::string& fname,
    const int a_doAverage,
    const int a_verbose,
    PMF::PmfData::DataContainer& h_pmf_data);

  static void host_initialize(PeleParams<PMF::PmfData::DataContainer>* parm_in)
  {
    amrex::ParmParse pp("pmf");
    std::string datafile;
    int verbose;
    pp.query("v", verbose);
    int do_average = 1;
    pp.query("do_cellAverage", do_average);
    if (!pp.contains("datafile")) {
      amrex::Abort("pmf.datafile is required when using pmf");
    }
    pp.get("datafile", datafile);
    read_pmf(datafile, do_average, verbose, parm_in->m_h_parm);
  };

  static void
  host_deallocate(PeleParams<PMF::PmfData::DataContainer>* /*parm_in*/)
  {
  }
};

} // namespace pele::physics
#endif
