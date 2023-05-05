#ifndef _DIAGBASE_H_
#define _DIAGBASE_H_

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "DiagFilter.H"
#include "Factory.H"

class DiagBase : public pele::physics::Factory<DiagBase>
{
public:
  static std::string base_identifier() { return "DiagBase"; }

  virtual void init(const std::string& a_prefix, std::string_view a_diagName);

  virtual void close() = 0;

  bool needUpdate() const { return need_update; }

  virtual bool doDiag(const amrex::Real& a_time, int a_nstep);

  virtual void prepare(
    int a_nlevels,
    const amrex::Vector<amrex::Geometry>& a_geoms,
    const amrex::Vector<amrex::BoxArray>& a_grids,
    const amrex::Vector<amrex::DistributionMapping>& a_dmap,
    const amrex::Vector<std::string>& a_varNames);

  virtual void processDiag(
    int a_nstep,
    const amrex::Real& a_time,
    const amrex::Vector<const amrex::MultiFab*>& a_state,
    const amrex::Vector<std::string>& a_varNames) = 0;

  virtual void addVars(amrex::Vector<std::string>& a_varList);

  static int getFieldIndex(
    const std::string& a_field, const amrex::Vector<std::string>& a_varList);

protected:
  std::string m_diagfile;
  int m_verbose{0};
  amrex::Real m_per{-1.0};
  int m_interval{-1};
  bool need_update{true};
  bool first_time{true};
  amrex::Vector<DiagFilter> m_filters{};
  amrex::Gpu::DeviceVector<DiagFilterData> m_filterData;
};

#endif
