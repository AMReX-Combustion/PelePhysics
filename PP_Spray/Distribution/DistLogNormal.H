#ifndef DISTLOGNORMAL_H
#define DISTLOGNORMAL_H

#include "DistBase.H"
#include "AMReX_Random.H"

class DistLogNormal : public DistBase::Register<DistLogNormal>
{
public:
  static std::string identifier() { return "DistLogNormal"; }

  void init(const std::string &a_prefix) override;

  amrex::Real get_dia() override;

private:
  amrex::Real m_log_mean;
  amrex::Real m_log_std;
};
#endif    
