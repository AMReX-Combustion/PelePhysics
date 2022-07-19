#ifndef DISTNORMAL_H
#define DISTNORMAL_H

#include "DistBase.H"
#include "AMReX_Random.H"

class DistNormal : public DistBase::Register<DistNormal>
{
public:
  static std::string identifier() { return "DistNormal"; }

  void init(const std::string &a_prefix) override;

  amrex::Real get_dia() override;

private:
  amrex::Real m_mean;
  amrex::Real m_std;
};
#endif    
