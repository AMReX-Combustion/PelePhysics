#ifndef DISTUNIFORM_H
#define DISTUNIFORM_H

#include "DistBase.H"

class DistUniform : public DistBase::Register<DistUniform>
{
public:
  static std::string identifier() { return "DistUniform"; }

  void init(const std::string &a_prefix) override;

  amrex::Real get_dia() override;

private:
  amrex::Real m_diam;
};
#endif    
