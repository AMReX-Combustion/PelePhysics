#ifndef DISTBASE_H
#define DISTBASE_H

#include "Factory.H"

class DistBase : public pele::physics::Factory<DistBase>
{
public:
  static std::string base_identifier() { return "DistBase"; }

  virtual void init(const std::string& a_prefix) = 0;

  virtual amrex::Real get_dia() = 0;
  virtual amrex::Real get_avg_dia() = 0;

protected:
  int m_verbose = 0;
};
#endif
