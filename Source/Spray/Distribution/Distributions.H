#ifndef DISTRIBS_H
#define DISTRIBS_H

#include "DistBase.H"

class Uniform : public DistBase::Register<Uniform>
{
public:
  static std::string identifier() { return "Uniform"; }

  void init(const std::string& a_prefix) override;

  void init(const amrex::Real& diam);

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_diam = 0.;
};

class Normal : public DistBase::Register<Normal>
{
public:
  static std::string identifier() { return "Normal"; }

  void init(const std::string& a_prefix) override;

  void init(const amrex::Real& mean_d, const amrex::Real& std);

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_mean = 0.;
  amrex::Real m_std = 0.;
};

class LogNormal : public DistBase::Register<LogNormal>
{
public:
  static std::string identifier() { return "LogNormal"; }

  void init(const std::string& a_prefix) override;

  void init(const amrex::Real& mean_d, const amrex::Real& std);

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_log_mean = 0.;
  amrex::Real m_log_std = 0.;
  amrex::Real m_mean = 0.;
};

class Weibull : public DistBase::Register<Weibull>
{
public:
  static std::string identifier() { return "Weibull"; }

  void init(const std::string& a_prefix) override;

  void init(const amrex::Real& mean_d, const amrex::Real& k);

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_mean = 0.;
  amrex::Real m_k = 0.;
};

class ChiSquared : public DistBase::Register<ChiSquared>
{
public:
  static std::string identifier() { return "ChiSquared"; }

  void init(const std::string& a_prefix) override;

  void init(const amrex::Real& d32);

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::GpuArray<amrex::Real, 100> rvals = {{0.0}};
  amrex::Real m_d32 = 0.;
};
#endif
