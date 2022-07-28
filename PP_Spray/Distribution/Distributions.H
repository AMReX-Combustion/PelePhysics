#ifndef DISTRIBS_H
#define DISTRIBS_H

#include "DistBase.H"

class Uniform : public DistBase::Register<Uniform>
{
public:
  static std::string identifier() { return "Uniform"; }

  void init(const std::string& a_prefix) override;

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_diam;
};

class Normal : public DistBase::Register<Normal>
{
public:
  static std::string identifier() { return "Normal"; }

  void init(const std::string& a_prefix) override;

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_mean;
  amrex::Real m_std;
};

class LogNormal : public DistBase::Register<LogNormal>
{
public:
  static std::string identifier() { return "LogNormal"; }

  void init(const std::string& a_prefix) override;

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_log_mean;
  amrex::Real m_log_std;
};

class Weibull : public DistBase::Register<Weibull>
{
public:
  static std::string identifier() { return "Weibull"; }

  void init(const std::string& a_prefix) override;

  amrex::Real get_dia() override;
  amrex::Real get_avg_dia() override;

private:
  amrex::Real m_mean;
  amrex::Real m_k;
};
#endif
