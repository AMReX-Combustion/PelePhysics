#include "Distributions.H"
#include "AMReX_Random.H"
#include "AMReX_ParmParse.H"

void
Uniform::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  pp.get("diameter", m_diam);
}

amrex::Real
Uniform::get_dia()
{
  return m_diam;
}

amrex::Real
Uniform::get_avg_dia()
{
  return get_dia();
}

void
Normal::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  pp.get("mean_diam", m_mean);
  pp.get("std_dev", m_std);
}

amrex::Real
Normal::get_dia()
{
  return amrex::RandomNormal(m_mean, m_std);
}

amrex::Real
Normal::get_avg_dia()
{
  return m_mean;
}

void
LogNormal::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  amrex::Real mean;
  amrex::Real std;
  pp.get("mean_diam", mean);
  pp.get("std_dev", std);
  m_log_mean = 2. * std::log(mean) - 0.5 * std::log(std * std + mean * mean);
  m_log_std = std::sqrt(
    amrex::max(-2. * std::log(mean) + std::log(std * std + mean * mean), 0.));
}

amrex::Real
LogNormal::get_dia()
{
  return std::exp(amrex::RandomNormal(m_log_mean, m_log_std));
}

amrex::Real
LogNormal::get_avg_dia()
{
  return std::exp(m_log_mean);
}

void
Weibull::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  pp.get("mean_diam", m_mean);
  pp.get("k", m_k);
}

amrex::Real
Weibull::get_dia()
{
  amrex::Real fact =
    -std::log(0.5 * (1. - std::erf(amrex::Random() / std::sqrt(2.))));
  return m_mean * std::pow(fact, 1. / m_k);
}

amrex::Real
Weibull::get_avg_dia()
{
  return m_mean;
}
