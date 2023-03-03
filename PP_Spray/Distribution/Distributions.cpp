#include "Distributions.H"
#include "AMReX_Random.H"
#include "AMReX_ParmParse.H"

void
Uniform::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  pp.get("diameter", m_diam);
}

void
Uniform::init(const amrex::Real& mean, const amrex::Real& /*std*/)
{
  m_diam = mean;
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
  pp.get("mean_dia", m_mean);
  pp.get("std_dev", m_std);
}

void
Normal::init(const amrex::Real& mean, const amrex::Real& std)
{
  m_mean = mean;
  m_std = std;
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
LogNormal::init(const amrex::Real& mean, const amrex::Real& std)
{
  m_mean = mean;
  m_log_mean = 2. * std::log(mean) - 0.5 * std::log(std * std + mean * mean);
  m_log_std = std::sqrt(
    amrex::max(-2. * std::log(mean) + std::log(std * std + mean * mean), 0.));
}
void
LogNormal::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  amrex::Real mean;
  amrex::Real std;
  pp.get("mean_dia", mean);
  pp.get("std_dev", std);
  init(mean, std);
}

amrex::Real
LogNormal::get_dia()
{
  return std::exp(amrex::RandomNormal(m_log_mean, m_log_std));
}

amrex::Real
LogNormal::get_avg_dia()
{
  return m_mean;
}

void
Weibull::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  pp.get("mean_dia", m_mean);
  pp.get("k", m_k);
}

void
Weibull::init(const amrex::Real& mean, const amrex::Real& k)
{
  m_mean = mean;
  m_k = k;
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

void
ChiSquared::init(const std::string& a_prefix)
{
  amrex::ParmParse pp(a_prefix);
  amrex::Real d32, dummy;
  pp.get("d32", d32);
  init(d32, dummy);
}

void
ChiSquared::init(const amrex::Real& d32, const amrex::Real& /*dummy*/)
{
  m_d32 = d32;
  amrex::Real xiend = 12.;
  amrex::Real dxi = xiend / 100.;
  amrex::Real rend =
    1. - std::exp(-xiend) * (1. + xiend * (1. + xiend * (0.5 + xiend / 6.)));
  for (int i = 0; i < 100; i++) {
    amrex::Real xi = dxi * (i + 1);
    amrex::Real rval =
      1. - std::exp(-xi) * (1. + xi * (1. + xi * (0.5 + xi / 6.)));
    rvals[i] = rval / rend;
  }
}

amrex::Real
ChiSquared::get_dia()
{
  amrex::Real dmean = m_d32 / 3.;
  amrex::Real dxi = 12. / 100.;
  amrex::Real fact = amrex::Random();
  int curn = 0;
  amrex::Real curr = rvals[0];
  amrex::Real curxi = 0.;
  while (fact > curr) {
    curn++;
    curr = rvals[curn];
    curxi += dxi;
  }
  return curxi * dmean;
}

amrex::Real
ChiSquared::get_avg_dia()
{
  // Rough estimate of mean
  amrex::Real dmean = m_d32 / 3.;
  amrex::Real dxi = 12. / 100.;
  amrex::Real fact = 0.5;
  int curn = 0;
  amrex::Real curr = rvals[0];
  amrex::Real curxi = 0.;
  while (fact > curr) {
    curn++;
    curr = rvals[curn];
    curxi += dxi;
  }
  return curxi * dmean;
}
