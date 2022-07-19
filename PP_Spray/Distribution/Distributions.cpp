#include "Distributions.H"
#include "AMReX_Random.H"
#include "AMReX_ParmParse.H"

void
DistUniform::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);
    pp.get("diameter", m_diam);
}

amrex::Real
DistUniform::get_dia()
{
    return m_diam;
}

void
DistNormal::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);
    pp.get("mean_diam", m_mean);
    pp.get("std_dev", m_std);
}

amrex::Real
DistNormal::get_dia()
{
    return amrex::RandomNormal(m_mean, m_std);
}

void
DistLogNormal::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);
    amrex::Real mean;
    amrex::Real std;
    pp.get("mean_diam", mean);
    pp.get("std_dev", std);
    m_log_mean = 2. * std::log(mean) - 0.5 * std::log(std*std + mean*mean);
    m_log_std = std::sqrt(amrex::max(-2. * std::log(mean) + std::log(std*std + mean*mean), 0.));
}

amrex::Real
DistLogNormal::get_dia()
{
    return std::exp(amrex::RandomNormal(m_log_mean, m_log_std));
}

void
DistWeibull::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);
    pp.get("mean_diam", m_mean);
    pp.get("k", m_k);
}

amrex::Real
DistWeibull::get_dia()
{
    amrex::Real fact =
      -std::log(0.5 * (1. - std::erf(amrex::Random() / std::sqrt(2.))));
    return m_mean * std::pow(fact, 1. / m_k);
}
