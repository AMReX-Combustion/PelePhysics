#include "DistLogNormal.H"
#include "AMReX_ParmParse.H"

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
