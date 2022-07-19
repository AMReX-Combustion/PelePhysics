#include "DistWeibull.H"
#include "AMReX_ParmParse.H"

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
