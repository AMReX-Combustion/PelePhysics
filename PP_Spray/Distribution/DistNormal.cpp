#include "DistNormal.H"
#include "AMReX_ParmParse.H"

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
