#include "DistUniform.H"
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
