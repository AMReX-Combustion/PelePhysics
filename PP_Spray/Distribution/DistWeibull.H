#ifndef DISTWEIBULL_H
#define DISTWEIBULL_H

#include "DistBase.H"
#include "AMReX_Random.H"

class DistWeibull : public DistBase::Register<DistWeibull>
{
public:
    static std::string identifier() { return "DistWeibull"; }

    void init(const std::string &a_prefix) override;

    amrex::Real get_dia() override;

private:
    amrex::Real m_mean;
    amrex::Real m_k;
};
#endif    
