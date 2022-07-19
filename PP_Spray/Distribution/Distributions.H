#ifndef DISTRIBS_H
#define DISTRIBS_H

#include "DistBase.H"

class DistUniform : public DistBase::Register<DistUniform>
{
public:
    static std::string identifier() { return "DistUniform"; }

    void init(const std::string &a_prefix) override;

    amrex::Real get_dia() override;

private:
    amrex::Real m_diam;
};

class DistNormal : public DistBase::Register<DistNormal>
{
public:
    static std::string identifier() { return "DistNormal"; }

    void init(const std::string &a_prefix) override;

    amrex::Real get_dia() override;

private:
    amrex::Real m_mean;
    amrex::Real m_std;
};

class DistLogNormal : public DistBase::Register<DistLogNormal>
{
public:
    static std::string identifier() { return "DistLogNormal"; }

    void init(const std::string &a_prefix) override;

    amrex::Real get_dia() override;

private:
    amrex::Real m_log_mean;
    amrex::Real m_log_std;
};

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
