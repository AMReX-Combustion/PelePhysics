// --------------------------------------------------------------------
// WritePlotFile.H
// --------------------------------------------------------------------
#ifndef BL_WRITEPLOTFILE_H
#define BL_WRITEPLOTFILE_H

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>


void WritePlotfile(const std::string         &pfversion,
                   const amrex::Vector<const amrex::MultiFab*> &data,
                   const amrex::Real                 time,
                   const amrex::Vector<amrex::Real>         &probLo,
                   const amrex::Vector<amrex::Real>         &probHi,
                   const amrex::Vector<int>          &refRatio,
                   const amrex::Vector<amrex::Box>          &probDomain,
                   const amrex::Vector<amrex::Vector<amrex::Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const amrex::Vector<std::string>  &names,
                   const bool                 verbose,
		   const bool                 isCartGrid = false,
		   const amrex::Real                *vfeps = 0,
		   const int                 *levelSteps = 0);


#endif
// --------------------------------------------------------------------
