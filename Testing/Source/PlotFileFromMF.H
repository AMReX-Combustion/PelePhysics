#ifndef _PlotFileFromMF_
#define _PlotFileFromMF_
#include <string>
#include <AMReX_MultiFab.H>

void PlotFileFromMF(const amrex::MultiFab& mf,
		    const std::string& oFile);

#endif // _PlotFileFromMF_
