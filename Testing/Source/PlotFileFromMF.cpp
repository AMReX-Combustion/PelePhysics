
#include <AMReX_Utility.H>
#include <WritePlotFile.H>

using namespace amrex;

void PlotFileFromMF(const amrex::MultiFab& mf,
		    const std::string& oFile)
{
    const std::string pfversion = "HyperCLaw-V1.1";
    Vector<const MultiFab*> data(1);
    data[0] = &mf;
    Real time = 0;
    const Box dbox = mf.boxArray().minimalBox();
    Vector<Box> probDomain(1,dbox);
    Vector<Real> probLo(BL_SPACEDIM,0), probHi(BL_SPACEDIM,0);
    for (int d=0; d<BL_SPACEDIM; ++d) {
	probHi[d] = dbox.length(d);
    }
    Vector<int> refRatio(0);
    Vector<Vector<Real> > dxLevel(1, Vector<Real>(BL_SPACEDIM, 1));
    int coordSys = 0;
    Vector<std::string> names(mf.nComp());
    for (int i=0; i<mf.nComp(); ++i) {
	names[i] = amrex::Concatenate("MF_",i,5);
    }
    
    WritePlotfile(pfversion,data,time,probLo,probHi,refRatio,probDomain,dxLevel,coordSys,oFile,names,false);
}
