#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

void
CKAWT(amrex::Real* /*awt*/)
{
}

void
CKNCF(int* /*ncf*/)
{
}

void
CKSYME_STR(amrex::Vector<std::string>& /*ename*/)
{
}
void
CKSYMS_STR(amrex::Vector<std::string>& /*kname*/)
{
}

#ifdef COMPILE_JACOBIAN
// compute sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* /*nJdata*/, int* /*consP/*, int /*NCELLS*/)
{
}

// compute sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* /*nJdata*/, int* /*consP/*, int /*NCELLS*/)
{
}

// compute sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(
  int* /*nJdata*/,
  int* /*consP/*)
{
}

// compute sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void
SPARSITY_PREPROC_CSC(int* /*rowVals*/
  ,
  int* /*colPtrs*/,
  int* /*consP/*, int /*NCELLS*/)
{
}

// compute sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void
SPARSITY_PREPROC_CSR(
  int* /*colVals*/,
  int* /*rowPtrs*/,
  int* /*consP/*, int /*NCELLS*/,
  int /*base*/)
{
}

// compute sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* /*colVals*/,
  int* /*rowPtr*/,
  int* /*consP/*, int /*NCELLS*/,
  int /*base*/)
{
}

// compute sparsity pattern of the simplified (for precond) system Jacobian on
// CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* /*rowVals*/,
  int* /*colPtrs*/,
  int* /*indx*/,
  int* /*consP/*)
{
}

// compute sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
int* /*colVals*/
  ,
  int* /*rowPtr*/,
  int* /*consP/*, int /*base*/)
{
}
#endif
#endif
