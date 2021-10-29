#ifndef _REACTORCVODE_H_
#define _REACTORCVODE_H_

#include "ReactorBase.H"
#include "ReactorCvodeUtils.H"
#include "ReactorCvodePreconditioner.H"
#include "ReactorCvodeJacobian.H"
#include "ReactorCvodeCustomLinSolver.H"

namespace pele {
namespace physics {
namespace reactions {

struct CVODEUserData
{
  amrex::Real dt_save;     // Internal cvode dt holder
  int ncells;              // Number of cells in solve
  int verbose;             // Verbose
  int maxOrder;            // CVODE max order
  int reactor_type;        // Either HP (LM) or UV (C)
  int analytical_jacobian; // Analytical Jacobian 'On' flag
  int solve_type;          // Type of linear solve for Newton direction
  int precond_type;        // Type of preconditioner (if iterative solve)
  int NNZ; // Number of non-zero entry for sparse representations
  amrex::Real* rhoe_init = nullptr; // Initial energy (rhoE for C, rhoH for LM)
  amrex::Real* rhoesrc_ext = nullptr; // External energy forcing
  amrex::Real* rYsrc_ext = nullptr;   // External species forcing
  int* mask =
    nullptr; // Masking tagging cells where integration should not be perfomed
  int* FCunt = nullptr; // Number of RHS evaluations (not used on GPU)
  amrex::Real gamma;    // System Chem. jacobian coefficient
  int nbBlocks;         // GPU kernel launch parameter
  int nbThreads;        // GPU kernel launch parameter

#ifdef AMREX_USE_GPU
  // Matrix data
  int* csr_row_count_h;   // Host CSR sparse mat
  int* csr_col_index_h;   // Host CSR sparse mat
  int* csr_row_count_d;   // Device CSR sparse mat
  int* csr_col_index_d;   // Device CSR sparse mat
  amrex::Real* csr_val_d; // Syst. Jac. CSR sparse mat values
  amrex::Real* csr_jac_d; // Chem. Jac. CSR sparse mat values

#ifdef AMREX_USE_CUDA
  // Batched QR preconditioner data
  void* buffer_qr = NULL;    // cuSolver internal work arrays
  csrqrInfo_t info;          // QR solve info
  cusparseMatDescr_t descrA; // Sparse mat. descriptor
  cusolverSpHandle_t cusolverHandle;
  cusparseHandle_t cuSPHandle;
#endif
  amrex::gpuStream_t stream;

#else
  bool FirstTimePrecond;
  amrex::Real** Jdata = NULL;
  // Sparse direct solve data
  // Ptrs to CSR/CSC matrix for each submatrices (cells)
  SUNMatrix* PS;
  int** rowVals = NULL;
  int** rowPtrs = NULL;
  int** colPtrs = NULL;
  int** colVals = NULL;
  int* indx = NULL;
  amrex::Real** JSPSmat = NULL;
#ifdef PELE_USE_KLU
  klu_common* Common;
  klu_symbolic** Symbolic;
  klu_numeric** Numeric;
#endif
  // Matrix data for denseSimpleAjac preconditioner
  amrex::Real**** Jbd;
  amrex::Real**** P;
  sunindextype*** pivot;
  // Custom direct linear solve: a CSR SUNMatrix and ptrs
  SUNMatrix PSc;
  int* colVals_c;
  int* rowPtrs_c;
#endif
};

class ReactorCvode : public ReactorBase::Register<ReactorCvode>
{
public:
  static std::string identifier() { return "ReactorCvode"; }

  using Ordering = utils::YCOrder;

  int init(int reactor_type, int ncells) override;

  void close() override;

  int react(
    const amrex::Box& box,
    amrex::Array4<amrex::Real> const& rY_in,
    amrex::Array4<amrex::Real> const& rYsrc_in,
    amrex::Array4<amrex::Real> const& T_in,
    amrex::Array4<amrex::Real> const& rEner_in,
    amrex::Array4<amrex::Real> const& rEner_src_in,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Array4<int> const& mask,
    amrex::Real& dt_react,
    amrex::Real& time
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  int react(
    amrex::Real* rY_in,
    amrex::Real* rYsrc_in,
    amrex::Real* rX_in,
    amrex::Real* rX_src_in,
    amrex::Real& dt_react,
    amrex::Real& time,
    int ncells
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  static int
  cF_RHS(amrex::Real t, N_Vector y_in, N_Vector ydot, void* user_data);

  void flatten(
    const amrex::Box box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, rhoY, frcExt, temperature, rhoE, frcEExt,
      y_vect, src_vect, vect_energy, src_vect_energy);
  }

  void unflatten(
    const amrex::Box box,
    const int ncells,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
    long int* FCunt,
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, rhoY, temperature, rhoE, frcEExt, FC_in,
      y_vect, vect_energy, FCunt, dt);
  }

private:
  void checkCvodeOptions() const;
  void allocUserData(
    CVODEUserData* udata,
    int a_ncells
#ifdef AMREX_USE_GPU
    ,
    SUNMatrix& a_A,
    amrex::gpuStream_t stream
#endif
  ) const;
  void freeUserData(CVODEUserData* data_wk);

  amrex::Real relTol{1e-6};
  amrex::Real absTol{1e-10};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};

#ifndef AMREX_USE_GPU
  N_Vector y = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A = NULL;
  void* cvode_mem = NULL;
  CVODEUserData* udata_g = NULL;
#endif
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
