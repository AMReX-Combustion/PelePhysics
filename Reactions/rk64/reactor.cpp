#include "reactor.H"

UserData data = NULL;
int eint_rho = 1;
int enth_rho = 2;
amrex::Real* rhoX_init = NULL;
amrex::Real* rhoXsrc_ext = NULL;
amrex::Real* rYsrc = NULL;
const int nstages_rk64 = 6;
const amrex::Real alpha_rk64[6] = {
  0.218150805229859,  // 3296351145737.0/15110423921029.0,
  0.256702469801519,  // 1879360555526.0/ 7321162733569.0,
  0.527402592007520,  // 10797097731880.0/20472212111779.0,
  0.0484864267224467, // 754636544611.0/15563872110659.0,
  1.24517071533530,   // 3260218886217.0/ 2618290685819.0,
  0.412366034843237,  // 5069185909380.0/12292927838509.0
};

const amrex::Real beta_rk64[6] = {
  -0.113554138044166,  //-1204558336989.0/10607789004752.0,
  -0.215118587818400,  //-3028468927040.0/14078136890693.0,
  -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
  -1.07992686223881,   //-17275898420483.0/15997285579755.0,
  -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
  0.0};

const amrex::Real err_rk64[6] = {
  -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
  0.158481845574980,   // 473021958881.0/ 2984707536468.0,
  -0.0905918835751907, //-947229622805.0/10456009803779.0,
  -0.219084567203338,  //-2921473878215.0/13334914072261.0,
  0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
  0.0426421977505659   // 167623581683.0/ 3930932046784.0
};

#ifdef _OPENMP
#pragma omp threadprivate(data)
#pragma omp threadprivate(rhoX_init, rhoXsrc_ext, rYsrc)
#endif

int
reactor_init(
  int reactor_type,
  int Ncells,
  amrex::Real rk64_errtol,
  int rk64_nsubsteps_guess,
  int rk64_nsubsteps_min,
  int rk64_nsubsteps_max)
{

#ifdef _OPENMP
  int omp_thread;
  omp_thread = omp_get_thread_num();
#endif
  data = AllocUserData(
    reactor_type, Ncells, rk64_errtol, rk64_nsubsteps_guess, rk64_nsubsteps_min,
    rk64_nsubsteps_max);

  rhoX_init = (amrex::Real*)malloc(data->ncells * sizeof(amrex::Real));
  rhoXsrc_ext = (amrex::Real*)malloc(data->ncells * sizeof(amrex::Real));
  rYsrc =
    (amrex::Real*)malloc((data->ncells * NUM_SPECIES) * sizeof(amrex::Real));

  return (0);
}

void
SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals)
{
  // WIP
}

void
SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol)
{
}

UserData
AllocUserData(
  int reactor_type,
  int num_cells,
  amrex::Real rk64_errtol,
  int rk64_nsubsteps_guess,
  int rk64_nsubsteps_min,
  int rk64_nsubsteps_max)
{
  UserData data_wk;
  data_wk = (UserData)malloc(sizeof *data_wk);
#ifdef _OPENMP
  int omp_thread;
  omp_thread = omp_get_thread_num();
#endif

  (data_wk->ncells) = num_cells;
  (data_wk->iverbose) = 1;
  (data_wk->ireactor_type) = reactor_type;

  (data_wk->errtol) = rk64_errtol;
  (data_wk->nsubsteps_guess) = rk64_nsubsteps_guess;
  (data_wk->nsubsteps_min) = rk64_nsubsteps_min;
  (data_wk->nsubsteps_max) = rk64_nsubsteps_max;

  return (data_wk);
}

int
react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time)
{

  amrex::Real time_init = time;
  amrex::Real time_out = time + dt_react;
  amrex::Real current_time = time_init;
  amrex::Real *soln_reg, *carryover_reg, *error_reg, *zero_reg, *rhs;
  amrex::Real dt_rk, dt_rk_min, dt_rk_max, change_factor;
  const amrex::Real exp1 = 0.25;
  const amrex::Real exp2 = 0.2;
  const amrex::Real beta = 1.0;
  const amrex::Real tinyval = 1e-50;
  int neq_tot = (NUM_SPECIES + 1) * (data->ncells);

#ifdef _OPENMP
  int omp_thread;
  omp_thread = omp_get_thread_num();
#endif

#ifdef _OPENMP
  if ((data->iverbose > 1) && (omp_thread == 0)) {
#else
  if (data->iverbose > 1) {
#endif
    amrex::Print() << "\n -------------------------------------\n";
  }

#ifdef _OPENMP
  if ((data->iverbose > 3) && (omp_thread == 0)) {
#else
  if (data->iverbose > 3) {
#endif
    amrex::Print() << "BEG : time curr is " << time_init << " and dt_react is "
                   << dt_react << " and final time should be " << time_out
                   << "\n";
  }

  soln_reg = (amrex::Real*)calloc(neq_tot, sizeof(amrex::Real));
  carryover_reg = (amrex::Real*)calloc(neq_tot, sizeof(amrex::Real));
  error_reg = (amrex::Real*)calloc(neq_tot, sizeof(amrex::Real));
  zero_reg = (amrex::Real*)calloc(neq_tot, sizeof(amrex::Real));
  rhs = (amrex::Real*)calloc(neq_tot, sizeof(amrex::Real));

  std::memcpy(soln_reg, rY_in, sizeof(amrex::Real) * neq_tot);
  std::memcpy(carryover_reg, rY_in, sizeof(amrex::Real) * neq_tot);
  std::memcpy(
    rYsrc, rY_src_in, (NUM_SPECIES * data->ncells) * sizeof(amrex::Real));
  std::memcpy(rhoX_init, rX_in, sizeof(amrex::Real) * data->ncells);
  std::memcpy(rhoXsrc_ext, rX_src_in, sizeof(amrex::Real) * data->ncells);

  dt_rk = dt_react / amrex::Real(data->nsubsteps_guess);
  dt_rk_min = dt_react / amrex::Real(data->nsubsteps_max);
  dt_rk_max = dt_react / amrex::Real(data->nsubsteps_min);

  int nsteps = 0;

  while (current_time < time_out) {
    current_time += dt_rk;
    nsteps++;
    std::memcpy(carryover_reg, soln_reg, sizeof(amrex::Real) * neq_tot);
    std::memcpy(error_reg, zero_reg, sizeof(amrex::Real) * neq_tot);

    for (int stage = 0; stage < nstages_rk64; stage++) {
      std::memcpy(rhs, zero_reg, sizeof(amrex::Real) * neq_tot);
      fKernelSpec(&current_time, soln_reg, rhs, data);

      for (int i = 0; i < neq_tot; i++) {
        error_reg[i] += err_rk64[stage] * dt_rk * rhs[i];
        soln_reg[i] = carryover_reg[i] + alpha_rk64[stage] * dt_rk * rhs[i];
        carryover_reg[i] = soln_reg[i] + beta_rk64[stage] * dt_rk * rhs[i];
      }
    }

    amrex::Real max_err = tinyval;
    for (int i = 0; i < neq_tot; i++) {
      if (fabs(error_reg[i]) > max_err) {
        max_err = fabs(error_reg[i]);
      }
    }

    if (max_err < data->errtol) {
      change_factor = beta * pow((data->errtol / max_err), exp1);
      dt_rk = std::min(dt_rk_max, dt_rk * change_factor);
    } else {
      change_factor = beta * pow((data->errtol / max_err), exp2);
      dt_rk = std::max(dt_rk_min, dt_rk * change_factor);
    }
  }

  // update guess from the current update
  (data->nsubsteps_guess) = nsteps;

#ifdef MOD_REACTOR
  time = time_init + dt_react;
#endif

#ifdef _OPENMP
  if ((data->iverbose > 3) && (omp_thread == 0)) {
#else
  if (data->iverbose > 3) {
#endif
    amrex::Print() << "END : time curr is " << current_time
                   << " and actual dt_react is " << dt_react << "\n";
  }

  std::memcpy(rY_in, soln_reg, neq_tot * sizeof(amrex::Real));
  for (int i = 0; i < data->ncells; i++) {
    rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
  }

  return nsteps;
}

void
fKernelSpec(
  amrex::Real* dt, amrex::Real* yvec_d, amrex::Real* ydot_d, void* user_data)
{
  UserData data_wk;
  data_wk = (UserData)user_data;
  int tid;

  for (tid = 0; tid < data_wk->ncells; tid++) {
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real Xi[NUM_SPECIES];
    amrex::Real cdot[NUM_SPECIES], molecular_weight[NUM_SPECIES];
    amrex::Real cX;
    amrex::Real temp, energy;

    int offset = tid * (NUM_SPECIES + 1);

    CKWT(molecular_weight);

    amrex::Real rho = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      rho = rho + yvec_d[offset + i];
    }

    temp = yvec_d[offset + NUM_SPECIES];

    for (int i = 0; i < NUM_SPECIES; i++) {
      massfrac[i] = yvec_d[offset + i] / rho;
    }

    energy = (rhoX_init[tid] + rhoXsrc_ext[tid] * (*dt)) / rho;

    auto eos = pele::physics::PhysicsType::eos();
    if (data_wk->ireactor_type == eint_rho) {
      eos.EY2T(energy, massfrac, temp);
      eos.TY2Cv(temp, massfrac, cX);
      eos.T2Ei(temp, Xi);
    } else {
      eos.HY2T(energy, massfrac, temp);
      eos.TY2Cp(temp, massfrac, cX);
      eos.T2Hi(temp, Xi);
    }
    eos.RTY2WDOT(rho, temp, massfrac, cdot);

    ydot_d[offset + NUM_SPECIES] = rhoXsrc_ext[tid];
    for (int i = 0; i < NUM_SPECIES; i++) {
      ydot_d[offset + i] = cdot[i] + rYsrc[tid * (NUM_SPECIES) + i];
      ydot_d[offset + NUM_SPECIES] =
        ydot_d[offset + NUM_SPECIES] - ydot_d[offset + i] * Xi[i];
    }
    ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] / (rho * cX);
  }
}

void
reactor_close()
{
  FreeUserData(data);
  free(rhoX_init);
  free(rhoXsrc_ext);
  free(rYsrc);
}

void
FreeUserData(UserData data_wk)
{
  free(data_wk);
}
