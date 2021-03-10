
// AMReX include statements
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

// PeleC-Physics include statements
#include <EOS.H>

// PeleC-MP include statements
#include <SootModel.H>

using namespace amrex;

void
SootModel::initializeReactData()
{
  SootConst sc;
  m_sootReact->SootDensityC = sc.SootDensityC;
  m_sootReact->SootChi = sc.SootChi;
  // Be sure that all other member data has been filled
  AMREX_ASSERT(m_memberDataDefined);
  if (m_sootVerbosity) {
    Print() << "SootModel::initializeReactData(): Filling reaction data"
            << std::endl;
  }
  // Units are CGS
  // Demonstration of reaction indices using a fake reaction
  // Soot-H + 2OH + C2H2 <=> Soot-* + 4H + 2CO
  // For the i-th reaction
  // rNum[i] = 2; // 2 gas phase reactants
  // nIndx_f[3*i+0] = GasSpecIndx::indxOH;
  // nu_f[3*i+0] = 2.; // Since there are 2 moles of OH
  // nIndx_f[3*i+1] = GasSpecIndx::indxC2H4;
  // nu_f[3*i+1] = 1; // Since there is 1 mole of C2H2
  // sIndx_f[i] = SootIndx::indxSootH;

  // pNum[i] = 2; // 2 gas phase products
  // nIndx_b[3*i+0] = GasSpecIndx::indxH;
  // nu_b[3*i+0] = 4.; // Since there are 4 moles of H
  // nIndx_b[3*i+1] = GasSpecIndx::indxCO;
  // nu_b[3*i+1] = 2.; // Since there are 2 moles of CO
  // sIndx_b[i] = SootIndx::indxSootS;

  // 1. Soot-H + OH <=> Soot-* + H2O
  m_sootReact->A_f[0] = 6.72E1;
  m_sootReact->n_f[0] = 3.33;
  m_sootReact->ER_f[0] = 6.09E10 / pele::physics::Constants::RU;
  m_sootReact->rNum[0] = 1;
  m_sootReact->nIndx_f[3 * 0 + 0] = SootGasSpecIndx::indxOH;
  m_sootReact->nu_f[0 * 3 + 0] = 1.;
  m_sootReact->sIndx_f[0] = SootIndx::indxSootH;

  m_sootReact->A_b[0] = 6.44E-1;
  m_sootReact->n_b[0] = 3.79;
  m_sootReact->ER_b[0] = 27.96E10 / pele::physics::Constants::RU;
  m_sootReact->pNum[0] = 1;
  m_sootReact->nIndx_b[0 * 3 + 0] = SootGasSpecIndx::indxH2O;
  m_sootReact->nu_b[0 * 3 + 0] = 1.;
  m_sootReact->sIndx_b[0] = SootIndx::indxSootS;

  // 2. Soot-H + H <=> Soot-* + H2
  m_sootReact->A_f[1] = 1.0E8;
  m_sootReact->n_f[1] = 1.80;
  m_sootReact->ER_f[1] = 68.42E10 / pele::physics::Constants::RU;
  m_sootReact->rNum[1] = 1;
  m_sootReact->nIndx_f[1 * 3 + 0] = SootGasSpecIndx::indxH;
  m_sootReact->nu_f[1 * 3 + 0] = 1.;
  m_sootReact->sIndx_f[1] = SootIndx::indxSootH;

  m_sootReact->A_b[1] = 8.68E4;
  m_sootReact->n_b[1] = 2.36;
  m_sootReact->ER_b[1] = 25.46E10 / pele::physics::Constants::RU;
  m_sootReact->pNum[1] = 1;
  m_sootReact->nIndx_b[1 * 3 + 0] = SootGasSpecIndx::indxH2;
  m_sootReact->nu_b[1 * 3 + 0] = 1.;
  m_sootReact->sIndx_b[1] = SootIndx::indxSootS;

  // 3. Soot-H <=> Soot-* + H
  m_sootReact->A_f[2] = 1.13E16;
  m_sootReact->n_f[2] = -0.06;
  m_sootReact->ER_f[2] = 476.05E10 / pele::physics::Constants::RU;
  m_sootReact->rNum[2] = 0;
  m_sootReact->sIndx_f[2] = SootIndx::indxSootH;

  m_sootReact->A_b[2] = 4.17E13;
  m_sootReact->n_b[2] = 0.15;
  m_sootReact->ER_b[2] = 0.;
  m_sootReact->pNum[2] = 1;
  m_sootReact->nIndx_b[2 * 3 + 0] = SootGasSpecIndx::indxH;
  m_sootReact->nu_b[2 * 3 + 0] = 1.;
  m_sootReact->sIndx_b[2] = SootIndx::indxSootS;

  // 4. Soot-* + C2H2 => Soot-H
  m_sootReact->A_f[3] = 2.52E9;
  m_sootReact->n_f[3] = 1.10;
  m_sootReact->ER_f[3] = 17.13E10 / pele::physics::Constants::RU;
  m_sootReact->rNum[3] = 1;
  m_sootReact->nIndx_f[3 * 3 + 0] = SootGasSpecIndx::indxC2H2;
  m_sootReact->nu_f[3 * 3 + 0] = 1.;
  m_sootReact->sIndx_f[3] = SootIndx::indxSootS;

  m_sootReact->sIndx_b[3] = SootIndx::indxSootH;

  // 5. Soot-* + O2 => Soot-* + 2CO
  m_sootReact->A_f[4] = 2.20E12;
  m_sootReact->n_f[4] = 0.;
  m_sootReact->ER_f[4] = 31.38E10 / pele::physics::Constants::RU;
  m_sootReact->rNum[4] = 1;
  m_sootReact->nIndx_f[4 * 3 + 0] = SootGasSpecIndx::indxO2;
  m_sootReact->nu_f[4 * 3 + 0] = 1.;
  m_sootReact->sIndx_f[4] = SootIndx::indxSootS;

  m_sootReact->pNum[4] = 1;
  m_sootReact->nIndx_b[4 * 3 + 0] = SootGasSpecIndx::indxCO;
  m_sootReact->nu_b[4 * 3 + 0] = 2.;
  m_sootReact->sIndx_b[4] = SootIndx::indxSootS;

  // 6. Soot-H + OH => Soot-H + CO
  // TODO: This transforms the Arrhenius formulation to be
  // reaction probability, 8.94*sqrt(T)*probGamma*A
  Real probGamma = 0.13;
  // FIXME: Find out what the units for this are
  m_sootReact->A_f[5] =
    8.94 * probGamma * pele::physics::Constants::Avna * 100.;
  m_sootReact->n_f[5] = 0.5;
  m_sootReact->ER_f[5] = 0.;
  m_sootReact->rNum[5] = 1;
  m_sootReact->nIndx_f[5 * 3 + 0] = SootGasSpecIndx::indxOH;
  m_sootReact->nu_f[5 * 3 + 0] = 1.;
  m_sootReact->sIndx_f[5] = SootIndx::indxSootH;

  m_sootReact->A_b[5] = 0.;
  m_sootReact->n_b[5] = 0.;
  m_sootReact->ER_b[5] = 0.;
  m_sootReact->pNum[5] = 1;
  m_sootReact->nIndx_b[5 * 3 + 0] = SootGasSpecIndx::indxCO;
  m_sootReact->nu_b[5 * 3 + 0] = 1.;
  m_sootReact->sIndx_b[5] = SootIndx::indxSootH;

  // Last reaction MUST be the dimerization reaction
  // 7. A# + A# => DIMER
  // TODO: Makes use of Arrhenius form similar to last reaction
  m_sootReact->A_f[6] =
    m_betaDimerFact * std::sqrt(sc.colFact) * 0.5 * m_gammaStick;
  m_sootReact->n_f[6] = 0.5;
  m_sootReact->ER_f[6] = 0.;
  m_sootReact->rNum[6] = 1;
  m_sootReact->nIndx_f[6 * 3 + 0] = SootGasSpecIndx::indxPAH;
  m_sootReact->nu_f[6 * 3 + 0] = 2.;
  m_sootReact->sIndx_f[6] = -1; // No soot in reactants

  m_reactDataFilled = true;
}
