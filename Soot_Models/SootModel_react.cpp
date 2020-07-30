
// AMReX include statements
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

// PeleC include statements
#include <PeleC.H>

// PeleC-Physics include statements
#include <EOS.H>

// PeleC-MP include statements
#include <SootModel.H>

using namespace amrex;

void
SootModel::initializeReactData()
{
  using HVReal = Gpu::HostVector<Real>;
  using HVInt = Gpu::HostVector<int>;
  // Be sure that all other member data has been filled
  AMREX_ASSERT(m_memberDataDefined);
  if (m_sootVerbosity) {
    Print() << "SootModel::initializeReactData(): Filling reaction data"
            << std::endl;
  }
  const int nsr = NUM_SOOT_REACT;
  // Vectors of the forward and backward per-exponential rates,
  // temperature exponents, and activation energies/R
  HVReal A_f(nsr, 0.);
  HVReal n_f(nsr, 0.);
  HVReal ER_f(nsr, 0.);
  HVReal A_b(nsr, 0.);
  HVReal n_b(nsr, 0.);
  HVReal ER_b(nsr, 0.);
  // Vector of number of reactants and products
  HVInt rNum(nsr, 0);
  HVInt pNum(nsr, 0);
  // Vector of species indices
  // Maximum 3 species per reaction side
  HVInt nIndx_f(3*nsr, -1);
  HVInt nIndx_b(3*nsr, -1);
  // Vector of stoichiometric coefficients
  HVReal nu_f(3*nsr, 0.);
  HVReal nu_b(3*nsr, 0.);
  /// Soot surface reaction reference indexing
  /// Forward and backward reaction soot indices
  // TODO: Assumes only 1 mole of soot is on each side
  HVInt sIndx_f(nsr, -1);
  HVInt sIndx_b(nsr, -1);
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
  A_f[0] = 6.72E1;
  n_f[0] = 3.33;
  ER_f[0] = 6.09E10/EOS::RU;
  rNum[0] = 1;
  nIndx_f[3*0+0] = SootConst::GasSpecIndx::indxOH;
  nu_f[0*3+0] = 1.;
  sIndx_f[0] = SootConst::SootIndx::indxSootH;

  A_b[0] = 6.44E-1;
  n_b[0] = 3.79;
  ER_b[0] = 27.96E10/EOS::RU;
  pNum[0] = 1;
  nIndx_b[0*3+0] = SootConst::GasSpecIndx::indxH2O;
  nu_b[0*3+0] = 1.;
  sIndx_b[0] = SootConst::SootIndx::indxSootS;

  // 2. Soot-H + H <=> Soot-* + H2
  A_f[1] = 1.0E8;
  n_f[1] = 1.80;
  ER_f[1] = 68.42E10/EOS::RU;
  rNum[1] = 1;
  nIndx_f[1*3+0] = SootConst::GasSpecIndx::indxH;
  nu_f[1*3+0] = 1.;
  sIndx_f[1] = SootConst::SootIndx::indxSootH;

  A_b[1] = 8.68E4;
  n_b[1] = 2.36;
  ER_b[1] = 25.46E10/EOS::RU;
  pNum[1] = 1;
  nIndx_b[1*3+0] = SootConst::GasSpecIndx::indxH2;
  nu_b[1*3+0] = 1.;
  sIndx_b[1] = SootConst::SootIndx::indxSootS;

  // 3. Soot-H <=> Soot-* + H
  A_f[2] = 1.13E16;
  n_f[2] = -0.06;
  ER_f[2] = 476.05E10/EOS::RU;
  rNum[2] = 0;
  sIndx_f[2] = SootConst::SootIndx::indxSootH;

  A_b[2] = 4.17E13;
  n_b[2] = 0.15;
  ER_b[2] = 0.;
  pNum[2] = 1;
  nIndx_b[2*3+0] = SootConst::GasSpecIndx::indxH;
  nu_b[2*3+0] = 1.;
  sIndx_b[2] = SootConst::SootIndx::indxSootS;

  // 4. Soot-* + C2H2 => Soot-H
  A_f[3] = 2.52E9;
  n_f[3] = 1.10;
  ER_f[3] = 17.13E10/EOS::RU;
  rNum[3] = 1;
  nIndx_f[3*3+0] = SootConst::GasSpecIndx::indxC2H2;
  nu_f[3*3+0] = 1.;
  sIndx_f[3] = SootConst::SootIndx::indxSootS;

  sIndx_b[3] = SootConst::SootIndx::indxSootH;

  // 5. Soot-* + O2 => Soot-* + 2CO
  A_f[4] = 2.20E12;
  n_f[4] = 0.;
  ER_f[4] = 31.38E10/EOS::RU;
  rNum[4] = 1;
  nIndx_f[4*3+0] = SootConst::GasSpecIndx::indxO2;
  nu_f[4*3+0] = 1.;
  sIndx_f[4] = SootConst::SootIndx::indxSootS;

  pNum[4] = 1;
  nIndx_b[4*3+0] = SootConst::GasSpecIndx::indxCO;
  nu_b[4*3+0] = 2.;
  sIndx_b[4] = SootConst::SootIndx::indxSootS;

  // 6. Soot-H + OH => Soot-H + CO
  // TODO: This transforms the Arrhenius formulation to be
  // reaction probability, 8.94*sqrt(T)*probGamma*A
  Real probGamma = 0.13;
  // FIXME: Find out what the units for this are
  A_f[5] = 8.94*probGamma*SootConst::avogadros*100.;
  n_f[5] = 0.5;
  ER_f[5] = 0.;
  rNum[5] = 1;
  nIndx_f[5*3+0] = SootConst::GasSpecIndx::indxOH;
  nu_f[5*3+0] = 1.;
  sIndx_f[5] = SootConst::SootIndx::indxSootH;

  A_b[5] = 0.;
  n_b[5] = 0.;
  ER_b[5] = 0.;
  pNum[5] = 1;
  nIndx_b[5*3+0] = SootConst::GasSpecIndx::indxCO;
  nu_b[5*3+0] = 1.;
  sIndx_b[5] = SootConst::SootIndx::indxSootH;

  // Last reaction MUST be the dimerization reaction
  // 7. A# + A# => DIMER
  // TODO: Makes use of Arrhenius form similar to last reaction
  A_f[6] = m_betaDimerFact*std::sqrt(SootConst::colFact)*0.5*m_gammaStick;
  n_f[6] = 0.5;
  ER_f[6] = 0.;
  rNum[6] = 1;
  nIndx_f[6*3+0] = SootConst::GasSpecIndx::indxPAH;
  nu_f[6*3+0] = 2.;
  sIndx_f[6] = -1; // No soot in reactants
  m_SootReactionContainer.build(A_f, n_f, ER_f, A_b, n_b, ER_b, rNum, pNum,
                                sIndx_f, sIndx_b, nIndx_f, nIndx_b, nu_f, nu_b);

  m_reactDataFilled = true;
}
