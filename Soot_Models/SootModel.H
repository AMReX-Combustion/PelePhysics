
#ifndef _SOOTMODEL_H_
#define _SOOTMODEL_H_

// Standard libray includes
#include <array>

// AMReX include statements
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>
#include <AMReX_Derive.H>

// PeleC and PelePhysics include statements
#include "EOS.H"
#include "IndexDefines.H"

#include "Constants_Soot.H"
#include "SootData.H"
#include "SootReactions.H"

class SootModel
{
public:
  //
  // Default constructor
  //
  SootModel();

  //
  // Destructor
  //
  ~SootModel()
  {
    m_sootData.reset();
    m_sootReact.reset();
  };

  //
  // Define instance of class
  //
  void define();

  //
  // Read soot parameters from input file
  //
  void readSootParams();

  //
  // Define member data
  //
  void defineMemberData(const Real dimerVol);

  //
  // Define the derived variable list (SootModel_setup.cpp)
  //
  void
  addSootDerivePlotVars(DeriveList& derive_lst, const DescriptorList& desc_lst);

  //
  // Compute HMOM source term
  //
  void addSootSourceTerm(
    const Box& vbox,
    Array4<const Real> const& Qstate,
    Array4<const Real> const& coeff_state,
    Array4<Real> const& soot_state,
    const Real time,
    const Real dt) const;

  /***********************************************************************
    Reaction functions (SootModel_react.cpp)
  ***********************************************************************/
  //
  // Initialize reaction and species data
  //
  void initializeReactData();

  /***********************************************************************
    Inline functions
  ***********************************************************************/
  //
  // Assign the member data for the time step
  //
  inline void setTimeStep(const Real soot_dt)
  {
    local_soot_dt = amrex::min(local_soot_dt, m_maxDtRate * soot_dt);
  }

  //
  // Return the estimated soot time step and reset the member variable
  //
  inline Real estSootTimeStep()
  {
    Real cur_dt = local_soot_dt;
    // Reset time step
    local_soot_dt = 1.E18;
    return cur_dt;
  }

  //
  // Return string of the names of each variable
  //
  inline std::string sootVariableName(const int indx) const
  {
    AMREX_ASSERT(indx >= 0 && indx < NUM_SOOT_VARS);
    return m_sootVarName[indx];
  }

  //
  // Return pointer to soot data
  //
  inline SootData* getSootData() { return m_sootData.get(); }

  /***********************************************************************
    Member data
  ***********************************************************************/

public:
  /// Units for moments are L^3(x + 2/3y)/L^3
  /// Units for the weight of the delta function are 1/L^3
  /// Where L is cm for PeleC and m for PeleLM

  /// Verbosity for soot model
  int m_sootVerbosity;
  /// Pointer to soot data and soot specific functions
  std::unique_ptr<SootData> m_sootData;
  /// Pointer for soot reaction data
  std::unique_ptr<SootReaction> m_sootReact;
  /// Bool for checking if parameters have been read from input file
  bool m_readSootParams;
  /// Bool for checking if member data has all been defined
  bool m_memberDataDefined;
  /// Mass is not conserved during OH oxidation
  /// Mass can be conserved again by adding remaining lost mass to H2
  bool m_conserveMass;

  // TODO: Make these vectors in case multiple inception species are used
  /// Species index of PAH
  int m_PAHindx;
  /// Name of PAH used for inception
  std::string m_PAHname;
  Real m_gammaStick;

  // Name of soot variables to output (moments and weight of delta function)
  Vector<std::string> m_sootVarName;

  // Time step estimation
  Real local_soot_dt;
  // Maximum allowable change for mass fractions, density, or moments
  // to adjust the time step
  Real m_maxDtRate;

  /***********************************************************************
    Reaction member data
  ***********************************************************************/

  /// Bool to ensure reaction data has been read and filled
  bool m_reactDataFilled;
  /// Number and names of species for surface reactions
  Vector<std::string> m_gasSpecNames;

  /************************************************************************
    Factor member data (to reduce computational costs)
  ************************************************************************/

  /// Frequently used constants
  /// Factors for beta
  Real m_betaDimerFact;
  Real m_betaNuclFact;
};

#endif
