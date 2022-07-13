
#ifndef _SOOTMODEL_H_
#define _SOOTMODEL_H_

// Standard libray includes
#include <array>

// AMReX include statements
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>

// Pele include statements
#include "PelePhysics.H"

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
    delete m_sootData;
    delete m_sootReact;
    amrex::The_Arena()->free(d_sootData);
    amrex::The_Arena()->free(d_sootReact);
  }

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
  void defineMemberData(const amrex::Real dimerVol);

  //
  // Compute HMOM source term
  //
  void computeSootSourceTerm(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate,
    amrex::Array4<const amrex::Real> const& coeff_mu,
    amrex::Array4<amrex::Real> const& soot_state,
    const amrex::Real time,
    const amrex::Real dt,
    const bool pres_term = true) const;

  //
  // Estimate the soot time step
  //
  amrex::Real estSootDt(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& Qstate) const;

  //
  // Access gas phase species name
  //
  std::string gasSpeciesName(const int sootGasIndx)
  {
    return m_gasSpecNames[sootGasIndx];
  }

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
  // Set the indices for primitive and conservative variables
  //
  inline void setIndices(const SootComps& sootIndx) { m_sootIndx = sootIndx; }

  //
  // Return string of the names of each variable
  //
  inline std::string sootVariableName(const int indx) const
  {
    AMREX_ASSERT(indx >= 0 && indx < NUM_SOOT_MOMENTS + 1);
    return m_sootVarName[indx];
  }

  //
  // Return pointer to soot data
  //
  inline SootData* getSootData() { return m_sootData; }
  inline SootData* getSootData_d() { return d_sootData; }

  /***********************************************************************
    Member data
  ***********************************************************************/

public:
  /// Units for moments are L^3(x + 2/3y)/L^3
  /// Units for the weight of the delta function are 1/L^3
  /// Where L is cm for PeleC and m for PeleLM

  /// Verbosity for soot model
  int m_sootVerbosity;
  /// Soot component indices
  SootComps m_sootIndx;
  bool m_setIndx;
  /// Pointer to soot data and soot specific functions
  SootData* m_sootData;
  /// Pointer for soot reaction data
  SootReaction* m_sootReact;
  /// Pointers for data on device
  SootData* d_sootData;
  SootReaction* d_sootReact;
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
  /// PAH inception path to be used, must be A2, A3, or A4
  std::string m_inceptPAH;
  /// Name of PAH used for inception
  std::string m_PAHname;
  amrex::Real m_gammaStick;

  // Name of soot variables to output (moments and weight of delta function)
  amrex::Vector<std::string> m_sootVarName;

  // Maximum allowable concentration for adjusting time step
  amrex::Real m_Xcutoff;
  // Temperature cutoff, any temparatures below this will not be solved
  amrex::Real m_Tcutoff;
  // Maximum possible subcycles
  int m_maxSubcycles;
  // Number of subcycles to use during source calculations
  int m_numSubcycles;

  /***********************************************************************
    Reaction member data
  ***********************************************************************/

  /// Bool to ensure reaction data has been read and filled
  bool m_reactDataFilled;
  /// Number and names of species for surface reactions
  amrex::Vector<std::string> m_gasSpecNames;

  /************************************************************************
    Factor member data (to reduce computational costs)
  ************************************************************************/

  /// Frequently used constants
  /// Factors for beta
  amrex::Real m_betaDimerFact;
  amrex::Real m_betaNuclFact;
};

#endif
