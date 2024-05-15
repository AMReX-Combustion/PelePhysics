#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>
#include "NeuralNetHomerolled.H"
#include "Table.H"

namespace pele::physics {
  // Need to instantiate these so Factory picks them up
  NNFuncParams2 dummy_1;
  TabFuncParams2 dummy_2;
}
