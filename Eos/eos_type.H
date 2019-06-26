#ifndef _EOS_TYPE_H_
#define _EOS_TYPE_H_

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

extern "C" {

void eos_build(EOSTYPE* eos_state);

void eos_destroy(EOSTYPE* eos_state);

}
#endif 
