//
// AMReX Interface to Mersenne Twistor
//

/* A C-program for MT19937: Real number version (1999/10/28)       */
/*   genrand() generates one pseudorandom real number (double)     */
/* which is uniformly distributed on [0,1]-interval, for each      */
/* call. sgenrand(seed) sets initial values to the working area    */
/* of 624 words. Before genrand(), sgenrand(seed) must be          */
/* called once. (seed is any 32-bit integer.)                      */
/* Integer generator is obtained by modifying two lines.           */
/*   Coded by Takuji Nishimura, considering the suggestions by     */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.               */

/* This library is free software under the Artistic license:       */
/* see the file COPYING distributed together with this code.       */
/* For the verification of the code, its output sequence file      */
/* mt19937-1.out is attached (2001/4/2)                            */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp                                        */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

#ifndef MERSENNETWISTER_H_
#define MERSENNETWISTER_H_

#include <AMReX_Vector.H>

namespace pele::physics::turbforcing::mersenne_twister {
class mt19937
{
public:
  using seed_type = unsigned long;
  explicit mt19937(seed_type seed = 4357UL);
  mt19937(seed_type array[], int array_len);
  void rewind();
  void reset(unsigned long seed);

  double d_value();  // [0,1] random numbers
  double d1_value(); // [0,1) random numbers
  double d2_value(); // (0,1) random numbers

  long l_value();          // [0,2^31-1] random numbers
  unsigned long u_value(); // [0,2^32-1] random numbers

  void save(amrex::Vector<unsigned long>& state) const;
  int RNGstatesize() const;
  void restore(const amrex::Vector<unsigned long>& state);

private:
  enum { N = 624 };
  static unsigned long init_seed;
  static unsigned long mt[N]; // the array for the state vector
  static int mti;             // mti==N+1 means mt[N] is not initialized
#ifdef _OPENMP
#pragma omp threadprivate(init_seed, mt, mti)
#endif
  void sgenrand(unsigned long seed);
  unsigned long igenrand();
  void reload();
};

/*
  Generates one pseudorandom real number (double) which is
  uniformly distributed on [0,1]-interval for each call.

  Accepts any 32-bit integer as a seed -- uses 4357 as the default.

  Has a period of 2**19937.

  Mersenne Twister Home Page: http://www.math.keio.ac.jp/matumoto/emt.html

  There is also an entry point for Fortran callable as:

  REAL_T rn
  call blutilrand(rn)

  Internally, blutilrand() calls a static Mersenne Twister object (the
  same one used by AMReX::Random()) to get a value in [0,1] and then
  sets "rn" to that value.
*/
double Random();                           // [0,1]
double Random1();                          // [0,1)
double Random2();                          // (0,1)
unsigned long Random_int(unsigned long n); // [0,n-1], where n<=2^32-1
/* Set the seed of the random number generator.

  There is also an entry point for Fortran callable as:

  INTEGER seed
  call blutilinitrand(seed)

  or

  INTEGER seed
  call blinitrand(seed)
*/
void InitRandom(unsigned long seed);

void ResetRandomSeed(unsigned long seed);
//
// Save and restore random state.
//
// state.size() == 626 on return from Save & on entry to Restore.
//
void SaveRandomState(amrex::Vector<unsigned long>& state);

int sizeofRandomState();

void RestoreRandomState(const amrex::Vector<unsigned long>& state);
//
// Create a unique subset of random numbers from a pool
//   of integers in the range [0, poolSize - 1]
//   the set will be in the order they are found
//   setSize must be <= poolSize
//   uSet will be resized to setSize
//   if you want all processors to have the same set,
//   call this on one processor and broadcast the array
//
void UniqueRandomSubset(
  amrex::Vector<int>& uSet, int setSize, int poolSize, bool printSet = false);
} // namespace pele::physics::turbforcing::mersenne_twister

#endif
