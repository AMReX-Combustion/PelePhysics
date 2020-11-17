#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>
#include <AMReX_Print.H>
#include <EOS.H>
/**********************************/
using namespace amrex;
typedef struct 
{
    /* Checks */
    /* Base items */
    int ncells;
    int iverbose;
    int ireactor_type;
    Real errtol;
    int nsubsteps_guess;
    int nsubsteps_min;
    int nsubsteps_max;
} *UserData;

const int nstages_rk64=6;
const Real alpha_rk64[6] = 
{
    0.218150805229859,  //            3296351145737.0/15110423921029.0,
    0.256702469801519,  //            1879360555526.0/ 7321162733569.0,
    0.527402592007520,  //            10797097731880.0/20472212111779.0,
    0.0484864267224467, //            754636544611.0/15563872110659.0,
    1.24517071533530,   //            3260218886217.0/ 2618290685819.0,
    0.412366034843237,  //            5069185909380.0/12292927838509.0
};

const Real beta_rk64[6] = 
{
    -0.113554138044166,  //-1204558336989.0/10607789004752.0,
    -0.215118587818400,  //-3028468927040.0/14078136890693.0,
    -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
    -1.07992686223881,   //-17275898420483.0/15997285579755.0,
    -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
    0.0
};

const Real err_rk64[6] = 
{
    -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
    0.158481845574980,   // 473021958881.0/ 2984707536468.0,
    -0.0905918835751907, //-947229622805.0/10456009803779.0,
    -0.219084567203338,  //-2921473878215.0/13334914072261.0,
    0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
    0.0426421977505659   // 167623581683.0/ 3930932046784.0
};

/**********************************/
/* Functions Called by the Program */
extern "C"
{
    int reactor_init(const int* reactor_type, const int* Ncells,Real rk64_errtol=1e-16,
            int rk64_nsubsteps_guess=10,int rk64_nsusbteps_min=5,int rk64_nsubsteps_max=500);

    int react(Real *rY_in, Real *rY_src_in, 
            Real *rX_in, Real *rX_src_in, 
            Real &dt_react, Real &time);

    void reactor_close();
}

void FreeUserData(UserData data);

void fKernelSpec(Real *dt, Real *yvec_d, Real *ydot_d,
        void *user_data);

UserData AllocUserData(int reactor_type, int num_cells,Real rk64_errtol,
        int rk64_nsubsteps_guess,int rk64_nsubsteps_min,int rk64_nsubsteps_max);
