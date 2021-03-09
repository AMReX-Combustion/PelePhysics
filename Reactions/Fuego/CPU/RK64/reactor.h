#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>
#include <AMReX_Print.H>
#include <PelePhysics.H>
/**********************************/

typedef struct 
{
    /* Checks */
    /* Base items */
    int ncells;
    int iverbose;
    int ireactor_type;
    double errtol;
    int nsubsteps_guess;
    int nsubsteps_min;
    int nsubsteps_max;
} *UserData;

/**********************************/
/* Functions Called by the Program */
int reactor_init(int reactor_type, int Ncells,double rk64_errtol=1e-16,
        int rk64_nsubsteps_guess=10,int rk64_nsusbteps_min=5,int rk64_nsubsteps_max=500);

int react(double *rY_in, double *rY_src_in, 
        double *rX_in, double *rX_src_in, 
        double &dt_react, double &time);

void reactor_close();

void FreeUserData(UserData data);

void fKernelSpec(double *dt, double *yvec_d, double *ydot_d,
        void *user_data);

UserData AllocUserData(int reactor_type, int num_cells,double rk64_errtol,
        int rk64_nsubsteps_guess,int rk64_nsubsteps_min,int rk64_nsubsteps_max);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);
