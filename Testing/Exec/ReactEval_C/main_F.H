#ifndef _main_F_H_
#define _main_F_H_

#include <AMReX_BLFort.H>

extern "C"
{
    void extern_init
    (const int* name, const int* namlen, 
     const int* fuelID,
     const int* oxyID,
     const int* bathID,
     const int* cvode_iE);

    void extern_close();
}

#endif // _main_F_H_
