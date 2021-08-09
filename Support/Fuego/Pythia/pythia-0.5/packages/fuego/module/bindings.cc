// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Michael A.G. Aivazis
//                        California Institute of Technology
//                        (C) 1998-2003 All Rights Reserved
// 
//  <LicenseText>
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include <Python.h>

#include "bindings.h"

#include "misc.h"          // miscellaneous methods

// the method table

struct PyMethodDef pyfuego_methods[] = {

    // dummy entry for testing
    {pyfuego_hello__name__, pyfuego_hello,
     METH_VARARGS, pyfuego_hello__doc__},

    {pyfuego_copyright__name__, pyfuego_copyright,
     METH_VARARGS, pyfuego_copyright__doc__},


// Sentinel
    {0, 0}
};

// version
// $Id$

// End of file
