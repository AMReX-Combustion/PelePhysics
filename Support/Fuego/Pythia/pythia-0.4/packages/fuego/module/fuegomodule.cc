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

#include "exceptions.h"
#include "bindings.h"


char pyfuego_module__doc__[] = "";

// Initialization function for the module (*must* be called initfuego)
extern "C"
void
initfuego()
{
    // create the module and add the functions
    PyObject * m = Py_InitModule4(
        "fuego", pyfuego_methods,
        pyfuego_module__doc__, 0, PYTHON_API_VERSION);

    // get its dictionary
    PyObject * d = PyModule_GetDict(m);

    // check for errors
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module fuego");
    }

    // install the module exceptions
    pyfuego_runtimeError = PyErr_NewException("fuego.runtime", 0, 0);
    PyDict_SetItemString(d, "RuntimeException", pyfuego_runtimeError);

    return;
}

// version
// $Id$

// End of file
