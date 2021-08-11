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


char pyjournal_module__doc__[] = "";

// Initialization function for the module (*must* be called initjournal)
extern "C"
void
init_journal()
{
    // create the module and add the functions
    PyObject * m = Py_InitModule4(
        "_journal", pyjournal_methods,
        pyjournal_module__doc__, 0, PYTHON_API_VERSION);

    // get its dictionary
    PyObject * d = PyModule_GetDict(m);

    // check for errors
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module journal");
    }

    // install the module exceptions
    pyjournal_runtimeError = PyErr_NewException("journal.runtime", 0, 0);
    PyDict_SetItemString(d, "RuntimeException", pyjournal_runtimeError);

    return;
}

// version
// $Id$

// End of file
