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

#include "journal.h"       // journal initialization
#include "misc.h"          // miscellaneous methods

// the method table

struct PyMethodDef pyjournal_methods[] = {

    // journal
    {pyjournal_initialize__name__, pyjournal_initialize,
     METH_VARARGS, pyjournal_initialize__doc__},

    // copyright
    {pyjournal_copyright__name__, pyjournal_copyright,
     METH_VARARGS, pyjournal_copyright__doc__},


// Sentinel
    {0, 0}
};

// version
// $Id$

// End of file
