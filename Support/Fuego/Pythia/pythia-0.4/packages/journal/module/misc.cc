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

#include "misc.h"


// copyright

char pyjournal_copyright__doc__[] = "";
char pyjournal_copyright__name__[] = "copyright";

static char pyjournal_copyright_note[] = 
    "journal python module: Copyright (c) 1998-2003 California Institute of Technology";


PyObject * pyjournal_copyright(PyObject *, PyObject *)
{
    return Py_BuildValue("s", pyjournal_copyright_note);
}
    
// version
// $Id$

// End of file
