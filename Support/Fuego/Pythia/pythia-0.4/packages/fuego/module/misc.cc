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
#include "src/hello.h"


// copyright

char pyfuego_copyright__doc__[] = "";
char pyfuego_copyright__name__[] = "copyright";

static char pyfuego_copyright_note[] = 
    "fuego python module: Copyright (c) 1998-2003 California Institute of Technology";


PyObject * pyfuego_copyright(PyObject *, PyObject *)
{
    return Py_BuildValue("s", pyfuego_copyright_note);
}
    
// hello

char pyfuego_hello__doc__[] = "";
char pyfuego_hello__name__[] = "hello";

PyObject * pyfuego_hello(PyObject *, PyObject *)
{
    return Py_BuildValue("s", hello());
}
    
// version
// $Id$

// End of file
