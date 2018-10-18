// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
// 
//  <LicenseText>
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include <Python.h>

#include <map>
#include <string>
#include <sstream>

#include "journal.h"

#include "../lib/Category.h"
#include "../lib/Diagnostic.h"

#include "../lib/Firewall.h"
#include "../lib/Debug.h"
#include "../lib/Info.h"
#include "../lib/Warning.h"
#include "../lib/Error.h"
#include "../lib/Device.h"
#include "../lib/Journal.h"
#include "../lib/Index.h"

#include "ProxyIndex.h"
#include "ProxyDevice.h"

#if defined(JOURNAL_TRACE)
#include <iostream>
#endif

using namespace journal;

// initialize
char pyjournal_initialize__doc__[] = "";
char pyjournal_initialize__name__[] = "initialize";

PyObject * pyjournal_initialize(PyObject *, PyObject * args)
{
    PyObject * py_journal;
    PyObject * py_firewall;
    PyObject * py_debug;
    PyObject * py_info;
    PyObject * py_warning;
    PyObject * py_error;

    int ok = PyArg_ParseTuple(
        args, "OOOOOO:initialize", 
        &py_journal, &py_firewall, &py_debug, &py_info, &py_warning, &py_error);

    if (!ok) {
        return 0;
    }

    Py_INCREF(py_journal);
    Py_INCREF(py_firewall);
    Py_INCREF(py_debug);
    Py_INCREF(py_info);
    Py_INCREF(py_warning);
    Py_INCREF(py_error);


    Device * device = new ProxyDevice(py_journal);
    Diagnostic::journal().device(device);

    // firewall
    Index * firewallIndex = new ProxyIndex(py_firewall);
    Firewall::newIndex(firewallIndex);
    
    // debug
#if defined(JOURNAL_TRACE)
    std::cout << " ** journal: initializing debug proxy index" << std::endl;
#endif
    Index * debugIndex = new ProxyIndex(py_debug);
    Debug::newIndex(debugIndex);
    
    // info
    Index * infoIndex = new ProxyIndex(py_info);
    Info::newIndex(infoIndex);

    // warning
    Index * warningIndex = new ProxyIndex(py_warning);
    Warning::newIndex(warningIndex);

    // error
    Index * errorIndex = new ProxyIndex(py_error);
    Error::newIndex(errorIndex);

    // return
    Py_INCREF(Py_None);
    return Py_None;
}

// $Id$

// End of file
