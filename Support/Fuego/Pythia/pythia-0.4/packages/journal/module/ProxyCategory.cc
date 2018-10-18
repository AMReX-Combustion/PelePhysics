// -*- C++ -*-
//
//--------------------------------------------------------------------------------
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
//
// <LicenseText>
//
//--------------------------------------------------------------------------------
//

#include <portinfo>

#include <Python.h>

#include "../lib/Category.h"
#include "ProxyCategory.h"

using namespace journal;

#if defined(JOURNAL_TRACE)
#include <iostream>
#endif

// interface
void journal::ProxyCategory::activate() {
    PyObject_SetAttrString(_category, "_state", Py_True);
    return;
}

void journal::ProxyCategory::deactivate() {
    PyObject_SetAttrString(_category, "_state", Py_False);
    return;
}

void journal::ProxyCategory::state(bool flag) {
    flag ? activate() : deactivate();
    return;
}

bool journal::ProxyCategory::state() const {
    PyObject * py_state = PyObject_GetAttrString(_category, "_state");
    bool flag = PyInt_AsLong(py_state);

#if defined(JOURNAL_TRACE)
    std::cout << " ** ProxyCategory: state=" << flag << std::endl;
#endif

    return flag;
}

// meta-methods
ProxyCategory::ProxyCategory(PyObject * category) :
    Category(),
    _category(category)
{
}

ProxyCategory::~ProxyCategory() {}

// version
// $Id$

// End of file
