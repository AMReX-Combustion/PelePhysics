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
#include <vector>
#include <string>

#include "../lib/Device.h"
#include "../lib/Entry.h"

#include "ProxyDevice.h"

using namespace journal;

// #include <iostream>

// interface
void ProxyDevice::record(const Entry & entry) {
    // std::cout << " ## ProxyDevice::record" << std::endl;

    // size_t lines = entry.lines();
    // std::cout << " ## ProxyDevice::record: entry with " << lines << " lines" << std::endl;

    // std::cout << " ## ProxyDevice::record: calling Journal.entry" << std::endl;
    PyObject * py_entry = PyObject_CallMethod(_journal, "entry", "");

    // check for errors
    if (PyErr_Occurred()) {
        PyErr_Print();
        Py_FatalError("FIREWALL: the journal interface has changed");
    }

    for (Entry::page_t::const_iterator line = entry.lineBegin(); line != entry.lineEnd(); ++line) {
        // std::cout 
            // << " ## ProxyDevice::record: converting line: {"
            // << (*line) << "}"
            // << std::endl;
        PyObject * py_line = PyString_FromStringAndSize((*line).c_str(), (*line).size());
        PyObject * py_none = PyObject_CallMethod(py_entry, "line", "O", py_line);
        Py_DECREF(py_none);
        Py_DECREF(py_line);
    }

    // std::cout << " ## ProxyDevice::record: extracting the attribute 'meta'" << std::endl;
    PyObject * py_meta = PyObject_GetAttrString(py_entry, "meta");
    // check for errors
    if (PyErr_Occurred()) {
        PyErr_Print();
        Py_FatalError("FIREWALL: the journal interface has changed");
    }

    // std::cout << " ## ProxyDevice::record: recording the entry attributes:" << std::endl;
    // std::cout << "    ";
    for (Entry::meta_t::const_iterator meta = entry.metaBegin(); meta != entry.metaEnd(); ++meta) {
        Entry::string_t key = (*meta).first;
        Entry::string_t value = (*meta).second;

        // std::cout << "{" << key << "=" << value << "}";
        PyObject * py_key = PyString_FromStringAndSize(key.c_str(), key.size());
        PyObject * py_value = PyString_FromStringAndSize(value.c_str(), value.size());

        PyDict_SetItem(py_meta, py_key, py_value);
    }
    // std::cout << std::endl;

    // std::cout << " ## ProxyDevice::record: calling Journal.record" << std::endl;
    PyObject * py_value = PyObject_CallMethod(_journal, "record", "O", py_entry);

    // check for errors
    if (PyErr_Occurred()) {
        PyErr_Print();
        Py_FatalError("FIREWALL: the journal interface has changed");
    }

    // clean up
    Py_DECREF(py_value);
    Py_DECREF(py_meta);
    Py_DECREF(py_entry);

    return;
}

// meta-methods
ProxyDevice::ProxyDevice(PyObject * journal) :
    Device(),
    _journal(journal)
{}


ProxyDevice::~ProxyDevice()
{}

// version
// $Id$

// End of file
