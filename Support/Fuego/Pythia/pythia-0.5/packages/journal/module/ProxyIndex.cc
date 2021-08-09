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
#include <map>
#include <string>

#include <Python.h>

#include "../lib/Index.h"
#include "ProxyIndex.h"

#include "../lib/Category.h"
#include "ProxyCategory.h"

#if defined(JOURNAL_TRACE)
#include <iostream>
#endif

using namespace journal;

// interface
Index::category_t & ProxyIndex::category(string_t name) {

#if defined(JOURNAL_TRACE)
    std::cout << " ** ProxyIndex: request for category '" << name << "'" << std::endl;
#endif

    index_t::const_iterator target = _cache.find(name);

    // is a proxy in the local cache?
    if (target != _cache.end()) {
#if defined(JOURNAL_TRACE)
        std::cout << " ** ProxyIndex: category '" << name << "' already in the cache" << std::endl;
#endif
        category_t * category = (*target).second;
        return *category;
    }

    PyObject * py_cat = PyObject_CallMethod(_index, "diagnostic", "s", name.c_str());
    // check for errors
    if (PyErr_Occurred()) {
        PyErr_Print();
        Py_FatalError("FIREWALL: the journal interface has changed");
    }

#if defined(JOURNAL_TRACE)
    std::cout << " ** ProxyIndex: inserting category '" << name << "' in the cache" << std::endl;
#endif
    category_t * category = new ProxyCategory(py_cat);
    _cache.insert(entry_t(name, category));

#if defined(JOURNAL_TRACE)
    std::cout << " ** ProxyIndex: category '" << name << "': " << category->state() << std::endl;
#endif

    return *category;
}

// meta-methods
ProxyIndex::ProxyIndex(PyObject * index) :
    Index(),
    _cache(),
    _index(index)
{}

ProxyIndex::~ProxyIndex() {}

// version
// $Id$

// End of file
