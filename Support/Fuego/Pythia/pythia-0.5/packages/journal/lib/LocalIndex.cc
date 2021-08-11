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

#include "Index.h"
#include "LocalIndex.h"
#include "Category.h"
#include "LocalCategory.h"

using namespace journal;

// interface
Index::category_t & LocalIndex::category(string_t name) {
    category_t * category;
    index_t::const_iterator target = _index.find(name);

    if (target == _index.end()) {
        category = new LocalCategory(_defaultState);
        _index.insert(entry_t(name, category));
    } else {
        category = (*target).second;
    }

    return *category;
}

// meta-methods
LocalIndex::LocalIndex(bool defaultState) :
    Index(),
    _index(),
    _defaultState(defaultState)
{}

LocalIndex::~LocalIndex() {}

// version
// $Id$

// End of file
