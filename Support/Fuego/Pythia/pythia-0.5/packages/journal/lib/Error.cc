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
#include <sstream>

#include "Category.h"
#include "Diagnostic.h"
#include "Error.h"
#include "Index.h"
#include "LocalIndex.h"

using namespace journal;

// statics
static Error::index_t * _index = new LocalIndex(true);

// interface
void Error::newIndex(index_t * newIndex) {
    delete _index;
    _index = newIndex;
    return;
}

// meta-methods
Error::~Error() {}

Error::Error(string_t name):
    Diagnostic("error", name, _index->category(name))
{}

// version
// $Id$

// End of file
