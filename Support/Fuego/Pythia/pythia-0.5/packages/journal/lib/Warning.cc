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
#include "Warning.h"
#include "Index.h"
#include "LocalIndex.h"

using namespace journal;

// statics
static Warning::index_t * _index = new LocalIndex(true);

// interface
void Warning::newIndex(index_t * newIndex) {
    delete _index;
    _index = newIndex;
    return;
}

// meta-methods
Warning::~Warning() {}

Warning::Warning(string_t name):
    Diagnostic("warning", name, _index->category(name))
{}

// version
// $Id$

// End of file
