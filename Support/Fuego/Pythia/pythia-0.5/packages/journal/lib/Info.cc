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
#include "Info.h"
#include "Index.h"
#include "LocalIndex.h"

using namespace journal;

// statics
static Info::index_t * _index = new LocalIndex(false);

// interface
void Info::newIndex(index_t * newIndex) {
    delete _index;
    _index = newIndex;
    return;
}

// meta-methods
Info::~Info() {}

Info::Info(string_t name):
    Diagnostic("info", name, _index->category(name))
{}

// version
// $Id$

// End of file
