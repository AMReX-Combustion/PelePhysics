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
#include <vector>
#include <string>

#include "Entry.h"

using namespace journal;

// helpers

static Entry::meta_t initializeDefaults();

// interface
// meta-methods

Entry::~Entry() {}

// static data

Entry::meta_t Entry::_defaults = initializeDefaults();

// helpers
Entry::meta_t initializeDefaults() {
    Entry::meta_t defaults;

    defaults["category"] = "<unknown>";
    defaults["facility"] = "<unknown>";
    defaults["filename"] = "<unknown>";
    defaults["function"] = "<unknown>";
    defaults["line"] = "<unknown>";

    return defaults;
}

// version
// $Id$

// End of file
