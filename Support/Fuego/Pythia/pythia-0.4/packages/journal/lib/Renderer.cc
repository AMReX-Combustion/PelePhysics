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

#include "Renderer.h"
#include "Entry.h"

using namespace journal;

// interface
Renderer::string_t Renderer::render(const Entry & entry) {
    string_t text = header(entry) + body(entry) + footer(entry);
    return text;
}

// meta-methods
Renderer::~Renderer() {}

// version
// $Id$

// End of file
