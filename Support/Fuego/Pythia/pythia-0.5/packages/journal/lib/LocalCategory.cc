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

#include "Category.h"
#include "LocalCategory.h"

using namespace journal;

// interface
void journal::LocalCategory::activate() {
    _state = true;
    return;
}

void journal::LocalCategory::deactivate() {
    _state = false;
    return;
}

void journal::LocalCategory::state(bool flag) {
    _state = flag;
    return;
}

bool journal::LocalCategory::state() const {
    return _state;
}

// meta-methods
LocalCategory::~LocalCategory() {}

// version
// $Id$

// End of file
