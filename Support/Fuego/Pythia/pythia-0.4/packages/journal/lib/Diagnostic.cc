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
#include <sstream>

#include "Category.h"
#include "Diagnostic.h"
#include "Entry.h"
#include "Device.h"
#include "Journal.h"

using namespace journal;

// interface
Diagnostic::journal_t & Diagnostic::journal() {
    // statics
    static Diagnostic::journal_t * _journal = new Journal();

    return *_journal;
}

void Diagnostic::newline() {

    if (!_category.state()) {
        return;
    }

    _newline();

    return;
}

void Diagnostic::record() {
    if (!_category.state()) {
        return;
    }

    // save the current buffer contents
    _newline();

    // record the journal entry
    journal().record(*_entry);

    // reset the entry
    delete _entry;
    _entry = new Entry;
    (*_entry)["category"] = _name;
    (*_entry)["facility"] = _facility;

    return;
}

void Diagnostic::attribute(Diagnostic::string_t key, Diagnostic::string_t value) {

    (*_entry)[key] = value;
    return;
}


Diagnostic::string_t Diagnostic::str() const {
    return _buffer.str();
}

// meta-methods
Diagnostic::~Diagnostic() {}

Diagnostic::Diagnostic(string_t facility, string_t name, category_t & category):
    _name(name),
    _facility(facility),
    _category(category),
    _buffer(),
    _entry(new Entry)
{
    (*_entry)["category"] = _name;
    (*_entry)["facility"] = _facility;
}

// implementation
void Diagnostic::_newline() {
    _entry->newline(_buffer.str());
    _buffer.str(string_t());

    return;
}

// version
// $Id$

// End of file
