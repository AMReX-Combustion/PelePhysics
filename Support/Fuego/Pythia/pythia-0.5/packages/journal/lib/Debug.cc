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
#include <cstdlib>

#include <map>
#include <string>
#include <sstream>

#include "Category.h"
#include "Diagnostic.h"
#include "Debug.h"
#include "Index.h"
#include "LocalIndex.h"

using namespace journal;

#if defined(JOURNAL_TRACE)
#include <iostream>
#endif

// helpers
static Debug::index_t * _initialize();

// statics
static Debug::index_t * _index = _initialize();

// interface
void Debug::newIndex(index_t * newIndex) {
#if defined(JOURNAL_TRACE)
    std::cout 
        << " ** Debug: replacing index@" << _index << " with index@" << newIndex
        << std::endl;
#endif

    delete _index;
    _index = newIndex;
    return;
}

// meta-methods
Debug::~Debug() {}

Debug::Debug(string_t name):
    Diagnostic("debug", name, _index->category(name))
{
#if defined(JOURNAL_TRACE)
    std::cout 
        << " ** Debug: request for debug:" << name 
        << " from index@" << _index
        << std::endl;
#endif
}

// helpers
Debug::index_t * _initialize() {

    Debug::index_t * index = new LocalIndex(false);

    char * env = std::getenv("DEBUG_OPT");
    if (!env) {
        return index;
    }

    std::string text = env;
    std::string::size_type pos = 0;
    std::string::size_type len = text.size();

    while (pos < len) {
        
        // find the next ':'
        std::string::size_type end = pos + 1;
        while ((end < len) && (text[end] != ':')) {
            ++end;
        }

        // activate the named category
        std::string name = text.substr(pos, end-pos);
        index->category(name).activate();

        pos = end + 1;
    }

    return index;
}


// version
// $Id$

// End of file
