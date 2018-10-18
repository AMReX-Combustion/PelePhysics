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
#include "DefaultRenderer.h"
#include "Entry.h"

using namespace journal;

// interface
DefaultRenderer::string_t DefaultRenderer::header(const Entry & entry) {
    string_t text = " >> ";

    string_t category = entry["category"];
    if (!category.empty()) {
        text += category;
        text += '(';

        string_t facility = entry["facility"];
        if (!facility.empty()) {
            text += facility;
        }

        text += ')';
    }

    text += "\n >> ";

    const size_t maxlen = 60;
    string_t filename = entry["filename"];
    if (!filename.empty()) {
        if (filename.size() > maxlen) {
            text += filename.substr(0, maxlen/2);
            text += "...";
            text += filename.substr(filename.size() - maxlen/2);
        } else {
            text += filename;
        }
    }

    string_t line = entry["line"];
    if (!line.empty()) {
        text += ":";
        text += line;
    }

    string_t function = entry["function"];
    if (!function.empty()) {
        text += ":";
        text += function;
    }

    text += '\n';
    
    return text;
}

DefaultRenderer::string_t DefaultRenderer::body(const Entry & entry) {
    string_t text;

    for (Entry::page_t::const_iterator l = entry.lineBegin(); l != entry.lineEnd(); ++l) {
        string_t line = *l;
        if (!line.empty()) {
            text += " -- " + line +"\n";
        }
    }

    return text;
}

DefaultRenderer::string_t DefaultRenderer::footer(const Entry & entry) {
    return string_t();
}

// meta-methods
DefaultRenderer::~DefaultRenderer() {}

// version
// $Id$

// End of file
