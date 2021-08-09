// -*- C++ -*-
//
//---------------------------------------------------------------------------
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
//
// <LicenseText>
//
//---------------------------------------------------------------------------
//

#if !defined(journal_manipulators_h)
#define journal_manipulators_h

// manipulators without any arguments

namespace journal {
    inline Diagnostic & end(Diagnostic &);
    inline Diagnostic & newline(Diagnostic &);

    typedef manip_2<const char *, long> loc2_t;
    typedef manip_3<const char *, long, const char *> loc3_t;

    typedef manip_2<const char *, const char *> manip2_t;
    typedef manip_3<const char *, const char *, const char *> manip3_t;

    inline manip2_t set(const char *, const char *);
    inline Diagnostic & __diagmanip_set(Diagnostic &, const char *, const char *);

    inline loc2_t loc(const char *, long);
    inline Diagnostic & __diagmanip_loc(Diagnostic &, const char *, long);

    inline loc3_t loc(const char *, long, const char *);
    inline Diagnostic & __diagmanip_loc(Diagnostic &, const char *, long, const char *);
}


// add a line
journal::Diagnostic & journal::newline(journal::Diagnostic & diag) {
    diag.newline();
    return diag;
}

// end of entry
journal::Diagnostic & journal::end(journal::Diagnostic & diag) {
    diag.record();
    return diag;
}

// set metadata key to value
journal::Diagnostic &
journal::__diagmanip_set(journal::Diagnostic & s, const char * key, const char * value) {
    s.attribute(key, value);
    return s;
}

journal::manip2_t journal::set(const char * key, const char * value) {
    return journal::manip2_t(journal::__diagmanip_set, key, value);
}

// location information
journal::Diagnostic &
journal::__diagmanip_loc(journal::Diagnostic & s, const char * filename, long line) {
    s.attribute("filename", filename);

    std::stringstream tmp;
    tmp << line;

    s.attribute("line", tmp.str());

    return s;
}

journal::loc2_t journal::loc(const char * file, long line) {
    return journal::loc2_t(journal::__diagmanip_loc, file, line);
}

journal::Diagnostic &
journal::__diagmanip_loc(
    journal::Diagnostic & s, const char * filename, long line, const char * function) 
{
    s.attribute("filename", filename);
    s.attribute("function", function);

    std::stringstream tmp;
    tmp << line;

    s.attribute("line", tmp.str());

    return s;
}

journal::loc3_t journal::loc(const char * file, long line, const char * function) {
    return journal::loc3_t(journal::__diagmanip_loc, file, line, function);
}

// get definition of __HERE__ macros
#include "macros.h"

#endif

// version
// $Id$

// End of file
