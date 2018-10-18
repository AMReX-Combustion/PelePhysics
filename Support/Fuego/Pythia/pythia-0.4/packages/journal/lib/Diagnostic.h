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

#if !defined(journal_Diagnostic_h)
#define journal_Diagnostic_h

namespace journal {

    // forward declarations
    class Entry;
    class Journal;
    class Category;
    class Diagnostic;

    template <typename datum_t>
    inline Diagnostic & operator<< (Diagnostic &, datum_t);
}

// 
class journal::Diagnostic {
// types
public:
    typedef Entry entry_t;
    typedef Journal journal_t;
    typedef Category category_t;

    typedef std::string string_t;
    typedef std::stringstream buffer_t;

// interface
public:
    // state management
    inline void activate();
    inline void deactivate();
    inline void state(bool);

    inline bool state() const;

    // entry manipulation
    void record();
    void newline();
    void attribute(string_t key, string_t value);

    string_t str() const;

    // builtin data type injection
    template <typename datum_t> 
    inline Diagnostic & inject(datum_t datum);

    // access to the singleton
    static journal_t & journal();

// meta-methods
public:
    Diagnostic(string_t facility, string_t name, category_t & category);
    virtual ~Diagnostic();

// implementation
    void _newline();


// data
private:
    string_t _name;
    string_t _facility;
    category_t & _category;

    buffer_t _buffer;
    entry_t * _entry;

// disable these
private:
    Diagnostic(const Diagnostic &);
    const Diagnostic & operator=(const Diagnostic &);
};

// get the inline definitions
#include "Diagnostic.icc"

// get the manipulators
#include "manip-decls.h"
#include "manipulators.h"

#endif

// version
// $Id$

// End of file
