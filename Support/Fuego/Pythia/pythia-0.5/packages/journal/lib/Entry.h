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

#if !defined(journal_Entry_h)
#define journal_Entry_h

// forward declarations
namespace journal {
    class Entry;
}

// 
class journal::Entry
{
// types
public:
    typedef std::string string_t;
    typedef std::vector<string_t> page_t;
    typedef std::map<string_t, string_t> meta_t;

// interface
public:

    inline void newline(string_t text);
    inline size_t lines() const;
    inline page_t::const_iterator lineEnd() const;
    inline page_t::const_iterator lineBegin() const;

    inline string_t & operator[](string_t key); 
    inline string_t  operator[](string_t key) const; 
    inline meta_t::const_iterator metaEnd() const;
    inline meta_t::const_iterator metaBegin() const;

    static inline void defaultAttributes(const meta_t & settings);

// meta-methods
public:
    virtual ~Entry();
    inline Entry();

// disable these
private:
    Entry(const Entry &);
    const Entry & operator=(const Entry &);

// data
private:
    meta_t _meta;
    page_t _text;

    static meta_t _defaults;
};

#endif

// get the inline definitions
#include "Entry.icc"

// version
// $Id$

// End of file
