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

#if !defined(journal_Index_h)
#define journal_Index_h

// forward declarations
namespace journal {

    class Index;
    class Category;

}

// 
class journal::Index {
// types
public:
    typedef Category category_t;
    typedef std::string string_t;

// interface
public:
    virtual category_t & category(string_t name) = 0;

// meta-methods
public:
    inline Index();
    virtual ~Index();

// disable these
private:
    Index(const Index &);
    const Index & operator=(const Index &);
};

// include the inline definitions
#include "Index.icc"

#endif

// version
// $Id$

// End of file
