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

#if !defined(journal_Category_h)
#define journal_Category_h

// forward declarations
namespace journal {

    class Category;

}

// 
class journal::Category {
// interface
public:

    virtual void activate() = 0;
    virtual void deactivate() = 0;
    virtual void state(bool) = 0;

    virtual bool state() const = 0;

// meta-methods
public:
    inline Category();
    virtual ~Category();

// disable these
private:
    Category(const Category &);
    const Category & operator=(const Category &);
};

// get the inline definitions
#include "Category.icc"

#endif

// version
// $Id$

// End of file
