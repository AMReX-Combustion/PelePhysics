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

#if !defined(journal_LocalCategory_h)
#define journal_LocalCategory_h

// forward declarations
namespace journal {

    class Category;
    class LocalCategory;

}

// 
class journal::LocalCategory : public journal::Category {
// types
public:
    typedef bool state_t;

// interface
public:

    virtual void activate();
    virtual void deactivate();
    virtual void state(bool);

    virtual bool state() const;

// meta-methods
public:
    inline LocalCategory(bool state);
    virtual ~LocalCategory();

// disable these
private:
    LocalCategory(const LocalCategory &);
    const LocalCategory & operator=(const LocalCategory &);

// data
private:
    state_t _state;
};


// get the inline definitions
#include "LocalCategory.icc"

#endif

// version
// $Id$

// End of file
