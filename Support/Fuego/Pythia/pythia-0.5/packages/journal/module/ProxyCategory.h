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

#if !defined(journal_ProxyCategory_h)
#define journal_ProxyCategory_h

// forward declarations
namespace journal {

    class Category;
    class ProxyCategory;

}

// 
class journal::ProxyCategory : public journal::Category {
// interface
public:

    virtual void activate();
    virtual void deactivate();
    virtual void state(bool);

    virtual bool state() const;

// meta-methods
public:
    ProxyCategory(PyObject * category);
    virtual ~ProxyCategory();

// data
private:
    PyObject * _category;

// disable these
private:
    ProxyCategory(const ProxyCategory &);
    const ProxyCategory & operator=(const ProxyCategory &);

};

#endif

// version
// $Id$

// End of file
