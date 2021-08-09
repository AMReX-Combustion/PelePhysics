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

#if !defined(journal_ProxyIndex_h)
#define journal_ProxyIndex_h

// forward declarations
namespace journal {

    class Index;
    class ProxyIndex;

}

// 
class journal::ProxyIndex : public journal::Index {
// types
public:
    typedef std::string string_t;
    typedef std::map<string_t, category_t *> index_t;
    typedef std::pair<string_t, category_t *> entry_t;

// interface
public:
    virtual category_t & category(string_t name);

// meta-methods
public:
    ProxyIndex(PyObject * index);
    virtual ~ProxyIndex();

// disable these
private:
    ProxyIndex(const ProxyIndex &);
    const ProxyIndex & operator=(const ProxyIndex &);

// data
private:
    index_t _cache;
    PyObject * _index;
};

#endif

// version
// $Id$

// End of file
