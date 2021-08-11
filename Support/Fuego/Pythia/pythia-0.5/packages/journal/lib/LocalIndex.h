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

#if !defined(journal_LocalIndex_h)
#define journal_LocalIndex_h

// forward declarations
namespace journal {

    class Index;
    class LocalIndex;

}

// 
class journal::LocalIndex : public journal::Index {
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
    LocalIndex(bool defaultState);
    virtual ~LocalIndex();

// disable these
private:
    LocalIndex(const LocalIndex &);
    const LocalIndex & operator=(const LocalIndex &);

// data
private:
    index_t _index;
    bool _defaultState;
};

#endif

// version
// $Id$

// End of file
