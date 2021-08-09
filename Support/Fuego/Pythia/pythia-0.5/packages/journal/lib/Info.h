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

#if !defined(journal_Info_h)
#define journal_Info_h

// forward declarations
namespace journal {

    class Info;
    class Index;
    class Diagnostic;

}

// 
class journal::Info : public journal::Diagnostic {
// types
public:
    typedef Index index_t;

// interface
public:
    static void newIndex(index_t *);

// meta-methods
public:
    Info(string_t name);
    virtual ~Info();

// disable these
private:
    Info(const Info &);
    const Info & operator=(const Info &);
};

#endif

// version
// $Id$

// End of file
