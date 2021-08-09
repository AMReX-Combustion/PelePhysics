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

#if !defined(journal_Error_h)
#define journal_Error_h

// forward declarations
namespace journal {

    class Error;
    class Index;
    class Diagnostic;

}

// 
class journal::Error : public journal::Diagnostic {
// types
public:
    typedef Index index_t;

// interface
public:
    static void newIndex(index_t *);

// meta-methods
public:
    Error(string_t name);
    virtual ~Error();

// disable these
private:
    Error(const Error &);
    const Error & operator=(const Error &);
};

#endif

// version
// $Id$

// End of file
