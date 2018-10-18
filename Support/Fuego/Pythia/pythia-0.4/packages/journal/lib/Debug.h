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

#if !defined(journal_Debug_h)
#define journal_Debug_h

// forward declarations
namespace journal {

    class Debug;
    class Index;
    class Diagnostic;

}

// 
class journal::Debug : public journal::Diagnostic {
// types
public:
    typedef Index index_t;

// interface
public:
    static void newIndex(index_t *);

// meta-methods
public:
    Debug(string_t name);
    virtual ~Debug();

// disable these
private:
    Debug(const Debug &);
    const Debug & operator=(const Debug &);
};

#endif

// version
// $Id$

// End of file
