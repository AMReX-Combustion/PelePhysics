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

#if !defined(journal_Warning_h)
#define journal_Warning_h

// forward declarations
namespace journal {

    class Warning;
    class Index;
    class Diagnostic;

}

// 
class journal::Warning : public journal::Diagnostic {
// types
public:
    typedef Index index_t;

// interface
public:
    static void newIndex(index_t *);

// meta-methods
public:
    Warning(string_t name);
    virtual ~Warning();

// disable these
private:
    Warning(const Warning &);
    const Warning & operator=(const Warning &);
};

#endif

// version
// $Id$

// End of file
