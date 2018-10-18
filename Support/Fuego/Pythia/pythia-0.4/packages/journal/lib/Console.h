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

#if !defined(journal_Console_h)
#define journal_Console_h

// forward declarations
namespace journal {
    class Console;
    class StreamDevice;
}

// 
class journal::Console : public journal::StreamDevice
{
// meta-methods
public:
    Console();
    virtual ~Console();

// disable these
private:
    Console(const Console &);
    const Console & operator=(const Console &);
};

#endif

// version
// $Id$

// End of file
