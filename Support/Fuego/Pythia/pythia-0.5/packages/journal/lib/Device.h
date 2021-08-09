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

#if !defined(journal_Device_h)
#define journal_Device_h

// forward declarations
namespace journal {
    class Entry;
    class Device;
}

// 
class journal::Device
{
//types
public:
    typedef Entry entry_t;
    typedef std::string string_t;

// interface
public:
    virtual void record(const entry_t &) = 0;

// meta-methods
public:
    virtual ~Device();
    inline Device();

// disable these
private:
    Device(const Device &);
    const Device & operator=(const Device &);
};

#endif

// get the inline definitions
#include "Device.icc"

// version
// $Id$

// End of file
