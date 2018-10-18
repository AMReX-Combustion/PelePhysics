// -*- C++ -*-
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
//
//  <LicenseText>
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#if !defined(journalmodule_ProxyDevice_h_)
#define journalmodule_Device_h_

namespace journal {
    class Device;
    class ProxyDevice;
}

// 
class journal::ProxyDevice : public journal::Device
{
// interface
public:
    virtual void record(const entry_t &);

// meta-methods
public:
    virtual ~ProxyDevice();
    ProxyDevice(PyObject * journal);

// data
private:
    PyObject * _journal;

// disable these
private:
    ProxyDevice(const ProxyDevice &);
    const ProxyDevice & operator=(const ProxyDevice &);
};

#endif


// version
// $Id$

// End of file
