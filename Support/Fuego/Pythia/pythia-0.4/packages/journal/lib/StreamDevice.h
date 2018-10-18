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

#if !defined(journal_StreamDevice_h)
#define journal_StreamDevice_h

// forward declarations
namespace journal {
    class Device;
    class Entry;
    class StreamDevice;
    class Renderer;
}

// 
class journal::StreamDevice : public journal::Device
{
// types
public:
    typedef Renderer renderer_t;
    typedef std::ostream stream_t;

// interface
public:
    virtual void record(const Entry &);

// meta-methods
public:
    virtual ~StreamDevice();
    StreamDevice(stream_t & stream);

// implementation
protected:
    virtual void _write(const string_t);

// data
private:
    stream_t & _stream;
    renderer_t * _renderer;

// disable these
private:
    StreamDevice(const StreamDevice &);
    const StreamDevice & operator=(const StreamDevice &);
};

#endif

// version
// $Id$

// End of file
