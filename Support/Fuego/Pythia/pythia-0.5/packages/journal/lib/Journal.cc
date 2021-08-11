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

#include <portinfo>
#include <map>
#include <string>
#include <vector>
#include <ostream>

#include "Device.h"
#include "Journal.h"
#include "Index.h"
#include "Entry.h"

#include "StreamDevice.h"
#include "Console.h"


using namespace journal;

// meta-methods

Journal::Journal() :
    _device(new Console)
{}

Journal::~Journal() {
    delete _device;
}

// version
// $Id$

// End of file
