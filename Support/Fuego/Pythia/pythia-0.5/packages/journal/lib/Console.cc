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
#include <ostream>
#include <string>

#include "Device.h"
#include "StreamDevice.h"
#include "Console.h"

#include <iostream>

using namespace journal;

// meta-methods
journal::Console::Console() :
    StreamDevice(std::cout)
{}

Console::~Console() {}

// version
// $Id$

// End of file
