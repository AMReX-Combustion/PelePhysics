// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Michael A.G. Aivazis
//                        California Institute of Technology
//                        (C) 1998-2003 All Rights Reserved
// 
//  <LicenseText>
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include <iostream>

#include "journal/journal.h"

// the tests
void test_basic();


// driver
int main() {
    test_basic();
    return 0;
}

// where

#if !defined(HAVE__FUNC__)
#define __FUNCTION__ "<unknown>"
#endif
#define TEST__HERE__ __FILE__ << ":" << __FUNCTION__ 

// test cases

void test_basic() {

    std::cout << " ** " << TEST__HERE__ << std::endl;

    journal::info_t info("test");

    info.activate();

    info 
        << journal::loc(__HERE__) 
        << journal::set("id", "0") 
        << "this is the first line" << journal::newline
        << "this is the second line" << journal::end;

    journal::info_t info2("test 1");
    info2.deactivate();

    return;
}

// $Id: info.cc,v 1.1.1.1 2003/02/16 04:23:02 aivazis Exp $

// End of file
