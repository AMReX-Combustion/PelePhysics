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
void test_inject();


// driver
int main() {
    test_inject();
    return 0;
}

// where

#if !defined(HAVE__FUNC__)
#define __FUNCTION__ "<unknown>"
#endif
#define TEST__HERE__ __FILE__ << ":" << __FUNCTION__ 

// test cases

void test_inject() {

    std::cout << " ** " << TEST__HERE__ << std::endl;

    journal::info_t info("test");

    info.activate();

    info.inject(__FILE__);
    info.inject(":");
    info.inject(__LINE__);
    info.inject(":");
    info.inject(__FUNCTION__);
    info.inject(": Hello world!"); 
    info.record(); 

    return;
}

// $Id: inject.cc,v 1.1 2003/04/14 23:13:21 aivazis Exp $

// End of file
