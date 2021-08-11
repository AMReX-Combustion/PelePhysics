/* -*- C -*-
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 *                               Julian C. Cummings
 *                        California Institute of Technology
 *                        (C) 1998-2003 All Rights Reserved
 * 
 *  <LicenseText>
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 */

#include <portinfo.h>
#include <stdio.h>

#include "journal/debuginfo.h"

/* the tests */
void test_basic();


/* driver */
int main() {
    test_basic();
    return 0;
}

/* where */

#if !defined(HAVE__FUNC__)
#define __FUNCTION__ "<unknown>"
#endif

/* test cases */

void test_basic() {

    printf(" ** %s:%s\n",__FILE__,__FUNCTION__);

    debuginfo_out("test",__HERE__,"this is the first line");
    debuginfo_out("test",__HERE__,"this is the second line");

    return;
}

/* $Id: c_info.c,v 1.1 2003/03/07 21:47:43 cummings Exp $ */

/* End of file */
