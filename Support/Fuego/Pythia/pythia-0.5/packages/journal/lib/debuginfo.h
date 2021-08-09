/*
 * -*- C++ -*-
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 *                              Michael A.G. Aivazis
 *                       California Institute of Technology
 *                       (C) 1998-2003  All Rights Reserved
 * 
 *  <LicenseText>
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */

#if !defined(journal_debuginfo_h)
#define journal_debuginfo_h


/* get definition of __HERE__ macros */
#include "macros.h"


#if defined(__cplusplus)
extern "C" {
#endif
    
    int debuginfo_active(const char * category);
    void debuginfo_out(const char * category, __HERE_DECL__, const char * fmt, ...);

#if defined(__cplusplus)
}
#endif

#endif

/* version
 * $Id$
 */

/* End of file */

