#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from __future__ import absolute_import
def session(name=None):
    from .Session import Session
    return Session(name)


def sessionManager(name=None):
    from .SessionManager import SessionManager
    return SessionManager(name)


def userManager(name=None):
    from .UserManager import UserManager
    return UserManager(name)


def daemon(name=None):
    from .Daemon import Daemon
    return Daemon(name)


# version
__id__ = "$Id$"

#  End of file 
