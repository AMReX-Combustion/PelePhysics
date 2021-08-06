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


from builtins import object
class AuthenticationRequest(object):


    def encode(self):
        import pickle
        return pickle.dumps(self)


    def __init__(self, username, password="", ticket=""):
        self.username = username
        self.password = password
        self.ticket = ticket
        return


# version
__id__ = "$Id$"

#  End of file 
