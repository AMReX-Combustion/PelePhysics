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
class Toggle(object):


    def activate(self):
        self._state = True
        return self


    def deactivate(self):
        self._state = False
        return self


    def flip(self):
        self._state ^= True
        return self


    def __init__(self, state=False):
        self._state = state
        return


    def _setState(self, flag):
        self._state = flag
        return


    def _getState(self):
        return self._state


    state = property(_getState, _setState, None, "")


    __slots__ = ("_state")


# version
__id__ = "$Id: Toggle.py,v 1.1.1.1 2003/02/16 04:22:52 aivazis Exp $"

#  End of file 
