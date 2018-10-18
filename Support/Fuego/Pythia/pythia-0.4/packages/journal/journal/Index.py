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


from Diagnostic import Diagnostic


class Index(object):


    def channel(self):
        return self._channel


    def categories(self):
        return self._index.keys()


    def diagnostic(self, name): 
        try:
            return self._index[name]
        except KeyError:
            diagnostic = Diagnostic(name, self._channel, self._defaultState, self._fatal)
            self._index[name] = diagnostic
            return diagnostic

        raise "Unknown error"


    def init(self, channel, defaultState, fatal=False):
        self._index = {}
        self._channel = channel
        self._defaultState = defaultState
        self._fatal = fatal

        import journal
        journal.journal().channel(channel, self)
        
        return


    def __new__(cls):
        index = cls.__dict__.get("__index__")
        if index is not None:
            return index

        cls.__index__ = index = object.__new__(cls)
        index.init()

        return index


# version
__id__ = "$Id$"

#  End of file 
