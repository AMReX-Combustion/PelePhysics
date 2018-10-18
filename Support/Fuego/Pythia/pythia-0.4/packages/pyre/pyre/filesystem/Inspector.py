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


class Inspector(object):


    def onCharacterDevice(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onCharacterDevice'" % self.__class__.__name__)


    def onBlockDevice(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onBlockDevice'" % self.__class__.__name__)


    def onDirectory(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onDirectory'" % self.__class__.__name__)


    def onFile(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onFile'" % self.__class__.__name__)


    def onLink(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onLink'" % self.__class__.__name__)


    def onNamedPipe(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onNamedPipe'" % self.__class__.__name__)


    def onSocket(self, node):
        raise NotImplementedError(
            "class '%s' must override 'onSocket'" % self.__class__.__name__)


# version
__id__ = "$Id$"

#  End of file 
