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


from AbstractNode import AbstractNode
from pyre.inventory.Registry import Registry


class Mill(AbstractNode):


    tag = "mill"


    def onAuthor(self, name):
        self._options.set("author", name)
        return


    def onOrganization(self, name):
        self._options.set("organization", name)
        return


    def onCopyright(self, name):
        self._options.set("copyright", name)
        return


    def onBannerCharacter(self, char):
        self._options.set("bannerCharacter", char)
        return


    def onBannerWidth(self, width):
        self._options.set("bannerWidth", width)
        return


    def notify(self, parent):
        parent.onMill(self._options)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        name = attributes["class"]
        self._options = Registry(name)

        return


# version
__id__ = "$Id$"

#  End of file 
