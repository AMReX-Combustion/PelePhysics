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


class Copyright(AbstractNode):


    tag = "copyright"


    def content(self, text):
        self._copyright = text
        return


    def notify(self, parent):
        parent.onCopyright(self._copyright)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._copyright = ""
        return


# version
__id__ = "$Id$"

#  End of file 
