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


class Entry(object):


    def fullname(self):
        parts = []
        self._fullname(parts)
        return os.path.join(*parts)


    def id(self, inspector):
        raise NotImplementedError("class '%s' must override 'id'" % self.__class__.__name__)


    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent
        return


    def _fullname(self, parts):
        if self.parent:
            self.parent._fullname(parts)

        parts.append(self.name)
        return


# version
__id__ = "$Id$"

#  End of file 
