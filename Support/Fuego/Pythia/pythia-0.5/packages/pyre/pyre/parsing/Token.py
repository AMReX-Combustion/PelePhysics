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


class Token(object):
    def __init__(self, match, groups):
        self.match = match
        self.lexeme = groups[self.__class__.__name__]
        self.size = len(self.lexeme)
        return

    def __str__(self):
        return "{token: %s}" % self.lexeme

    __slots__ = ("lexeme", "match", "size")


# version
__id__ = "$Id$"

#  End of file
