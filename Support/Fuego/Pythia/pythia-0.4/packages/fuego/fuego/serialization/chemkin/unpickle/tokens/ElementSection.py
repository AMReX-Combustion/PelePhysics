#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from Token import Token


class ElementSection(Token):

    pattern = r"[Ee][Ll][Ee][Mm]([Ee][Nn][Tt][Ss])?"


    def identify(self, auth):
        return auth.anElementSection(self)


    def __str__(self):
        return "{Element section}"


# version
__id__ = "$Id$"

#
# End of file
