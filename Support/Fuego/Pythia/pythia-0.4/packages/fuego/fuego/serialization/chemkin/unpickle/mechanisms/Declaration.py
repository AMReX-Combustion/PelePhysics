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


class Declaration:


    def __init__(self, locator):
        self.locator = locator
        return


    def __str__(self):
        if self.locator:
            return "source='%s', line=%d, column=%d" % (
                self.locator.filename, self.locator.line, self.locator.column)

        return "source=<unknown>"

# version
__id__ = "$Id$"

# End of file
