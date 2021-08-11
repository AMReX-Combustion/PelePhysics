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

from __future__ import absolute_import

from .Trans import Trans


class TransDatabaseParser(Trans):

    # the main parsing loop

    def parse(self, mechanism, file):
        import journal
        import pyre

        self._mechanism = mechanism

        # prepare the parsing machinery
        tokenizer = pyre.parsing.tokenizer(file)
        tokenizer._info = journal.debug("fuego")

        # enter the parsing loop
        return Trans.parse(self, self._scanner, tokenizer)

    # transitions

    def aTransSection(self, token):
        self._transAll = 1
        return 0

    def anEndSection(self, token):
        return 1

    def onEndOfFile(self):
        return 1

    def onWarning(self, msg, locator):
        return

    def __init__(self):
        Trans.__init__(self, None, None)
        return


# version
__id__ = "$Id$"

#
# End of file
