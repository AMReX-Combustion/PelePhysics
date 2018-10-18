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

from BaseParser import BaseParser


class Elements(BaseParser):


    # the interesting tokens

    def anElementName(self, token):
        try:
            element = self._mechanism.newElement(token.name, token.weight, self.locator())

        except self._mechanism.DuplicateElement, msg:
            self.onWarning(str(msg), self.locator())
        
        return 0


    # transitions

    def anElementSection(self, token):
        self._info.log("element parser: section start")
        self._parse(self._scanner, self._tokenizer)
        return 0


    # other methods

    def __init__(self, mechanism, tokenizer):
        import pyre
        BaseParser.__init__(self, mechanism)

        self._tokenizer = tokenizer

        import fuego
        self._scanner = fuego.serialization.chemkin.unpickle.scanners.elements()

        return


# version
__id__ = "$Id$"

#
# End of file
