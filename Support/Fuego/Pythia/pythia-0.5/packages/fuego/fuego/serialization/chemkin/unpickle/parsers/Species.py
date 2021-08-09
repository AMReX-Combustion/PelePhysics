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

from builtins import str

from .BaseParser import BaseParser


class Species(BaseParser):

    # the interesting tokens

    def aSpeciesName(self, token):
        try:
            species = self._mechanism.newSpecies(token.name, self.locator())

        except self._mechanism.DuplicateSpecies as msg:
            self.onWarning(str(msg), self.locator())

        return 0

    # transitions

    def aSpeciesSection(self, token):
        self._info.log("species parser: section start")
        self._parse(self._scanner, self._tokenizer)
        return 0

    # other methods

    def __init__(self, mechanism, tokenizer):
        import pyre

        BaseParser.__init__(self, mechanism)

        self._tokenizer = tokenizer

        import fuego

        self._scanner = fuego.serialization.chemkin.unpickle.scanners.species()

        return


# version
__id__ = "$Id$"

#
# End of file
