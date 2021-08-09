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

from __future__ import absolute_import, print_function

from builtins import str

from .BaseParser import BaseParser


class Trans(BaseParser):

    # the interesting tokens

    def aTransLine(self, token):

        # dispatch to appropriate line info parser
        self._lineParsers[0](token)

        return 0

    # transitions

    def aTransSection(self, token):
        self._info.log("trans parser: section start")

        # self._transAll = token._all
        # self._mechanism.transAll(self._transAll)
        # self._range = self._mechanism.thermoRange()
        self._parse(self._scanner, self._tokenizer)

        return 0

    # def anEndSection(self, token):
    # BaseParser.anEndSection(self, token)

    # if self._transAll:
    #    species = self._mechanism.species()

    #    candidates = []
    #    for s in species:
    #        if not s.trans:
    #            candidates.append(s.symbol)

    #    if candidates:
    #        msg = "no species information for the folowing species: %s" % candidates
    #        self.onWarning(msg, self.locator())

    # return 1

    # def onEndOfFile(self):
    #    self._mechanism.transDone()
    #
    #    return 1

    # other methods

    def __init__(self, mechanism, tokenizer):
        import pyre

        BaseParser.__init__(self, mechanism)

        self._tokenizer = tokenizer

        import fuego

        self._scanner = fuego.serialization.chemkin.unpickle.scanners.trans()

        # Private data
        # self._range = ()
        self._transAll = 0
        self._parameters = []
        # self._currentRange = None

        import re

        from fuego.serialization.chemkin.unpickle.tokens.RegularExpressions import (
            species,
        )

        self._speciesScanner = re.compile(species)

        self._nextId = 0  # line ids are zero-based

        self._currentSpecies = None

        self._lineParsers = [self._parseLine1]

        # self._thermoAllWarned = 0

        return

    def _parseLine1(self, token):
        spec_tmp = token.text
        lin = token.id
        params_tmp = token.text_2

        spec = spec_tmp.strip()
        params = params_tmp.strip()

        match = self._speciesScanner.match(spec)
        if not match:
            msg = "Could not match a valid species name in '%s'" % spec
            self.onError(msg, self.locator())

        speciesName = match.group()
        print(speciesName)

        species = self._mechanism.species(speciesName)
        if not species:
            msg = "trans section: undeclared species '%s'" % speciesName
            self.onWarning(msg, self.locator())
            species = self._mechanism.newSpecies(speciesName)

        species.locator(self.locator())
        # Save this information
        self._currentSpecies = species

        # Get the spec transport params
        self._parameters = [lin]
        # for i in range(0, 5):
        #    self._parameters.append(params.split()[i])
        self.EPS = params.split()[0]
        self.SIG = params.split()[1]
        self.DIP = params.split()[2]
        self.POL = params.split()[3]
        self.ZROT = params.split()[4]
        # store in the species
        # transParametrization(self, type, EPS, SIG, DIP, POL, ZROT, locator, parameters)
        self._currentSpecies.transParametrization(
            str(lin).strip(),
            self.EPS,
            self.SIG,
            self.DIP,
            self.POL,
            self.ZROT,
            self.locator(),
            self._parameters,
        )
