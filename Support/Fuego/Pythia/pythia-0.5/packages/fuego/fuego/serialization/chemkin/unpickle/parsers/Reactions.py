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

from .BaseParser import BaseParser


class Reactions(BaseParser):

    # the interesting tokens

    def aReaction(self, token):

        if self._currentReaction:
            self._finishReaction()

        self._reactionIndex += 1

        record = self._mechanism.newReaction(
            self._reactionIndex, self.locator()
        )
        # record.token = token

        if token._externalProduct != token._externalReactant:
            str = "The third body specification must be identical on both sides of a reaction"
            self.onWarning(str, self.locator())
            record.thirdBody = ("<mixture>", 1)
        else:
            record.thirdBody = token._externalReactant

        record.products = token._products
        record.reactants = token._reactants

        record.reversible = token._reversibleReaction
        record.falloff = token._falloff

        # Unit conversions happen in mechanism.reactionFactory
        record.arrhenius = token.arrhenius

        # Store the record
        self._currentReaction = record

        return 0

    def aReactionUnitsA(self, token):
        self._units["prefactor"] = token.units_A
        return 0

    def aReactionUnitsE(self, token):
        self._units["activation"] = token.units_E
        return 0

    def someArrheniusCoefficients(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.arrhenius = token.parameters
        return 0

    def aReactionDuplicate(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.duplicate = 1

        return 0

    def someReactionEfficiencies(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.efficiencies += token.efficiencies
        return 0

    def aReactionFORD(self, token):
        print("ford !!")
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        (record.ford).append(self._parameterParser.extractRaw("FORD", 2))
        return 0

    def aReactionHV(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        if self._currentReaction.radiation:
            msg = "Duplicate HV line"
            self.onWarning(msg, self.locator())
            return 0

        record.radiation = self._parameterParser.extract("HV", 1)

        return 0

    def aReactionLT(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.lt = self._parameterParser.extract("LT", 2)
        return 0

    def aReactionRLT(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.rlt = self._parameterParser.extract("RLT", 2)
        return 0

    def aReactionReverse(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        if not record.reversible:
            msg = "REV specification provided for an irreversible reaction"
            self.onWarning(msg, self.locator())
            return 0

        if record.low or record.sri or record.troe:
            msg = "ignoring REV specification for this falloff reaction"
            self.onWarning(msg, self.locator())
            return 0

        record.rev = self._parameterParser.extract("REV", 3)
        return 0

    def aReactionLOW(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        if record.rev:
            msg = "ignoring LOW parameters for this reaction with REV specification"
            self.onWarning(msg, self.locator())
            return 0

        if not record.falloff:
            msg = "LOW parameters provided for a reaction with no parenthesized species"
            self.onWarning(msg, self.locator())
            return 0

        record.low = self._parameterParser.extract("LOW", 3)
        return 0

    def aReactionSRI(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        if record.rev:
            msg = "ignoring SRI parameters for this reaction with REV specification"
            self.onWarning(msg, self.locator())
            return 0

        if not record.falloff:
            msg = "SRI parameters provided for a reaction with no parenthesized species"
            self.onWarning(msg, self.locator())
            return 0

        if record.troe:
            msg = "ignoring SRI parameters for this reaction with TROE falloff"
            self.onWarning(msg, self.locator())
            return 0

        record.sri = self._parameterParser.extractOpt("SRI", 3, 5)
        return 0

    def aReactionTROE(self, token):
        record = self._currentReaction
        if not record:
            msg = "no current reaction"
            self.onWarning(msg, self.locator())
            return 0

        if record.rev:
            msg = "ignoring TROE parameters for this reaction with REV specification"
            self.onWarning(msg, self.locator())
            return 0

        if not record.falloff:
            msg = "TROE parameters provided for a reaction with no parenthesized species"
            self.onWarning(msg, self.locator())
            return 0

        if record.sri:
            msg = "ignoring TROE parameters for this reaction with SRI falloff"
            self.onWarning(msg, self.locator())
            return 0

        record.troe = self._parameterParser.extractOpt("TROE", 3, 4)
        return 0

    # transitions

    def aReactionSection(self, token):
        self._info.log("reaction parser: section start")

        self._units["activation"] = "cal/mole"
        self._units["prefactor"] = "mole/cm**3"

        self._parse(self._scanner, self._tokenizer)
        return 0

    def anEndSection(self, token):
        if self._currentReaction:
            self._finishReaction()
        self._currentReaction = None
        return 1

    def onEndOfFile(self):
        if self._currentReaction:
            self._finishReaction()
        return 1

    # other methods

    def __init__(self, mechanism, tokenizer):
        import pyre

        BaseParser.__init__(self, mechanism)

        self._tokenizer = tokenizer
        import fuego

        self._scanner = (
            fuego.serialization.chemkin.unpickle.scanners.reactions()
        )

        from .Parameters import Parameters

        self._parameterParser = Parameters(tokenizer)

        self._reactionIndex = 0
        self._units = {}
        self._currentReaction = None

        return

    def _finishReaction(self):
        record = self._currentReaction

        if not record:
            import journal

            pyre.firewall("fuego").log("Null reaction record!")

        # Check for valid combinations of options
        # SRI and TROW require LOW
        if record.falloff and not record.low:
            msg = "Falloff reactions require LOW parameters"
            self.onWarning(msg, self.locator())
            return None

        # RLT requires LT
        if record.rlt and not record.lt:
            msg = "RLT requires LT"
            self.onWarning(msg, self.locator())
            return None

        # REV and LT requires RLT
        if record.rev and record.lt and not record.rlt:
            msg = (
                "Must specify RLT for reactions with REV and LT specifications"
            )
            self.onWarning(msg, self.locator())
            return None

        # Store the unit conversion factors
        record.units = self._units

        self._currentReaction = None
        return

    def _resolveParticipants(self, token, participants):

        resolved = []
        for participant in participants:
            name = participant[0]
            coefficient = participant[1]
            try:
                species = self._mechanism.species[name]
                resolved.append((species, coefficient))
            except KeyError:
                msg = "Unknown species '%s' in reaction" % name
                self.onWarning(msg, self.locator())

        return resolved

    def _unexpectedStatement(self, token):
        str = ""
        return str


# version
__id__ = "$Id$"

#
# End of file
