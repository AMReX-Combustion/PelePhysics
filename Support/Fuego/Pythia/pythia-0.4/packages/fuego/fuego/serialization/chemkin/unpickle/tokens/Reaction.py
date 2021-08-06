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
from builtins import map
from builtins import range
import re

from .Token import Token

from .RegularExpressions import eol, whitespace, whitespaceOpt
from .RegularExpressions import species, coeff, namedNumbers_3

_maxParticipants = 6 # Maximum number of species on a reaction side


def _makeParticipant(id):
    namedCoeffOpt = r"(?P<coeff_%d>" % id + coeff + r")?"
    namedSpecies = r"(?P<species_%d>" % id + species + r")"
    participant = namedCoeffOpt + whitespaceOpt + namedSpecies
    return participant


def _makeParticipantOpt(id):

    pattern = r"(" + whitespaceOpt \
             + r"[+]" + whitespaceOpt + _makeParticipant(id) \
             + r")?"
    
    return pattern


def _makeParticipants():
    participants = [_makeParticipant(0)] + list(map(_makeParticipantOpt, list(range(1, _maxParticipants))))
    return "".join(participants) + whitespaceOpt + eol


class Reaction(Token):

    arrow = r"(?P<reaction_arrow><=>|=>|=)"
    notArrow = r"[^<=>]"
    
    _paramNames = ("reaction_A", "reaction_beta", "reaction_E")

    parameters = namedNumbers_3 % _paramNames

    pattern = \
           r"(?P<reactants>" + notArrow + r"+)" \
           + whitespaceOpt + arrow + whitespaceOpt \
           + r"(?P<products>" + notArrow + r"+)" \
           + whitespace + parameters
           
    parenthesizedParticipant = \
                             r"[(][+]" \
                             + _makeParticipant(0) \
                             + r"[)]"

    pressureDependent = re.compile(parenthesizedParticipant)

    participantPattern = re.compile(_makeParticipants())


    def identify(self, auth): return auth.aReaction(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        # extract reaction participants
        self._thirdBody = None
        self._falloff = 0

        products = groups["products"].strip()
        self._products, self._externalProduct = self._extractParticipants(products)

        reactants = groups["reactants"].strip()
        self._reactants, self._externalReactant = self._extractParticipants(reactants)

        # extract the reaction type
        arrow = groups["reaction_arrow"]
        if arrow == "=>":
            self._reversibleReaction = 0
        elif arrow == "=" or arrow == "<=>":
            self._reversibleReaction = 1
        else:
            import journal
            journal.firewall("fuego").log("Reaction: Unknown arrow '%s'" % arrow)

        paramList = list(map(groups.get, self._paramNames))
        try:
            self.arrhenius = list(map(float, paramList))
        except ValueError:
            # this can't happen because the regexp requires three floats here
            import journal
            str = "Could not convert '%s' into a list of numbers" % paramList
            journal.firewall("fuego").log(str)
            return

        return


    def _extractParticipants(self, text):
        thirdBody = None
        # check for pressure dependent reaction
        match = self.pressureDependent.search(text)
        if match:
            self._falloff = 1
            species, coeff = match.group("species_0", "coeff_0")

            if species == "M" or species == "m":
                species = "<mixture>"

            coeff = self._extractCoefficient(coeff)
            if coeff != 1:
                import journal
                msg = "Third body '%s' with coefficient=%d" % (species, coeff)
                journal.firewall("fuego").log(msg)

            thirdBody = (species, coeff)

            text = text[:match.start()] + text[match.end():]

        # extract the participants
        match = self.participantPattern.match(text)
        if not match:
            str = "syntax: '%s' is not well formed" % text
            raise self.TokenizationException(str)

        groups = match.groupdict()

        participants = []
        for i in range(_maxParticipants):
            species = groups["species_%d" % i]
            if not species: break
            coeff = self._extractCoefficient(groups["coeff_%d" % i])

            if species == "M" or species == "m":
                if thirdBody:
                    str = "more than one species acting as a third body"
                    raise self.TokenizationException(str)
                    
                thirdBody = ("<mixture>", coeff)
            else:
                participants.append((species, coeff))

        return (participants, thirdBody)


    def _extractCoefficient(self, coefficient):

        # Convert coefficient to an integer
        try:
            c = float(coefficient)
        except TypeError: # attempt to int(None); set to 1
            c = 1
        except ValueError:
            # this can't happen since t he regexp requires a number
            import journal
            msg = "Could not convert '%s' into a number" % coefficient
            journal.firewall("fuego").log(str)
            return

        return c


    def __str__(self):
        str = "{reaction: arrhenius=%s" % self.arrhenius
        str = str + "}"
        return str


# version
__id__ = "$Id$"

#
# End of file
