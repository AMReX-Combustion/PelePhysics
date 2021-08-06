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
from .BaseParser import BaseParser


class Parameters(BaseParser):


    def extract(self, id, count):
        
        plist = self.parse()

        if plist.count != count:
            self._wrongCount(plist, id, count)
            
        return tuple(self._convert(plist))


    def extractRaw(self, id, count):
        plist = self.parse()

        if plist.count != count:
            self._wrongCount(plist, id, count)
            
        return tuple(plist.parameters)


    def extractOpt(self, id, lowCount, highCount):
        plist = self.parse()

        if plist.count < lowCount or plist.count > highCount:
            self._wrongVarCount(plist, id, lowCount, highCount)

        parameters = self._convert(plist)

        return tuple(parameters)


    def parse(self):
        self._parameters = None
        self._parse(self._scanner, self._tokenizer)
        return self._parameters


    # token handlers

    def aParameterList(self, token):
        self._parameters = token
        return 1
    
    
    def aComment(self, token):
        return 0


    def aWhitespace(self, token):
        return 0


    # other methods

    def __init__(self, tokenizer):
        import pyre
        BaseParser.__init__(self)

        self._tokenizer = tokenizer
        import fuego
        self._scanner = fuego.serialization.chemkin.unpickle.scanners.parameters()

        self._parameters = None

        return

            
    def _convert(self, token):
        try:
            parameters = list(map(float, token.parameters))
        except ValueError:
            msg = "Could not convert /%s/ into a list of numbers" % token.text
            self.onError(msg, self.locator())

        return parameters
            

    def _wrongCount(self, token, id, expectedCount):
        str = "The %s reaction specification requires exactly %d parameters." % (id, expectedCount)
        str = str + " '/%s/' apparently contains %d numbers" % (token.text, token.count)
        self.onError(str, self.locator())
        return


    def _wrongVarCount(self, token, id, lowCount, highCount):
        str = "The %s reaction specification requires" % id
        str = str + " %d-%d parameters." % (lowCount, highCount)
        str = str + " '/%s/' apparently contains %d numbers" % (token.text, token.count)
        self.onError(msg, self.locator())
        return


# version
__id__ = "$Id$"

#
# End of file
