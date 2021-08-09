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


import journal
from pyre.parsing.Parser import Parser


class BaseParser(Parser):


    # transitions

    def anElementSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)

    def aSpeciesSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)

    def aQssSpeciesSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)

    def aThermoSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)
        
    def aTransSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)

    def aReactionSection(self, token):
        self._tokenizer.unfetch(token)
        return self.anEndSection(token)

    def anEndSection(self, token):
        return 1


    # callback for the end-of-file

    def onEndOfFile(self):
        return 1
    

    # handlers for the boring tokens

    def aComment(self, token):
        return 0


    def aWhitespace(self, token):
        return 0


    # error handlers

    def onWarning(self, msg, locator=None): 
        channel = journal.warning("fuego")
        self._diagnostic(channel, msg, locator)
        return


    def onError(self, msg, locator=None):
        channel = journal.error("fuego")
        self._diagnostic(channel, msg, locator)
        return


    # others

    def __init__(self, mechanism=None):
        Parser.__init__(self)
        self._mechanism = mechanism
        return


    def _diagnostic(self, channel, message, locator):
        
        file = locator.filename
        line = locator.line
        column = locator.column
        
        channel.line("file '%s'(%d:%d): %s:" % (file, line, column, type))
        channel.line(message)
        channel.line(self._tokenizer.text)
        channel.log(" "*column + "^")

        return


    _info = journal.debug("pyre.chemistry.chemkin-parser")


# version
__id__ = "$Id$"

#
# End of file
