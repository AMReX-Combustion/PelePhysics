#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from __future__ import absolute_import
from pyre.xml.Parser import Parser as BaseParser


class Parser(BaseParser):


    def parse(self, stream, parserFactory=None):
        from .parser.Document import Document
        return BaseParser.parse(self, stream, Document(stream.name), parserFactory)


    def __init__(self):
        BaseParser.__init__(self)
        return
           

# version
__id__ = "$Id$"

# End of file 
