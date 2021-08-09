#!/usr/bin/env python
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
#  <LicenseText>
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import xml.sax

import journal


class Parser(xml.sax.ContentHandler):
    def document(self):
        return self._documentNode

    # parsing
    def parse(self, file, documentNode, parserFactory=None):

        self._currentNode = documentNode
        self._documentNode = documentNode

        if parserFactory:
            parser = parserFactory()
        else:
            parser = xml.sax.make_parser()

        parser.setContentHandler(self)
        parser.parse(file)

        # break a circular reference introduced above
        parser.setContentHandler(None)

        # clean up
        self._nodeStack = []
        self._currentNode = None

        return

    # content demultiplexing

    def startDocument(self):
        # give document nodes a way to get at line info
        self._documentNode.setLocator(self._locator)
        return

    def endDocument(self):
        # break a circular reference introduced above
        self._documentNode.setLocator(None)
        return

    def startElement(self, name, attributes):

        line, column = (
            self._locator.getLineNumber(),
            self._locator.getColumnNumber(),
        )
        self._info.log(
            "startElement: '%s', at (%d, %d)" % (name, line, column)
        )

        self._nodeStack.append(self._currentNode)
        self._currentNode = self._documentNode.node(name, attributes)

        return

    def characters(self, content):
        content = content.strip()
        if content:
            line, column = (
                self._locator.getLineNumber(),
                self._locator.getColumnNumber(),
            )
            self._info.log(
                "characters: '%s', at (%d, %d)" % (content, line, column)
            )
            self._currentNode.content(content)

        return

    def endElement(self, name):
        line, column = (
            self._locator.getLineNumber(),
            self._locator.getColumnNumber(),
        )
        self._info.log("endElement: '%s', at (%d, %d)" % (name, line, column))

        node = self._currentNode
        self._currentNode = self._nodeStack.pop()

        node.notify(self._currentNode)

        return

    # constructor
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self._nodeStack = []
        self._currentNode = None
        self._documentNode = None
        return

    _info = journal.debug("pyre.xml.parsing")


# version
__id__ = "$Id$"

#  End of file
