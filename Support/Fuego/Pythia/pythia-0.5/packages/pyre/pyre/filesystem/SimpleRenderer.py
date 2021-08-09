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


from __future__ import absolute_import, print_function

from .Inspector import Inspector


class SimpleRenderer(Inspector):

    _INDENT = " " * 4

    def onCharacterDevice(self, node):
        self._render(node, "c")
        return

    def onBlockDevice(self, node):
        self._render(node, "b")
        return

    def onDirectory(self, node):
        self._render(node, "d")
        self._indent += 1

        for entry in node.children():
            entry.id(self)

        self._indent -= 1

        return

    def onFile(self, node):
        self._render(node, "f")
        return

    def onLink(self, node):
        self._render(node, "l")
        return

    def onNamedPipe(self, node):
        self._render(node, "p")
        return

    def onSocket(self, node):
        self._render(node, "s")
        return

    def render(self, node):
        node.id(self)
        return

    def __init__(self):
        self._indent = 0
        return

    def _render(self, node, code):
        print("%s(%s) %s" % (self._INDENT * self._indent, code, node.name))
        return


# version
__id__ = "$Id$"

#  End of file
