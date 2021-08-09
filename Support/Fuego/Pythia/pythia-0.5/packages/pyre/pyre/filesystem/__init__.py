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


from __future__ import absolute_import

from .fastfind import fastfind


def root(directoryName):
    from .FileSystem import FileSystem

    return FileSystem(directoryName)


def listing(fs):
    from .SimpleRenderer import SimpleRenderer

    renderer = SimpleRenderer()
    renderer.render(fs.root())
    return


def explore(fs):
    from .Explorer import Explorer

    renderer = Explorer()
    renderer.render(fs.root())

    gtk.mainloop()

    return


def renderTree(fs):
    from .TreeRenderer import TreeRenderer

    renderer = TreeRenderer()
    renderer.render(fs.root())

    return


def find(fs, name):
    from .Finder import Finder

    root = fs.root()
    finder = Finder()

    return finder.find(root, name)


# version
__id__ = "$Id$"

#  End of file
