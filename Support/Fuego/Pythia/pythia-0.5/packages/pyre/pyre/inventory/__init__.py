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


def renderer(mode="pml"):
    if mode == "pml":
        from .pml.Renderer import Renderer

        renderer = Renderer()
    else:
        renderer = None
        import journal

        journal.error.log("'%s': unknown registry rendering mode" % mode)

    return renderer


def parser(mode="pml"):
    if mode == "pml":
        from .pml.Parser import Parser

        parser = Parser()
    else:
        parser = None
        import journal

        journal.error.log("'%s': unknown registry parsing mode" % mode)

    return parser


# version
__id__ = "$Id$"

# End of file
