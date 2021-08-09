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


from builtins import object


class Renderer(object):
    def render(self, entry):

        str = []
        meta = entry.meta

        if self.header:
            filename = meta["filename"]
            if len(filename) > 53:
                filename = filename[0:20] + "..." + filename[-30:]
                meta["filename"] = filename
            str.append(self.header % meta)

        for line in entry.text:
            str.append(self.format % line)

        if self.footer:
            str.append(self.footer % meta)

        return str

    def __init__(self, header=None, format=None, footer=None):
        if header is None:
            header = " >> %(category)s(%(facility)s) -- %(filename)s:%(line)s:%(function)s"

        if format is None:
            format = " -- %s"

        if footer is None:
            footer = ""

        self.header = header
        self.format = format
        self.footer = footer

        return


# version
__id__ = "$Id$"

#  End of file
