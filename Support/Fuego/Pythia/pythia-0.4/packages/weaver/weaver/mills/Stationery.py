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

import journal
from pyre.inventory.Configurable import Configurable


class Stationery(Configurable):


    def header(self):

        options = self.inventory
        separator = self.separator()

        h = ['', separator, '']
        h += self.copyright()
        h += ['']
        h += options.licenseText
        h += ['', separator]

        return [self.firstLine] + self.commentBlock(h)


    def footer(self):
        f = ['', self.line(self.inventory.lastLine) ]
        return f
            

    def copyright(self):
        c = []
        options = self.inventory

        # required
        width = options.bannerWidth

        # optional
        author = options.author
        copyright = options.copyright
        organization = options.organization

        if author:
            c.append(author.center(width).rstrip())

        if organization:
            c.append(organization.center(width).rstrip())

        if copyright:
            c.append((options.copyrightLine % copyright).center(width).rstrip())

        return c


    def separator(self):
        options = self.inventory
        banner = options.bannerCharacter
        cycles = options.bannerWidth/len(banner)
        separator = ' ' + banner * cycles
        return separator


    def blankLine(self):
        return ''


    def __init__(self, name):
        Configurable.__init__(self, name)
        self._debug = journal.debug(name)
        return


    # properties
    class Inventory(Configurable.Inventory):


        import pyre.properties


        inventory = (
            pyre.properties.property("author", default="-*- author -*-"),
            pyre.properties.property("organization", default="-*- organization -*-"),
            pyre.properties.property("copyright", default="-*- years -*-"),
            pyre.properties.int("bannerWidth", default=78),
            pyre.properties.property("bannerCharacter", default='~'),
            pyre.properties.property("lastLine", default=" End of file "),
            pyre.properties.property("copyrightLine", default="(C) %s  All Rights Reserved"),
            pyre.properties.property("licenseText", default=[ " <LicenseText>" ]),
            )


# version
__id__ = "$Id$"

#  End of file 
