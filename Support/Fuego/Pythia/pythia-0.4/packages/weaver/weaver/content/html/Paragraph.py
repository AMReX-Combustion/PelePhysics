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
from .Weaver import Weaver


class Paragraph(Weaver):


    def add(self, lines):
        self._lines += lines
        return


    def __init__(self):
        Weaver.__init__(self)
        self._lines = []
        return


    def content(self):
        p = self.properties()
        content = [
            '',
            '<!-- paragraph -->',
            '<p class="%s">' % (p.style),
            ]

        content += self._lines

        content += [
             '</p>',
            '<!-- end of: paragraph -->',
            ]

        return content


    class Properties(Weaver.Properties):


        def set(self, options):
            style = options.get("class")
            if style:
                self.style = style

            return


        def __init__(self):
            self.style = ""
            return


# version
__id__ = "$Id$"

#  End of file 
