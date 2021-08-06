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


class TextBlock(Weaver):


    def __init__(self, message):
        Weaver.__init__(self)
        self._text = message
        return


    def content(self):
        p = self.properties()
        content = [
            '',
            '<!-- text block -->',
            '<h2 class="%s">' % p.style,
            self._text,
            '',
            '</h2>',
            '<!-- end of: text block -->',
             ''
            ]

        return content


    class Properties(Weaver.Properties):


        def set(self, options):
            style = options.get("class")
            if style:
                self.style = style

            return


        def __init__(self):
            self.style = "blank"
            return


# version
__id__ = "$Id$"

#  End of file 
