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


class BodyWeaver(Weaver):


    def add(self, contentGenerator):
        self._contents.append(contentGenerator)
        return


    def __init__(self):
        Weaver.__init__(self)
        self._contents = []
        return


    def content(self):
        body = self._renderContents()
        if not body:
            body = ["<center><h2>This page was left blank accidentally</h2></center>"]

        p = self.properties()
        content = (
            ['', '<!-- body -->', '<body class="%s">' % p.style ]
            + body
            + ['', '</body>', '<!-- end of body -->', ''])

        return content


    def _renderContents(self):

        body = []
        for generator in self._contents:
            body += [''] + generator.content()

        return body

        
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
