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


class Section(Weaver):
    def paragraph(self, paragraph=None):
        if paragraph is None:
            from .Paragraph import Paragraph

            paragraph = Paragraph()
            return paragraph

        self._entries.append(paragraph)
        return

    def __init__(self):
        Weaver.__init__(self)
        self._entries = []
        return

    def content(self):
        p = self.properties()
        content = [
            "",
            "<!-- section -->",
            '<h%d class="%s">%s</h%d>' % (p.level, p.style, p.title, p.level),
            "",
        ]

        for entry in self._entries:
            content += entry.content()

        content += ["<!-- end of: section -->", ""]

        return content

    class Properties(Weaver.Properties):
        def set(self, options):
            level = options.get("level")
            if level:
                self.style = int(level)

            style = options.get("class")
            if style:
                self.style = style

            title = options.get("title")
            if title:
                self.title = title

            return

        def __init__(self):
            self.level = 2
            self.style = "blank"
            self.title = ""
            return


# version
__id__ = "$Id$"

#  End of file
