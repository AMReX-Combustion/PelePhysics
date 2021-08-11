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

from .FormField import FormField
from .Weaver import Weaver


class FormWeaver(Weaver):
    def register(self, name, state="ok", value=""):
        field = FormField(name, state, value)
        self._fields[name] = field
        return field

    def field(self, name):
        return self._fields.get(name)

    def content(self):
        p = self.properties()

        content = [
            "",
            "<!-- form: %s -->" % p.name,
            "",
            '<table class="form" id="%s" cellpadding="5" width="100%%">'
            % p.id,
            '  <tr class="form_header">',
            "    <td>%s</td>" % p.title,
            "  </tr>",
        ]

        if p.info:
            content += [
                '  <tr class="form_info">',
                "    <td>%s</td>" % p.info,
                "  </tr>",
            ]

        content += ['  <tr class="form_body">', "    <td>"]

        content += self.body()

        content += [
            "    </td>",
            "  </tr>",
            '  <tr class="form_footer">',
            "    <td>",
        ]

        content += self.footer()

        content += [
            "    </td>",
            "  </tr>",
            "</table>",
            "",
            "<!-- end of form: %s -->" % p.name,
            "",
        ]

        return content

    def __init__(self):
        Weaver.__init__(self)
        self._fields = {}
        return

    class Properties(Weaver.Properties):
        def set(self, options):
            home = options.get("home")
            if home:
                self.home = home

            project = options.get("project")
            if project:
                self.project = project

            name = options.get("name")
            if name:
                self.name = name

            title = options.get("title")
            if title:
                self.title = title

            fid = options.get("id")
            if fid:
                self.id = fid

            info = options.get("info")
            if info:
                self.info = info

            return

        def __init__(self):
            Weaver.Properties.__init__(self)
            self.home = ""
            self.name = ""
            self.id = ""
            self.title = ""
            self.project = ""
            self.info = ""
            return

        __slots__ = ("home", "name", "id", "title", "project", "info")


# version
__id__ = "$Id$"

#  End of file
