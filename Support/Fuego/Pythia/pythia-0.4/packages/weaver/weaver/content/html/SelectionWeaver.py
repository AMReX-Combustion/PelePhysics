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


from Weaver import Weaver


class SelectionWeaver(Weaver):


    def options(self, ops):
        self._options = ops
        return


    def select(self, value):
        self._selection = value
        return


    def content(self):
        p = self.properties()
        
        content = [
            '',
            '<!-- selection: %s -->' % p.name,
            '',
            '<select name="%s" class="%s">' % (p.name, p.style),
            ]

        for value, text in self._options:
            if value == self._selection:
                selected = " selected"
            else:
                selected = ""
            option = '  <option%s value="%s">%s' % (selected, value, text)
            content.append(option)

        content += [
            '</select>',
            '<!-- end of selection: %s -->' % p.name,
            ''
            ]

        return content


    def __init__(self):
        Weaver.__init__(self)
        self._options = []
        self._selection = ""
        return


    class Properties(Weaver.Properties):


        def set(self, options):
            name = options.get("name")
            if name:
                self.name = name

            style = options.get("style")
            if style:
                self.style = style

            return


        def __init__(self):
            Weaver.Properties.__init__(self)
            self.name = ""
            self.style = ""
            return


        __slots__ = ("name", "style")


# version
__id__ = "$Id$"

#  End of file 
