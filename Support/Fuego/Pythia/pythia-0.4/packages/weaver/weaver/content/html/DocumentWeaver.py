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
from weaver.mills.HTMLMill import HTMLMill


class DocumentWeaver(Weaver, HTMLMill):


    def body(self, body=None):
        if body:
            self._body = body

        return self._body


    def javascript(self, library=None):
        if library:
            self._javascript = library

        return self._javascript


    def weave(self):
        for line in self.pickle():
            print line
        return


    def pickle(self):
        self._begin()
        self._rep += self.content()
        self._end()
        return self._rep


    def content(self):
        p = self.properties()
        
        content = [
            '',
            '<html>',
            '<head>',
            '<title>%s</title>' % (p.title),
            '<link href="%s" title="stylesheet" rel="stylesheet" type="text/css">' % p.stylesheet,
            '</head>',
            '',
            ]

        if self._javascript:
            content += self._javascript.content()

        if self._body:
            self._body.properties().footer = self._timestamp()
            content += self._body.content()
        else:
            content += [
                '<body class="blank">',
                '<center><h2>This page was left blank accidentally</h2></center>',
                '</body>',
                ]
            
        content += [
            '',
            '</html>',
            ''
            ]

        return content


    def __init__(self):
        Weaver.__init__(self)
        HTMLMill.__init__(self)
        self.initialize()

        self._body = None
        self._javascript = None

        return


    class Properties(Weaver.Properties, HTMLMill.Properties):


        def set(self, options):
            HTMLMill.Properties.set(self, options)

            title = options.get("title")
            if title:
                self.title = title

            stylesheet = options.get("stylesheet")
            if stylesheet:
                self.stylesheet = stylesheet

            return


        def __init__(self):
            HTMLMill.Properties.__init__(self)

            self.title = ""
            self.stylesheet = ""
            
            return


# version
__id__ = "$Id$"

#  End of file 
