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
from .BodyWeaver import BodyWeaver


class StandardBody(BodyWeaver):


    def __init__(self):
        BodyWeaver.__init__(self)
        return


    def _renderContents(self):

        body = [
            '<table class="page" width="100%">',
            '  <tr class="page_header">',
            '    <td colspan="2">',
            '    <!-- page header -->',
            ]
        
        body += self._header()

        body += [
            '    <!-- end of: page header -->',
            '    </td>',
            '  </tr>',
            '  <tr>',
            '    <td class="page_logo">',
            '    <!-- page logo -->',
            ]

        body += self._logo()

        body += [
            '    <!-- end of: page logo -->',
            '    </td>',
            '    <td class="page_banner">',
            '    <!-- page banner -->',
            ]

        body += self._banner()

        body += [
            '    <!-- end of: page banner -->',
            '    </td>',
            '  </tr>',  
            '  <tr>',
            '    <td class="page_navigation">',
            '    <!-- page navigation bar -->',
            ]

        body += self._navigation()

        body += [
            '    <!-- end of: page navigation bar -->',
            '    </td>',
            '    <td class="page_body">',
            '    <!-- page body -->',
            ]
        
        body += BodyWeaver._renderContents(self)

        body += [
            '    <!-- end of: page body -->',
            '    </td>',
            '  </tr>',
            '  <tr class="page_footer">',
            '    <td colspan="2">',
            '    <!-- page footer -->',
            ]

        body += self._footer()

        body += [
            '    <!-- end of: page footer -->',
            '    </td>',
            '  </tr>',
            '</table>',
            ]
        
        return body


    def _header(self):
        return ['']


    def _logo(self):
        return ['']


    def _banner(self):
        return ['']


    def _navigation(self):
        return ['']


    def _footer(self):
        return [self.properties().footer]


    class Properties(BodyWeaver.Properties):


        def set(self, options):

            footer = options.get(footer)
            if footer:
                self.footer = footer

            return


        def __init__(self):
            BodyWeaver.Properties.__init__(self)
            self.footer = ""
            return


# version
__id__ = "$Id$"

#  End of file 
