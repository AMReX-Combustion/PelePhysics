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


class JavaScript(Weaver):


    def scriptLibrary(self):
        return []


    def content(self):
        content = [
            '',
            '<script language="javascript">',
            '<!-- hide javascript from older browsers',
            ]
        content += self.scriptLibrary()
        content += [
            '// end of javascript -->',
            '</script>',
            ''
            ]

        return content


    def __init__(self):
        Weaver.__init__(self)
        self._fields = {}
        return


# version
__id__ = "$Id$"

#  End of file 
