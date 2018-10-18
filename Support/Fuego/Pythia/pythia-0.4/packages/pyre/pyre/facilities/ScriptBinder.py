#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from Binder import Binder


class ScriptBinder(Binder):


    def bind(self, facility, value):
        script = file(self._resolve(value))
        context = {}
        exec script in context
        component = eval("%s()" % facility, context)
        return component


    def _resolve(self, name):
        import os
        base, ext = os.path.splitext(name)
        if not ext:
            name += ".py"

        return name

    
# version
__id__ = "$Id$"

# End of file 
