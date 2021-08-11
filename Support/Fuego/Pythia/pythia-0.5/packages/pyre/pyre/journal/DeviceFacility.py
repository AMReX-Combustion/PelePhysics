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

from __future__ import absolute_import

from pyre.facilities.Facility import Facility
from pyre.facilities.ScriptBinder import ScriptBinder


class DeviceFacility(Facility):
    def __init__(self, default=None):
        if default is None:
            from .Console import Console

            default = Console()

        Facility.__init__(
            self, name="device", default=default, binder=self.DeviceBinder()
        )
        return

    class DeviceBinder(ScriptBinder):
        def bind(self, facility, value):
            try:
                return self._builtins[value]()
            except KeyError:
                pass

            return ScriptBinder.bind(self, facility, value)

        def __init__(self):
            ScriptBinder.__init__(self)

            from .Console import Console
            from .File import File
            from .Remote import Remote

            self._builtins = {
                "console": Console,
                "file": File,
                "remote": Remote,
            }

            return


# version
__id__ = "$Id$"

# End of file
