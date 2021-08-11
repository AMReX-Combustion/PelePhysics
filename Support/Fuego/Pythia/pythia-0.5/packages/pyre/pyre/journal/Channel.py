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

from pyre.components.Component import Component

from . import journal

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


class Channel(Component):
    def configure(self, registry):
        from pyre.util.bool import bool

        listing = self._listing(registry)

        for category, state in listing:
            journal.journal().channel(self.name).diagnostic(
                category
            ).state = bool(state)

        return []

    def __init__(self, name):
        Component.__init__(self, name, name)
        return

    def _listing(self, registry):
        listing = [
            (name, value) for name, value in registry.properties.items()
        ]

        listing += [
            ("%s.%s" % (nodename, name), value)
            for nodename, node in registry.facilities.items()
            for name, value in self._listing(node)
        ]

        return listing


# version
__id__ = "$Id$"

# End of file
