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


from pyre.inventory.Configurable import Configurable


class Component(Configurable):
    def stage(self):
        return

    def execute(self):
        raise NotImplementedError(
            "class '%s' must override 'execute'" % self.__class__.__name__
        )

    def __init__(self, name, facility):
        Configurable.__init__(self, name)
        self.facility = facility

        import journal

        self._info = journal.info(self.name)
        self._debug = journal.debug(self.name)
        return


# version
__id__ = "$Id$"

# End of file
