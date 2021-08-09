#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from builtins import object


class EntitySet(object):
    def insert(self, key, entity):
        self._entities.append(entity)
        self._index[key] = entity
        return

    def size(self):
        return len(self._entities)

    def find(self, key=None):
        if key:
            return self._index.get(key)

        return self._entities

    def __init__(self):
        self._index = {}
        self._entities = []
        return


# version
__id__ = "$Id$"

# End of file
