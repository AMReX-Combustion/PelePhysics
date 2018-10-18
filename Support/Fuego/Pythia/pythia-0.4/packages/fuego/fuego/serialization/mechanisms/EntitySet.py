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


class EntitySet:


    def insert(self, key, entity):
        self._entities.append(entity)
        self._index[key] = entity
        return


    def replace(self, key, oldEntity, newEntity):
        newEntity.id = oldEntity.id
        self._entities[self._entities.index(oldEntity)] = newEntity
        self._index[key] = newEntity
        return

    def replace2(self, key, i, newEntity):
        self._entities[i] = newEntity
        self._index[key] = newEntity

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
