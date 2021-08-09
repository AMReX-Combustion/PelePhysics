#!/usr/bin/env python
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
#  <LicenseText>
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from builtins import map, object


class Nodal(object):
    def read(self, points, fields, input):

        for zone in fields:
            self.fields[zone] = []

        count = 0
        self._info.log("scanning for %d points" % points)
        while 1:
            line = input.readline()
            self._dump.log("processing point %d" % count)
            tokens = line.split()

            for zone, value in map(None, fields, tokens):
                self.fields[zone].append(float(value))

            count += 1
            if count == points:
                break

        self._info.log("read %d nodes" % count)

        return

    def __init__(self):
        self.fields = {}
        return

    import journal

    _info = journal.debug("tecplot")
    _dump = journal.debug("tecplot.scanning")


# version
__id__ = "$Id$"

#  End of file
