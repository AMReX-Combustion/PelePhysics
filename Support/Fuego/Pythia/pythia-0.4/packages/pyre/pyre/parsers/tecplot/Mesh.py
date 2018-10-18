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


class Mesh(object):


    def read(self, simplices, file):

        count = 0
        self._info.log("scanning for simplices")
        while count < simplices:
            self._dump.log("processing simplex %d" % count)

            simplex = [ (int(node) - 1) for node in file.readline().split() ]
            self.simplices.append(simplex)
            count += 1

        self._info.log("read %d simplices" % count)

        return


    def __init__(self):
        self.simplices = []
        return


    import journal
    _info = journal.debug("tecplot")
    _dump = journal.debug("tecplot.scanning")


# version
__id__ = "$Id$"

#  End of file 
