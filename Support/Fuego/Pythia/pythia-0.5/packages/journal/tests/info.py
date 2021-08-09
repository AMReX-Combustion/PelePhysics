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

from __future__ import print_function

if __name__ == "__main__":

    import journal

    print(journal.copyright())

    info = journal.info("info-1")
    debug = journal.debug("debug")
    info = journal.info("info-2")

    print("state of %s(%s): %s" % (info.name, info.facility, info.state))
    info.activate()
    print("state of %s(%s): %s" % (info.name, info.facility, info.state))

    info.log("hello")

    print("info categories:", journal.infoIndex().categories())
    print("debug categories:", journal.debugIndex().categories())


# version
__id__ = "$Id: info.py,v 1.1.1.1 2003/02/16 04:23:02 aivazis Exp $"

#  End of file
