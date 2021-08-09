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


def ipad():
    import pyre.ipa

    app = pyre.ipa.daemon()
    app.main()
    return


# main
if __name__ == "__main__":
    import journal

    journal.info("ipad").activate()
    journal.debug("ipad").activate()
    journal.info("user-manager").activate()
    journal.debug("user-manager").activate()
    journal.info("session-manager").activate()
    journal.debug("session-manager").activate()

    ipad()


# version
__id__ = "$Id$"

# End of file
