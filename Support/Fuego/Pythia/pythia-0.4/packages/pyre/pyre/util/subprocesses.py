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


def spawn(onParent, onChild):
    import os

    pid = os.fork()
    if pid > 0:
        return onParent(pid)

    import os
    pid = os.getpid()

    return onChild(pid)


def spawn_pty(onParent, onChild):
    import pty

    pid, fd = pty.fork()
    if pid > 0:
        return onParent(pid, fd)

    import os
    pid = os.getpid()

    return onChild(pid)


# version
__id__ = "$Id$"

# End of file 
