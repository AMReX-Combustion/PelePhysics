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


from __future__ import absolute_import
def tmp():
    from .tmpdir import tmp
    return tmp()


def spawn(onParent, onChild):
    from .subprocesses import spawn
    return spawn(onParent, onChild)


def spawn_pty(onParent, onChild):
    from .subprocesses import spawn_pty
    return spawn_pty(onParent, onChild)

# version
__id__ = "$Id$"

#  End of file 
