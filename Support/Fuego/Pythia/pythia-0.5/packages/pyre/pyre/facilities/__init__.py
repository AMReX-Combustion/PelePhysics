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


# bultin facilities

from __future__ import absolute_import


def facility(
    name,
    default=None,
    public=None,
    binder=None,
    requirements=None,
    tip="",
    doc="",
):

    from .Facility import Facility

    return Facility(name, default, public, binder, requirements)


# version
__id__ = "$Id$"

# End of file
