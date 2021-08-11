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

TRUE = true = on = yes = 1
FALSE = false = off = no = 0


# convert strings to bool


def bool(response):
    return _stringToBool[response.lower()]


_stringToBool = {
    "y": 1,
    "yes": 1,
    "on": 1,
    "t": 1,
    "true": 1,
    "n": 0,
    "no": 0,
    "off": 0,
    "f": 0,
    "false": 0,
}


# version
__id__ = "$Id$"

#  End of file
