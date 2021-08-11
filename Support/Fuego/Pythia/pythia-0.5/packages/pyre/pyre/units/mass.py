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

from __future__ import absolute_import, division

from .SI import centi, kilo, kilogram, milli

#
# Definitions of common mass units
# Data taken from Appendix F of Halliday, Resnick, Walker, "Fundamentals of Physics",
#     fourth edition, John Willey and Sons, 1993

gram = kilogram / kilo
centigram = centi * gram
milligram = milli * gram


# aliases

kg = kilogram
g = gram
cg = centigram
mg = milligram


# other

metric_ton = 1000 * kilogram

ounce = 28.35 * gram
pound = 16 * ounce
ton = 2000 * pound


# version
__id__ = "$Id$"

#
# End of file
