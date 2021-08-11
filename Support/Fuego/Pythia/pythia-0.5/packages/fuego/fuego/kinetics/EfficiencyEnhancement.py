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


from builtins import object


class EfficiencyEnhancement(object):
    def __init__(self, mixture, efficiencies):
        self._mixture = mixture
        self._efficiencies = efficiencies
        return

    def __call__(self):
        alpha = self._mixture.concentration()

        for species, factor in self._efficiencies:
            alpha += factor * species.concentration()

        return alpha


# version
__id__ = "$Id$"

#
# End of file
