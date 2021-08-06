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
class Reaction(object):


    def update(self, T=None):
        self.forwardRate.update(T)
        self.reverseRate.update(T)
        return


    def id(self):
        return self.proxy.id


    def rate(self):
        return (self.forwardRate.rate(), self.reverseRate.rate())


    def progressRate(self):
        return (self.forwardRate(), self.reverseRate())


    def equilibriumConstant(self, T):
        return self.equillibrium(T)


    def __init__(self, proxy, equilibrium, forwardRate, reverseRate):
        self.proxy = proxy
        self.equlibrium = equilibrium
        self.forwardRate = forwardRate
        self.reverseRate = reverseRate

        return


    def __str__(self):
        return self.proxy.equation()


# version
__id__ = "$Id$"

#
# End of file
