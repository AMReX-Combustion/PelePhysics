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

from .AbstractNode import AbstractNode


class NASA(AbstractNode):

    tag = "NASA"

    def notify(self, parent):
        parameters = (
            self._parameters["a1"],
            self._parameters["a2"],
            self._parameters["a3"],
            self._parameters["a4"],
            self._parameters["a5"],
            self._parameters["a6"],
            self._parameters["a7"],
        )
        parent.onNASA(self._lowT, self._highT, parameters)
        return

    def onNASA_a1(self, a1):
        self._parameters["a1"] = a1
        return

    def onNASA_a2(self, a2):
        self._parameters["a2"] = a2
        return

    def onNASA_a3(self, a3):
        self._parameters["a3"] = a3
        return

    def onNASA_a4(self, a4):
        self._parameters["a4"] = a4
        return

    def onNASA_a5(self, a5):
        self._parameters["a5"] = a5
        return

    def onNASA_a6(self, a6):
        self._parameters["a6"] = a6
        return

    def onNASA_a7(self, a7):
        self._parameters["a7"] = a7
        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._lowT = float(attributes["lowT"])
        self._highT = float(attributes["highT"])
        self._parameters = {}

        return


# version
__id__ = "$Id$"

#  End of file
