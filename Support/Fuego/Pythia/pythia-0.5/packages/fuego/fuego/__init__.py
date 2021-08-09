#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import

import os

from . import serialization


def mixture(speciesSet, type="pyre"):
    factory = registrar().retrieve(type)
    if not factory:
        import journal

        journal.error("fuego").log("unknown chemistry solver mode '%s'" % type)
        factory = registrar().retrieve("pyre")

    return factory.mixture(speciesSet)


def mechanism(mix, proxy, type="pyre"):
    factory = registrar().retrieve(type)
    if not factory:
        import journal

        journal.error("fuego").log("unknown chemistry solver mode '%s'" % type)
        factory = registrar().retrieve("pyre")

    return factory.mechanism(mix, proxy)


def registrar():
    global _registrar
    if not _registrar:
        from pyre.support.Registrar import Registrar

        _registrar = Registrar()

        # register the always available python factory
        from pyre.chemistry import mechanisms

        _registrar.register(mechanisms, "pyre", "native")

        # try to load fuego
        try:
            from pyre.solvers import fuego

            _registrar.register(fuego, "fuego")
        except:
            import journal

            journal.warning("fuego").log("could not register 'fuego'")

    return _registrar


# the registrar singleton
_registrar = None


# access to the chemistry related files


def mechanismPath():
    import pyre

    pathlist = FUEGO_MECHANISM_PATH

    if not pathlist:
        import os

        pathlist = [defaultMechanismDirectory()]

    return pathlist


def defaultMechanismDirectory():

    import os

    dir = os.path.abspath(
        os.path.join(home(), FUEGO_ETC_DIR, FUEGO_MECHANISM_DIR)
    )

    return dir


def chemkinMechanismFile(filename):
    import os

    return os.path.join(chemkinMechanismDirectory(), filename)


def home():
    return __path__[0]


def copyright():
    return "fuego pyre module: Copyright (c) 1998-2003 California Institute of Technology"


FUEGO_ETC_DIR = "../../etc/fuego"
FUEGO_MECHANISM_DIR = "mechanisms"
FUEGO_MECHANISM_PATH = []


# version
__id__ = "$Id$"

# End of file
