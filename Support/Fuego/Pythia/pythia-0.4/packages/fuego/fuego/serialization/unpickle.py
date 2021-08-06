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


# extract a mechanism from a file

from __future__ import print_function
from __future__ import absolute_import
def load(filename, format=None, mechanism=None):

    import journal
    journal.debug("fuego").log("loading file '%s', format='%s')" % (filename, format))

    if not mechanism:
        from .mechanisms.Mechanism import Mechanism
        mechanism = Mechanism(filename)

    if not format:
        format = guessMechanismType(filename)
 
    print("***** Using format", format)

    factory = registrar().retrieve(format)

    if not factory:
        journal.error("fuego").log("unknown mechanism file format '%s'" % format)
        return None

    # GO TO packages/fuego/fuego/serialization/chemkin/unpickle/parsers
    parser = factory.parser()

    file = locate(filename)
    if not file:
        return None

    # chemkin/unpickle/parsers/ChemkinParser.py
    parser.parse(mechanism, file)

    return mechanism


def loadThermoDatabase(filename, format="chemkin", mechanism=None):

    if format != "chemkin":
        import journal
        journal.firewall("fuego").log(
            "cannot import thermodynamic databases in '%s' format" % format)

        return None

    file = locate(filename)
    if not file:
        return None

    if not mechanism:
        from .mechanisms.Mechanism import Mechanism
        mechanism = Mechanism(filename)



    from .chemkin.unpickle.parsers.ThermoDatabaseParser import ThermoDatabaseParser
    parser = ThermoDatabaseParser()
    parser.parse(mechanism, file)
    
    return mechanism


# find a mechanism file in the mechanism path
def locate(filename):
    import os
    import fuego
    import journal
    import pyre.util.locate

    pathlist = ["."] + fuego.mechanismPath()
    candidate = pyre.util.locate.locate(filename, pathlist)
    if candidate:
        journal.debug("fuego").log("resolved name '%s' as file '%s'" % (filename, candidate))
        return open(candidate, "r")

    import journal
    journal.error("fuego").log("could not locate '%s' in %s" % (filename, pathlist))

    return None


# use the filename to guess the mechanism format
def guessMechanismType(filename):
    import os
    name, ext = os.path.splitext(filename)

    factory = registrar().retrieve(ext)
    if factory:
        return factory.format()
    
    return "chemkin"
        

# factory method

def unpickler(format="chemkin"):
    factory = registrar().retrieve(format)
    if factory:
        return factory()
    
    return None


# the file format registrar
def registrar():
    global _registrar
    if not _registrar:
        from .Registrar import Registrar
        _registrar = Registrar()

        from . import native
        _registrar.register(native, native.format(), native.extensions())

        from . import chemkin
        _registrar.register(chemkin, chemkin.format(), chemkin.extensions())

        from . import ckml
        _registrar.register(ckml, ckml.format(), ckml.extensions())

    return _registrar


# access  to the registrar
def unpicklers():
    return registrar()


# the registrar singleton
_registrar = None


# version
__id__ = "$Id$"

#  End of file 
