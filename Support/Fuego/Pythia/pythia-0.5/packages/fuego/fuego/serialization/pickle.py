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
def save(mechanism, format="chemkin"):

    import journal
    journal.debug("fuego").log("pickling mechanism, format='%s')" % format)
    
    factory = registrar().retrieve(format)
    if not factory:
        journal.error("fuego").log("unknown mechanism file format '%s'" % format)
        return []
        
    pickler = factory.pickler()
    pickler.initialize()

    return pickler.pickle(mechanism)


# factory methods for the serializers

def pickler(format="chemkin"):
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
        _registrar.register(native, native.format())

        from . import chemkin
        _registrar.register(chemkin, chemkin.format())

        from . import ckml
        _registrar.register(ckml, ckml.format())

        from . import c
        _registrar.register(c, c.format())

        from . import f
        _registrar.register(f, f.format())

        from . import html
        _registrar.register(html, html.format())

        from . import python
        _registrar.register(python, python.format())

    return _registrar


# access  to the registrar
def picklers():
    return registrar()


# the registrar singleton
_registrar = None


# version
__id__ = "$Id$"

#  End of file 
