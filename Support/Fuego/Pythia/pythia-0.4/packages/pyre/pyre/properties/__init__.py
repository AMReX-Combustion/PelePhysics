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


# builtin property types

def property(name, default=None, public=None, validator=None, tip="", doc=""):
    from Property import Property
    return Property(name, default, public, validator, tip, doc)


def str(name, default="", public=None, validator=None, tip="", doc=""):
    from String import String
    return String(name, default, public, validator, tip, doc)


def bool(name, default=False, public=None, validator=None, tip="", doc=""):
    from Bool import Bool
    return Bool(name, default, public, validator, tip, doc)


def int(name, default=0, public=None, validator=None, tip="", doc=""):
    from Integer import Integer
    return Integer(name, default, public, validator, tip, doc)


def float(name, default=0.0, public=None, validator=None, tip="", doc=""):
    from Float import Float
    return Float(name, default, public, validator, tip, doc)


def list(name, default=[], public=None, validator=None, tip="", doc=""):
    from List import List
    return List(name, default, public, validator, tip, doc)


def sequence(name, default=[], public=None, validator=None, tip="", doc=""):
    from Sequence import Sequence
    return Sequence(name, default, public, validator, tip, doc)


def dimensional(name, default=0.0, public=None, validator=None, tip="", doc=""):
    from Dimensional import Dimensional
    return Dimensional(name, default, public, validator, tip, doc)


# bultin validators

def less(value):
    from Less import Less
    return Less(value)


def greater(value):
    from Greater import Greater
    return Greater(value)


def range(low, high):
    from Range import Range
    return Range(low, high)


def choice(set):
    from Choice import Choice
    return Choice(set)


# version
__id__ = "$Id$"

# End of file 
