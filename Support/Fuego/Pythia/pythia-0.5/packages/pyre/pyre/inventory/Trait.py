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


from __future__ import print_function

from builtins import object


class Trait(object):

    prefix = "_pRoPeRtY_"

    def __init__(self, name, default=None, public=None, tip="", doc=""):
        self.doc = doc
        self.tip = tip
        self.name = name
        self.type = "generic trait"
        self.default = default
        self.mangled = self.prefix + name

        if public is None:
            public = name
        self.public = public

        # print "Trait.__init__: created name='%s', public='%s'" % (name, public)

        return

    def __del__(self):
        # print("Trait.__del__: deleted '%s'" % self.name)
        return

    def __get__(self, instance, cls):
        try:
            value = instance.__dict__[self.mangled]
        except KeyError:
            value = self.default
            instance.__dict__[self.mangled] = value
        # instance is None when accessed as a class attribute
        except AttributeError:
            value = cls.__dict__[self.name]

        return value


# version
__id__ = "$Id$"

# End of file
