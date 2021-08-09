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


from pyre.applications.Registrar import Registrar


class MillRegistrar(Registrar):
    def mill(self, language):
        key = language.lower()
        factory = self._registry.get(key)
        if factory:
            return factory()
        return None

    def mills(self):
        return tuple(self._mills)

    # the singleton initializer
    def init(self):
        self._mills, self._registry = _buildRegistry()

        import journal

        self._debug = journal.debug("weaver.registry")
        return


# builtins

_builtinModules = [
    "CMill",
    "CshMill",
    "CxxMill",
    "Fortran77Mill",
    "Fortran90Mill",
    "HTMLMill",
    "MakeMill",
    "PerlMill",
    "PythonMill",
    "ShMill",
    "TeXMill",
    "XMLMill",
]


def _buildRegistry():

    mills = []
    registry = {}

    for module in _builtinModules:
        exec("from %s import %s as mill" % (module, module))
        names = mill.names
        mills.append(names[0])
        for name in names:
            registry[name] = mill

    return mills, registry


# version
__id__ = "$Id$"

#  End of file
