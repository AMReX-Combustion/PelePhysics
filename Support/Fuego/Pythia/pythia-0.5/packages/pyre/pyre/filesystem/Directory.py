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

from .File import File


class Directory(File):
    def fullname(self):
        return self._path

    def id(self, inspector):
        return inspector.onDirectory(self)

    def children(self):
        return tuple(self._children.values())

    def subdirectories(self):
        return tuple(self._subdirectories)

    def expand(self):

        import stat

        from .BlockDevice import BlockDevice
        from .CharacterDevice import CharacterDevice
        from .File import File
        from .Link import Link
        from .NamedPipe import NamedPipe
        from .Socket import Socket

        subdirectories = []

        import os

        children = os.listdir(self.fullname())

        for name in children:

            if name in self._children:
                continue

            pathname = os.path.join(self.fullname(), name)
            # PORTABILITY: lstat is unix only
            mode = os.lstat(pathname)[stat.ST_MODE]

            if stat.S_ISDIR(mode):
                node = Directory(name, self)
                subdirectories.append(node)
            elif stat.S_ISREG(mode):
                node = File(name, self)
            elif stat.S_ISLNK(mode):
                node = Link(name, self)
            elif stat.S_ISSOCK(mode):
                node = Socket(name, self)
            elif stat.S_ISFIFO(mode):
                node = NamedPipe(name, self)
            elif stat.S_ISCHR(mode):
                node = CharacterDevice(name, self)
            elif stat.S_ISBLK(mode):
                node = BlockDevice(name, self)
            else:
                Firewall.hit("unknown file type: mode=%x" % mode)

            self._children[node.name] = node

        self._subdirectories = subdirectories
        return subdirectories

    def __init__(self, name, parent):
        File.__init__(self, name, parent)

        if not parent:
            self._path = name
        else:
            import os

            self._path = os.path.join(parent.fullname(), name)

        self._children = {}
        self._subdirectories = []
        return


# version
__id__ = "$Id$"

#  End of file
