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
from builtins import object
class FileSystem(object):

    def root(self):
        return self._root


    def expand(self, levels=0):
        level = 0
        todo = [self._root]

        while todo:
            working = todo
            todo = []
            for directory in working:
                todo += directory.expand()

            level += 1
            if levels and level >= levels: break

        return
        

    def __init__(self, root):
        import os
        from .Directory import Directory

        directory = os.path.abspath(root)
        self._root = Directory(directory, None) 
        return


# version
__id__ = "$Id$"

#  End of file 
