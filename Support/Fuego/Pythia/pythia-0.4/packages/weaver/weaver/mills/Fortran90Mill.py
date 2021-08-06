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
from .LineMill import LineMill

class Fortran90Mill(LineMill):


    names = ["fortran90", "fortran95", "f90", "f95"]


    def __init__(self):
        LineMill.__init__(self, "!", "! -*- F90 -*-")
        return


# version
__id__ = "$Id$"

#  End of file 
