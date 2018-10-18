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


from Mill import Mill
from BlockComments import BlockComments

class BlockMill(Mill, BlockComments):


    def __init__(self, begin, line, end, firstline):
        Mill.__init__(self)
        BlockComments.__init__(self)

        self.commentBlockLine = line
        self.commentBeginBlock = begin
        self.commentEndBlock = end
        self.firstLine = firstline

        return


# version
__id__ = "$Id$"

#  End of file 
