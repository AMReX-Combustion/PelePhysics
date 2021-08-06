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
from .CommentingStrategy import CommentingStrategy


class BlockComments(CommentingStrategy):


    def line(self, line=""):
        return self.commentBeginBlock + line + self.commentEndBlock


    def _beginCommentBlock(self, text=''):
        return self.commentBeginBlock + text


    def _commentLineInBlock(self, line=''):
        return self.commentBlockLine + line


    def _endCommentBlock(self, text=''):
        return self.commentEndBlock + text


# version
__id__ = "$Id$"

#  End of file 
