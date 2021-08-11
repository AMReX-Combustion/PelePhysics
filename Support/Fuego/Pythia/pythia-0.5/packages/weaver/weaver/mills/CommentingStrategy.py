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


from builtins import object


class CommentingStrategy(object):
    def commentBlock(self, lines):
        block = []

        if not lines:
            return block

        block.append(self._beginCommentBlock(lines[0]))

        for line in lines[1:]:
            block.append(self._commentLineInBlock(line))
        block.append(self._endCommentBlock())

        return block

    def __init__(self):
        return

    def _beginCommentBlock(self, text=""):
        raise NotImplementedError(
            "class '%s' should override 'beginCommentBlock'"
            % self.__class__.__name__
        )

    def _commentLineInBlock(self, line=""):
        raise NotImplementedError(
            "class '%s' should override 'commentLineInBlock'"
            % self.__class__.__name__
        )

    def _endCommentBlock(self, text=""):
        raise NotImplementedError(
            "class '%s' should override 'beginCommentBlock'"
            % self.__class__.__name__
        )


# version
__id__ = "$Id$"

#  End of file
