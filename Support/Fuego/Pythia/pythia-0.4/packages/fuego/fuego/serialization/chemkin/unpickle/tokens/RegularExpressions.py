#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Common regular expressions
eol = r"$"
whitespace = r"\s+"
whitespaceOpt = r"\s*"

element = r"[A-Za-z][\w+-]?"
namedElement = r"(?P<%s>" + element + r")"
species = r"[A-Za-z(][\w()=*,-]{0,15}[+]*"

coeff = r"\d+[.]?\d*"
coeffOpt = r"\d+[.]?\d*"

number = r"[+-]?(\d+[.]\d*|[.]\d+|\d+)([eE][-+]?\d{1,3})?"
numberOpt = r"(" + number + r")*"
namedNumber = r"(?P<%s>" + number + r")"

namedNumbers_3 = namedNumber + whitespace + namedNumber + whitespace + namedNumber

inlineNumber = r"/" + whitespaceOpt + number + whitespaceOpt + "/"
inlineNumbers = r"/" + whitespaceOpt \
                + number + r"(" + whitespace + number + ")*" \
                + whitespaceOpt + "/"

namedInlineNumber = r"/" + whitespaceOpt + namedNumber + whitespaceOpt + "/"

namedInlineNumbers = r"/" + whitespaceOpt + r"(?P<%s>" \
                     + number + r"(" + whitespace + number + r")*" \
                     + r")" + whitespaceOpt + "/"

namedInlineParameters = r"/" + whitespaceOpt\
                        + r"(?P<%s>" \
                        + r"(%s|%s)" % (species, number) \
                        + r"(" + whitespace + r"(%s|%s)" % (species, number) + r")*" \
                        + r")" + whitespaceOpt + "/"

# Simpler version of the above
#namedInlineParameters = r"/" + whitespaceOpt + "(?P<%s>[^/]+)" + whitespaceOpt + r"/"

namedInlineNumbers_2 = r"/" + whitespaceOpt \
                       + namedNumber + whitespace + namedNumber \
                       + whitespaceOpt + "/"

namedInlineNumbers_3 = r"/" + whitespaceOpt \
                       + namedNumber + whitespace + namedNumber + whitespace + namedNumber \
                       + whitespaceOpt + "/"


# version
__id__ = "$Id$"

#
# End of file
