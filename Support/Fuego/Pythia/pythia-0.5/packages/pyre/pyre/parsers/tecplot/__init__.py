#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from __future__ import absolute_import
def parse(input, mode="ascii"):

    if mode not in _supportedModes:
        import journal
        journal.error("tecplot").log("error: tecplot mode '%s' not supported" % mode)
        return [], {}


    if mode == "ascii":
        from .ReaderASCII import ReaderASCII as Reader
    else:
        import journal
        journal.firewall("tceplot").log("mode '%s' known but not supported" % mode)
        return [], {}
        
    reader = Reader()
    reader.read(input)

    return reader.mesh, reader.fields


# statics
_supportedModes = ("ascii", )

# version
__id__ = "$Id$"

#  End of file 
