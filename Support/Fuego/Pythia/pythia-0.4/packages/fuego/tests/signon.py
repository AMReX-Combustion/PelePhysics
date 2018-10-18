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

if __name__ == "__main__":

    import fuego
    from fuego import fuego as fuegomodule

    print "copyright information:"
    print "   ", fuego.copyright()
    print "   ", fuegomodule.copyright()

    print
    print "module information:"
    print "    file:", fuegomodule.__file__
    print "    doc:", fuegomodule.__doc__
    print "    contents:", dir(fuegomodule)

    print
    print fuegomodule.hello()

# version
__id__ = "$Id$"

#  End of file 
