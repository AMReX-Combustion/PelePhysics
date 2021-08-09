#!/usr/bin/env python
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                           (C) 2001 All Rights Reserved
# 
#  <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import print_function
def test():

    import socket
    import select

    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    while 1:
        reads, writes, excepts = select.select([s], [], [], 180.00)

        # if reads is empty, we ran out of time
        if not reads:
            break

        print("read ready:", reads)
        print("write ready:", writes)
        print("exception ready:", excepts)

    return

# main

if __name__ == "__main__":
    nodelist = test()


# version
__id__ = "$Id$"

# End of file
