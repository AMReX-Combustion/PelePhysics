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

from builtins import range


def serve(port=8091):
    return


def machinefile(nodes, filename=None):
    if not filename:
        import os

        filename = os.path.join(os.environ["HOME"], ".nodes")

    file = open(filename, "w")

    for node, load in nodes:
        file.write(node + "\n")

    file.close()

    return


def refresh(port=8092, maxnode=100):

    import select
    import socket

    nodelist = []

    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    for node in range(maxnode):
        host = "a%03d" % (node + 1)
        s.sendto(host, (host, port))

    while 1:
        reads, writes, excepts = select.select([s], [], [], 1.00)

        # if reads is empty, we ran out of time
        if not reads:
            break

        msg = s.recv(8092)

        fields = msg.split(":")
        host = fields[0]
        load = fields[1]

        nodelist.append((host, load))

    nodelist.sort()

    return nodelist


# main

if __name__ == "__main__":
    nodelist = refresh()
    machinefile(nodelist)

    print("%3d" % len(nodelist), "nodes:", end=" ")
    for node, load in nodelist:
        print(node, end=" ")
    print()


# version
__id__ = "$Id: beomon.py,v 1.2 2001/06/06 23:56:59 aivazis Exp $"

# End of file
