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


from __future__ import print_function


def run(options):

    import journal

    daemon = journal.daemon()
    daemon.properties().set(options)
    daemon.serve()

    return


# usage


def usage(program):

    print("Usage: %s [options ...]" % program)
    print("Options: (default values in brackets)")

    print("    --home=<directory> [%s]" % defaults["home"])
    print("    --logfile=<filename> [%s]" % defaults["logfile"])
    print("    --pidfile=<filename> [%s]" % defaults["pidfile"])
    print("    --port=<int> [%s]" % defaults["port"])
    print("    --timeout=<double> [%s]" % defaults["timeout"])

    return


defaults = {
    "home": "/",
    "port": "50000",
    "timeout": "30.0",
    "logfile": "/tmp/journald",
    "pidfile": "/tmp/journald.pid",
}

# main

if __name__ == "__main__":

    import pyre.applications

    options = pyre.applications.options(defaults, usage)
    run(options)


# version
__id__ = "$Id: journald.py,v 1.1.1.1 2003/02/16 04:23:02 aivazis Exp $"

#  End of file
