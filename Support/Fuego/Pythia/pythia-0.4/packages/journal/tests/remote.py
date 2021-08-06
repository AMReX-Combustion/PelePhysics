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
from builtins import range
import pyre
import journal


def run(argv):
    host = argv.get("host")
    port = int(argv.get("port"))
    sleep = float(argv.get("sleep"))
    timeout = float(argv.get("timeout"))

    from pyre.util.bool import bool
    if bool(argv.get("client")):
        client(host, port, sleep)
        return

    if bool(argv.get("server")):
        server(port, timeout)
        return

    import os
    clientId = os.fork()

    if not clientId:
        client(host, port, sleep)
        return

    server(port, timeout)
    os.waitpid(clientId, 0)

    return


def server(port, timeout):
    import select
    import socket
    import pickle
    import journal
    from journal.Entry import Entry
    import pyre.network

    # create the monitor
    monitor = pyre.network.monitor("udp")

    # install the monitor
    # the returned port may not be the same as we requested
    port, s = monitor.install(port, "")
    file = s.makefile("rb")

    while 1:
        reads, writes, excepts = select.select([s], [], [], timeout)
        if not reads:
            print(" -- timeout")
            continue
        
        entry = pickle.load(file)
        journal.journal().record(entry)

    return


def client(host, port, sleep):
    # sleep for the requested amount of time
    import select
    select.select([], [], [], sleep)

    import journal
    journal.remote(port, host, mode="udp")

    for id in range(5):
        info = journal.debug("test-%02d" % id).activate()
        info.log("test message on channel %d" % id)

    return


# usage

def usage(program):

    print("Usage: %s [options ...]" % program)
    print("Options: (default values in brackets)")

    # connection
    print("  connection:")
    print("    --host=<address> [%s]" % defaults["host"])
    print("    --port=<int> [%s]" % defaults["port"])

    # control
    print("  control:")
    print("    --sleep=<double> [%s]" % defaults["sleep"])
    print("    --timeout=<double> [%s]" % defaults["timeout"])

    # advanced
    print("  advanced:")
    print("    --client=<bool> [%s]" % defaults["client"])
    print("    --server=<bool> [%s]" % defaults["server"])

    return
        

defaults = {
    "host": "localhost",
    "port": "50000",
    "sleep": "2.0",
    "timeout": "5.0",
    "client": "false",
    "server": "false",
    }

# main

if __name__ == "__main__":

    import pyre.applications
    options = pyre.applications.options(defaults, usage)
    run(options)


# version
__id__ = "$Id: remote.py,v 1.1.1.1 2003/02/16 04:23:02 aivazis Exp $"

#  End of file 
