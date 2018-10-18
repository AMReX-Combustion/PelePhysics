#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import pty
import pyre.util.subprocesses
from pyre.applications.Application import Application


class ExpectApp(Application):


    def run(self):
        pyre.util.subprocesses.spawn_pty(self._read, self._print)
        return


    def __init__(self):
        Application.__init__(self, "expect")
        return


    def _read(self, childPid, childFd):
        import os
        import sys
        import select
        import signal

        pid = os.getpid()
        
        #print "%d: child pid: %d, child fd: %d" % (pid, childPid, childFd)
        index = 10
        while index:
            index -= 1
            r,w,e = select.select([childFd], [], [], 10)
            if not r:
                raise "timeout"

            if childFd in r:
                s = os.read(childFd, 10000)
                print "%d: got {%s}" % (pid, s.strip())
                os.write(childFd, "print 1+1%s" % os.linesep)


        os.kill(childPid, signal.SIGHUP)
        
        return


    def _print (self, pid):
        import os
        os.execvp("python", ("-i",))
        return
        raise "UNREACHABLE"
    

    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = [
            ]


# main
if __name__ == "__main__":
    app = ExpectApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
