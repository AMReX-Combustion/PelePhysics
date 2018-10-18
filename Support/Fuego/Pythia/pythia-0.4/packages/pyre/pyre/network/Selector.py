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

import journal


class Selector(object):


    def notifyOnReadReady(self, fd, handler):
        self._input.setdefault(fd, []).append(handler)
        return


    def notifyOnWriteReady(self, fd, handler):
        self._output.setdefault(fd, []).append(handler)
        return


    def notifyOnException(self, fd, handler):
        self._exception.setdefault(fd, []).append(handler)
        return


    def notifyOnInterrupt(self, handler):
        self._interrupt.append(handler)
        return


    def notifyWhenIdle(self, handler):
        self._idle.append(handler)
        return


    def watch(self, timeout=None):
        if timeout:
            self._timeout = timeout

        try:
            self._watch()
        except ValueError:
            self._cleanup()
        except KeyboardInterrupt:
            self._cleanup()

        return


    def __init__(self):
        self.state = True
        self._timeout = self._TIMEOUT

        # the fd activity clients
        self._input = {}
        self._output = {}
        self._exception = {}

        # clients to notify when there is nothing else to do
        self._idle = []
        self._interrupt = []
        
        return


    def _watch(self):
        import select

        while self.state:

            self._info.log("constructing list of watchers")
            iwtd = self._input.keys()
            owtd = self._output.keys()
            ewtd = self._exception.keys()

            self._info.log("input: %s" % iwtd)
            self._info.log("output: %s" % owtd)
            self._info.log("exception: %s" % ewtd)

            self._info.log("checking for indefinite block")
            if not iwtd and not owtd and not ewtd and not self._idle:
                self._info.log("no registered handlers left; exiting")
                return

            self._info.log("calling select")
            try:
                reads, writes, excepts = select.select(iwtd, owtd, ewtd, self._timeout)
            except select.error:
                self._info.log("signal received")
                continue
                
            self._info.log("returned from select")

            if not reads and not writes and not excepts:
                self._info.log("no activity; dispatching to idle handlers")
                for handler in self._idle:
                    if not handler(self):
                        self._idle.remove(handler)

            self._info.log("dispatching to exception handlers")
            self._dispatch(self._exception, excepts)
            self._info.log("dispatching to input handlers")
            self._dispatch(self._input, reads)
            self._info.log("dispatching to output handlers")
            self._dispatch(self._output, writes)

        return


    def _dispatch(self, handlers, entities):

        for fd in entities:
            for handler in handlers[fd]:
                if not handler(self, fd):
                    handlers[fd].remove(handler)
                    if not handlers[fd]:
                        del handlers[fd]

        return


    def _cleanup(self):
        self._info.log("cleaning up")
        for fd in self._input:
            fd.close()
        for fd in self._output:
            fd.close()
        for fd in self._exception:
            fd.close()

        for handler in self._interrupt:
            handler(self)

        return
        

    # static members
    _info = journal.debug("pyre.network.selector")

    # constants
    _TIMEOUT = .5


# version
__id__ = "$Id$"

#  End of file 
