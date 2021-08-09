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

# journal

from __future__ import absolute_import
def journal():
    global _theJournal
    if _theJournal is not None:
        return _theJournal
    
    from .Journal import Journal
    _theJournal = Journal("journal")
    return _theJournal

# channels

def firewall(name):
    from .IndexFirewall import IndexFirewall
    return IndexFirewall().diagnostic(name)


def debug(name):
    from .IndexDebug import IndexDebug
    return IndexDebug().diagnostic(name)


def info(name):
    from .IndexInfo import IndexInfo
    return IndexInfo().diagnostic(name)


def warning(name):
    from .IndexWarning import IndexWarning
    return IndexWarning().diagnostic(name)


def error(name):
    from .IndexError import IndexError
    return IndexError().diagnostic(name)


# indices

def firewallIndex():
    from .IndexFirewall import IndexFirewall
    return IndexFirewall()


def debugIndex():
    from .IndexDebug import IndexDebug
    return IndexDebug()


def infoIndex():
    from .IndexInfo import IndexInfo
    return IndexInfo()


def warningIndex():
    from .IndexWarning import IndexWarning
    return IndexWarning()


def errorIndex():
    from .IndexError import IndexError
    return IndexError()


# special setups

def remote(port, host="localhost", mode="udp"):

    if mode == "udp":
        from .UDPDevice import UDPDevice
        device = UDPDevice(port, host)
    else:
        from .Console import Console
        device = Console()
        
    journal().device = device
    return device


def daemon():
    from JournalDaemon import JournalDaemon
    return JournalDaemon()
    

# misc

def copyright():
    return "journal pyre module: Copyright (c) 1998-2003 Michael A.G. Aivazis";


# statics
_theJournal = None


# initialize
try:
    #print " ** __init__.py: importing _journal"
    import _journal
except ImportError:
    err = error("journal")
    err.line("could not import the C++ bindings for journal")
    err.log("control of diagnostics from extension modules is unavailable")
else:
    #print " ** __init__.py: initializing C++ bindings"
    _journal.initialize(
        journal(),
        firewallIndex(), debugIndex(),
        infoIndex(), warningIndex(), errorIndex()
        )
        
# version
__id__ = "$Id$"

#  End of file 
