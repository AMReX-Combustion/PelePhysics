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


from pyre.applications.Application import Application


class SessionApp(Application):


    def run(self):
        session = self.inventory.session

        ticket = session.login("aivazis", "mga5demo")
        print "login: aivazis, ticket:", ticket
 
        ticket = session.login("aivazis", "mga4demo")
        print "login: aivazis, ticket:", ticket

        return
    

    def __init__(self):
        Application.__init__(self, "session")
        self._session = None
        return


    class Inventory(Application.Inventory):

        import pyre.ipa
        import pyre.facilities

        inventory = [
            pyre.facilities.facility("session", default=pyre.ipa.session())
            ]


# main
if __name__ == "__main__":
    import journal
    journal.info("session").activate()
    
    app = SessionApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
