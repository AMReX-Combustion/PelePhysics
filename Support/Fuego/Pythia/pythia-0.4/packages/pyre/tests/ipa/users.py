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


class UsersApp(Application):


    def run(self):
        manager = self._userManager
        manager.authenticate("aivazis", "mga5demo")
        manager.authenticate("aivazis", "mga4demo")
        print "\n".join(manager.pickle())
        return
    

    def __init__(self):
        Application.__init__(self, "users")
        self._userManager = None
        return


    def _init(self, parent):
        import pyre.ipa
        self._userManager = pyre.ipa.userManager("user-manager")
        self._userManager.db = self.inventory.db
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = [
            pyre.properties.str("db", default="userdb.crypt"),
            ]


# main
if __name__ == "__main__":
    import journal
    journal.info("user-manager").activate()
    
    app = UsersApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
