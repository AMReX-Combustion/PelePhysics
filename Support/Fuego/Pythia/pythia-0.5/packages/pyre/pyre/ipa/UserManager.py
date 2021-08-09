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


from builtins import str

from past.builtins import execfile
from pyre.components.Component import Component


class UserManager(Component):
    def authenticate(self, username, cleartext):
        if self._reload:
            self._load()

        if not username or not cleartext:
            self._info.log("empty username or password")
            return ""

        cryptotext = self._users.get(username)
        if not cryptotext:
            self._info.log("bad username '%s'" % username)
            return ""

        if not self._verifyPassword(cleartext, cryptotext):
            self._info.log("bad password for '%s'" % username)
            return ""

        self._info.log("accepted password for '%s'" % username)
        return cryptotext

    def pickle(self):
        text = ["users = ["]
        text += [
            '    ("%s", "%s"),' % record
            for record in list(self._users.items())
        ]
        text += ["    ]"]

        return text

    def reload(self, *unused):
        self._reload = True
        return

    def __init__(self, name=None):
        if name is None:
            name = "user-manager"

        Component.__init__(self, name, facility="userManager")

        # public data
        self.db = None

        # private data
        self._users = {}
        self._reload = True
        self._verifyPassword = None

        self._encryptionMethods = {
            "md5": self._md5,
            "sha": self._sha,
            "crypt": self._crypt,
        }

        import journal

        self._info = journal.info(name)

        return

    def _load(self):
        # read the user database
        context = {}

        try:
            self._info.log("reading user records from '%s'" % self.db)
            execfile(self.db, context)
            users = context["users"]
            method = context.get("method", "crypt")

        except EnvironmentError as error:
            self._info.log("EnvironmentError: {%s}" % str(error))
            return
        except LookupError as error:
            self._info.log("LookupError: {%s}" % str(error))
            return
        except Exception as error:
            self._info.log("StandardError: {%s}" % str(error))
            return

        self._info.log("encryption method is '%s'" % method)

        count = len(users)
        if count == 1:
            suffix = ""
        else:
            suffix = "s"
        self._info.log("found %d user record%s" % (count, suffix))

        self._users = users
        self._reload = False
        self._verifyPassword = self._encryptionMethods[method]

        return

    def _crypt(self, cleartext, cryptotext):
        import crypt

        return crypt.crypt(cleartext, cryptotext[:2]) == cryptotext

    def _md5(self, cleartext, cryptotext):
        import md5

        return md5.new(cleartext).hexdigest() == cryptotext

    def _sha(self, cleartext, cryptotext):
        import sha

        return sha.new(cleartext).hexdigest() == cryptotext


# version
__id__ = "$Id$"

#  End of file
