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


from builtins import object

import weaver.content


class LoginPage(object):
    def form(self):
        return self._form

    def render(self):
        self._weaver.weave()
        return

    def __init__(self, message=""):

        body = self._newBody()

        if message:
            msg = self._newMessage(message)
            body.add(msg)

        form = self._newForm()
        body.add(form)

        weaver = self._newWeaver()
        weaver.body(body)

        self._form = form
        self._weaver = weaver

        return

    def _newWeaver(self):
        w = weaver.content.weaver()
        p = w.properties()
        p.title = "Testing user authorization"
        p.stylesheet = "/weaver/test-style.html"
        return w

    def _newBody(self):
        body = weaver.content.body()
        p = body.properties()
        p.style = "test"
        return body

    def _newMessage(self, text):
        msg = weaver.content.textBlock(text)
        return msg

    def _newForm(self):
        from LoginForm import LoginForm

        form = LoginForm()
        p = form.properties()
        p.project = "test"
        p.mailto = "aivazis@caltech.edu"
        return form


# version
__id__ = "$Id$"

#  End of file
