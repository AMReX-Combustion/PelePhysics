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

import os
import cgi

class AuthorizationRequest:

    def __init__(self):
        self.username = ""
        self.cleartext = ""
        self.origin = ""
        self.attempts = 0
        self.ticket = ""
        return


def authorize():
    form = cgi.FieldStorage()

    request = AuthorizationRequest()
    request.username = form.getvalue("username")
    request.cleartext = form.getvalue("password")
    request.origin = os.environ.get("REMOTE_ADDR")

    attempts = form.getvalue("attempts")
    if attempts:
        request.attempts = int(attempts) + 1

    if request.username is None or request.cleartext is None:
        missingFields(request, "Please supply the required information")
        return

    import weaver.ipa
    session = weaver.ipa.session("auth")

    props = session.properties()
    props.port = 51000
    props.host = "localhost"
    props.mailto = "aivazis@caltech.edu"
    
    request.ticket = session.login(request.username, request.cleartext)
    if not request.ticket:
        missingFields(request, "Bad username or password: %d attempts" % request.attempts)
        return

    reportSuccess(request)

    return


def missingFields(request, message=""):
    from LoginPage import LoginPage

    page = LoginPage(message)
    form = page.form()

    if request.username is None:
        form.field("username").state = "error"
    else:
        form.field("username").value = request.username
    
    if request.cleartext is None:
        form.field("password").state = "error"
    else:
        form.field("password").value = request.cleartext

    form.field("attempts").value = request.attempts

    page.render()

    return


def reportSuccess(request):
    print '<html>'

    print '<head>'
    print '<title>Welcome</title>'
    print '</head>'

    print '<body>'
    print '<h2>Login information:</h2>'
    print '<table cols="2">'
    print '<tr><td align="right">username:</td><td>%s</td></tr>' % request.username
    print '<tr><td align="right">password:</td><td>%s</td></tr>' % request.cleartext
    print '<tr><td align="right">host:</td><td>%s</td></tr>' % request.origin
    print '<tr><td align="right">ticket:</td><td>%s</td></tr>' % request.ticket
    print '<tr><td align="right">attempts:</td><td>%d</td></tr>' % request.attempts
    print '</table>'
    print '</body>'
    print '</html>'

    return


# main
if __name__ == "__main__":
    print 'Content-type: text/html'
    print

    import settings
    settings.configure()

    authorize()



# version
__id__ = "$Id$"

#  End of file 
