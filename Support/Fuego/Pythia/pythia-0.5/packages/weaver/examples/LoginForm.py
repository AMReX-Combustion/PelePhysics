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


from weaver.content.FormWeaver import FormWeaver


class LoginForm(FormWeaver):


    def action(self):
       return "/cgi-bin/weaver/authorize.py"


    def body(self):
        p = self.properties()
        
        usr = self._fields["username"]
        pswd = self._fields["password"]
        attempts = self._fields["attempts"]
        
        content = [
            '<form name="%s" method="post" action="%s">' % (p.name, self.action()),
            '  <input name="attempts" type="hidden" value="%d">' % attempts.value,
            '  <table border="0" cols="2" cellspacing="10">',
            '    <tr>',
            '      <td align="right">',
            '        <table border="0" cellpadding="4" cellspacing="5">',
            '          <tr>',
            '            <td align="right"><b>Username:</b></td>',
            '            <td width="250" class="field_%s" align="left">' % usr.state,
            '              <input name="username" type="text" value="%s">' % usr.value,
            '            </td>',
            '          </tr>',
            '          <tr>',
            '            <td align="right"><b>Password:</b></td>',
            '            <td width="250" class="field_%s" align="left">' % pswd.state,
            '              <input name="password" type="password" value="%s">' % pswd.value,
            '            </td>',
            '          </tr>',
            '        </table>',
            '      </td>',
            '      <td align="left"><input type="submit" value="Login"></td>',
            '    </tr>',
            '  </table>',
            '</form>'
            ]
        
        return content


    def footer(self):
        p = self.properties()
        
        content = [
            '<a href="mailto:%s">questions?</a>' % p.mailto,
            ]

        return content


    def __init__(self):
        FormWeaver.__init__(self)

        self.register("attempts", "ok", 0)
        self.register("username", "required")
        self.register("password", "required")

        return


    class Properties(FormWeaver.Properties):


        def set(self, options):
            FormWeaver.set(self, options)

            mailto = options.get("mailto")
            if mailto:
                self.mailto = mailto

            return


        def __init__(self):
            FormWeaver.Properties.__init__(self)
            self.name = "login"
            self.title = "Access authorization"
            self.info = "Please enter your username and password"
            self.mailto = ""
            return


# version
__id__ = "$Id$"

#  End of file 
