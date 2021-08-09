#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from weaver.content.FormWeaver import FormWeaver


class RegistrationForm(FormWeaver):


    def action(self):
       return "/cgi-bin/weaver/authorize.py"


    def body(self):
        p = self.properties()

        first = self._fields["firstname"]
        last = self._fields["lastname"]
        affiliation = self._fields["affiliation"]

        address1 = self._fields["address1"]
        city = self._fields["city"]
        state = self._fields["state"]
        zip = self._fields["zip"]

        email = self._fields["email"]
        phone = self._fields["phone"]
        fax = self._fields["fax"]
        
        first.value = "Michael"
        last.value = "Aivazis"
        email.value = "aivazis@caltech.edu"
        phone.value = "(626) 395-3424"

        content = [
            '<form name="%s" method="post" action="%s">' % (p.name, self.action()),
            '  <table border="0" cols="3" cellspacing="10">',
            '    <tr>',
            '      <td class="field_%s">' % first.state,
            '        <div>First Name:</div>',
            '        <input name="firstName" type="text" value="%s">' % first.value,
            '      </td>',
            '      <td class="field_%s" colspan="2">' % last.state,
            '        <div>Last Name:</div>',
            '        <input name="lastName" type="text" value="%s">' % last.value,
            '      </td>',
            '    </tr>',
            '    <tr>',
            '      <td class="field_%s" colspan="3">' % affiliation.state,
            '        <div>Company:</div>',
            '        <input name="affiliation" type="text" value="%s">' % affiliation.value,
            '      </td>',
            '    </tr>',
            '    <tr>',
            '      <td class="field_%s" colspan="3">' % address1.state,
            '        <div>Address:</div>',
            '        <input name="affiliation" type="text" value="%s">' % address1.value,
            '      </td>',
            '    </tr>',
            '    <tr>',
            '      <td class="field_%s">' % city.state,
            '        <div>City:</div>',
            '        <input name="city" type="text" value="%s">' % city.value,
            '      </td>',
            '      <td class="field_%s">' % state.state,
            '        <div>State:</div>',
            '        <select name="state">',
            '          <option value="" selected></option>',
            '          <option value="AL">Alabama</option>',
            '          <option value="AK">Alaska</option>',
            '          <option value="AZ">Arizona</option>',
            '          <option value="AR">Arkansas</option>',
            '          <option value="CA">California</option>',
            '          <option value="CO">Colorado</option>',
            '          <option value="CT">Connecticut</option>',
            '          <option value="DE">Delaware</option>',
            '          <option value="DC">District of Columbia</option>',
            '          <option value="FL">Florida</option>',
            '          <option value="GA">Georgia</option>',
            '          <option value="HI">Hawaii</option>',
            '          <option value="ID">Idaho</option>',
            '          <option value="IL">Illinois</option>',
            '          <option value="IN">Indiana</option>',
            '          <option value="IA">Iowa</option>',
            '          <option value="KS">Kansas</option>',
            '          <option value="KY">Kentucky</option>',
            '          <option value="LA">Louisiana</option>',
            '          <option value="ME">Maine</option>',
            '          <option value="MD">Maryland</option>',
            '          <option value="MA">Massachusetts</option>',
            '          <option value="MI">Michigan</option>',
            '          <option value="MN">Minnesota</option>',
            '          <option value="MS">Mississippi</option>',
            '          <option value="MO">Missouri</option>',
            '          <option value="MT">Montana</option>',
            '          <option value="NE">Nebraska</option>',
            '          <option value="NV">Nevada</option>',
            '          <option value="NH">New Hampshire</option>',
            '          <option value="NJ">New Jersey</option>',
            '          <option value="NM">New Mexico</option>',
            '          <option value="NY">New York</option>',
            '          <option value="NC">North Carolina</option>',
            '          <option value="ND">North Dakota</option>',
            '          <option value="OH">Ohio</option>',
            '          <option value="OK">Oklahoma</option>',
            '          <option value="OR">Oregon</option>',
            '          <option value="PA">Pennsylvania</option>',
            '          <option value="RI">Rhode Island</option>',
            '          <option value="SC">South Carolina</option>',
            '          <option value="SD">South Dakota</option>',
            '          <option value="TN">Tennessee</option>',
            '          <option value="TX">Texas</option>',
            '          <option value="UT">Utah</option>',
            '          <option value="VT">Vermont</option>',
            '          <option value="VA">Virginia</option>',
            '          <option value="WA">Washington</option>',
            '          <option value="WV">West Virginia</option>',
            '          <option value="WI">Wisconsin</option>',
            '          <option value="WY">Wyoming</option>',
            '        </select>',
            '      </td>',
            '      <td class="field_%s">' % zip.state,
            '        <div>Zip Code:</div>',
            '        <input name="zip" type="text" value="%s">' % zip.value,
            '      </td>',
            '    </tr>',
            '    <tr>',
            '      <td class="field_%s">' % email.state,
            '        <div>Email:</div>',
            '        <input name="email" type="text" value="%s">' % email.value,
            '      </td>',
            '      <td class="field_%s">' % phone.state,
            '        <div>Phone:</div>',
            '        <input name="phone" type="text" value="%s">' % phone.value,
            '      </td>',
            '      <td class="field_%s">' % fax.state,
            '        <div>Fax:</div>',
            '        <input name="fax" type="text" value="%s">' % fax.value,
            '      </td>',
            '    </tr>',
            '    <tr>',
            '      <td></td>',
            '      <td></td>',
            '      <td align="center"><input type="submit" value="Update"></td>',
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

        return []


    def __init__(self):
        FormWeaver.__init__(self)

        self.register("firstname", "required", "")
        self.register("lastname", "required", "")
        self.register("affiliation", "required", "")
        self.register("address1", "required", "")
        self.register("city", "required", "")
        self.register("state", "required", "")
        self.register("zip", "required", "")
        self.register("email", "required", "")
        self.register("phone", "required", "")
        self.register("fax", "ok", "")

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
            self.name = "registration"
            self.title = "Registration form"
            self.info = "Please fill out the following form to register</br>Registration confirmation will be sent by email"
            self.mailto = ""
            return


# version
__id__ = "$Id$"

#  End of file 
