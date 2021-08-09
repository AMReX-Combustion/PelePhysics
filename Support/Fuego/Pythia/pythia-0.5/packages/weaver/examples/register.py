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


from __future__ import print_function


def run():
    from RegistrationPage import RegistrationPage

    page = RegistrationPage()
    page.render()

    return


# main

if __name__ == "__main__":
    print("Content-type: text/html")
    print()

    import settings

    settings.configure()

    run()


# version
__id__ = "$Id"

#  End of file
