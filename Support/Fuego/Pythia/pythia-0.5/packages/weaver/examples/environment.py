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
    header()
    body()
    footer()
    return


def header():
    print("Content-type: text/html")
    print("")
    print("<html>")

    print("<head>")
    print("<title>Welcome</title>")
    print("</head>")

    return


def body():
    print("<body>")
    print("<h1>Welcome</h1>")

    print("<br>")
    runtime()
    print("<br>")
    environment()
    print("<br>")
    modules()

    print("</body>")
    return


def runtime():
    import os
    import sys

    print("<h2>Runtime</h2>")
    print("<ul>")
    print(
        '<li>python: "%s" from "%s"</li>' % (sys.version_info, sys.executable)
    )
    print('<li>current working directory: "%s"</li>' % (os.getcwd()))
    import journal

    print('<li>journal: "%s"</li>' % (journal.__path__[0]))
    import weaver

    print('<li>weaver: "%s"</li>' % (weaver.__path__[0]))
    import srcmill

    print('<li>srcmill: "%s"</li>' % (srcmill.__path__[0]))
    from srcmill.html.HTMLMill import HTMLMill

    print('<li>HTMLMill: "%s"</li>' % (HTMLMill))
    print("</ul>")

    return


def environment():
    import os

    print("<h2>Environment variables</h2>")
    print("<ul>")

    for var, value in list(os.environ.items()):
        print('<li><b>%s</b> = "%s"</li>' % (var, value))

    print("</ul>")
    return


def modules():
    print("<h2>Modules</h2>")
    print("<blockquote>")
    print("<pre>")
    import journal

    info = journal.info("hello").activate()
    info.log("help")
    print("</pre>")
    print("</blockquote>")
    return


def footer():
    print("</html>")
    return


# main
if __name__ == "__main__":
    import settings

    settings.configure()

    run()


# version
__id__ = "$Id$"

#  End of file
