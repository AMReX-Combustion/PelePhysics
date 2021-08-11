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

import journal
from pyre.components.Component import Component


class CommandlineParser(Component):
    def parse(self, argv=None, root=None):
        if argv is None:
            import sys

            argv = sys.argv[1:]

        if root is None:
            from Configuration import Configuration

            root = Configuration("root")

        return self._parse(argv, root)

    def __init__(self):
        Component.__init__(self, "commandline-parser", "commandlineParser")

        self.help = ["?", "h", "help"]
        self.assignment = "="
        self.prefixes = ["--", "-"]
        self.separator = "."

        return

    def _parse(self, argv, root):
        unprocessed = []

        for arg in argv:
            self._debug.line("processing '%s'" % arg)

            # is this an option
            for prefix in self.prefixes:
                if arg.startswith(prefix):
                    self._debug.line(
                        "    prefix: '%s starts with '%s'" % (arg, prefix)
                    )
                    candidate = arg[len(prefix) :]
                    break
            else:
                # prefix matching failed; leave this argument alone
                self._debug.line("    prefix: '%s' is not an option" % arg)
                unprocessed.append(arg)
                continue

            self._debug.line(
                "    prefix: arg='%s' after prefix stripping" % candidate
            )

            # check for assignment
            tokens = candidate.split(self.assignment)
            self._debug.line("    tokens: %s" % repr(candidate))

            # dangling =
            if len(tokens) > 1 and not tokens[1]:
                self._debug.log("tokens: bad expression: %s" % arg)
                raise CommandlineParser.CommandlineException(
                    "bad expression: '%s': no rhs" % arg
                )

            # lhs, rhs
            lhs = tokens[0]
            if len(tokens) > 1:
                rhs = tokens[1]
            else:
                rhs = "true"
            self._debug.line("    tokens: key={%s}, value={%s}" % (lhs, rhs))

            if lhs in self.help:
                root.help = True

            # store this option
            self._processOption(lhs, rhs, root)

        self._debug.log()

        return unprocessed

    def _processOption(self, key, value, root):
        separator = self.separator
        fields = key.split(separator)

        self._debug.line("    sub: fields=%s" % repr(fields))
        node = root
        for field in fields[:-1]:
            self._debug.line("    tokens: looking for node '%s'" % field)
            node = node.open(field)

        # store the attribute
        self._debug.line("    option: setting '%s'='%s'" % (fields[-1], value))
        node.set(fields[-1], value)

        return

    # the class used for reporting errors
    class CommandlineException(Exception):
        def __init__(self, msg):
            self._msg = msg
            return

        def __str__(self):
            return self._msg


# version
__id__ = "$Id$"

#  End of file
