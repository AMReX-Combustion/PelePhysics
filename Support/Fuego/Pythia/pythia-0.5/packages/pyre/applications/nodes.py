#!/usr/bin/env python
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2001 All Rights Reserved
#
#  <LicenseText>
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from __future__ import print_function

from builtins import range


def nodes(args):
    # extract user options
    spec = args["--range"]
    itemSep = args["--item-separator"]
    rangeSep = args["--range-separator"]

    format = args["--format"]

    maxload = float(args["--max-load"])
    port = int(args["--port"])
    maxnode = int(args["--max-node"])

    filename = args["--file"]

    # disable for now until the mond response format stabilizes
    # alive = query(maxload, format, port, maxnode)

    if spec:
        candidates = userRange(spec, format, itemSep, rangeSep)
        # NYI: validate against the live nodes
    else:
        candidates = []

    if not candidates:
        import sys

        sys.stderr.write(
            " ### no valid nodes specified: file '%s' will not be written\n"
            % filename
        )
        return []

    if filename:
        write(filename, candidates)

    return candidates


def query(maxload, format, port, maxnode):

    import select
    import socket

    nodelist = []

    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    for node in range(maxnode):
        host = format % (node + 1)
        s.sendto(host, (host, port))

    while 1:
        reads, writes, excepts = select.select([s], [], [], 1.00)

        # if reads is empty, we ran out of time
        if not reads:
            break

        msg = s.recv(port)
        print(msg)
        continue

        fields = msg.split(":")
        host = fields[0]
        # print "host =",host, "load = {", fields[1], "}"
        continue
        load = float(fields[1])
        if load < maxload:
            nodelist.append(host)

    nodelist.sort()

    return nodelist


def userRange(spec, format, itemSep, rangeSep):
    try:
        ids = expand(spec, itemSep, rangeSep)
    except ValueError:
        usage(program)
        return []

    candidates = [format % x for x in ids]

    return candidates


def expand(spec, itemSep, rangeSep):

    candidates = []

    initial = spec.split(itemSep)
    for entry in initial:
        token = entry.split(rangeSep)

        if len(token) == 1:
            candidates.append(int(token[0]))
        else:
            lower = int(token[0])
            upper = int(token[1])
            if upper > lower:
                candidates += list(range(lower, upper + 1))
            else:
                candidates += list(range(lower, upper - 1, -1))

    return candidates


def write(filename, candidates):

    file = open(filename, "w")
    for line in candidates:
        file.write(line + "\n")
    file.close()

    return


# main


def usage(program):

    print("Usage: %s [options ...]" % program)
    print("Options: (default values in brackets)")
    print("    --range=<range spec> [%s]" % defaults["--range"])
    print("    --file=<filename> [%s]" % defaults["--file"])
    print("    --format=<format spec> [%s]" % defaults["--format"])
    print(
        "    --item-separator=<character> [%s]" % defaults["--item-separator"]
    )
    print(
        "    --range-separator=<character> [%s]"
        % defaults["--range-separator"]
    )
    print("    --port=<daemon port> [%s]" % defaults["--port"])
    print(
        "    --max-load=<maximum acceptable load> [%s]"
        % defaults["--max-load"]
    )
    print(
        "    --max-node=<max node id to query> [%s]" % defaults["--max-node"]
    )
    print()
    print("Example: (with default values for the separators):")
    print("    %s --range=1-20,23,54,30-33 --file=beonodes" % program)
    return


defaults = {
    "--range": None,
    "--format": "a%03d",
    "--file": None,
    "--item-separator": ",",
    "--range-separator": "-",
    "--port": "8092",
    "--max-load": ".5",
    "--max-node": "100",
}


if __name__ == "__main__":

    import sys

    from pythlets.applications import main
    from pythlets.debug import DebugCenter

    program = sys.argv[0]

    args = main(defaults, usage)

    nodelist = nodes(args)
    for node in nodelist:
        print(node)


# version
__id__ = "$Id: nodes.py,v 1.5 2001/08/23 22:01:46 aivazis Exp $"

#  End of file
