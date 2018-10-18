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


def sequence(spec, rangeSep='-', stepSep=':', itemSep=','):

    if not spec:
        return []

    if spec[0] == '[':
        spec = spec[1:]
    if spec[-1] == ']':
        spec = spec[:-1]
    
    if not spec:
        return []

    candidates = []

    initial = spec.split(itemSep)
    for entry in initial:
        token = entry.split(rangeSep)

        if len(token) == 1:
            candidates.append(int(token[0]))
        else:
            lower = int(token[0])

            stepspec = token[1].split(stepSep)
            upper = int(stepspec[0])
            if len(stepspec) == 1:
                step = 1
            else:
                step = int(stepspec[1])
                
            if upper > lower:
                candidates += range(lower,upper+1,step)
            else:
                candidates += range(lower,upper-1,-step)
            
    return candidates


# version
__id__ = "$Id$"

#  End of file 
