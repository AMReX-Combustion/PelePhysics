#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2001  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

import pyre

# convenient definitions for interactive use
from pyre.handbook.units.SI import meter, second, mole, kelvin
from pyre.handbook.units.length import cm
from pyre.handbook.units.energy import cal, kcal
from pyre.handbook.units.pressure import atm

T = 1000.0 * kelvin
R = pyre.handbook.constants.fundamental.gas_constant


timer = pyre.timers.timer("meshanism")


def test(options):
    mechanism = load(options)
    mixture = prepareMixture(mechanism)
    calc = calculator(mixture, mechanism)

    T = 1000.0 * kelvin
    progressRate(calc, T)

    nasty = mismatches(calc)
    for reaction in nasty:
        reactionDetails(reaction)

    validateThermo(mixture, T)
    validateRates(calc)

    return mechanism, calc


def progressRate(calc, T):
    print " --- computing the rate of progress at T = %s" % str(T),

    timer.reset()
    timer.start()

    i = 0
    for reaction in calc.reaction():
        i += 1
        reaction.update(T)

    timer.stop()
    print "... done in %g sec" % timer.read()
    #print

    return


def reactionDetails(reaction):

    header = "Reaction %d: " % reaction.id() + str(reaction)
    line = "-"*len(header)

    print
    print header
    print line


    print
    print "phase space dimension:", reaction.equlibrium.dimPhaseSpace
    print "equilibrium constant:", reaction.equlibrium(T)
    print

    f = reaction.forwardRate
    print "forward:"
    print "        rate calculator:", f
    print "              rate type:", f.rateCalculator
    print
    print "          progress rate:", f()
    print
    print "            enhancement:", f.enhancement
    print "           phase factor:", f.phaseSpace
    print "          reaction rate:", f.rate
    #print "               Preduced:", f.rateCalculator.P_red
    print

    r = reaction.reverseRate
    print "reverse:"
    print "        rate calculator:", r
    print "              rate type:", r.rateCalculator
    print
    print "          progress rate:", r()
    print
    print "            enhancement:", r.enhancement
    print "           phase factor:", r.phaseSpace
    print "   equilibrium constant:", reaction.equlibrium(T)
    print "          reaction rate:", r.rate
    print

    q = f() - r()
    from chemkin_rates import chemkin_rates

    n = reaction.id()
    cq, cK_c = chemkin_rates[n-1]
    cK_c *= (cm**3 / mole)
    cq *= (mole/cm**3)/second

    print "chemkin:"
    print "    K_c = %s" % cK_c
    print "      q = %s" % cq
    print "  ratio = %g, error = %g" % (q/cq, abs(1.0 - q/cq))

    return


def validateRates(calc):
    
    width = 14
    precision = width - 7

    header1 = "".center(4) \
             + '   ' + "".center(40) \
             + ' ' + 'progress rate'.center(width+2) \
             + ' ' + 'chemkin value'.center(width+2)
             
    header2 = "Id".center(4) \
              + ' | ' + "Reaction".center(40) \
              + ('|' + 'mole/cm^3/s'.center(width+2))*2 \
              + '|' + 'relative error'.center(width+2)

    line = "-----+-" + "-"*40 + ("+" + "-"*(width+2))*3

    print header1
    print header2
    print line

    i = 0
    format = "%4d : %-40s" + ("| %%%d.%de " % (width, precision))*3
    #format = "%4d : %-40s" + "| %s "*2

    from chemkin_rates import chemkin_rates

    for reaction in calc.reaction():
        i += 1
        q_f, q_r = reaction.progressRate()

        if q_f.derivation != (mole/meter**3/second).derivation:
            pyre.debug.Firewall.hit("bad units in forward progress rate for reaction %d" % i)

        if q_r.derivation != (mole/meter**3/second).derivation:
            pyre.debug.Firewall.hit("bad units in reverse progress rate for reaction %d" % i)

        progress = 1e-6 * (q_f - q_r).value
        chemkin, k_c = chemkin_rates[i-1]
        err = abs((progress-chemkin)/chemkin)

        print format % (i, str(reaction), progress, chemkin, err)

    return


def mismatches(calc):
    from chemkin_rates import chemkin_rates

    candidates = []

    for reaction in calc.reaction():
        i = reaction.id()
        q_f, q_r = reaction.progressRate()

        if q_f.derivation != (mole/meter**3/second).derivation:
            pyre.debug.Firewall.hit("bad units in forward progress rate for reaction %d" % i)

        if q_r.derivation != (mole/meter**3/second).derivation:
            pyre.debug.Firewall.hit("bad units in reverse progress rate for reaction %d" % i)

        progress = 1e-6 * (q_f - q_r).value
        chemkin, k_c = chemkin_rates[i-1]
        err = abs((progress-chemkin)/chemkin)

        if err > 1e-2:
            print " ### reaction %d: mismatch = %g" % (i, err)
            candidates.append(reaction)

    return candidates



def prepareMixture(mechanism):
    print " --- preparing a mixture",

    timer.reset()
    timer.start()
    mix = pyre.chemistry.mixture(mechanism)

    from pyre.handbook.units.SI import *
    c = mole / (centi*meter)**3
    
    for species in mix.find():
        species.updateConcentration(c)
    mix.find("<mixture>").updateConcentration((mix.size() - 1) * c)

    #mix.find("O").updateConcentration(c)
    #mix.find("CO").updateConcentration(c)
    #mix.find("<mixture>").updateConcentration(2*c)
    
    timer.stop()
    print "... done in %g sec" % timer.read()
    print " +++ mixture: concentration", mix.find("<mixture>").concentration()
    #print

    return mix


def calculator(mix, mechanism):

    print " --- coverting mechanism to a rate calculator",
    timer.reset()
    timer.start()
    calc = pyre.chemistry.mechanism(mix, mechanism)
    timer.stop()
    print "... done in %g sec" % timer.read()

    return calc



def validateThermo(mixture, T):

    msg = "Thermal properties at T=%s" % T
    print msg
    print "-" * len(msg)

    print "Tabulating thermal properties"
    table = []
    timer.reset()
    timer.start()
    for species in mixture.find():
        symbol = species.symbol()
        if symbol == "<mixture>":
            continue
        cp_R = species.cp_R(T)
        h_RT = species.h_RT(T)
        s_R = species.s_R(T)
        g_RT = species.g_RT(T)
        table.append((symbol, cp_R, h_RT, s_R, g_RT))

    timer.stop()
    print "Tabulation done in %g sec" % timer.read()
    print

    width = 18
    precision = width - 7
    
    header = "Species".center(12) \
             + '|' + "cp_R".center(width+2) \
             + '|' + "h_RT".center(width+2) \
             + '|' + "s_R".center(width+2) \
             + '|' + "g_RT".center(width+2)
    line = "-"*12 + ("+" + "-"*(width+2))*4

    print header
    print line
    format = "  %-10s" + ("| %%%d.%df " % (width, precision))*4
    for symbol, cp_R, h_RT, s_R, g_RT in table:
        print  format % (symbol, cp_R, h_RT, s_R, g_RT)

    print

    print "Verifying against cantera data"
    print

    print header
    print line
    from cantera_thermo import cantera
    format = "  %-10s" + ("| %%%d.%de " % (width, precision))*4
    for symbol, mcp_R, mh_RT, ms_R, mg_RT in table:
        ccp_R, ch_RT, cs_R, cg_RT = cantera[symbol]
        
        dcp_R = abs((mcp_R - ccp_R)/mcp_R)
        dh_RT = abs((mh_RT - ch_RT)/mh_RT)
        ds_R = abs((ms_R - cs_R)/ms_R)
        dg_RT = abs((mg_RT - cg_RT)/mg_RT)

        print format % (symbol, dcp_R, dh_RT, ds_R, dg_RT)

    print

    return


def load(options):

    file = options["--file"]
    format = options["--format"]

    print " --- loading mechanism '%s'" % file,

    timer.start()
    mechanism = pyre.chemistry.serialization.load(file, format)
    timer.stop()

    print "... done in %g sec" % timer.read()
    #mechanism.printStatistics()

    return mechanism


# usage

def usage(program):
    print "Usage: %s [options ...]" % program
    print "Options: (default values in brackets)"
    print "    --file=<mechanism filename> [%s]" % defaults["--file"]
    print "    --format=<chemkin|ckml> [%s]" % defaults["--file"]
    print
    return
        

defaults = {
    "--file": "GRIMech-3.0.ck2",
    "--format": "chemkin"
    }

# main

if __name__ == "__main__":
        
    # pyre.support.debug.activate("pyre.chemistry.serialization")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-scanner")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-parser")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-tokenizer")
    
    import pyre

    options = pyre.applications.main(defaults, usage)
    mechanism, calc = test(options)
    

# version
__id__ = "$Id$"

#
# End of file
