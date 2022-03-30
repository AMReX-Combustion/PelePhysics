#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import

from .unit import dimensionless, unit

# basic SI units

meter = unit(1.0, (1, 0, 0, 0, 0, 0, 0))
kilogram = unit(1.0, (0, 1, 0, 0, 0, 0, 0))
second = unit(1.0, (0, 0, 1, 0, 0, 0, 0))
ampere = unit(1.0, (0, 0, 0, 1, 0, 0, 0))
kelvin = unit(1.0, (0, 0, 0, 0, 1, 0, 0))
mole = unit(1.0, (0, 0, 0, 0, 0, 1, 0))
candela = unit(1.0, (0, 0, 0, 0, 0, 0, 1))

# the 22 derived SI units with special names

radian = dimensionless  # plane angle
steradian = dimensionless  # solid angle

hertz = 1 / second  # frequency

newton = meter * kilogram / second**2  # force
pascal = newton / meter**2  # pressure
joule = newton * meter  # work, heat
watt = joule / second  # power, radiant flux

coulomb = ampere * second  # electric charge
volt = watt / ampere  # electric potential difference
farad = coulomb / volt  # capacitance
ohm = volt / ampere  # electric resistance
siemens = ampere / volt  # electric conductance
weber = volt * second  # magnetic flux
tesla = weber / meter**2  # magnetic flux density
henry = weber / ampere  # inductance

celcius = kelvin  # Celcius temperature

lumen = candela * steradian  # luminus flux
lux = lumen / meter**2  # illuminance

becquerel = 1 / second  # radioactivity
gray = joule / kilogram  # absorbed dose
sievert = joule / kilogram  # dose equivalent
katal = mole / second  # catalytic activity

# prefixes

yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1000
hecto = 100
deka = 10
deci = 0.1
centi = 0.01
milli = 0.001
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21
yocto = 1e-24


# version
__id__ = "$Id$"

#
# End of file
