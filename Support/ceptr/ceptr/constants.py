import cantera as ct
from pint import UnitRegistry

smallnum = 1e-100
ureg = UnitRegistry()
kb = (ct.boltzmann * ureg.joule / ureg.kelvin).to(ureg.erg / ureg.kelvin).m
R = 8.31446261815324e7 * (ureg.erg / (ureg.mole / ureg.kelvin))
Rc = 1.98721558317399615845 * ((ureg.cal / ureg.mole) / ureg.kelvin)
Patm = 1013250.0
Patm_pa = Patm / 10.0
sym = ""
