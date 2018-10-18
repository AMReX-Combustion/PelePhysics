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

from weaver.mills.PythonMill import PythonMill

from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ, erg
from pyre.handbook.constants.fundamental import avogadro
from pyre.handbook.constants.fundamental import gas_constant as R

smallnum = 1e-100
R = 8.314e7 * erg/mole/kelvin
Rc = 1.987 * cal/mole/kelvin
Patm = 1013250.0

class speciesDb:
    def __init__(self, id, name, mwt):
        self.id = id
        self.symbol = name
        self.weight = mwt
        return

class PythonPickler(PythonMill):


    def __init__(self):
        PythonMill.__init__(self)
        self.species = [] 
        return
 
    def _setSpecies(self, mechanism):
        """ For internal use """
        import pyre
        periodic = pyre.handbook.periodicTable()
        
        nSpecies = len(mechanism.species())
        self.species = [ 0.0 for x in range(nSpecies) ]
        
        for species in mechanism.species():
            weight = 0.0 
            for elem, coef in species.composition:
                aw = mechanism.element(elem).weight
                if not aw:
                    aw = periodic.symbol(elem.capitalize()).atomicWeight
                weight += coef * aw

            tempsp = speciesDb(species.id, species.symbol, weight)
            self.species[species.id] = tempsp
        return

    def _renderDocument(self, mechanism, options=None):

        self._setSpecies(mechanism)
        self._rep += ['']
	self._write(self.line('Mechanism is summarized as follows'))
        #self._analyzeThermodynamics(mechanism)
        self._rendersummary(mechanism)
        self._renderspeciesDict(mechanism)
        self._renderzeroVect(mechanism)
        self._rendernormalizeVect()
        self._rendernormalizeVectConst()
        
        self._renderpx(mechanism)
        self._renderpy(mechanism)
        self._renderpc(mechanism)
        self._renderrhox(mechanism)
        self._renderrhoy(mechanism)
        self._renderrhoc(mechanism)
        self._renderwt(mechanism)
        self._rendermmwx(mechanism)
        self._rendermmwy(mechanism)
        self._rendermmwc(mechanism)
        
        self._renderxty(mechanism)
        self._renderxtcp(mechanism)
        self._renderxtcr(mechanism)
        self._renderytx(mechanism)
        self._renderytphi(mechanism)
        self._renderytcp(mechanism)
        self._renderytcr(mechanism)
        self._renderctx(mechanism)
        self._rendercty(mechanism)
        self._renderthermo(mechanism)
        
        self._rendercvml(mechanism)
        self._rendercpml(mechanism)
        self._renderuml(mechanism)
        self._renderhml(mechanism)
        self._rendergml(mechanism)
        self._renderaml(mechanism)
        self._rendersml(mechanism)
        
        self._rendercvms(mechanism)
        self._rendercpms(mechanism)
        self._renderums(mechanism)
        self._renderhms(mechanism)
        self._rendergms(mechanism)
        self._renderams(mechanism)
        self._rendersms(mechanism)

        self._rendercpbl(mechanism)
        self._rendercpbs(mechanism)
        self._rendercvbl(mechanism)
        self._rendercvbs(mechanism)
        
        self._renderhbml(mechanism)
        self._renderhbms(mechanism)
        self._renderubml(mechanism)
        self._renderubms(mechanism)
        self._rendersbml(mechanism)
        self._rendersbms(mechanism)
        self._rendergbml(mechanism)
        self._rendergbms(mechanism)
        self._renderabml(mechanism)
        self._renderabms(mechanism)

        self._renderprogressRate(mechanism)
        self._renderproductionRate(mechanism)
        
        self._renderwc(mechanism)
        self._renderqc(mechanism)
        
        self._render_phity(mechanism)
        self._render_eytt(mechanism)
        self._render_hytt(mechanism)
        
        self._rendermain()
        self._rep += ['']

        return

    def _rendermain(self):
        self._write()
        self._write()
        self._write('if __name__ == "__main__": ')
        self._indent()
        self._write('x = zeroVect(1.0)')
        self._write('y = xty(x)')
        self._write('x = ytx(y)')
        self._write('print "after a ytx and xty"')
        self._write('for i in range(min(6,len(x))): ')
        self._indent()
        self._write('print "x[%d] = %20.15e " % (i,x[i]) ')
        self._outdent()
        self._write('phi = ytphi(y)')
        self._write('y = phity(phi)')
        self._write('x = ytx(y)')
        self._write("print 'after a ytphi, phity and then ytx'")
        self._write('for i in range(min(6,len(x))): ')
        self._indent()
        self._write('print "x[%d] = %20.15e " % (i,x[i]) ')
        self._outdent()
        self._write('c = ytcr(1.0,2.0,y)')
        self._write('y = cty(c)')
        self._write('x = ytx(y)')
        self._write("print 'after a ytcr, then cty, then ytx'")
        self._write('for i in range(min(6,len(x))): ')
        self._indent()
        self._write('print "x[%d] = %20.15e " % (i,x[i]) ')
        self._outdent()
        self._write('x = zeroVect()')
        self._write('x[species["H2"]] = 2')
        self._write('x[species["O2"]] = 1')
        self._write('x[species["N2"]] = 3.72')
        self._write('y=xty(x)')
        self._write('x=ytx(y)')
        self._write('ene = 1.28e10')
        self._write('t=eytt(ene,y)')
        self._write('print "Temperature from eytt is: %g (target ~ 1541.4 for GRImech)" % t')
        self._write('h=hbms(t,y)')
        self._write('t=hytt(h,y)')
        self._write('print "Temperature from hytt is: %g" % t')
        self._write('c = ytcr(1.0,2.0,y)')
        self._write('print "mean mwt from x: %g" % mmwx(x)')
        self._write('print "mean mwt from y: %g" % mmwy(y)')
        self._write('print "mean mwt from c: %g" % mmwc(c)')
        
        self._outdent()
        return

    def _rendersummary(self,mechanism):
        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        nElements = len(mechanism.element())
        self._write("source = '%s'" % mechanism.name())
        self._write('ns = %d' % nSpecies)
        self._write('nr = %d' % nReactions)
        self._write('ne = %d' % nElements)
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        self._write('tmin = %g' % lowT)
        self._write('tmax = %g' % highT)
        
        return
        
    
    def _renderspeciesDict(self,mechanism):
        self._write('species = {}')
	knt = 0
	for species in mechanism.species():
		self._write("species['%s'] = %d" % (species.symbol, knt) )
		self._write("species[%d] = '%s'" % (knt, species.symbol) )
		knt = knt + 1
        self._write()
        return

    def _renderzeroVect(self,mechanism):
        nSpecies = len(mechanism.species())
        self._write("def zeroVect( val = 0.0):")
        self._indent()
        self._write("return [ val for x in range(%d) ]" % nSpecies)
        self._outdent()
        self._write()
        return

    def _rendernormalizeVect(self):
        self._write("def normalizeVect(vect):")
        self._indent()
        self._write('nn = len(vect)')
        self._write('xx = 0.0')
        self._write('for yy in vect:')
        self._indent()
        self._write('xx += yy')
        self._outdent()
        self._write('for id in range(nn):')
        self._indent()
        self._write('vect[id] = vect[id] / xx')
        self._outdent()
        self._write('return vect')
        self._outdent()
        self._write()
        return

    def _rendernormalizeVectConst(self):
        self._write("def normalizeVectConst(vect):")
        self._indent()
        self._write('nn = len(vect)')
        self._write('xx = 0.0')
        self._write('out = [ 0.0 for x in range(nn) ]')
        self._write('for yy in vect:')
        self._indent()
        self._write('xx += yy')
        self._outdent()
        self._write('for id in range(nn):')
        self._indent()
        self._write('out[id] = vect[id] / xx')
        self._outdent()
        self._write('return out')
        self._outdent()
        self._write()
        return

    def _renderwt(self,mechanism):
        self._write()
        docstring = 'takes no argument, returns vector of molecular wts'
        self._write(self.line(docstring))
        self._write("def wt():")
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('mwt = zeroVect()')
        for species in self.species:
            self._write('mwt[%d] = %f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol) )

        self._write('return mwt')
        self._outdent()
        return
   
    def _rendermmwy(self, mechanism):
        self._write()
        self._write()
        docstring = 'wtm = mmwy(y): returns mean molecular weight (gm/mole)'
        self._write(self.line(docstring))
        self._write('def mmwy(y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # molecular weights of all species
        self._write('YOW = 0.0')
        for species in self.species:
            self._write('YOW += y[%d]/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('wtm = 1.0 / YOW')
        
        self._write()
        self._write('return wtm')
        self._outdent()

        return 
 
    def _rendermmwx(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'wtm = mmwx(x): returns mean molecular weight (gm/mole)'
        self._write(self.line(docstring))
        self._write('def mmwx(x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # molecular weights of all species
        self._write('wtm = 0.0')
        for species in self.species:
            self._write('wtm += x[%d]*%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('return wtm')
        self._outdent()

        return 
 
    def _rendermmwc(self, mechanism):
        self._write()
        self._write()
        docstring = 'wtm = mmwc(c): returns mean molecular weight (gm/mole)'
        self._write('def mmwc(c):')
        self._indent()
        self._write('"""%s"""' % docstring)

        
        # molecular weights of all species
        self._write('W = 0.0')
        for species in self.species:
            self._write('W += c[%d]*%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        nSpecies = len(mechanism.species())
        self._write('sumC = 0.0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('sumC += c[id]')
        self._outdent()

        self._write(self.line('CK provides no guard against divison by zero'))
        self._write('wtm = W/sumC')
        
        self._write()
        self._write('return wtm')
        self._outdent()

        return 
 
    def _renderxty(self,mechanism):
        self._write()
        docstring = 'convert x[species] (mole fracs) to y[species] (mass fracs)'
        self._write(self.line(docstring))
        self._write('def xty(x):')
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('y = zeroVect()')
        
        # now compute conversion
        self._write(self.line('compute conversion (needs normalization)'))
        for species in self.species:
            self._write('y[%d] = x[%d]*%f ' % ( species.id, species.id, species.weight) +
                        self.line('%s' % species.symbol) )

        self._write()
        self._write('return normalizeVect(y)')
        self._outdent()

        return

    def _renderytx(self, mechanism):
        self._write()
        self._write()
        docstring = 'convert y[species] (mass fracs) to x[species] (mole fracs)'
        self._write(self.line(docstring))
        self._write('def ytx(y):')
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('x = zeroVect()')
        
        self._write(self.line('compute conversion (needs normalization)'))
        for species in self.species:
            self._write('x[%d] = y[%d]/%f ' % (
                species.id, species.id, species.weight) +
                        self.line('%s' % species.symbol) )

        self._write()
        self._write('return normalizeVect(x)')
        self._outdent()

        return 
    
    def _renderytphi(self, mechanism):
        self._write()
        self._write()
        docstring = 'convert y[species] (mass fracs) to phi[species] (specific mole num)'
        self._write(self.line(docstring))
        self._write('def ytphi(y):')
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('phi = zeroVect()')
        
        for species in self.species:
            self._write('phi[%d] = y[%d]*1000.0/%f  ' % (
                species.id, species.id, species.weight) +
                        self.line('%s (wt in kg)' % species.symbol))
 
        self._write()
        self._write('return phi')
        self._outdent()

        return
 
    def _render_phity(self, mechanism):
        
        self._write()
        self._write()
        docstring = 'convert phi[species] (specific mole nums) to y[species] (mass fracs)'
        self._write(self.line(docstring))
        self._write('def phity(phi):')
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('y   = zeroVect()')
        
        for species in self.species:
            self._write('y[%d] = phi[%d]*%f ' % (species.id, species.id, species.weight) +
                        self.line('%s' % species.symbol ))
            
        self._write()
        self._write('return normalizeVect(y)')
        self._outdent()

        return 
 
    def _renderxtcp(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring =  'c = xtcp(P,T,x): convert x[species] (mole fracs) to c[species] (molar conc)'
        self._write(self.line(docstring))
        self._write('def xtcp(P,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)
        self._write('c = zeroVect()')
        self._write('PORT = P/(%g * T) ' % (R*kelvin*mole/erg) +
                    self.line('P/RT'))
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT')
        self._outdent()

        self._write()
        self._write('return c')
        self._outdent()

        return 
 
    def _renderxtcr(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'c = xtcr(rho,T,x): convert x[species] (mole fracs) to c[species] (molar conc)'
        self._write(self.line(docstring))
        self._write('def xtcr(rho, T, x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('c = zeroVect()')
        self._write(self.line('See Eq 4, 11 in CK Manual'))
        self._write('XW = 0')
        for species in self.species:
            self._write('XW += x[%d]*%f ' %
                        (species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write('ROW = rho / XW')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW')
        self._outdent()

        self._write()
        self._write('return c')
        self._outdent()

        return
 
    def _renderytcp(self, mechanism):

        self._write()
        self._write()
        docstring = 'c = ytcp(P,T,y): convert y[species] (mass fracs) to c[species] (molar conc)'
        self._write(self.line(docstring))
        self._write('def ytcp(P,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('c = zeroVect')
        self._write('YOW = 0')
        self._write('PWORT')
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
 
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = P/(YOW * %g * T) ' % (R*kelvin*mole/erg) )

        # now compute conversion
        self._write(self.line('Now compute conversion'))
        for species in self.species:
            self._write('c[%d] = PWORT * y[%d]/%f ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return c')
        self._outdent()

        return 
 
    def _renderytcr(self, mechanism):

        self._write()
        self._write()
        docstring = 'c = ytcr(rho,T,y): convert y[species] (mass fracs) to c[species] (molar conc)'
        self._write(self.line(docstring))
        self._write('def ytcr(rho, T, y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('c = zeroVect()')
        self._write(self.line('See Eq 8 (Temperature not used)'))
        for species in self.species:
            self._write('c[%d] = rho * y[%d]/%f  ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return c')
        self._outdent()

        return 
 
    def _renderctx(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'convert c[species] (molar conc) to x[species] (mole fracs)'
        self._write(self.line(docstring))
        self._write('def ctx(c):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write()
        self._write('return normalizeVectConst(c)')
        self._outdent()

        return 
 
    def _rendercty(self, mechanism):

        self._write()
        self._write()
        docstring = 'convert c[species] (molar conc) to y[species] (mass fracs)'
        self._write(self.line(docstring))
        self._write('def cty(c):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('y = zeroVect()')
        
        for species in self.species:
            self._write('y[%d] = c[%d]*%f ' % ( species.id, species.id, species.weight) )

        self._write()
        self._write('return normalizeVect(y)')
        self._outdent()

        return 

        
    def _renderpx(self, mechanism):
        self._write()
        self._write()
        docstring = 'P = px(rho, T, x): Compute P = rhoRT/W(x)'
        self._write(self.line(docstring))
        self._write('def px(rho,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # molecular weights of all species
        self._write('XW = 0.0')
        for species in self.species:
            self._write('XW += x[%d]*%f  ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('P = rho * %g * T / XW ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return P')
        self._outdent()
        
        return

    def _renderpy(self, mechanism):
        self._write()
        self._write()
        docstring = 'P = py(rho,T,y): Compute P = rhoRT/W(y)'
        self._write(self.line(docstring))
        self._write('def py(rho,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # molecular weights of all species
        self._write('YOW = 0.0')
        for species in self.species:
            self._write('YOW += y[%d]/%f  ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self.line('YOW holds the reciprocal of the mean molecular wt')
        self._write('P = rho * %g * T * YOW ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return P')
        self._outdent()

        return 
 
    def _renderpc(self, mechanism):
        self._write()
        self._write()
        docstring = 'P = pc(rho,T,c): Compute P = rhoRT/W(c)'
        self._write(self.line(docstring))
        self._write('def pc(rho,T,c):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write(self.line('See Eq 5 in CK Manual'))
        
        # molecular weights of all species
        self._write('W = 0.0')
        for species in self.species:
            self._write('W += c[%d]*%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        nSpecies = len(mechanism.species())
        self._write('sumC = 0.0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('sumC += c[id]')
        self._outdent()

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            'P = rho * %g * T * sumC / W ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return P')
        self._outdent()

        return 

    def _renderrhox(self, mechanism):
        self._write()
        self._write()
        docstring = 'rho = rhox(P,T,x): Compute rho = PW(x)/RT'
        self._write(self.line(docstring))
        self._write('def rhox(P,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # molecular weights of all species
        self._write('XW = 0.0')
        for species in self.species:
            self._write('XW += x[%d]*%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(
            'rho = P * XW / (%g * T) ' % (R*kelvin*mole/erg)
            + self.line('rho = P*W/(R*T)'))
        
        self._write()
        self._write('return rho')
        self._outdent()
        
        return

    def _renderrhoy(self, mechanism):
        self._write()
        self._write()
        docstring = 'rho = rhoy(P,T,y): Compute rho = P*W(y)/RT'
        self._write(self.line(docstring))
        self._write('def rhoy(P,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        
        self._write('YOW = 0.0')
        # molecular weights of all species
        for species in self.species:
            self._write('YOW += y[%d]/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self.line('YOW holds the reciprocal of the mean molecular wt')
        self._write(
            'rho = P / (%g * T * YOW) ' % (R*kelvin*mole/erg)
            + self.line('rho = P*W/(R*T)'))
        
        
        self._write()
        self._write('return rho')
        self._outdent()

        return 
 
    def _renderrhoc(self, mechanism):
        self._write()
        self._write()
        docstring = 'rho = rhoc(P,T,c): Compute rho = P*W(c)/(R*T)'
        self._write(self.line(docstring))
        self._write('def rhoc(P,T,c):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write(self.line('See Eq 5 in CK Manual'))
        
        # molecular weights of all species
        self._write('W = 0.0')
        for species in self.species:
            self._write('W += c[%d]*%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        nSpecies = len(mechanism.species())
        self._write('sumC = 0.0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('sumC += c[id]')
        self._outdent()

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            'rho = P * W / (sumC * T * %g) ' % (R*kelvin*mole/erg)
            + self.line('rho = PW/(R*T)'))
        
        self._write()
        self._write('return rho')
        self._outdent()

        return

    
    def _generateThermoRoutine(self, name, expressionGenerator, speciesInfo):

        lowT, highT, midpoints = speciesInfo
        
        self._write('def %s(T):' % name)
        
        self._indent()

        # declarations
        self._write()
        self._write('import math')
        self._write(self.line('temperature'))
        self._write('T = T*1.0')
        self._write('tc = [ math.log(T), T, T*T, T*T*T, T*T*T*T ]')
        self._write('species = zeroVect()')

        # temperature check
        # self._write()
        # self._write(self.line('check the temperature value'))
        # self._write('if (T < %g || T > %g) :' % (lowT, highT))
        # self._indent()
        # self._write(
        #     'print('temperature %%g is outside the range (%g, %g)", T) '
        #     % (lowT, highT))
        # self._write('return ')
        # self._outdent()
                    
        for midT, speciesList in midpoints.items():

            self._write('')
            self._write(self.line('species with midpoint at T=%g kelvin' % midT))
            self._write('if (T < %g): ' % midT)
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] = (' % species.id)
                self._indent()
                expressionGenerator(lowRange.parameters)
                self._outdent()
                self._write(')')

            self._outdent()
            self._write('else :')
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] = (' % species.id)
                self._indent()
                expressionGenerator(highRange.parameters)
                self._outdent()
                self._write(')')

            self._outdent()
            
        self._write('return species')
        self._outdent()

        return

    def _analyzeThermodynamics(self, mechanism):
        # Copied from CPickler
        lowT = 0.0
        highT = 1000000.0

        midpoints = {}

        for species in mechanism.species():
            models = species.thermo
            if len(models) > 2:
                import pyre
                pyre.debug.Firewall.hit("unsupported configuration in species.thermo")
                return
            
            m1 = models[0]
            m2 = models[1]

            if m1.lowT < m2.lowT:
                lowRange = m1
                highRange = m2
            else:
                lowRange = m2
                highRange = m1

            low = lowRange.lowT
            mid = lowRange.highT
            high = highRange.highT

            if low > lowT:
                lowT = low
            if high < highT:
                highT = high

            midpoints.setdefault(mid, []).append((species, lowRange, highRange))
            
        return lowT, highT, midpoints


    def _Kc(self, mechanism, reaction):

        dim = 0
        dG = ""

        terms = []
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
                    
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim -= coefficient
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in reaction.products:
            if coefficient == 1:
                factor = ""
            else:
                factor = "%d * " % coefficient
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim += coefficient
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'math.exp(' + dG + ')'

        if dim == 0:
            conversion = ""
        elif dim > 0:
            conversion = "*".join(["refC"] * dim) + ' * '
        else:
            conversion = "1.0 / (" + "*".join(["refC"] * abs(dim)) + ') * '

        K_c = conversion + K_p

        return K_c


    def _cpNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4] ' % parameters[4])
        return


    def _cvNASA(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4] ' % parameters[4])
        return


    def _enthalpyNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e / tc[1] ' % (parameters[5]))
        return


    def _internalEnergy(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e / tc[1] ' % (parameters[5]))
        return

    
    def _gibbsNASA(self, parameters):
        self._write('%+15.8e / tc[1]' % parameters[5])
        self._write('%+15.8e' % (parameters[0] - parameters[6]))
        self._write('%+15.8e * tc[0]' % (-parameters[0]))
        self._write('%+15.8e * tc[1]' % (-parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (-parameters[2]/6))
        self._write('%+15.8e * tc[3]' % (-parameters[3]/12))
        self._write('%+15.8e * tc[4] ' % (-parameters[4]/20))
        return
    
    def _helmholtzNASA(self, parameters):
        self._write('%+15.8e / tc[1]' % parameters[5])
        self._write('%+15.8e' % (parameters[0] - parameters[6] - 1.0))
        self._write('%+15.8e * tc[0]' % (-parameters[0]))
        self._write('%+15.8e * tc[1]' % (-parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (-parameters[2]/6))
        self._write('%+15.8e * tc[3]' % (-parameters[3]/12))
        self._write('%+15.8e * tc[4] ' % (-parameters[4]/20))
        return

    def _entropyNASA(self, parameters):
        self._write('%+15.8e * tc[0]' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]))
        self._write('%+15.8e * tc[2]' % (parameters[2]/2))
        self._write('%+15.8e * tc[3]' % (parameters[3]/3))
        self._write('%+15.8e * tc[4]' % (parameters[4]/4))
        self._write('%+15.8e ' % (parameters[6]))
        return


    def _cv(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cv/R at the given temperature'))
        self._generateThermoRoutine("cv_R", self._cvNASA, speciesInfo)

        return
    
    def _cp(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cp/R at the given temperature'))
        self._generateThermoRoutine("cp_R", self._cpNASA, speciesInfo)

        return


    def _gibbs(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the g/(RT) at the given temperature'))
        self._generateThermoRoutine("gibbs", self._gibbsNASA, speciesInfo)

        return

    def _helmholtz(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the a/(RT) at the given temperature'))
        self._generateThermoRoutine("helmholtz", self._helmholtzNASA, speciesInfo)

        return

    def _speciesEntropy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the S/R at the given temperature (Eq 21)'))
        self._generateThermoRoutine("speciesEntropy", self._entropyNASA, speciesInfo)

        return

    def _speciesInternalEnergy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the e/(RT) at the given temperature'))
        self._generateThermoRoutine("speciesInternalEnergy", self._internalEnergy, speciesInfo)

        return

    def _speciesEnthalpy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the h/(RT) at the given temperature (Eq 20)'))
        self._generateThermoRoutine("speciesEnthalpy", self._enthalpyNASA, speciesInfo)

        return

    def _renderthermo(self, mechanism):
        speciesInfo = self._analyzeThermodynamics(mechanism)

        self._gibbs(speciesInfo)
        self._helmholtz(speciesInfo)
        self._cv(speciesInfo)
        self._cp(speciesInfo)
        self._speciesInternalEnergy(speciesInfo)
        self._speciesEnthalpy(speciesInfo)
        self._speciesEntropy(speciesInfo)

        return
      
    def _rendercvml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get specific heat at constant volume as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def cvml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cvor = cv_R(T)')
        
        # convert cv/R to cv
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('cvor[id] *= %g' % (R*kelvin*mole/erg) )
        self._outdent()

        self._write('return cvor')
        self._outdent()

        return
       
    def _rendercpml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get specific heat at constant pressure as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def cpml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cpor = cp_R(T)')
        
        # convert cp/R to cp
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('cpor[id] *= %g' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('return cpor') 
        self._outdent()

        return
     
    def _renderuml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get internal energy as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def uml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write( 'RT = %15.8e*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('uml = speciesInternalEnergy(T)')
        
        # convert e/RT to e with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('uml[id] *= RT')
        self._outdent()
        self._write('return uml') 
        self._outdent()

        return
      
    def _renderhml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get enthalpy as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def hml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write( 'RT = %15.8e*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('hml = speciesEnthalpy(T)')
        
        # convert h/RT to h with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('hml[id] *= RT')
        self._outdent()

        self._write('return hml')
        self._outdent()

        return
    
    def _rendergml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get standard-state Gibbs energy as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def gml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %15.8e*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('gml = gibbs(T)')
        
        # convert g/RT to g with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('gml[id] *= RT')
        self._outdent()
        self._write('return gml')
        self._outdent()

        return
    
    def _renderaml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'get standard-state Helmholtz free energy as a function of T for all species (molar units)'
        self._write(self.line(docstring))
        self._write('def aml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %15.8e*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('aml = helmholtz(T)')
        
        # convert A/RT to A with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('aml[id] *= RT')
        self._outdent()
        self._write('return aml')
        self._outdent()

        return
   
    def _rendersml(self, mechanism):

        nSpecies = len(mechanism.species())
        
        self._write()
        self._write()
        docstring = 'Returns the standard-state entropies in molar units'
        self._write(self.line(docstring))
        self._write('def sml(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('sml = speciesEntropy(T)')
        
        # convert s/R to s
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('sml[id] *= %g' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('return sml')
       
        self._outdent()

        return
 
    def _renderproductionRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write()
        docstring = 'wdot = productionRate( concentration, T ): compute species production rates'
        self._write(self.line(docstring))
        self._write('def productionRate(sc, T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # declarations
        self._initializeRateCalculation(mechanism)

        # clear out wdot
        self._write()
        self._write('wdot = zeroVect()')
        
        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            # compute the rates
            self._forwardRate(mechanism, reaction)
            self._reverseRate(mechanism, reaction)

            # store the progress rate
            self._write("qdot = q_f - q_r")

            for symbol, coefficient in reaction.reactants:
                self._write(
                    "wdot[%d] -= %d * qdot" % (mechanism.species(symbol).id, coefficient))

            for symbol, coefficient in reaction.products:
                self._write(
                    "wdot[%d] += %d * qdot" % (mechanism.species(symbol).id, coefficient))

        self._write()
        self._write('return wdot')
        self._outdent()

        return


    def _renderprogressRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write()
        docstring = 'qdot = progressRate( conc, T): compute the progress rate for each reaction'
        self._write(self.line(docstring))
        self._write('def progressRate(sc, T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # declarations
        self._initializeRateCalculation(mechanism)

        self._write(self.line('qdot has dimension: number of reactions'))
        self._write('qdot = [ 0.0 for x in range(nr) ]')
        
        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            # compute the rates
            self._forwardRate(mechanism, reaction)
            self._reverseRate(mechanism, reaction)

            # store the progress rate
            self._write("qdot[%d] = q_f - q_r" % (reaction.id - 1))

        self._write()
        self._write('return qdot')
        self._outdent()

        return


    def _initializeRateCalculation(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # compute the reference concentration
        self._write()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('import math')
        self._write('T =  T*1.0')
        self._write('tc = [ math.log(T), T, T*T, T*T*T, T*T*T*T ]')
        self._write('refC = %g / %g / T' % (atm.value, R.value))

        # compute the mixture concentration
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('mixture = 0.0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('mixture += sc[id]')
        self._outdent()

        # compute the Gibbs free energies
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('g_RT = gibbs(T)')
        
        return


    def _forwardRate(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._phaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s" % phi_f)

        arrhenius = self._arrhenius(reaction, reaction.arrhenius)

        thirdBody = reaction.thirdBody
        if not thirdBody:
            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim)
            self._write("k_f = %g * %s" % (uc.value, arrhenius))
            self._write("q_f = phi_f * k_f")
            return
            
        alpha = self._enhancement(mechanism, reaction)
        self._write("alpha = %s" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            uc = self._prefactorUnits(reaction.units["prefactor"], -dim)
            self._write("k_f = %g * alpha * %s" % (uc.value, arrhenius))
            self._write("q_f = phi_f * k_f")
            return

        uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim)
        self._write("k_f = %g * %s" % (uc.value, arrhenius))
        k_0 = self._arrhenius(reaction, reaction.low)
        redP = "alpha / k_f * " + k_0
        self._write("redP = 1e-12 * %s" % redP)
        self._write("F = redP / (1 + redP)")

        if sri:
            self._write("logPred = math.log10(redP)")

            self._write("X = 1.0 / (1.0 + logPred*logPred)")

            SRI = "math.exp(X * math.log(%g*math.exp(%g/T) + math.exp(T/%g))" % (sri[0], -sri[1], -sri[2])
            if len(sri) > 3:
                SRI += " * %g * math.exp(%g*tc[0])" % (sri[3], sri[4])

            self._write("F_sri = %s" % SRI)
            self._write("F *= Ftroe")

        elif troe:
            self._write("logPred = math.log10(redP)")

            logF_cent = "logFcent = math.log10("
            logF_cent += "(%g*math.exp(T/%g))" % (1-troe[0], -troe[1])
            logF_cent += "+ (%g*math.exp(T/%g))" % (troe[0], -troe[2])
            if len(troe) == 4:
                logF_cent += "+ (math.exp(%g/T))" % (-troe[3])
            logF_cent += ')'
            self._write(logF_cent)
            
            d = .14
            self._write("troe_c = -.4 - .67 * logFcent")
            self._write("troe_n = .75 - 1.27 * logFcent")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred))")
            self._write("F_troe = math.pow(10, logFcent / (1.0 + troe*troe))")
            self._write("F *= F_troe")

        self._write("k_f *= F")
        self._write("q_f = phi_f * k_f")
        return
        

    def _reverseRate(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r = 0.0")
            return

        phi_r = self._phaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:
            arrhenius = self._arrhenius(reaction, reaction.rev)
            thirdBody = reaction.thirdBody
            if thirdBody:
                uc = self._prefactorUnits(reaction.units["prefactor"], -dim)
                self._write("k_r = %g * alpha * %s" % (uc.value, arrhenius))
            else:
                uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim)
                self._write("k_r = %g * %s" % (uc.value, arrhenius))

            self._write("q_f = phi_r * k_r")
            return
        
        K_c = self._Kc(mechanism, reaction)
        self._write("Kc = %s" % K_c)

        self._write("k_r = k_f / Kc")
        self._write("q_r = phi_r * k_r")

        return


    def _arrhenius(self, reaction, parameters):
        A, beta, E = parameters
        if A == 0:
            return "0.0"
        expr = "%g" % A
        if beta == 0 and E == 0:
            return expr
        expr +="*math.exp("
        if beta != 0:
            expr += "%g*tc[0]" % beta
        if E != 0:
            uc = self._activationEnergyUnits(reaction.units["activation"])
            expr += "%+g/tc[1]" % (- uc * E / Rc / kelvin) # catch unit conversion errors!
        expr += ')'
        
        return expr


    def _prefactorUnits(self, code, exponent):

        if code == "mole/cm**3":
            units = mole / cm**3
        elif code == "moles":
            units = mole / cm**3
        elif code == "molecules":
            import pyre
            units = 1.0 / avogadro / cm**3
        else:
            import pyre
            pyre.debug.Firewall.hit("unknown prefactor units '%s'" % code)
            return 1

        return units ** exponent / second


    def _activationEnergyUnits(self, code):
        if code == "cal/mole":
            units = cal / mole
        elif code == "kcal/mole":
            units = kcal /mole
        elif code == "joules/mole":
            units = J / mole
        elif code == "kjoules/mole":
            units = kJ / mole
        elif code == "kelvins":
            units = Rc * kelvin
        else:
            pyre.debug.Firewall.hit("unknown activation energy units '%s'" % code)
            return 1

        return units


    def _phaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in reagents:
            conc = "sc[%d]" % mechanism.species(symbol).id
            phi += [conc] * coefficient

        return "*".join(phi)


    def _phaseSpaceUnits(self, reagents):
        dim = 0
        for symbol, coefficient in reagents:
            dim += coefficient

        return dim


    def _enhancement(self, mechanism, reaction):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                return "mixture"
            return "sc[%d]" % mechanism.species(species).id

        alpha = ["mixture"]
        for symbol, efficiency in efficiencies:
            factor = efficiency - 1
            conc = "sc[%d]" % mechanism.species(symbol).id
            if factor == 1:
                alpha.append(conc)
            else:
                alpha.append("%g*%s" % (factor, conc))

        return " + ".join(alpha)


    def _renderwc(self, mechanism):

        nSpecies = len(mechanism.species())

        self._write()
        self._write()
        docstring = 'wdot = wc(T, conc): compute the production rate for each species'
        self._write('def wc(T, C):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6')
        self._outdent()
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('wdot = productionRate(C, T)')

        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6')
        self._write('wdot[id] *= 1.0e-6')
        self._outdent()
        self._write('return wdot')
        
        self._outdent()
        return


    def _renderqc(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write()
        docstring = 'qdot = qc(T,C): Returns the rate of progress for each reaction'
        self._write('def qc(T,C):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6')
        self._outdent()
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('qdot = progressRate(C, T)')

        # convert C to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6')
        self._outdent()

        # convert qdot to chemkin units
        self._write()
        self._write('for id in range(%d):' % nReactions)
        self._indent()
        self._write('qdot[id] *= 1.0e-6')
        self._outdent()
        self._write('return qdot')
        
        self._outdent()
        return

    def _render_eytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        docstring = 'T = eytt(e, y): get temperature given internal energy in mass units and mass fracs'
        self._write('def eytt(e, y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('maxiter = 50')
        self._write('tol  = 0.001')
        self._write('ein  = e')
        self._write('tmin = %g' % lowT + self.line('max lower bound for thermo def'))
        self._write('tmax = %g' % highT + self.line('min upper bound for thermo def'))
        self._write('emin = ubms(tmin, y)')
        self._write('emax = ubms(tmax, y)')
        self._write('if (ein < emin):')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('cv = cvbs(tmin, y)')
        self._write('t = tmin - (emin-ein)/cv')
        self._write('return t')
        self._outdent()
        
        self._write('if (ein > emax):')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('cv = cvbs(tmax, y)')
        self._write('t = tmax - (emax-ein)/cv')
        self._write('return t')
        self._outdent()

        self._write('t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin)')
        self._write('for i in range(maxiter):')
        self._indent()
        self._write('e1 = ubms(t1,y)')
        self._write('cv = cvbs(t1,y)')
        self._write('dt = (ein - e1) / cv')
        self._write('if (dt > 100):  dt = 100 ')
        self._write('elif (dt < -100):  dt = -100')
        self._write('elif (abs(dt) < tol): break')
        self._write('t1 += dt')
        self._outdent()
        
        self._write('t = t1')
        self._write('return t')
        
        self._outdent()
        return

 
    def _render_hytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        docstring = 'T = hytt(h, y): get temperature given enthalpy in mass units and mass fracs'
        self._write('def hytt(h, y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('maxiter = 50')
        self._write('tol  = 0.001')
        self._write('hin  = h')
        self._write('tmin = %g' % lowT + self.line('max lower bound for thermo def'))
        self._write('tmax = %g' % highT + self.line('min upper bound for thermo def'))
        self._write('hmin = hbms(tmin, y)')
        self._write('hmax = hbms(tmax, y)')
        self._write('if (hin < hmin):')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('cp = cpbs(tmin, y)')
        self._write('t = tmin - (hmin-hin)/cp')
        self._write('return t')
        self._outdent()
        
        self._write('if (hin > hmax):')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('cp = cpbs(tmax, y)')
        self._write('t = tmax - (hmax-hin)/cp')
        self._write('return t')
        self._outdent()

        self._write('t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin)')
        self._write('for i in range(maxiter):')
        self._indent()
        self._write('h1 = hbms(t1,y)')
        self._write('cp = cpbs(t1,y)')
        self._write('dt = (hin - h1) / cp')
        self._write('if (dt > 100):  dt = 100 ')
        self._write('elif (dt < -100):  dt = -100')
        self._write('elif (abs(dt) < tol): break')
        self._write('t1 += dt')
        self._outdent()
        
        self._write('t = t1')
        self._write('return t')
        
        self._outdent()
        return

 
    def _renderums(self, mechanism):
        self._write()
        self._write()
        docstring = 'ums = ums(T): Returns internal energy in mass units (Eq 30.)'
        self._write(self.line(docstring))
        self._write('def ums(T):')
        self._indent()
        self._write('"""%s"""' % docstring)
        
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('ums = speciesInternalEnergy(T)')
        
        # convert e/RT to e with mass units
        for species in self.species:
            self._write('ums[%d] *= RT/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('return ums')
        self._outdent()
        return

 
    def _renderhms(self, mechanism):
        self._write()
        self._write()
        docstring = 'hms = hms(T): Returns enthalpy in mass units (Eq 27.)'
        self._write(self.line('docstring'))
        self._write('def hmx(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        # call routine
        self._write('hms = speciesEnthalpy(T)')
        
        # convert h/RT to h with mass units
        for species in self.species:
            self._write('hms[%d] *= RT/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('return hms')
        self._outdent()

        return

    def _renderams(self, mechanism):
        self._write()
        self._write()
        docstring = 'ams = ams(T): Returns helmholtz in mass units (Eq 32.)'
        self._write(self.line(docstring))
        self._write('def ams(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('ams = helmholtz(T)')
        
        # convert A/RT to A with mass units
        for species in self.species:
            self._write('ams[%d] *= RT/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('return ams')
        self._outdent()
        
        return

    def _rendergms(self, mechanism):
        self._write()
        self._write()
        docstring = 'gms = gms(T): Returns gibbs in mass units (Eq 31.)'
        self._write(self.line(docstring))
        self._write('def gms(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('gms = gibbs(T)')

        # convert g/RT to g with mass units
        for species in self.species:
            self._write('gms[%d] *= RT/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write('return gms')
        self._outdent()

        return


    def _rendercvms(self, mechanism):
        self._write()
        self._write()
        docstring = 'cvms = cvms(T): Returns the specific heats at constant volume in mass units (Eq 29.'
        self._write(self.line(docstring))
        self._write('def cvms(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cvms = cv_R(T)')

        # convert cv/R to cv with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cvms[%d] *= %15.8e ' % (
                species.id, ROW) + self.line('%s' % species.symbol))
                    
        self._write('return cvms')
        self._outdent()

        return

    def _rendercpms(self, mechanism):
        self._write()
        self._write()
        docstring = 'cpms = cpms(T): Returns the specific heats at constant pressure in mass units (Eq. 26)'
        self._write(self.line(docstring))
        self._write('def cpms(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cpms = cp_R(T)')

        # convert cp/R to cp with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cpms[%d] *= %15.8e ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

        self._write('return cpms')
        self._outdent()

        return

    def _rendersms(self, mechanism):
        self._write()
        self._write()
        docstring = 'sms = sms(T): Returns the entropies in mass units (Eq 28.)'
        self._write(self.line(docstring))
        self._write('def sms(T):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('sms = speciesEntropy(T)')
        
        # convert s/R to s with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('sms[%d] *= %15.8e ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

       
        self._write('return sms')
        self._outdent()

        return
    
    
    def _rendercpbl(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'cpbl = cpbl(T,x): Returns the mean specific heat Cp (Eq. 33)'
        self._write(self.line(docstring))
        self._write('def cpbl(T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cpor = cp_R(T)')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*cpor[id]')
        self._outdent()

        self._write()
        self._write('cpbl = result * %g' % (R*kelvin*mole/erg) )
        
        self._write('return cpbl')
        self._outdent()

        return

 
    def _rendercpbs(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'cpbs = cpbs(T,y): Returns the mean specific heat Cp (Eq. 34)'
        self._write(self.line(docstring))
        self._write('def cpbs(T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cpor = cp_R(T)')
        
        # do dot product
        self._write(self.line('multiply by y/molecularweight'))
        self._write('result = 0')
        for species in self.species:
            self._write('result += cpor[%d]*y[%d]/%g ' % (
                species.id, species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('cpbs = result * %g' % (R*kelvin*mole/erg) )
        
        self._write('return cpbs')
        self._outdent()

        return

   
    def _rendercvbl(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'cvbl = cvbl(T,x): Returns the mean specific heat Cv (Eq. 35)'
        self._write(self.line(docstring))
        self._write('def cvbl(T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cvor = cv_R(T)')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*cvor[id]')
        self._outdent()

        self._write()
        self._write('cvbl = result * %g' % (R*kelvin*mole/erg) )
        
        self._write('return cvbl')
        self._outdent()

        return


    def _rendercvbs(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'cvbs = cvbs(T,y): Returns the mean specific heat Cv (Eq. 36)'
        self._write(self.line(docstring))
        self._write('def cvbs(T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # call routine
        self._write('cvor = cv_R(T)')
        
        # do dot product
        self._write(self.line('multiply by y/molecularweight'))
        self._write('result = 0')
        for species in self.species:
            self._write('result += cvor[%d]*y[%d]/%g  ' % (
                species.id, species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('cvbs = result * %g' % (R*kelvin*mole/erg) )
        
        self._write('return cvbs')
        self._outdent()

        return

      
    def _renderhbml(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'hbml = hbml(T,x): Returns the mean enthalpy of the mixture in molar units'
        self._write(self.line(docstring))
        self._write('def hbml(T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('hml = speciesEnthalpy(T)')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*hml[id]')
        self._outdent()

        self._write()
        self._write('hbml = result * RT')
        
        self._write('return hbml')
        self._outdent()

        return
 
 
    def _renderhbms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'hbms = hbms(T,y): Returns mean enthalpy of mixture in mass units'
        self._write(self.line(docstring))
        self._write('def hbms(T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('hml = speciesEnthalpy(T)')

        # convert e/RT to e with mass units
        self._write(self.line('perform dot product + scaling by wt'))
        self._write('result = 0')
        for species in self.species:
            self._write('result += y[%d]*hml[%d]/%f ' % (
                species.id, species.id, species.weight)
                        + self.line('%s' % species.symbol))

        self._write()
        # finally, multiply by RT
        self._write('hbms = result * RT')
        
        self._write('return hbms')
        self._outdent()
        
        return

    
    def _renderubml(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'ubml = ubml(T,x): get mean internal energy in molar units'
        self._write(self.line(docstring))
        self._write('def ubml(T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('uml = speciesInternalEnergy(T)')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*uml[id]')
        self._outdent()

        self._write()
        self._write('ubml = result * RT')
        
        self._write('return ubml')
        self._outdent()

        return

 
    def _renderubms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'ubms = ubms(T,y): get mean internal energy in mass units'
        self._write(self.line(docstring))
        self._write('def ubms(T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write('ums = speciesInternalEnergy(T)')

        # convert e/RT to e with mass units
        self._write(self.line('perform dot product + scaling by wt'))
        self._write('result = 0')
        for species in self.species:
            self._write('result += y[%d]*ums[%d]/%f ' % (
                species.id, species.id, species.weight)
                        + self.line('%s' % species.symbol))

        
        self._write()
        # finally, multiply by RT
        self._write('ubms = result * RT')
        
        self._write('return ubms')
        self._outdent()
        
        return

 
    def _rendersbml(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'sbml = sbml(P,T,x): get mixture entropy in molar units'
        self._write(self.line(docstring))
        self._write('def sbml(P,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( P / 1013250.0 )')
        
        # call routine
        self._write('sor = speciesEntropy(T)')
        
        # Equation 42
        self._write()
        self._write(self.line('Compute Eq 42'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*(sor[id]-math.log((x[id]+%g))-logPratio)' %
                    smallnum )
        self._outdent()

        self._write()
        
        self._write('sbml = result * %g' % (R*kelvin*mole/erg) )
        
        self._write('return sbml')
        self._outdent()

        return

    def _rendersbms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'sbms = sbms(P,T,y): get mixture entropy in mass units'
        self._write(self.line(docstring))
        self._write('def sbms(P,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( P / 1013250.0 ) ')
        self._write('x = ytx(y)')
            
        # call routine
        self._write('sor = speciesEntropy(T)')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 42 and 43'))
        self._write('result = 0')
        for species in self.species:
            self._write('result += x[%d]*(sor[%d]-math.log((x[%d]+%g))-logPratio) ' %
                        (species.id, species.id, species.id, smallnum) )
            
        self._write('YOW = 0.0')
        for species in self.species:
            self._write('YOW += y[%d]/%f  ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(self.line('Scale by R/W'))
        self._write('sbms = result * %g * YOW' % (R*kelvin*mole/erg) )
        
        self._write('return sbms')
        self._outdent()

        return


    def _rendergbml(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'gbml = gbml(P,T,x): Returns mean gibbs free energy in molar units'
        self._write(self.line(docstring))
        self._write('def gbml(P,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( P / 1013250.0 ) ')
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write(self.line('Compute g/RT'))
        self._write('gort = gibbs(T)')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('result = 0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*(gort[id]+math.log((x[id]+%g))+logPratio)' %
                    smallnum )
        self._outdent()

        self._write()
        self._write('gbml = result * RT')
        
        self._write('return gbml')
        self._outdent()

        return


    def _rendergbms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'gbms(P,T,y): Returns mixture gibbs free energy in mass units'
        self._write(self.line(docstring))
        self._write('def gbms(P,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( *P / 1013250.0 ) ')
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        self._write('YOW = 0.0')
        for species in self.species:
            self._write('YOW += y[%d]/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write('x = ytx(y)')
            
        # call routine
        self._write('gort = gibbs(T)')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        self._write('result = 0.0')
        for species in self.species:
            self._write('result += x[%d]*(gort[%d]+math.log((x[%d]+%g))+logPratio)' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('gbms = result * RT * YOW')
        
        self._write('return gbms')
        self._outdent()

        return
    

    def _renderabml(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'abml = abml(P,T,x): Returns mean helmholtz free energy in molar units'
        self._write(self.line(docstring))
        self._write('def abml(P,T,x):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( *P / 1013250.0 ) ')
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # call routine
        self._write(self.line('Compute a/RT'))
        self._write('aort = helmholtz(T)')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('result = 0.0')
        self._write('for id in range(%d):' % nSpecies)
        self._indent()
        self._write('result += x[id]*(aort[id]+math.log((x[id]+%g))+logPratio)' %
                    smallnum )
        self._outdent()

        self._write()
        self._write('abml = result * RT')
        
        self._write('return abml')
        self._outdent()

        return
    

    def _renderabms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        docstring = 'abms = abms(P,T,y): Returns mixture helmholtz free energy in mass units'
        self._write(self.line(docstring))
        self._write('def abms(P,T,y):')
        self._indent()
        self._write('"""%s"""' % docstring)

        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write('import math')
        self._write('logPratio = math.log ( P / 1013250.0 ) ')
        self._write('RT = %g*T ' % (R*kelvin*mole/erg) + self.line('R*T'))
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        self._write('YOW = 0.0')
        for species in self.species:
            self._write('YOW += y[%d]/%f ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
 
        # now do ytx
        self._write('x = ytx(y)')
            
        # call routine
        self._write('aort = helmholtz(T)')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        self._write('result = 0.0')
        for species in self.species:
            self._write('result += x[%d]*(aort[%d]+math.log((x[%d]+%g))+logPratio)' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('abms = result * RT * YOW')
        
        self._write('return abms')
        self._outdent()

        return
     
#  End of file 
