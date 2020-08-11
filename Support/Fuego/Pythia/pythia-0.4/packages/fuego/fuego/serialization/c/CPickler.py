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
from collections import defaultdict

from weaver.mills.CMill import CMill

from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ, erg
from pyre.handbook.constants.fundamental import boltzmann
from pyre.handbook.constants.fundamental import avogadro
from pyre.handbook.constants.fundamental import gas_constant as R

import sys
import numpy as np

smallnum = 1e-100
R = 8.31446261815324e7 * erg/mole/kelvin
Rc = 1.98721558317399615845 * cal/mole/kelvin
Patm = 1013250.0
sym  = ""
fsym = "_"

class speciesDb:
    def __init__(self, id, name, mwt):
        self.id = id
        self.symbol = name
        self.weight = mwt
        return


class CPickler(CMill):


    def __init__(self):
        CMill.__init__(self)
        self.species = []
        self.nSpecies = 0
        self.reactionIndex = []
        self.lowT = 100.0
        self.highT = 10000.0
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

        self.nSpecies = nSpecies
        return


    ##########################
    #This is the main routine
    #called in weaver/weaver/mills/Mill.py
    ##########################
    def _renderDocument_CHOP(self, mechanism, options=None):

        reorder_reactions=False

        if(reorder_reactions):

            plot_react_matrix = True
            use_tsp = True #traveling salesman reordering

            if(plot_react_matrix):
                import matplotlib.pyplot as mplt
                (fig,ax)=mplt.subplots(1,4,figsize=(20,5))
                rmat=mechanism._get_reaction_matrix()
                ax[0].matshow(rmat)

            #sort reactions by type
            self.reactionIndex = mechanism._sort_reactions()
            if(plot_react_matrix):
                rmat=mechanism._get_reaction_matrix()
                ax[1].matshow(rmat)

            #reorder reactions
            if(use_tsp):
                mechanism._sort_reactions_within_type_tsp(self.reactionIndex)
            else:
                mechanism._sort_reactions_within_type_random(self.reactionIndex)
            if(plot_react_matrix):
                rmat=mechanism._get_reaction_matrix()
                ax[2].matshow(rmat)

            #reorder species
            if(use_tsp):
                mechanism._sort_species_ids_tsp()
            else:
                mechanism._sort_species_ids_random()
            if(plot_react_matrix):
                rmat=mechanism._get_reaction_matrix()
                ax[3].matshow(rmat)
                mplt.savefig("rmat_all.pdf")

            #set species after reordering    
            self._setSpecies(mechanism)

        else:
            self._setSpecies(mechanism)
            self.reactionIndex = mechanism._sort_reactions()

        #chemistry_file.H
        self._includes(True)
        self._header_chop(mechanism)
        self._namespace(mechanism)
        #chemistry_file.H

        self._includes_chop()
        self._statics_chop(mechanism)
        self._ckinit_chop(mechanism)

        # chemkin wrappers
        self._ckindx(mechanism)
        self._ckxnum(mechanism)
        self._cksnum(mechanism)
        self._cksyme_str(mechanism)
        self._cksyme(mechanism)
        self._cksyms_str(mechanism)
        self._cksyms(mechanism)
        self._ckrp(mechanism)
        
        self._ckpx(mechanism)
        self._ckpy(mechanism)
        self._vckpy(mechanism)
        self._ckpc(mechanism)
        self._ckrhox(mechanism)
        self._ckrhoy(mechanism)
        self._ckrhoc(mechanism)
        self._ckwt(mechanism)
        self._ckawt(mechanism)
        self._ckmmwy(mechanism)
        self._ckmmwx(mechanism)
        self._ckmmwc(mechanism)
        self._ckytx(mechanism)
        self._vckytx(mechanism)
        self._ckytcp(mechanism)
        self._ckytcr(mechanism)
        self._ckxty(mechanism)
        self._ckxtcp(mechanism)
        self._ckxtcr(mechanism)
        self._ckctx(mechanism)
        self._ckcty(mechanism)
        
        self._ckcpor(mechanism)
        self._ckhort(mechanism)
        self._cksor(mechanism)
        
        self._ckcvml(mechanism)
        self._ckcpml(mechanism)
        self._ckuml(mechanism)
        self._ckhml(mechanism)
        self._ckgml(mechanism)
        self._ckaml(mechanism)
        self._cksml(mechanism)
        
        self._ckcvms(mechanism)
        self._ckcpms(mechanism)
        self._ckums(mechanism)
        self._ckhms(mechanism)
        self._vckhms(mechanism)
        self._ckgms(mechanism)
        self._ckams(mechanism)
        self._cksms(mechanism)

        self._ckcpbl(mechanism)
        self._ckcpbs(mechanism)
        self._ckcvbl(mechanism)
        self._ckcvbs(mechanism)
        
        self._ckhbml(mechanism)
        self._ckhbms(mechanism)
        self._ckubml(mechanism)
        self._ckubms(mechanism)
        self._cksbml(mechanism)
        self._cksbms(mechanism)
        self._ckgbml(mechanism)
        self._ckgbms(mechanism)
        self._ckabml(mechanism)
        self._ckabms(mechanism)

        self._ckwc(mechanism)
        self._ckwyp(mechanism)
        self._ckwxp(mechanism)
        self._ckwyr(mechanism)
        self._vckwyr(mechanism)
        self._ckwxr(mechanism)
        
        self._ckqc(mechanism)
        self._ckkfkr(mechanism)
        self._ckqyp(mechanism)
        self._ckqxp(mechanism)
        self._ckqyr(mechanism)
        self._ckqxr(mechanism)

        self._cknu(mechanism)
        self._ckinu(mechanism)
        self._ckncf(mechanism)
        
        self._ckabe(mechanism)
        
        self._ckeqc(mechanism)
        self._ckeqyp(mechanism)
        self._ckeqxp(mechanism)
        self._ckeqyr(mechanism)
        self._ckeqxr(mechanism)
        
        # Fuego Functions
        # GPU version
        self._productionRate_GPU(mechanism)
        # ORI CPU version
        self._productionRate(mechanism)
        # ORI CPU vectorized version
        self._vproductionRate(mechanism)
        self._DproductionRatePrecond(mechanism)
        #self._DproductionRateSPSPrecond(mechanism)
        self._DproductionRate(mechanism)
        self._sparsity(mechanism)
        # GPU version
        self._ajac_GPU(mechanism)
        # ORI CPU version
        self._ajac(mechanism)
        self._ajacPrecond(mechanism)
        self._dthermodT(mechanism)
        self._progressRate(mechanism)
        self._progressRateFR(mechanism)
        self._equilibriumConstants(mechanism)
        self._thermo_GPU(mechanism)
        self._atomicWeight(mechanism)
        self._T_given_ey(mechanism)
        self._T_given_hy(mechanism)
        self._getCriticalParameters(mechanism)
        #AF: add transport data
        self._trans_chop(mechanism)
        #AF: dummy gjs routines
        self._emptygjs(mechanism)

        ### MECH HEADER -- second file starts here
        self._print_mech_header(mechanism)
        ### MECH HEADER

        return



    #Pieces for the file chemistry_file.H#
    def _includes(self, header):
        self._rep += [
            '#include <math.h>',
            '#include <stdio.h>',
            '#include <string.h>'
        ]
        if header:
            self._rep += [
                '#include <stdlib.h>',
                '#include <vector>',
                '#include <AMReX_Gpu.H>'
            ]
        else:
            self._rep += [
                '#include <stdlib.h>'
            ]

        return


    def _header_chop(self, mechanism):
        self._rep += [
            '',
            'extern "C"',
            '{',
            'AMREX_GPU_HOST_DEVICE void get_imw(double imw_new[]);',
            'AMREX_GPU_HOST_DEVICE void get_mw(double mw_new[]);',
            'void egtransetEPS(double *  EPS);',
            'void egtransetSIG(double* SIG);',
            'void atomicWeight(double *  awt);',
            'void molecularWeight(double *  wt);',
            'AMREX_GPU_HOST_DEVICE void gibbs(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void helmholtz(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void speciesEnthalpy(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void speciesEntropy(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void cp_R(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void cv_R(double *  species, double *  tc);',
            'void equilibriumConstants(double *  kc, double *  g_RT, double T);',
            'AMREX_GPU_HOST_DEVICE void productionRate(double *  wdot, double *  sc, double T);',
            'AMREX_GPU_HOST_DEVICE void comp_qfqr(double *  q_f, double *  q_r, double *  sc, double *  tc, double invT);',
            '#ifndef AMREX_USE_CUDA',
            'void comp_k_f(double *  tc, double invT, double *  k_f);',
            'void comp_Kc(double *  tc, double invT, double *  Kc);',
            '#endif',
            'AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  speciesConc, double T);',
            'AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  speciesConc, double T);',
            ##'#ifndef AMREX_USE_CUDA',
            'AMREX_GPU_HOST_DEVICE void CKINIT'+sym+'();',
            'AMREX_GPU_HOST_DEVICE void CKFINALIZE'+sym+'();',
            '#ifndef AMREX_USE_CUDA',
            'void GET_REACTION_MAP(int *  rmap);',
            'void SetAllDefaults();',
            '#endif',
            'void CKINDX'+sym+'(int * mm, int * kk, int * ii, int * nfit );',
            'void CKXNUM'+sym+'(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline);',
            'void CKSNUM'+sym+'(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray);',
            'void CKSYME_STR(amrex::Vector<std::string>& ename);',
            'void CKSYME(int * kname, int * lenkname);',
            'void CKSYMS_STR(amrex::Vector<std::string>& kname);',
            'void CKSYMS(int * kname, int * lenkname);',
            'void CKRP'+sym+'(double *  ru, double *  ruc, double *  pa);',
            'void CKPX'+sym+'(double *  rho, double *  T, double *  x, double *  P);',
            'AMREX_GPU_HOST_DEVICE void CKPY'+sym+'(double *  rho, double *  T, double *  y, double *  P);',
            'void CKPC'+sym+'(double *  rho, double *  T, double *  c, double *  P);',
            'void CKRHOX'+sym+'(double *  P, double *  T, double *  x, double *  rho);',
            'AMREX_GPU_HOST_DEVICE void CKRHOY'+sym+'(double *  P, double *  T, double *  y, double *  rho);',
            'void CKRHOC'+sym+'(double *  P, double *  T, double *  c, double *  rho);',
            'void CKWT'+sym+'(double *  wt);',
            'void CKAWT'+sym+'(double *  awt);',
            'AMREX_GPU_HOST_DEVICE void CKMMWY'+sym+'(double *  y, double *  wtm);',
            'void CKMMWX'+sym+'(double *  x, double *  wtm);',
            'void CKMMWC'+sym+'(double *  c, double *  wtm);',
            'AMREX_GPU_HOST_DEVICE void CKYTX'+sym+'(double *  y, double *  x);',
            'void CKYTCP'+sym+'(double *  P, double *  T, double *  y, double *  c);',
            'AMREX_GPU_HOST_DEVICE void CKYTCR'+sym+'(double *  rho, double *  T, double *  y, double *  c);',
            'AMREX_GPU_HOST_DEVICE void CKXTY'+sym+'(double *  x, double *  y);',
            'void CKXTCP'+sym+'(double *  P, double *  T, double *  x, double *  c);',
            'void CKXTCR'+sym+'(double *  rho, double *  T, double *  x, double *  c);',
            'void CKCTX'+sym+'(double *  c, double *  x);',
            'void CKCTY'+sym+'(double *  c, double *  y);',
            'void CKCPOR'+sym+'(double *  T, double *  cpor);',
            'void CKHORT'+sym+'(double *  T, double *  hort);',
            'void CKSOR'+sym+'(double *  T, double *  sor);',
            
            'void CKCVML'+sym+'(double *  T, double *  cvml);',
            'void CKCPML'+sym+'(double *  T, double *  cvml);',
            'void CKUML'+sym+'(double *  T, double *  uml);',
            'void CKHML'+sym+'(double *  T, double *  uml);',
            'void CKGML'+sym+'(double *  T, double *  gml);',
            'void CKAML'+sym+'(double *  T, double *  aml);',
            'void CKSML'+sym+'(double *  T, double *  sml);',
            
            'AMREX_GPU_HOST_DEVICE void CKCVMS'+sym+'(double *  T, double *  cvms);',
            'AMREX_GPU_HOST_DEVICE void CKCPMS'+sym+'(double *  T, double *  cvms);',
            'AMREX_GPU_HOST_DEVICE void CKUMS'+sym+'(double *  T, double *  ums);',
            'AMREX_GPU_HOST_DEVICE void CKHMS'+sym+'(double *  T, double *  ums);',
            'void CKGMS'+sym+'(double *  T, double *  gms);',
            'void CKAMS'+sym+'(double *  T, double *  ams);',
            'void CKSMS'+sym+'(double *  T, double *  sms);',
            
            'void CKCPBL'+sym+'(double *  T, double *  x, double *  cpbl);',
            'AMREX_GPU_HOST_DEVICE void CKCPBS'+sym+'(double *  T, double *  y, double *  cpbs);',
            'void CKCVBL'+sym+'(double *  T, double *  x, double *  cpbl);',
            'AMREX_GPU_HOST_DEVICE void CKCVBS'+sym+'(double *  T, double *  y, double *  cpbs);',
            
            'void CKHBML'+sym+'(double *  T, double *  x, double *  hbml);',
            'AMREX_GPU_HOST_DEVICE void CKHBMS'+sym+'(double *  T, double *  y, double *  hbms);',
            'void CKUBML'+sym+'(double *  T, double *  x, double *  ubml);',
            'AMREX_GPU_HOST_DEVICE void CKUBMS'+sym+'(double *  T, double *  y, double *  ubms);',
            'void CKSBML'+sym+'(double *  P, double *  T, double *  x, double *  sbml);',
            'void CKSBMS'+sym+'(double *  P, double *  T, double *  y, double *  sbms);',
            'void CKGBML'+sym+'(double *  P, double *  T, double *  x, double *  gbml);',
            'void CKGBMS'+sym+'(double *  P, double *  T, double *  y, double *  gbms);',
            'void CKABML'+sym+'(double *  P, double *  T, double *  x, double *  abml);',
            'void CKABMS'+sym+'(double *  P, double *  T, double *  y, double *  abms);',

            
            'AMREX_GPU_HOST_DEVICE void CKWC'+sym+'(double *  T, double *  C, double *  wdot);',
            'void CKWYP'+sym+'(double *  P, double *  T, double *  y, double *  wdot);',
            'void CKWXP'+sym+'(double *  P, double *  T, double *  x, double *  wdot);',
            'AMREX_GPU_HOST_DEVICE void CKWYR'+sym+'(double *  rho, double *  T, double *  y, double *  wdot);',
            'void CKWXR'+sym+'(double *  rho, double *  T, double *  x, double *  wdot);',

            
            'void CKQC'+sym+'(double *  T, double *  C, double *  qdot);',
            'void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r);',
            'void CKQYP'+sym+'(double *  P, double *  T, double *  y, double *  qdot);',
            'void CKQXP'+sym+'(double *  P, double *  T, double *  x, double *  qdot);',
            'void CKQYR'+sym+'(double *  rho, double *  T, double *  y, double *  qdot);',
            'void CKQXR'+sym+'(double *  rho, double *  T, double *  x, double *  qdot);',
            
            'void CKNU'+sym+'(int * kdim, int * nuki);',
            '#ifndef AMREX_USE_CUDA',
            'void CKINU'+sym+'(int * i, int * nspec, int * ki, int * nu);',
            '#endif',
            'void CKNCF'+sym+'(int * ncf);',
            
            'void CKABE'+sym+'(double *  a, double *  b, double *  e );',
            'void CKEQC'+sym+'(double *  T, double *  C , double *  eqcon );',
            'void CKEQYP'+sym+'(double *  P, double *  T, double *  y, double *  eqcon);',
            'void CKEQXP'+sym+'(double *  P, double *  T, double *  x, double *  eqcon);',
            'void CKEQYR'+sym+'(double *  rho, double *  T, double *  y, double *  eqcon);',
            'void CKEQXR'+sym+'(double *  rho, double *  T, double *  x, double *  eqcon);',
            'AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  T, int * consP);',
            'AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(double *  J, double *  sc, double *  Tp, int * HP);',
            #'AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(double *  Jsps, int * indx, int * len, double * sc, double * Tp, int * HP, double * gamma);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);',
            'AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);',
            'AMREX_GPU_HOST_DEVICE void aJacobian(double *  J, double *  sc, double T, int consP);',
            'AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP);',
            'AMREX_GPU_HOST_DEVICE void dcvpRdT(double *  species, double *  tc);',
            'AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int *ierr);',
            'AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int *ierr);',
            'void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i);',
            self.line('vector version'),
            'void VCKYTX'+sym+'(int *  np, double *  y, double *  x);',
            'void VCKHMS'+sym+'(int *  np, double *  T, double *  ums);',
            'void VCKWYR'+sym+'(int *  np, double *  rho, double *  T,',
            '            double *  y,',
            '            double *  wdot);',
            '#ifndef AMREX_USE_CUDA',
            'void vproductionRate(int npt, double *  wdot, double *  c, double *  T);',
            'void VCKPY'+sym+'(int *  np, double *  rho, double *  T, double *  y, double *  P);',
            'void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT);',
            'void vcomp_gibbs(int npt, double *  g_RT, double *  tc);',
            'void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT);',
            
            ]
        nReactions = len(mechanism.reaction())
        if nReactions <= 50:
            self._rep += [
                'void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,',
                '                double *  k_f_s, double *  Kc_s,',
                '                double *  tc, double *  invT, double *  T);',
                ]
        else:
            for i in range(0,nReactions,50):
                self._rep += [
                    'void vcomp_wdot_%d_%d(int npt, double *  wdot, double *  mixture, double *  sc,' % (i+1,min(i+50,nReactions)),
                    '                double *  k_f_s, double *  Kc_s,',
                    '                double *  tc, double *  invT, double *  T);',
                    ]                

        self._rep += [
                '#endif',
                self.line('Transport function declarations'),
                'void egtransetLENIMC(int* LENIMC);',
                'void egtransetLENRMC(int* LENRMC);',
                'void egtransetNO(int* NO);',
                'void egtransetKK(int* KK);',
                'void egtransetNLITE(int* NLITE);',
                'void egtransetPATM(double* PATM);',
                'void egtransetWT(double* WT);',
                'void egtransetEPS(double* EPS);',
                'void egtransetSIG(double* SIG);',
                'void egtransetDIP(double* DIP);',
                'void egtransetPOL(double* POL);',
                'void egtransetZROT(double* ZROT);',
                'void egtransetNLIN(int* NLIN);',
                'void egtransetCOFETA(double* COFETA);',
                'void egtransetCOFLAM(double* COFLAM);',
                'void egtransetCOFD(double* COFD);',
                'void egtransetKTDIF(int* KTDIF);',
            ]

        self._rep += [
                self.line('gauss-jordan solver external routine'),
                'AMREX_GPU_HOST_DEVICE void sgjsolve(double* A, double* x, double* b);',
                'AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(double* A, double* x, double* b);',
                '}',
                ]

        return


    def _namespace(self,mechanism):
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write('namespace thermo')
        self._write('{')

        self._indent()

        nReactions = len(mechanism.reaction())
        self._write()
        self._write('extern double fwd_A[%d], fwd_beta[%d], fwd_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double low_A[%d], low_beta[%d], low_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double rev_A[%d], rev_beta[%d], rev_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double troe_a[%d],troe_Ts[%d], troe_Tss[%d], troe_Tsss[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('extern double sri_a[%d], sri_b[%d], sri_c[%d], sri_d[%d], sri_e[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('extern double activation_units[%d], prefactor_units[%d], phase_units[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('extern int is_PD[%d], troe_len[%d], sri_len[%d], nTB[%d], *TBid[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('extern double *TB[%d];' 
                    % (nReactions))

        self._write('extern std::vector<std::vector<double>> kiv; ')
        self._write('extern std::vector<std::vector<double>> nuv; ')

        self._write()
        self._write('extern double fwd_A_DEF[%d], fwd_beta_DEF[%d], fwd_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double low_A_DEF[%d], low_beta_DEF[%d], low_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double rev_A_DEF[%d], rev_beta_DEF[%d], rev_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('extern double troe_a_DEF[%d],troe_Ts_DEF[%d], troe_Tss_DEF[%d], troe_Tsss_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('extern double sri_a_DEF[%d], sri_b_DEF[%d], sri_c_DEF[%d], sri_d_DEF[%d], sri_e_DEF[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('extern double activation_units_DEF[%d], prefactor_units_DEF[%d], phase_units_DEF[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('extern int is_PD_DEF[%d], troe_len_DEF[%d], sri_len_DEF[%d], nTB_DEF[%d], *TBid_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('extern double *TB_DEF[%d];' 
                    % (nReactions))
        self._write('extern std::vector<int> rxn_map;')

        self._outdent()

        self._write('}')
        self._write('#endif')

        return
    #Pieces for the file chemistry_file.H#



    #Pieces for mechanism.cpp
    def _includes_chop(self):
        self._rep += [
            '#include "chemistry_file.H"'
            ]
        return


    def _statics_chop(self,mechanism):
        nReactions = len(mechanism.reaction())
        nSpecies = len(mechanism.species())

        ispecial   = self.reactionIndex[5:7]

        nspecial   = ispecial[1]   - ispecial[0]

        self._write()

        self._write('#ifndef AMREX_USE_CUDA')
        self._write('namespace thermo')
        self._write('{')
        self._indent()
        self._write('double fwd_A[%d], fwd_beta[%d], fwd_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double low_A[%d], low_beta[%d], low_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double rev_A[%d], rev_beta[%d], rev_Ea[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double troe_a[%d],troe_Ts[%d], troe_Tss[%d], troe_Tsss[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('double sri_a[%d], sri_b[%d], sri_c[%d], sri_d[%d], sri_e[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('double activation_units[%d], prefactor_units[%d], phase_units[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('int is_PD[%d], troe_len[%d], sri_len[%d], nTB[%d], *TBid[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('double *TB[%d];' 
                    % (nReactions))

        if nspecial > 0:  
                self._write('double prefactor_units_rev[%d], activation_units_rev[%d];' 
                            % (nReactions,nReactions))

        self._write('std::vector<std::vector<double>> kiv(%d); ' % (nReactions))
        self._write('std::vector<std::vector<double>> nuv(%d); ' % (nReactions))

        self._write()
        self._write('double fwd_A_DEF[%d], fwd_beta_DEF[%d], fwd_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double low_A_DEF[%d], low_beta_DEF[%d], low_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double rev_A_DEF[%d], rev_beta_DEF[%d], rev_Ea_DEF[%d];' 
                    % (nReactions,nReactions,nReactions))
        self._write('double troe_a_DEF[%d],troe_Ts_DEF[%d], troe_Tss_DEF[%d], troe_Tsss_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions))
        self._write('double sri_a_DEF[%d], sri_b_DEF[%d], sri_c_DEF[%d], sri_d_DEF[%d], sri_e_DEF[%d];'
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('double activation_units_DEF[%d], prefactor_units_DEF[%d], phase_units_DEF[%d];'
                    % (nReactions,nReactions,nReactions))
        self._write('int is_PD_DEF[%d], troe_len_DEF[%d], sri_len_DEF[%d], nTB_DEF[%d], *TBid_DEF[%d];' 
                    % (nReactions,nReactions,nReactions,nReactions,nReactions))
        self._write('double *TB_DEF[%d];' 
                    % (nReactions))
        self._write('std::vector<int> rxn_map;')

        self._outdent()

        self._write('};')

        self._write()

        self._write('using namespace thermo;')
        self._write('#endif')
        self._write()


        self._write(self.line(' Inverse molecular weights'))
        self._write(self.line(' TODO: check necessity on CPU'))
        self._write('static AMREX_GPU_DEVICE_MANAGED double imw[%d] = {' %nSpecies )
        self._indent()
        for i in range(0,self.nSpecies):
            species = self.species[i]
            text = '1.0 / %f' % (species.weight)
            if (i<self.nSpecies-1):
               text += ',  '
            else:
               text += '};  '
            self._write(text + self.line('%s' % species.symbol))
        self._outdent()
        self._write()

        self._write(self.line(' Inverse molecular weights'))
        self._write(self.line(' TODO: check necessity because redundant with molecularWeight'))
        self._write('static AMREX_GPU_DEVICE_MANAGED double molecular_weights[%d] = {' %nSpecies )
        self._indent()
        for i in range(0,self.nSpecies):
            species = self.species[i]
            text = '%f' % (species.weight)
            if (i<self.nSpecies-1):
               text += ',  '
            else:
               text += '};  '
            self._write(text + self.line('%s' % species.symbol))
        self._outdent()
        self._write()

        self._write('AMREX_GPU_HOST_DEVICE')
        self._write('void get_imw(double imw_new[]){')
        ##self._write('#pragma unroll')
        self._indent()
        self._write('for(int i = 0; i<%d; ++i) imw_new[i] = imw[i];' %nSpecies )
        self._outdent()
        self._write('}')
        self._write()

        self._write(self.line(' TODO: check necessity because redundant with CKWT'))
        self._write('AMREX_GPU_HOST_DEVICE')
        self._write('void get_mw(double mw_new[]){')
        ##self._write('#pragma unroll')
        self._indent()
        self._write('for(int i = 0; i<%d; ++i) mw_new[i] = molecular_weights[i];' %nSpecies )
        self._outdent()
        self._write('}')
        self._write()

        self._write()

        return


    def _ckinit_chop(self, mechanism):

        nElement = len(mechanism.element())
        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line(' Initializes parameter database'))
        self._write('void CKINIT'+sym+'()')
        self._write('{')
        self._write()

        self._indent()

        # build reverse reaction map
        rmap = {}
        for i, reaction in zip(range(nReactions), mechanism.reaction()):
            rmap[reaction.orig_id-1] = i
        
        self._write('rxn_map = {%s};' % (",".join(str(rmap[x]) for x in range(len(rmap)))))

        self._write()

        for j in range(nReactions):
            reaction = mechanism.reaction()[rmap[j]]
            id = reaction.id - 1

            ki = []
            nu = []
            for symbol, coefficient in reaction.reactants:
                ki.append(mechanism.species(symbol).id)
                nu.append(-coefficient)
            for symbol, coefficient in reaction.products:
                ki.append(mechanism.species(symbol).id)
                nu.append(coefficient)

            self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
            kistr = "{" + ','.join(str(x) for x in ki) + "}"
            nustr = "{" + ','.join(str(x) for x in nu) + "}"
            self._write("kiv[%d] = %s;" % (id,kistr))
            self._write("nuv[%d] = %s;" % (id,nustr))

            A, beta, E = reaction.arrhenius
            self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
            self._write("fwd_A[%d]     = %.17g;" % (id,A))
            self._write("fwd_beta[%d]  = %.17g;" % (id,beta))
            self._write("fwd_Ea[%d]    = %.17g;" % (id,E))

            thirdBody = reaction.thirdBody
            low = reaction.low

            if (reaction.rev):
                Ar, betar, Er = reaction.rev
                self._write("rev_A[%d]     = %.17g;" % (id,Ar))
                self._write("rev_beta[%d]  = %.17g;" % (id,betar))
                self._write("rev_Ea[%d]    = %.17g;" % (id,Er))
                dim_rev       = self._phaseSpaceUnits(reaction.products)
                if not thirdBody:
                    uc_rev = self._prefactorUnits(reaction.units["prefactor"], 1-dim_rev)
                elif not low:
                    uc_rev = self._prefactorUnits(reaction.units["prefactor"], -dim_rev)
                else:
                    uc_rev = self._prefactorUnits(reaction.units["prefactor"], 1-dim_rev)
                self._write("prefactor_units_rev[%d]  = %.17g;" % (id,uc_rev.value))
                aeuc_rev = self._activationEnergyUnits(reaction.units["activation"])
                self._write("activation_units_rev[%d] = %.17g;" % (id,aeuc_rev / Rc / kelvin))

            if (len(reaction.ford) > 0) :
                if (reaction.rev):
                    print '\n\n ***** WARNING: Reac is FORD and REV. Results might be wrong !\n'
                dim = self._phaseSpaceUnits(reaction.ford)
            else:
                dim = self._phaseSpaceUnits(reaction.reactants)

            if not thirdBody:
                uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB
            elif not low:
                uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
            else:
                uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
                low_A, low_beta, low_E = low
                self._write("low_A[%d]     = %.17g;" % (id,low_A))
                self._write("low_beta[%d]  = %.17g;" % (id,low_beta))
                self._write("low_Ea[%d]    = %.17g;" % (id,low_E))
                if reaction.troe:
                    troe = reaction.troe
                    ntroe = len(troe)
                    is_troe = True
                    self._write("troe_a[%d]    = %.17g;" % (id,troe[0]))
                    if ntroe>1:
                        self._write("troe_Tsss[%d] = %.17g;" % (id,troe[1]))
                    if ntroe>2:
                        self._write("troe_Ts[%d]   = %.17g;" % (id,troe[2]))
                    if ntroe>3:
                        self._write("troe_Tss[%d]  = %.17g;" % (id,troe[3]))
                    self._write("troe_len[%d]  = %d;" % (id,ntroe))
                if reaction.sri:
                    sri = reaction.sri
                    nsri = len(sri)
                    is_sri = True
                    self._write("sri_a[%d]     = %.17g;" % (id,sri[0]))
                    if nsri>1:
                        self._write("sri_b[%d]     = %.17g;" % (id,sri[1]))
                    if nsri>2:
                        self._write("sri_c[%d]     = %.17g;" % (id,sri[2]))
                    if nsri>3:
                        self._write("sri_d[%d]     = %.17g;" % (id,sri[3]))
                    if nsri>4:
                        self._write("sri_e[%d]     = %.17g;" % (id,sri[4]))
                    self._write("sri_len[%d]   = %d;" % (id,nsri))

            self._write("prefactor_units[%d]  = %.17g;" % (id,uc.value))
            aeuc = self._activationEnergyUnits(reaction.units["activation"])
            self._write("activation_units[%d] = %.17g;" % (id,aeuc / Rc / kelvin))
            self._write("phase_units[%d]      = pow(10,-%f);" % (id,dim*6))

            if low:
                self._write("is_PD[%d] = 1;" % (id) )
            else:
                self._write("is_PD[%d] = 0;" % (id) )


            if thirdBody:
                efficiencies = reaction.efficiencies
                if (len(efficiencies) > 0):
                    self._write("nTB[%d] = %d;" % (id, len(efficiencies)))
                    self._write("TB[%d] = (double *) malloc(%d * sizeof(double));" % (id, len(efficiencies)))
                    self._write("TBid[%d] = (int *) malloc(%d * sizeof(int));" % (id, len(efficiencies)))
                    for i, eff in enumerate(efficiencies):
                        symbol, efficiency = eff
                        self._write("TBid[%d][%d] = %.17g; TB[%d][%d] = %.17g; // %s"
                                    % (id, i, mechanism.species(symbol).id, id, i, efficiency, symbol ))
                else:
                    self._write("nTB[%d] = 0;" % (id))
            else:
                self._write("nTB[%d] = 0;" % (id))

            self._write()

        self._write("SetAllDefaults();")
        self._outdent()
        self._write("}")
        self._write()
            
        self._write('void GET_REACTION_MAP(int *rmap)')
        self._write('{')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()        
        self._write('rmap[i] = rxn_map[i] + 1;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write()

        self._write("#include <ReactionData.H>")
        self._write("double* GetParamPtr(int                reaction_id,")
        self._write("                    REACTION_PARAMETER param_id,")
        self._write("                    int                species_id,")
        self._write("                    int                get_default)")
        self._write("{")
        if nReactions == 0:
            self._write("  printf(\"No reactions in this model\");")
            self._write("  abort();")
            self._write("  return 0;")
        else:
            self._write("  double* ret = 0;")
            self._write("  if (reaction_id<0 || reaction_id>=%d) {" % (nReactions))
            self._write("    printf(\"Bad reaction id = %d\",reaction_id);")
            self._write("    abort();")
            self._write("  };")
            self._write("  int mrid = rxn_map[reaction_id];")
            self._write()
            self._write("  if (param_id == THIRD_BODY) {")
            self._write("    if (species_id<0 || species_id>=%d) {" % (self.nSpecies))
            self._write("      printf(\"GetParamPtr: Bad species id = %d\",species_id);")
            self._write("      abort();")
            self._write("    }")
            self._write("    if (get_default) {")
            self._write("      for (int i=0; i<nTB_DEF[mrid]; ++i) {")
            self._write("        if (species_id == TBid_DEF[mrid][i]) {")
            self._write("          ret = &(TB_DEF[mrid][i]);")
            self._write("        }")
            self._write("      }")
            self._write("    }")
            self._write("    else {")
            self._write("      for (int i=0; i<nTB[mrid]; ++i) {")
            self._write("        if (species_id == TBid[mrid][i]) {")
            self._write("          ret = &(TB[mrid][i]);")
            self._write("        }")
            self._write("      }")
            self._write("    }")
            self._write("    if (ret == 0) {")
            self._write("      printf(\"GetParamPtr: No TB for reaction id = %d\",reaction_id);")
            self._write("      abort();")
            self._write("    }")
            self._write("  }")
            self._write("  else {")
            self._write("    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}")
            self._write("      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}")
            self._write("      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}")
            self._write("      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}")
            self._write("      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}")
            self._write("      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}")
            self._write("      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}")
            self._write("      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}")
            self._write("      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}")
            self._write("      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}")
            self._write("      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}")
            self._write("      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}")
            self._write("      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}")
            self._write("      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}")
            self._write("      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}")
            self._write("      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}")
            self._write("      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}")
            self._write("      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}")
            self._write("    else {")
            self._write("      printf(\"GetParamPtr: Unknown parameter id\");")
            self._write("      abort();")
            self._write("    }")
            self._write("  }")
            self._write("  return ret;")
        self._write("}")
        self._write()

        self._write("void ResetAllParametersToDefault()")
        self._write("{")
        self._write("    for (int i=0; i<%d; i++) {" % (nReactions))
        self._write("        if (nTB[i] != 0) {")
        self._write("            nTB[i] = 0;")
        self._write("            free(TB[i]);")
        self._write("            free(TBid[i]);")
        self._write("        }")
        self._write("")
        self._write("        fwd_A[i]    = fwd_A_DEF[i];")
        self._write("        fwd_beta[i] = fwd_beta_DEF[i];")
        self._write("        fwd_Ea[i]   = fwd_Ea_DEF[i];")
        self._write("")
        self._write("        low_A[i]    = low_A_DEF[i];")
        self._write("        low_beta[i] = low_beta_DEF[i];")
        self._write("        low_Ea[i]   = low_Ea_DEF[i];")
        self._write("")
        self._write("        rev_A[i]    = rev_A_DEF[i];")
        self._write("        rev_beta[i] = rev_beta_DEF[i];")
        self._write("        rev_Ea[i]   = rev_Ea_DEF[i];")
        self._write("")
        self._write("        troe_a[i]    = troe_a_DEF[i];")
        self._write("        troe_Ts[i]   = troe_Ts_DEF[i];")
        self._write("        troe_Tss[i]  = troe_Tss_DEF[i];")
        self._write("        troe_Tsss[i] = troe_Tsss_DEF[i];")
        self._write("")
        self._write("        sri_a[i] = sri_a_DEF[i];")
        self._write("        sri_b[i] = sri_b_DEF[i];")
        self._write("        sri_c[i] = sri_c_DEF[i];")
        self._write("        sri_d[i] = sri_d_DEF[i];")
        self._write("        sri_e[i] = sri_e_DEF[i];")
        self._write("")
        self._write("        is_PD[i]    = is_PD_DEF[i];")
        self._write("        troe_len[i] = troe_len_DEF[i];")
        self._write("        sri_len[i]  = sri_len_DEF[i];")
        self._write("")
        self._write("        activation_units[i] = activation_units_DEF[i];")
        self._write("        prefactor_units[i]  = prefactor_units_DEF[i];")
        self._write("        phase_units[i]      = phase_units_DEF[i];")
        self._write("")
        self._write("        nTB[i]  = nTB_DEF[i];")
        self._write("        if (nTB[i] != 0) {")
        self._write("           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);")
        self._write("           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);")
        self._write("           for (int j=0; j<nTB[i]; j++) {")
        self._write("             TB[i][j] = TB_DEF[i][j];")
        self._write("             TBid[i][j] = TBid_DEF[i][j];")
        self._write("           }")
        self._write("        }")
        self._write("    }")
        self._write("}")
        self._write()
        self._write("void SetAllDefaults()")
        self._write("{")
        self._write("    for (int i=0; i<%d; i++) {" % (nReactions))
        self._write("        if (nTB_DEF[i] != 0) {")
        self._write("            nTB_DEF[i] = 0;")
        self._write("            free(TB_DEF[i]);")
        self._write("            free(TBid_DEF[i]);")
        self._write("        }")
        self._write("")
        self._write("        fwd_A_DEF[i]    = fwd_A[i];")
        self._write("        fwd_beta_DEF[i] = fwd_beta[i];")
        self._write("        fwd_Ea_DEF[i]   = fwd_Ea[i];")
        self._write("")
        self._write("        low_A_DEF[i]    = low_A[i];")
        self._write("        low_beta_DEF[i] = low_beta[i];")
        self._write("        low_Ea_DEF[i]   = low_Ea[i];")
        self._write("")
        self._write("        rev_A_DEF[i]    = rev_A[i];")
        self._write("        rev_beta_DEF[i] = rev_beta[i];")
        self._write("        rev_Ea_DEF[i]   = rev_Ea[i];")
        self._write("")
        self._write("        troe_a_DEF[i]    = troe_a[i];")
        self._write("        troe_Ts_DEF[i]   = troe_Ts[i];")
        self._write("        troe_Tss_DEF[i]  = troe_Tss[i];")
        self._write("        troe_Tsss_DEF[i] = troe_Tsss[i];")
        self._write("")
        self._write("        sri_a_DEF[i] = sri_a[i];")
        self._write("        sri_b_DEF[i] = sri_b[i];")
        self._write("        sri_c_DEF[i] = sri_c[i];")
        self._write("        sri_d_DEF[i] = sri_d[i];")
        self._write("        sri_e_DEF[i] = sri_e[i];")
        self._write("")
        self._write("        is_PD_DEF[i]    = is_PD[i];")
        self._write("        troe_len_DEF[i] = troe_len[i];")
        self._write("        sri_len_DEF[i]  = sri_len[i];")
        self._write("")
        self._write("        activation_units_DEF[i] = activation_units[i];")
        self._write("        prefactor_units_DEF[i]  = prefactor_units[i];")
        self._write("        phase_units_DEF[i]      = phase_units[i];")
        self._write("")
        self._write("        nTB_DEF[i]  = nTB[i];")
        self._write("        if (nTB_DEF[i] != 0) {")
        self._write("           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);")
        self._write("           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);")
        self._write("           for (int j=0; j<nTB_DEF[i]; j++) {")
        self._write("             TB_DEF[i][j] = TB[i][j];")
        self._write("             TBid_DEF[i][j] = TBid[i][j];")
        self._write("           }")
        self._write("        }")
        self._write("    }")
        self._write("}")

        self._write()
        self._write(self.line(' Finalizes parameter database'))
        self._write('void CKFINALIZE()')
        self._write('{')
        self._write('  for (int i=0; i<%d; ++i) {' % (nReactions))
        self._write('    free(TB[i]); TB[i] = 0; ')
        self._write('    free(TBid[i]); TBid[i] = 0;')
        self._write('    nTB[i] = 0;')
        self._write()
        self._write('    free(TB_DEF[i]); TB_DEF[i] = 0; ')
        self._write('    free(TBid_DEF[i]); TBid_DEF[i] = 0;')
        self._write('    nTB_DEF[i] = 0;')
        self._write('  }')
        self._write('}')
        self._write()

        self._write('#else')

        self._write(self.line(' TODO: Remove on GPU, right now needed by chemistry_module on FORTRAN'))
        self._write('AMREX_GPU_HOST_DEVICE void CKINIT'+sym+'()')
        self._write('{')
        self._write('}')
        self._write()
        self._write('AMREX_GPU_HOST_DEVICE void CKFINALIZE()')
        self._write('{')
        self._write('}')
        self._write()

        self._write('#endif')

        return



    ################################
    # CHEMKIN WRAPPERS
    ################################

    def _ckindx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('A few mechanism parameters'))
        self._write('void CKINDX'+sym+'(int * mm, int * kk, int * ii, int * nfit)')
        self._write('{')
        self._indent()
        self._write('*mm = %d;' % len(mechanism.element()))
        self._write('*kk = %d;' % len(mechanism.species()))
        self._write('*ii = %d;' % len(mechanism.reaction()))
        self._write('*nfit = -1; ' + self.line(
            'Why do you need this anyway ? '))
        
        # done
        self._outdent()
        self._write('}')
        return


    def _ckxnum(self, mechanism):
        self._write()
        self._write()
        self._write()
        self._write(self.line(' ckxnum... for parsing strings '))
        self._write('void CKXNUM'+sym+'(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline )')
        self._write('{')
        self._indent()
        self._write('int n,i; /*Loop Counters */')
        self._write('char cstr[1000];')
        self._write('char *saveptr;')
        self._write('char *p; /*String Tokens */')
        self._write(self.line(' Strip Comments '))
        self._write('for (i=0; i<lenline; ++i) {')
        self._indent()
        self._write('if (line[i]==\'!\') {')
        self._indent()
        self._write('break;')
        self._outdent()
        self._write('}')
        self._write('cstr[i] = line[i];')
        self._outdent()
        self._write('}')
        self._write('cstr[i] = \'\\0\';')
        self._write()
        self._write('p = strtok_r(cstr," ", &saveptr);')
        self._write('if (!p) {')
        self._indent()
        self._write('*nval = 0;')
        self._write('*kerr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('for (n=0; n<*nexp; ++n) {')
        self._indent()
        self._write('rval[n] = atof(p);')
        self._write('p = strtok_r(NULL, \" \", &saveptr);')
        self._write('if (!p) break;')
        self._outdent()
        self._write('}')
        self._write('*nval = n+1;')
        self._write('if (*nval < *nexp) *kerr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        return


    def _cksnum(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(' cksnum... for parsing strings '))
        self._write('void CKSNUM'+sym+'(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray)')
        self._write('{')
        self._indent()
        
        self._write(self.line('Not done yet ...'))
        
        # done
        self._outdent()
        self._write('}')
        return


    def _cksyme_str(self, mechanism):
        nElement = len(mechanism.element())
        self._write()
        self._write()
        self._write(
            self.line(' Returns the vector of strings of element names'))
        self._write('void CKSYME_STR'+sym+'(amrex::Vector<std::string>& ename)')
        self._write('{')
        self._indent()
        self._write('ename.resize(%d);' % nElement)
        for element in mechanism.element():
            self._write('ename[%d] = "%s";' %(element.id, element.symbol))
        # done
        self._outdent()
        self._write('}')
        return


    def _cksyme(self, mechanism):
        nElement = len(mechanism.element())
        self._write()
        self._write()
        self._write(
            self.line(' Returns the char strings of element names'))
        self._write('void CKSYME'+sym+'(int * kname, int * plenkname )')
        self._write('{')
        self._indent()

        self._write('int i; '+self.line('Loop Counter'))
        self._write('int lenkname = *plenkname;')
        self._write(self.line('clear kname'))
        self._write('for (i=0; i<lenkname*%d; i++) {' % nElement)
        self._indent()
        self._write('kname[i] = \' \';')
        self._outdent()
        self._write('}')
        self._write()
        for element in mechanism.element():
            self._write(self.line(' %s ' % element.symbol))
            ii = 0
            for char in element.symbol:
                self._write('kname[ %d*lenkname + %d ] = \'%s\';' %
                           (element.id, ii, char.capitalize()))
                ii = ii+1
            self._write('kname[ %d*lenkname + %d ] = \' \';' %
                           (element.id, ii))
            self._write()
            
        # done
        self._outdent()
        self._write('}')
        return
        

    def _cksyms_str(self, mechanism):
        nSpecies = len(mechanism.species())  
        self._write() 
        self._write()
        self._write(
            self.line(' Returns the vector of strings of species names'))
        self._write('void CKSYMS_STR'+sym+'(amrex::Vector<std::string>& kname)')
        self._write('{')
        self._indent()
        self._write('kname.resize(%d);' % nSpecies)
        for species in mechanism.species():
            self._write('kname[%d] = "%s";' %(species.id, species.symbol))

        self._outdent() 
        self._write('}') 
        return


    def _cksyms(self, mechanism):
        nSpecies = len(mechanism.species())
        self._write()
        self._write()
        self._write(
            self.line(' Returns the char strings of species names'))
        self._write('void CKSYMS'+sym+'(int * kname, int * plenkname )')
        self._write('{')
        self._indent()
        
        self._write('int i; '+self.line('Loop Counter'))
        self._write('int lenkname = *plenkname;')
        self._write(self.line('clear kname'))
        self._write('for (i=0; i<lenkname*%d; i++) {' % nSpecies)
        self._indent()
        self._write('kname[i] = \' \';')
        self._outdent()
        self._write('}')
        self._write()
        for species in mechanism.species():
            self._write(self.line(' %s ' % species.symbol))
            ii = 0
            for char in species.symbol:
                self._write('kname[ %d*lenkname + %d ] = \'%s\';' %
                           (species.id, ii, char.capitalize()))
                ii = ii+1
            self._write('kname[ %d*lenkname + %d ] = \' \';' %
                           (species.id, ii))
            self._write()

        # done
        self._outdent()
        self._write('}')
        return


    def _ckrp(self, mechanism):
        self._write()
        self._write()
        self._write(
            self.line(' Returns R, Rc, Patm' ))
        self._write('void CKRP'+sym+'(double *  ru, double *  ruc, double *  pa)')
        self._write('{')
        self._indent()
        
        self._write(' *ru  = %1.14e; ' % (R * mole * kelvin / erg))
        self._write(' *ruc = %.20f; ' % (Rc * mole * kelvin / cal))
        self._write(' *pa  = %g; ' % (Patm) )
        
        # done
        self._outdent()
        self._write('}')
        return
        

    def _ckpx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(x)'))
        self._write('void CKPX'+sym+'(double *  rho, double *  T, double *  x, double *  P)')
        self._write('{')
        self._indent()

        self._write('double XW = 0;'+
                    self.line(' To hold mean molecular wt'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(
            '*P = *rho * %g * (*T) / XW; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        return


    def _ckpy(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(y)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKPY'+sym+'(double *  rho, double *  T, double *  y,  double *  P)')
        self._write('{')
        self._indent()

        self._write('double YOW = 0;'+self.line(' for computing mean MW'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))

        self.line('YOW holds the reciprocal of the mean molecular wt')
        self._write(
            '*P = *rho * %g * (*T) * YOW; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 

    def _vckpy(self, mechanism):
        self._write()
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line('Compute P = rhoRT/W(y)'))
        self._write('void VCKPY'+sym+'(int *  np, double *  rho, double *  T, double *  y,  double *  P)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW[*np];')
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 0.0;')
        self._outdent()
        self._write('}')        
        self._write('')
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] += y[n*(*np)+i] * imw[n];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write('')

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write(
            'P[i] = rho[i] * %g * T[i] * YOW[i]; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write('#endif')

        return 


    def _ckpc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute P = rhoRT/W(c)'))
        self._write('void CKPC'+sym+'(double *  rho, double *  T, double *  c,  double *  P)')
        
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        nSpecies = len(mechanism.species())
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            '*P = *rho * %g * (*T) * sumC / W; ' % (R*kelvin*mole/erg)
            + self.line('P = rho*R*T/W'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 


    def _ckrhox(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute rho = PW(x)/RT'))
        self._write('void CKRHOX'+sym+'(double *  P, double *  T, double *  x,  double *  rho)')
        self._write('{')
        self._indent()

        self._write('double XW = 0;'+
                    self.line(' To hold mean molecular wt'))
        
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write(
            '*rho = *P * XW / (%g * (*T)); ' % (R*kelvin*mole/erg)
            + self.line('rho = P*W/(R*T)'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        return


    def _ckrhoy(self, mechanism):
        species = self.species
        nSpec = len(species)
        self._write()
        self._write()
        self._write(self.line('Compute rho = P*W(y)/RT'))
        self._write('AMREX_GPU_HOST_DEVICE void CKRHOY'+sym+'(double *  P, double *  T, double *  y,  double *  rho)')
        self._write('{')
        self._indent()
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('*rho = *P / (%g * (*T) * YOW);' % (R * mole * kelvin / erg) + self.line('rho = P*W/(R*T)'))
        self._write('return;')
        self._outdent()
        self._write('}')
        return 


    def _ckrhoc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Compute rho = P*W(c)/(R*T)'))
        self._write('void CKRHOC'+sym+'(double *  P, double *  T, double *  c,  double *  rho)')
        
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        self.line('W/sumC holds the mean molecular wt')
        self._write(
            '*rho = *P * W / (sumC * (*T) * %g); ' % (R*kelvin*mole/erg)
            + self.line('rho = PW/(R*T)'))
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 


    def _ckwt(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get molecular weight for all species'))
        self._write('void CKWT'+sym+'( double *  wt)')
        self._write('{')
        self._indent()
        # call molecularWeight
        self._write('get_mw(wt);')
        self._outdent()
        self._write('}')
        return

      
    def _ckawt(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get atomic weight for all elements'))
        self._write('void CKAWT'+sym+'( double *  awt)')
        self._write('{')
        self._indent()
        # call atomicWeight
        self._write('atomicWeight(awt);')
        self._outdent()
        self._write('}')
        self._write()
        return


    def _ckmmwy(self, mechanism):
        self._write()
        self._write(self.line('given y[species]: mass fractions'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKMMWY'+sym+'(double *  y,  double *  wtm)')
        self._write('{')
        self._indent()
        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('*wtm = 1.0 / YOW;')
        self._write('return;')
        self._outdent()
        self._write('}')

        return 
 

    def _ckmmwx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('given x[species]: mole fractions'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('void CKMMWX'+sym+'(double *  x,  double *  wtm)')
        self._write('{')
        self._indent()
        self._write('double XW = 0;'+self.line(' see Eq 4 in CK Manual'))
        # molecular weights of all species
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
        self._write('*wtm = XW;')
        self._write()
        self._write('return;')
        self._outdent()
        self._write('}')

        return 

 
    def _ckmmwc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('given c[species]: molar concentration'))
        self._write(self.line('returns mean molecular weight (gm/mole)'))
        self._write('void CKMMWC'+sym+'(double *  c,  double *  wtm)')
        self._write('{')
        self._indent()
        self._write('int id; ' + self.line('loop counter'))
        self._write(self.line('See Eq 5 in CK Manual'))
        self._write('double W = 0;')
        self._write('double sumC = 0;')
        # molecular weights of all species
        for species in self.species:
            self._write('W += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')
        self._write(self.line(' CK provides no guard against divison by zero'))
        self._write('*wtm = W/sumC;')
        self._write()
        self._write('return;')
        self._outdent()
        self._write('}')

        return 


    def _ckytx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to x[species] (mole fracs)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKYTX'+sym+'(double *  y,  double *  x)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double tmp[%d];' % (nSpec))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tmp[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += tmp[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write('double YOWINV = 1.0/YOW;')
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('x[i] = y[i]*imw[i]*YOWINV;')
        self._outdent()
        self._write('}')
        self._write('return;')
        self._outdent()
        self._write('}')

        return 

 
    def _vckytx(self, mechanism):
        self._write()
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line(
            'convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs)'))
        self._write('void VCKYTX'+sym+'(int *  np, double *  y,  double *  x)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW[*np];')
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 0.0;')
        self._outdent()
        self._write('}')        
        self._write('')
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('x[n*(*np)+i] = y[n*(*np)+i] * imw[n];')
        self._write('YOW[i] += x[n*(*np)+i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write('')

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('YOW[i] = 1.0/YOW[i];')
        self._outdent()
        self._write('}')

        self._write('')
        
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('x[n*(*np)+i] *=  YOW[i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('#else') 
        self._write(self.line('TODO: remove this on GPU'))
        self._write('void VCKYTX'+sym+'(int *  np, double *  y,  double *  x)')
        self._write('{')
        self._write('}')
        self._write('#endif') 

        return 


    def _ckytcp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to c[species] (molar conc)'))
        self._write('void CKYTCP'+sym+'(double *  P, double *  T, double *  y,  double *  c)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)
        self._write('double YOW = 0;')
        self._write('double PWORT;')
        self._write('')
        self._write(self.line('Compute inverse of mean molecular wt first'))
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = y[i]*imw[i];')
        self._outdent()
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('YOW += c[i];')
        self._outdent()
        self._write('}')
        self._write('')
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )

        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = PWORT * y[i] * imw[i];')
        self._outdent()
        self._write('}')
        self._write('return;')
        self._outdent()
        self._write('}')

        return 


    def _ckytcr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to c[species] (molar conc)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKYTCR'+sym+'(double *  rho, double *  T, double *  y,  double *  c)')
        self._write('{')
        self._indent()
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('c[i] = (*rho)  * y[i] * imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        return 



    #Pieces for mechanism.h#
    def _print_mech_header(self, mechanism):
        self._write()
        self._write("#ifndef MECHANISM_h")
        self._write("#define MECHANISM_h")
        self._write()
        self._write("#if 0")
        self._write("/* Elements")
        nb_elem = 0
        for element in mechanism.element():
            self._write('%d  %s' % (element.id, element.symbol) )
            nb_elem += 1
        self._write('*/')
        self._write("#endif")
        self._write()
        self._write('/* Species */')
        nb_spec = 0
        for species in mechanism.species():
            s = species.symbol.strip()
            # Ionic species
            if s[-1] == '-': s = s[:-1] + 'n'
            if s[-1] == '+': s = s[:-1] + 'p'
            # Excited species
            s = s.replace('*', 'D')
            # Remove other characters not allowed in preprocessor defines
            s = s.replace('-', '').replace('(','').replace(')','')
            self._write('#define %s_ID %d' % (s, species.id))
            nb_spec += 1
        self._write()
        self._write("#define NUM_ELEMENTS %d" % (nb_elem))
        self._write("#define NUM_SPECIES %d" % (nb_spec))
        self._write("#define NUM_REACTIONS %d" %(len(mechanism.reaction())))
        self._write()
        self._write("#define NUM_FIT 4")
        self._write("#endif")

        return

    #Pieces for mechanism.h#



    def _thermo_GPU(self, mechanism):
        speciesInfo = self._analyzeThermodynamics(mechanism)

        self._gibbs_GPU(speciesInfo)
        self._helmholtz_GPU(speciesInfo)
        ##self._dcvpRdT_GPU(speciesInfo)
        self._cv_GPU(speciesInfo)
        self._cp_GPU(speciesInfo)
        self._speciesInternalEnergy_GPU(speciesInfo)
        self._speciesEnthalpy_GPU(speciesInfo)
        self._speciesEntropy_GPU(speciesInfo)

        return


    def _trans_chop(self, mechanism):
        speciesTransport = self._analyzeTransport(mechanism)
        NLITE=0
        idxLightSpecs = []
        for spec in self.species:
            if spec.weight < 5.0:
                NLITE+=1
                idxLightSpecs.append(spec.id)
        self._miscTransInfo(KK=self.nSpecies, NLITE=NLITE, do_declarations=False)
        self._wt(False)
        self._eps(mechanism, speciesTransport, False)
        self._sig(mechanism, speciesTransport, False)
        self._dip(mechanism, speciesTransport, False)
        self._pol(mechanism, speciesTransport, False)
        self._zrot(mechanism, speciesTransport, False)
        self._nlin(mechanism, speciesTransport, False)

        self._viscosity(speciesTransport, False, NTFit=50)
        self._diffcoefs(speciesTransport, False, NTFit=50)
        self._lightSpecs(idxLightSpecs, False)
        self._thermaldiffratios(speciesTransport, idxLightSpecs, False, NTFit=50)

        return


    def _dthermodT(self, mechanism):
        speciesInfo = self._analyzeThermodynamics(mechanism)
        self._dcvpdT(speciesInfo)
        return

      
    def _ckcvml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get specific heat at constant volume as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKCVML'+sym+'(double *  T,  double *  cvml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cv_R(cvml, tc);')
        
        # convert cv/R to cv
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('cvml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return

       
    def _ckcpml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get specific heat at constant pressure as a '))
        self._write(self.line('function of T for all species (molar units)'))
        self._write('void CKCPML'+sym+'(double *  T,  double *  cpml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpml, tc);')
        
        # convert cp/R to cp
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('cpml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
     
    def _ckuml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get internal energy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKUML'+sym+'(double *  T,  double *  uml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(uml, tc);')
        
        # convert e/RT to e with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('uml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
      

    def _ckhml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get enthalpy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKHML'+sym+'(double *  T,  double *  hml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')
        
        # convert h/RT to h with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('hml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return
    

    def _ckgml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get standard-state Gibbs energy as a function '))
        self._write(self.line('of T for all species (molar units)'))
        self._write('void CKGML'+sym+'(double *  T,  double *  gml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('gibbs(gml, tc);')
        
        # convert g/RT to g with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('gml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return

    
    def _ckaml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get standard-state Helmholtz free energy as a '))
        self._write(self.line('function of T for all species (molar units)'))
        self._write('void CKAML'+sym+'(double *  T,  double *  aml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('helmholtz(aml, tc);')
        
        # convert A/RT to A with molar units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('aml[id] *= RT;')
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return

   
    def _cksml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the standard-state entropies in molar units'))
        self._write('void CKSML'+sym+'(double *  T,  double *  sml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sml, tc);')
        
        # convert s/R to s
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sml[id] *= %g;' % (R*kelvin*mole/erg) )
        self._outdent()
        self._write('}')
       
        self._outdent()

        self._write('}')

        return

 
    def _ckums(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns internal energy in mass units (Eq 30.)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKUMS'+sym+'(double *  T,  double *  ums)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(ums, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('ums[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return

    def _ckhms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns enthalpy in mass units (Eq 27.)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKHMS'+sym+'(double *  T,  double *  hms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hms, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('hms[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return


    def _vckhms(self, mechanism):
        self._write()
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line('Returns enthalpy in mass units (Eq 27.)'))
        self._write('void VCKHMS'+sym+'(int *  np, double *  T,  double *  hms)')
        self._write('{')
        self._indent()

        species = self.species
        nSpec = len(species)

        self._write('double tc[5], h[%d];' % nSpec)

        self._write()

        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('tc[0] = 0.0;')
        self._write('tc[1] = T[i];')
        self._write('tc[2] = T[i]*T[i];')
        self._write('tc[3] = T[i]*T[i]*T[i];')
        self._write('tc[4] = T[i]*T[i]*T[i]*T[i];')

        self._write()

        self._write('speciesEnthalpy(h, tc);')

        self._write()

        for ispec in range(nSpec):
            self._write('hms[%d*(*np)+i] = h[%d];' % (ispec, ispec))
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('for (int n=0; n<%d; n++) {' % (nSpec))
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('hms[n*(*np)+i] *= %g * T[i] * imw[n];' % (R*kelvin*mole/erg))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._outdent()
        self._write('}')
        self._write('#else')
        self._write(self.line('TODO: remove this on GPU'))
        self._write('void VCKHMS'+sym+'(int *  np, double *  T,  double *  hms)')
        self._write('{')
        self._write('}')
        self._write('#endif')

        return

    def _ckams(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns helmholtz in mass units (Eq 32.)'))
        self._write('void CKAMS'+sym+'(double *  T,  double *  ams)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('helmholtz(ams, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('ams[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return

    def _ckgms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns gibbs in mass units (Eq 31.)'))
        self._write('void CKGMS'+sym+'(double *  T,  double *  gms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('gibbs(gms, tc);')
        
        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('gms[i] *= RT*imw[i];')
        self._outdent()
        self._write('}')
        self._outdent()

        self._write('}')

        return


    def _ckcvms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the specific heats at constant volume'))
        self._write(self.line('in mass units (Eq. 29)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKCVMS'+sym+'(double *  T,  double *  cvms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cv_R(cvms, tc);')

        # convert cv/R to cv with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cvms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

        self._outdent()

        self._write('}')

        return


    def _ckcpms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the specific heats at constant pressure'))
        self._write(self.line('in mass units (Eq. 26)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKCPMS'+sym+'(double *  T,  double *  cpms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpms, tc);')
        

        # convert cp/R to cp with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('cpms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

        self._outdent()

        self._write('}')

        return


    def _cksms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the entropies in mass units (Eq 28.)'))
        self._write('void CKSMS'+sym+'(double *  T,  double *  sms)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sms, tc);')
        

        # convert s/R to s with mass units
        self._write(self.line('multiply by R/molecularweight'))
        for species in self.species:
            ROW = (R*kelvin*mole/erg) / species.weight
            self._write('sms[%d] *= %20.15e; ' % (
                species.id, ROW) + self.line('%s' % species.symbol))

        self._outdent()

        self._write('}')

        return

    
    def _ckcpbl(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CP (Eq. 33)'))
        self._write('void CKCPBL'+sym+'(double *  T, double *  x,  double *  cpbl)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cpor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*cpor[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cpbl = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

 
    def _ckcpbs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CP (Eq. 34)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKCPBS'+sym+'(double *  T, double *  y,  double *  cpbs)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cpor[%d], tresult[%d]; ' % (self.nSpecies,self.nSpecies) + self.line(' temporary storage'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        

        species = self.species
        nSpec = len(species)
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('tresult[i] = cpor[i]*y[i]*imw[i];')
        self._outdent()
        self._write('')
        self._write('}')
        self._write('for (int i = 0; i < %d; i++)' % (nSpec))
        self._write('{')
        self._indent()
        self._write('result += tresult[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cpbs = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return


    def _ckcvbl(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CV (Eq. 35)'))
        self._write('void CKCVBL'+sym+'(double *  T, double *  x,  double *  cvbl)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cvor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cv_R(cvor, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*cvor[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*cvbl = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _ckcvbs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean specific heat at CV (Eq. 36)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKCVBS'+sym+'(double *  T, double *  y,  double *  cvbs)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double cvor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write('cv_R(cvor, tc);')
        
        # do dot product
        self._write(self.line('multiply by y/molecularweight'))
        for species in self.species:
            self._write('result += cvor[%d]*y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) + self.line('%s' % species.symbol))

        self._write()
        self._write('*cvbs = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _ckhbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the mean enthalpy of the mixture in molar units'))
        self._write('void CKHBML'+sym+'(double *  T, double *  x,  double *  hbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double hml[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*hml[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*hbml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
 
 
    def _ckhbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean enthalpy of mixture in mass units'))
        self._write('AMREX_GPU_HOST_DEVICE void CKHBMS'+sym+'(double *  T, double *  y,  double *  hbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0;')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double hml[%d], tmp[%d]; ' % (self.nSpecies,self.nSpecies) + self.line(' temporary storage'))
        
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesEnthalpy(hml, tc);')

        self._write('int id;')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('tmp[id] = y[id]*hml[id]*imw[id];')
        self._outdent()
        self._write('}')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += tmp[id];')
        self._outdent()
        self._write('}')

        self._write()
        # finally, multiply by RT
        self._write('*hbms = result * RT;')
        
        self._outdent()

        self._write('}')
        
        return

    def _ckubml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mean internal energy in molar units'))
        self._write('void CKUBML'+sym+'(double *  T, double *  x,  double *  ubml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double uml[%d]; ' % self.nSpecies + self.line(' temporary energy array'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(uml, tc);')
        
        # dot product
        self._write()
        self._write(self.line('perform dot product'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*uml[id];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*ubml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
 
    def _ckubms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mean internal energy in mass units'))
        self._write('AMREX_GPU_HOST_DEVICE void CKUBMS'+sym+'(double *  T, double *  y,  double *  ubms)')
        self._write('{')
        self._indent()

        self._write('double result = 0;')
        
        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double ums[%d]; ' % self.nSpecies + self.line(' temporary energy array'))
        
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        
        # call routine
        self._write('speciesInternalEnergy(ums, tc);')

        # convert e/RT to e with mass units
        self._write(self.line('perform dot product + scaling by wt'))
        for species in self.species:
            self._write('result += y[%d]*ums[%d]*imw[%d]; ' % (
                species.id, species.id, species.id)
                        + self.line('%s' % species.symbol))

        
        self._write()
        # finally, multiply by RT
        self._write('*ubms = result * RT;')
        
        self._outdent()

        self._write('}')
        
        return

    def _cksbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mixture entropy in molar units'))
        self._write('void CKSBML'+sym+'(double *  P, double *  T, double *  x,  double *  sbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double sor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        # Equation 42
        self._write()
        self._write(self.line('Compute Eq 42'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(sor[id]-log((x[id]+%g))-logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*sbml = result * %g;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _cksbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get mixture entropy in mass units'))
        self._write('void CKSBMS'+sym+'(double *  P, double *  T, double *  y,  double *  sbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double sor[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write('double YOW = 0; '+self.line('See Eq 4, 6 in CK Manual'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 42 and 43'))
        for species in self.species:
            self._write('result += x[%d]*(sor[%d]-log((x[%d]+%g))-logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by R/W'))
        self._write('*sbms = result * %g * YOW;' % (R*kelvin*mole/erg) )
        
        self._outdent()

        self._write('}')

        return

    def _ckgbml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean gibbs free energy in molar units'))
        self._write('void CKGBML'+sym+'(double *  P, double *  T, double *  x,  double *  gbml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double gort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write(self.line('Compute g/RT'))
        self._write('gibbs(gort, tc);')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(gort[id]+log((x[id]+%g))+logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*gbml = result * RT;')
        
        self._outdent()

        self._write('}')

        return


    def _ckgbms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mixture gibbs free energy in mass units'))
        self._write('void CKGBMS'+sym+'(double *  P, double *  T, double *  y,  double *  gbms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double gort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write(
            'double YOW = 0; '
            + self.line('To hold 1/molecularweight'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('gibbs(gort, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        for species in self.species:
            self._write('result += x[%d]*(gort[%d]+log((x[%d]+%g))+logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('*gbms = result * RT * YOW;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckabml(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mean helmholtz free energy in molar units'))
        self._write('void CKABML'+sym+'(double *  P, double *  T, double *  x,  double *  abml)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double aort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        
        # call routine
        self._write(self.line('Compute g/RT'))
        self._write('helmholtz(aort, tc);')
        
        # Equation 44
        self._write()
        self._write(self.line('Compute Eq 44'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('result += x[id]*(aort[id]+log((x[id]+%g))+logPratio);' %
                    smallnum )
        self._outdent()
        self._write('}')

        self._write()
        
        self._write('*abml = result * RT;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckabms(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns mixture helmholtz free energy in mass units'))
        self._write('void CKABMS'+sym+'(double *  P, double *  T, double *  y,  double *  abms)')
        self._write('{')
        self._indent()

        self._write('double result = 0; ')
        
        # get temperature cache
        self._write(self.line('Log of normalized pressure in cgs units dynes/cm^2 by Patm'))
        self._write( 'double logPratio = log ( *P / 1013250.0 ); ')
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double RT = %g*tT; ' % (R*kelvin*mole/erg)
            + self.line('R*T'))
        self._write(
            'double aort[%d]; ' % self.nSpecies + self.line(' temporary storage'))
        self._write(
            'double x[%d]; ' % self.nSpecies + self.line(' need a ytx conversion'))

        self._write(
            'double YOW = 0; '
            + self.line('To hold 1/molecularweight'))
        
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        # now to ytx
        self._write(self.line('Now compute y to x conversion'))
        for species in self.species:
            self._write('x[%d] = y[%d]/(%f*YOW); ' % (
                species.id, species.id, species.weight) )
            
        # call routine
        self._write('helmholtz(aort, tc);')
        
        # Equation 42 and 43
        self._write(self.line('Perform computation in Eq 44'))
        for species in self.species:
            self._write('result += x[%d]*(aort[%d]+log((x[%d]+%g))+logPratio);' %
                        (species.id, species.id, species.id, smallnum) )

        self._write(self.line('Scale by RT/W'))
        self._write('*abms = result * RT * YOW;')
        
        self._outdent()

        self._write('}')

        return
    

    def _ckwc(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('compute the production rate for each species'))
        self._write('AMREX_GPU_HOST_DEVICE void CKWC'+sym+'(double *  T, double *  C,  double *  wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, C, *T);')

        # convert C and wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6;')
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    def _ckwyp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKWYP'+sym+'(double *  P, double *  T, double *  y,  double *  wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        self._write('double YOW = 0; ')
        self._write('double PWORT; ')
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )
        
        self._write(self.line('multiply by 1e6 so c goes to SI'))
        self._write('PWORT *= 1e6; ')

        # now compute conversion
        self._write(self.line('Now compute conversion (and go to SI)'))
        for species in self.species:
            self._write('c[%d] = PWORT * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )

        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckwxp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKWXP'+sym+'(double *  P, double *  T, double *  x,  double *  wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckwyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('AMREX_GPU_HOST_DEVICE void CKWYR'+sym+'(double *  rho, double *  T, double *  y,  double *  wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        for species in self.species:
            self._write('c[%d] = 1e6 * (*rho) * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )
            
        # call productionRate
        self._write()
        self._write(self.line('call productionRate'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _vckwyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void VCKWYR'+sym+'(int *  np, double *  rho, double *  T,')
        self._write('	    double *  y,')
        self._write('	    double *  wdot)')
        self._write('{')
        self._write('#ifndef AMREX_USE_CUDA')
        self._indent()

        self._write('double c[%d*(*np)]; ' % self.nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        self._write('for (int n=0; n<%d; n++) {' % self.nSpecies)
        self._indent()
        self._write('for (int i=0; i<(*np); i++) {')
        self._indent()
        self._write('c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        # call productionRate
        self._write()
        self._write(self.line('call productionRate'))
        self._write('vproductionRate(*np, wdot, c, T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (int i=0; i<%d*(*np); i++) {' % self.nSpecies)
        self._indent()
        self._write('wdot[i] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()
        self._write('#endif')

        self._write('}')

        return


    def _ckwxr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('Returns the molar production rate of species'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKWXR'+sym+'(double *  rho, double *  T, double *  x,  double *  wdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % self.nSpecies + self.line('temporary storage'))
        
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Extra 1e6 factor to take c to SI'))
        self._write('ROW = 1e6*(*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('productionRate(wdot, c, *T);')

        # convert wdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('wdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _cknu(self, mechanism):

        nSpecies  = len(mechanism.species())
        nReaction = len(mechanism.reaction())

        self._write()
        self._write()
        self._write(self.line('Returns the stoichiometric coefficients'))
        self._write(self.line('of the reaction mechanism. (Eq 50)'))
        self._write('void CKNU'+sym+'(int * kdim,  int * nuki)')
        self._write('{')
        self._indent()

 
        self._write('int id; ' + self.line('loop counter'))
        self._write('int kd = (*kdim); ')
        self._write(self.line('Zero nuki'))
        self._write('for (id = 0; id < %d * kd; ++ id) {' % (nSpecies) )
        self._indent()
        self._write(' nuki[id] = 0; ')
        self._outdent()
        self._write('}')
        
        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            for symbol, coefficient in reaction.reactants:
                self._write(
                    "nuki[ %d * kd + %d ] += -%f ;"
                    % (mechanism.species(symbol).id, reaction.id-1, coefficient))

            for symbol, coefficient in reaction.products:
                self._write(
                    "nuki[ %d * kd + %d ] += +%f ;"
                    % (mechanism.species(symbol).id, reaction.id-1, coefficient))
       
        # done
        self._outdent()
        self._write('}')

        return


    def _ckinu(self, mechanism):

        nSpecies  = len(mechanism.species())
        nReaction = len(mechanism.reaction())

        self._write()
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line('Returns a count of species in a reaction, and their indices'))
        self._write(self.line('and stoichiometric coefficients. (Eq 50)'))
        self._write('void CKINU'+sym+'(int * i, int * nspec, int * ki, int * nu)')
        self._write('{')
        self._indent()

        self._write("if (*i < 1) {")
        self._indent()

        maxsp = 0
        for reaction in mechanism.reaction():
            maxsp = max(maxsp,len(reaction.reactants) + len(reaction.products))

        self._write(self.line('Return max num species per reaction'))
        self._write("*nspec = %d;" % (maxsp))
        self._outdent()
        self._write("} else {")
        self._indent()
        self._write("if (*i > %d) {" % (nReaction))
        self._indent()
        self._write("*nspec = -1;")
        self._outdent()
        self._write("} else {")
        self._indent()
        self._write("*nspec = kiv[*i-1].size();")
        self._write("for (int j=0; j<*nspec; ++j) {")
        self._indent()
        self._write("ki[j] = kiv[*i-1][j] + 1;")
        self._write("nu[j] = nuv[*i-1][j];")
        self._outdent()
        self._write("}")
        self._outdent()
        self._write("}")
        self._outdent()
        self._write("}")
        # done
        self._outdent()
        self._write('}')
        self._write('#endif')

        return


    def _ckncf(self, mechanism):

        nSpecies  = len(mechanism.species())
        nElement  = len(mechanism.element())

        self._write()
        self._write()
        self._write(self.line('Returns the elemental composition '))
        self._write(self.line('of the speciesi (mdim is num of elements)'))
        self._write('void CKNCF'+sym+'(int * ncf)')
        self._write('{')
        self._indent()

 
        self._write('int id; ' + self.line('loop counter'))
        self._write('int kd = %d; ' % (nElement))
        self._write(self.line('Zero ncf'))
        self._write('for (id = 0; id < kd * %d; ++ id) {' % (self.nSpecies) )
        self._indent()
        self._write(' ncf[id] = 0; ')
        self._outdent()
        self._write('}')
        
        self._write()
        for species in mechanism.species():
           self._write(self.line('%s' % species.symbol))
           for elem, coef in species.composition:
               self._write('ncf[ %d * kd + %d ] = %d; ' % (
                   species.id, mechanism.element(elem).id, coef) +
                       self.line('%s' % elem) )
                           
           self._write()
                            
        # done
        self._outdent()

        self._write('}')

        return


    def _ckabe(self, mechanism):

        nElement  = len(mechanism.element())

        self._write()
        self._write()
        self._write(self.line('Returns the arrehenius coefficients '))
        self._write(self.line('for all reactions'))
        self._write('void CKABE'+sym+'( double *  a, double *  b, double *  e)')
        self._write('{')
        self._indent()

        nReactions = len(mechanism.reaction())
        for j in range(nReactions):
            reaction   = mechanism.reaction(id=j)
            A, beta, E = reaction.arrhenius
            self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
            self._write("a[%d] = %.17g;" %(j,A))
            self._write("b[%d] = %.17g;" %(j,beta))
            self._write("e[%d] = %.17g;" %(j,E))
            self._write()

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return
                            
 
    def _ckxty(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to y[species] (mass fracs)'))
        self._write('AMREX_GPU_HOST_DEVICE void CKXTY'+sym+'(double *  x,  double *  y)')
        self._write('{')
        self._indent()

        self._write('double XW = 0; '+self.line('See Eq 4, 9 in CK Manual'))
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))
 
        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('double XWinv = 1.0/XW;')
        for species in self.species:
            self._write('y[%d] = x[%d]*%f*XWinv; ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckxtcp(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to c[species] (molar conc)'))
        self._write('void CKXTCP'+sym+'(double *  P, double *  T, double *  x,  double *  c)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double PORT = (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('P/RT'))
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckxtcr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert x[species] (mole fracs) to c[species] (molar conc)'))
        self._write('void CKXTCR'+sym+'(double *  rho, double *  T, double *  x, double *  c)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write('ROW = (*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckctx(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert c[species] (molar conc) to x[species] (mole fracs)'))
        self._write('void CKCTX'+sym+'(double *  c, double *  x)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))
        self._write('double sumC = 0; ')

        self._write()
        self._write(self.line('compute sum of c '))
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('sumC += c[id];')
        self._outdent()
        self._write('}')

        # now compute conversion
        self._write()
        self._write(self.line(' See Eq 13 '))
        self._write('double sumCinv = 1.0/sumC;')
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('x[id] = c[id]*sumCinv;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckcty(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert c[species] (molar conc) to y[species] (mass fracs)'))
        self._write('void CKCTY'+sym+'(double *  c, double *  y)')
        self._write('{')
        self._indent()

        self._write('double CW = 0; '+self.line('See Eq 12 in CK Manual'))
        
        # compute denominator in eq 12
        self._write(self.line('compute denominator in eq 12 first'))
        for species in self.species:
            self._write('CW += c[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Now compute conversion'))
        self._write('double CWinv = 1.0/CW;')
        for species in self.species:
            self._write('y[%d] = c[%d]*%f*CWinv; ' % (
                species.id, species.id, species.weight) )

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
    def _ckcpor(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get Cp/R as a function of T '))
        self._write(self.line('for all species (Eq 19)'))
        self._write('void CKCPOR'+sym+'(double *  T, double *  cpor)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('cp_R(cpor, tc);')
        
        self._outdent()

        self._write('}')

        return
    
    def _ckhort(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get H/RT as a function of T '))
        self._write(self.line('for all species (Eq 20)'))
        self._write('void CKHORT'+sym+'(double *  T, double *  hort)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEnthalpy(hort, tc);')
        
        self._outdent()

        self._write('}')

        return
 
    def _cksor(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('get S/R as a function of T '))
        self._write(self.line('for all species (Eq 21)'))
        self._write('void CKSOR'+sym+'(double *  T, double *  sor)')
        self._write('{')
        self._indent()

        # get temperature cache
        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        
        # call routine
        self._write('speciesEntropy(sor, tc);')
        
        self._outdent()

        self._write('}')

        return


    def _ckqc(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write()
        self._write(self.line('Returns the rate of progress for each reaction'))
        self._write('void CKQC'+sym+'(double *  T, double *  C, double *  qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        # convert C to SI units
        self._write()
        self._write(self.line('convert to SI'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e6;')
        self._outdent()
        self._write('}')
        
        # call productionRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, C, *T);')

        # convert C to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('C[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')

        # convert qdot to chemkin units
        self._write()
        self._write('for (id = 0; id < %d; ++id) {' % nReactions)
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    
    def _ckkfkr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKKFKR'+sym+'(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call progressRateFR
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRateFR(q_f, q_r, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('q_f[id] *= 1.0e-6;')
        self._write('q_r[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqyp(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKQYP'+sym+'(double *  P, double *  T, double *  y, double *  qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        self._write('double YOW = 0; ')
        self._write('double PWORT; ')
        
        # compute inverse of mean molecular weight first (eq 3)
        self._write(self.line('Compute inverse of mean molecular wt first'))
        for species in self.species:
            self._write('YOW += y[%d]*imw[%d]; ' % (
                species.id, species.id) + self.line('%s' % species.symbol))
 
        self._write(self.line('PW/RT (see Eq. 7)'))
        self._write('PWORT = (*P)/(YOW * %g * (*T)); ' % (R*kelvin*mole/erg) )
        
        self._write(self.line('multiply by 1e6 so c goes to SI'))
        self._write('PWORT *= 1e6; ')

        # now compute conversion
        self._write(self.line('Now compute conversion (and go to SI)'))
        for species in self.species:
            self._write('c[%d] = PWORT * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )

        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqxp(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKQXP'+sym+'(double *  P, double *  T, double *  x, double *  qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double PORT = 1e6 * (*P)/(%g * (*T)); ' % (R*kelvin*mole/erg) +
                    self.line('1e6 * P/RT so c goes to SI units'))
        
        # now compute conversion
        self._write()
        self._write(self.line('Compute conversion, see Eq 10'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*PORT;')
        self._outdent()
        self._write('}')
        
        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqyr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void CKQYR'+sym+'(double *  rho, double *  T, double *  y, double *  qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))

        # now compute conversion
        self._write(self.line('See Eq 8 with an extra 1e6 so c goes to SI'))
        for species in self.species:
            self._write('c[%d] = 1e6 * (*rho) * y[%d]*imw[%d]; ' % (
                species.id, species.id, species.id) )
            
        # call progressRate
        self._write()
        self._write(self.line('call progressRate'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return


    def _ckqxr(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())
        
        self._write()
        self._write()
        self._write(self.line('Returns the progress rates of each reactions'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKQXR'+sym+'(double *  rho, double *  T, double *  x, double *  qdot)')
        self._write('{')
        self._indent()

        self._write('int id; ' + self.line('loop counter'))

        self._write('double c[%d]; ' % nSpecies + self.line('temporary storage'))
        
        self._write('double XW = 0; '+self.line('See Eq 4, 11 in CK Manual'))
        self._write('double ROW; ')
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('XW += x[%d]*%f; ' % (
                species.id, species.weight) + self.line('%s' % species.symbol))

        # now compute conversion
        self._write(self.line('Extra 1e6 factor to take c to SI'))
        self._write('ROW = 1e6*(*rho) / XW;')
        self._write()
        self._write(self.line('Compute conversion, see Eq 11'))
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('c[id] = x[id]*ROW;')
        self._outdent()
        self._write('}')
        
        # call progressRate
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('progressRate(qdot, c, *T);')

        # convert qdot to chemkin units
        self._write()
        self._write(self.line('convert to chemkin units'))
        self._write('for (id = 0; id < %d; ++id) {' % nReactions )
        self._indent()
        self._write('qdot[id] *= 1.0e-6;')
        self._outdent()
        self._write('}')
        
        self._outdent()

        self._write('}')

        return

    
    def __ckeqcontent(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write(
            'double tT = *T; '
            + self.line('temporary temperature'))
        self._write(
            'double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; '
            + self.line('temperature cache'))
        self._write(
            'double gort[%d]; ' % nSpecies + self.line(' temporary storage'))

        # compute the gibbs free energy
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(gort, tc);')

        # compute the equilibrium constants
        self._write()
        self._write(self.line('compute the equilibrium constants'))
        self._write('equilibriumConstants(eqcon, gort, tT);')

        for reaction in mechanism.reaction():

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            somepow = 0
            for symbol, coefficient in reaction.reactants:
                somepow = somepow - coefficient

            for symbol, coefficient in reaction.products:
                somepow = somepow + coefficient

            if somepow == 0:
                self._write(self.line(
                    'eqcon[%d] *= %g; ' % (reaction.id-1, (1e-6)**somepow) ) )
                
            else:
                self._write( 'eqcon[%d] *= %g; ' % (reaction.id-1, (1e-6)**somepow) ) 


    def _ckeqc(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write('void CKEQC'+sym+'(double *  T, double *  C, double *  eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
                
        self._outdent()

        self._write('}')

        return

    
    def _ckeqyp(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given P, T, and mass fractions'))
        self._write('void CKEQYP'+sym+'(double *  P, double *  T, double *  y, double *  eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqxp(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given P, T, and mole fractions'))
        self._write('void CKEQXP'+sym+'(double *  P, double *  T, double *  x, double *  eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqyr(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given rho, T, and mass fractions'))
        self._write('void CKEQYR'+sym+'(double *  rho, double *  T, double *  y, double *  eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


    def _ckeqxr(self, mechanism):

        import pyre
        periodic = pyre.handbook.periodicTable()

        self._write()
        self._write()
        self._write(self.line('Returns the equil constants for each reaction'))
        self._write(self.line('Given rho, T, and mole fractions'))
        self._write('void CKEQXR'+sym+'(double *  rho, double *  T, double *  x, double *  eqcon)')
        self._write('{')
        self._indent()

        self.__ckeqcontent(mechanism)
        
        self._outdent()

        self._write('}')

        return


# Fuego Extensions. All functions in this section has the fe prefix
# All fuctions in this section uses the standard fuego chemkin functions
    def _ck_eytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        self._write(self.line(
            'get temperature given internal energy in mass units and mass fracs'))
        self._write('int feeytt'+fsym+'(double *  e, double *  y, double *  t)')
        self._write('{')
        self._indent()

        self._write('const int maxiter = 50;')
        self._write('const double tol  = 0.001;')
        self._write('double ein  = *e;')
        self._write('double tmin = %g; // max lower bound for thermo def' % lowT)
        self._write('double tmax = %g; // min upper bound for thermo def' % highT)
        self._write('double e1,emin,emax,cv,t1,dt;')
        self._write('int i; // loop counter')
        self._write('CKUBMS'+sym+'(&tmin, y, &emin);')
        self._write('CKUBMS'+sym+'(&tmax, y, &emax);')
        self._write('if (ein < emin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCVBS'+sym+'(&tmin, y, &cv);')
        self._write('*t = tmin - (emin-ein)/cv;')
        self._write('return 1;')
        self._outdent()
        self._write('}')
        
        self._write('if (ein > emax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCVBS'+sym+'(&tmax, y, &cv);')
        self._write('*t = tmax - (emax-ein)/cv;')
        self._write('return 1;')
        self._outdent()
        self._write('}')

        self._write('t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKUBMS'+sym+'(&t1,y,&e1);')
        self._write('CKCVBS'+sym+'(&t1,y,&cv);')
        self._write('dt = (ein - e1) / cv;')
        self._write('if (dt > 100) { dt = 100; }')
        self._write('else if (dt < -100) { dt = -100; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        
        self._write('*t = t1;')
        self._write('return 0;')
        
        self._outdent()

        self._write('}')

        return

 
    def _ck_hytt(self, mechanism):

        nSpecies = len(mechanism.species())
        lowT,highT,dummy = self._analyzeThermodynamics(mechanism)
        
        self._write()
        self._write()
        self._write(self.line(
            'get temperature given enthalpy in mass units and mass fracs'))
        self._write('int fehytt'+fsym+'(double *  h, double *  y, double *  t)')
        self._write('{')
        self._indent()

        self._write('const int maxiter = 50;')
        self._write('const double tol  = 0.001;')
        self._write('double hin  = *h;')
        self._write('double tmin = %g; // max lower bound for thermo def' % lowT)
        self._write('double tmax = %g; // min upper bound for thermo def' % highT)
        self._write('double h1,hmin,hmax,cp,t1,dt;')
        self._write('int i; // loop counter')
        self._write('CKHBMS'+sym+'(&tmin, y, &hmin);')
        self._write('CKHBMS'+sym+'(&tmax, y, &hmax);')
        self._write('if (hin < hmin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCPBS'+sym+'(&tmin, y, &cp);')
        self._write('*t = tmin - (hmin-hin)/cp;')
        self._write('return 1;')
        self._outdent()
        self._write('}')
        
        self._write('if (hin > hmax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCPBS'+sym+'(&tmax, y, &cp);')
        self._write('*t = tmax - (hmax-hin)/cp;')
        self._write('return 1;')
        self._outdent()
        self._write('}')

        self._write('t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKHBMS'+sym+'(&t1,y,&h1);')
        self._write('CKCPBS'+sym+'(&t1,y,&cp);')
        self._write('dt = (hin - h1) / cp;')
        self._write('if (dt > 100) { dt = 100; }')
        self._write('else if (dt < -100) { dt = -100; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        
        self._write('*t = t1;')
        self._write('return 0;')
        
        self._outdent()

        self._write('}')

        return

 
    def _ck_phity(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert phi[species] (specific mole nums) to y[species] (mass fracs)'))
        self._write('void fephity'+fsym+'(double *  phi, double *  y)')
        self._write('{')
        self._indent()

        self._write('double XW  = 0; ')
        self._write('int id; ' + self.line('loop counter'))
        
        # compute mean molecular weight first (eq 3)
        self._write(self.line('Compute mean molecular wt first'))
        for species in self.species:
            self._write('y[%d] = phi[%d]*%f;   XW += y[%d]; ' % (
                species.id, species.id, species.weight, species.id) +
                        self.line('%s' % species.symbol))
 
        self._write('for (id = 0; id < %d; ++id) {' % self.nSpecies)
        self._indent()
        self._write('y[id] = y[id]/XW;')
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
 
    def _ck_ytphi(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'convert y[species] (mass fracs) to phi[species] (specific mole num)'))
        self._write('void feytphi'+fsym+'(double *  y, double *  phi)')
        self._write('{')
        self._indent()

        for species in self.species:
            self._write('phi[%d] = y[%d]/%15.8e; ' % (
                species.id, species.id, species.weight/1000.0) +
                        self.line('%s (wt in kg)' % species.symbol))
 
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _ck_ctyr(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'reverse of ytcr, useful for rate computations'))
        self._write('void fectyr'+fsym+'(double *  c, double *  rho, double *  y)')
        self._write('{')
        self._indent()

        # now compute conversion
        for species in self.species:
            self._write('y[%d] = c[%d] * %f / (*rho); ' % (
                species.id, species.id, species.weight) )
        
        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return 
 
 
    # CAREFULL : need to remove rwrk dependencies before using this one
    def _ck_cvrhs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'ddebdf compatible right hand side of CV burner'))
        self._write(self.line(
            'rwrk[0] and rwrk[1] should contain rho and ene respectively'))
        self._write(self.line(
            'working variable phi contains specific mole numbers'))
        self._write('void fecvrhs'+fsym+'(double *  time, double *  phi, double *  phidot)')

	self._write('{')
	self._indent()
	# main body
        self._write('double rho,ene; ' + self.line('CV Parameters'))
        self._write('double y[%s], wdot[%s]; ' % (self.nSpecies, self.nSpecies) +
                    self.line('temporary storage'))
        self._write('int i; ' + self.line('Loop counter'))
        self._write('double temperature,pressure; ' + self.line('temporary var'))
        self._write('rho = rwrk[0];')
        self._write('ene = rwrk[1];')
        self._write('fephity'+fsym+'(phi, y);')
        self._write('feeytt'+fsym+'(&ene, y, &temperature);')
        self._write('CKPY'+sym+'(&rho, &temperature,  y, &pressure);')
        self._write('CKWYP'+sym+'(&pressure, &temperature,  y, wdot);')
        self._write('for (i=0; i<%s; ++i) phidot[i] = wdot[i] / (rho/1000.0); ' % self.nSpecies)
        self._write()
        self._write('return;')

	self._outdent()
	self._write('}')
	return


    def _ck_cvdim(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the dimensionality of the cv burner (number of species)'))
        self._write('int fecvdim'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return %d;' % self.nSpecies)

	self._outdent()
	self._write('}')
	return

 
    # CAREFULL : need to remove rwrk dependencies before using this one
    def _ck_zndrhs(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'ddebdf compatible right hand side of ZND solver'))
        self._write(self.line( 'rwrk[0] : scaling factor for pressure'))
        self._write(self.line( 'rwrk[1] : preshock density (g/cc) '))
        self._write(self.line( 'rwrk[2] : detonation velocity (cm/s) '))
        self._write(self.line( 'solution vector: [P; rho; y0 ... ylast] '))
        self._write('void fezndrhs'+fsym+'(double *  time, double *  z, double *  zdot)')

	self._write('{')
	self._indent()
	# main body
        self._write('double psc,rho1,udet; ' + self.line('ZND Parameters'))
        self._write('double wt[%s], hms[%s], wdot[%s]; ' %
                    (self.nSpecies, self.nSpecies, self.nSpecies) +
                    self.line('temporary storage'))
        self._write('int i; ' + self.line('Loop counter'))
        self._write(self.line('temporary variables'))
        self._write('double ru, T, uvel, wtm, p, rho, gam, son, xm, sum, drdy, eta, cp, cv ;')
        self._write('double *  y; ' + self.line('mass frac pointer'))
        self._write()
        self._write('ru = %g;' % (R * mole * kelvin / erg))
        self._write()
        self._write('psc = rwrk[0];')
        self._write('rho1 = rwrk[1];')
        self._write('udet = rwrk[2];')
        self._write()
        self._write('p = z[0] * psc;')
        self._write('rho = z[1];')
        self._write()
        self._write('y = &z[3];')
        self._write()
        self._write('CKMMWY'+sym+'(y, 0, 0, &wtm);')
        self._write()
        self._write('T = p * wtm / rho / ru;')
        self._write()
        self._write('uvel = (rho1 * udet)/ rho;')
        self._write()
        self._write('CKCPBS'+sym+'(&T, y, 0, 0, &cp);')
        self._write('CKCVBS'+sym+'(&T, y, 0, 0, &cv);')
        self._write('gam = cp/cv;')
        self._write()
        self._write('son = sqrt(fabs(gam*ru*T/wtm));')
        self._write('xm = uvel/son;')
        self._write()
        self._write('CKHMS'+sym+'(&T, 0, 0, hms);')
        self._write('CKWT'+sym+'(0, 0, wt);')
        self._write('CKWYP'+sym+'(&p, &T, y, 0, 0, wdot);')
        self._write()
        self._write('sum = 0.0;')
        self._write('for (i=0; i<%s; ++i) {' % self.nSpecies)
        self._indent()
        self._write('zdot[i+3] = wdot[i] * wt[i] / rho;')
        self._write('drdy = -rho * wtm / wt[i];')
        self._write('sum += -( drdy + rho * hms[i]/ (cp*T) ) * zdot[i+3];')
        self._outdent()
        self._write('}')
        self._write()
        self._write('eta = 1.0 - xm*xm;')
        self._write('zdot[0] = -(uvel*uvel/eta/psc)*sum;')
        self._write('zdot[1] = -sum/eta;')
        self._write('zdot[2] = uvel;')
        self._write()
        self._write('return;')

	self._outdent()
	self._write('}')
	return


    def _ck_znddim(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the dimensionality of the ZND solver (3+number of species)'))
        self._write('int feznddim'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return %d;' % (self.nSpecies + 3) )

	self._outdent()
	self._write('}')
	return
    
    def _ck_mechfile(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the name of the source mechanism file '))
        self._write('char* femechfile'+fsym+'()')

	self._write('{')
	self._indent()
	# main body
        self._write('return "%s";' % mechanism.name())

	self._outdent()
	self._write('}')
	return

    def _ck_symnum(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the species number'))
        self._write('int fesymnum'+fsym+'(const char* s1)')

	self._write('{')
	self._indent()
        
        for species in self.species:
            self._write('if (strcmp(s1, "%s")==0) return %d; ' % (
                species.symbol, species.id))
 
        self._write(self.line( 'species name not found' ))
        self._write('return -1;')

	self._outdent()
	self._write('}')
	return
    
    def _ck_symname(self, mechanism):
        self._write()
        self._write()
        self._write(self.line(
            'returns the species name'))
        self._write('char* fesymname'+fsym+'(int sn)')

	self._write('{')
	self._indent()

        for species in self.species:
            self._write('if (sn==%d) return "%s"; ' % (
                species.id, species.symbol))
 
        self._write(self.line( 'species name not found' ))
        self._write('return "NOTFOUND";')

	self._outdent()
	self._write('}')
	return
    
# Fuego's core routines section begins here
    def _atomicWeight(self, mechanism):

        self._write()
        self._write()
        self._write(self.line('save atomic weights into array'))
        self._write('void atomicWeight(double *  awt)')
        self._write('{')
        self._indent()
        import pyre
        periodic = pyre.handbook.periodicTable()
        for element in mechanism.element():
            aw = mechanism.element(element.symbol).weight
            if not aw:
                aw = periodic.symbol(element.symbol.capitalize()).atomicWeight

            self._write('awt[%d] = %f; ' % (
                element.id, aw) + self.line('%s' % element.symbol))

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        return 


    ##def _initialize_chemistry_GPU(self, mechanism):

    ##    self._write()
    ##    self._write('void initialize_chemistry_device(UserData user_data)')
    ##    self._write('{')
    ##    self._indent()
    ##    

    ##    nElement = len(mechanism.element())
    ##    nSpecies = len(mechanism.species())
    ##    nReactions = len(mechanism.reaction())
    ##    
    ##    # build reverse reaction map
    ##    rmap = {}
    ##    for i, reaction in zip(range(nReactions), mechanism.reaction()):
    ##        rmap[reaction.orig_id-1] = i

    ##    for j in range(nReactions):
    ##        reaction = mechanism.reaction()[rmap[j]]
    ##        id = reaction.id - 1

    ##        A, beta, E = reaction.arrhenius
    ##        self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
    ##        self._write("user_data->fwd_A[%d]     = %.17g;" % (id,A))
    ##        self._write("user_data->fwd_beta[%d]  = %.17g;" % (id,beta))
    ##        self._write("user_data->fwd_Ea[%d]    = %.17g;" % (id,E))

    ##        dim = self._phaseSpaceUnits(reaction.reactants)
    ##        thirdBody = reaction.thirdBody
    ##        low = reaction.low
    ##        if not thirdBody:
    ##            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB
    ##        elif not low:
    ##            uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
    ##        else:
    ##            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
    ##            low_A, low_beta, low_E = low
    ##            self._write("user_data->low_A[%d]     = %.17g;" % (id,low_A))
    ##            self._write("user_data->low_beta[%d]  = %.17g;" % (id,low_beta))
    ##            self._write("user_data->low_Ea[%d]    = %.17g;" % (id,low_E))
    ##            if reaction.troe:
    ##                troe = reaction.troe
    ##                ntroe = len(troe)
    ##                is_troe = True
    ##                self._write("user_data->troe_a[%d]    = %.17g;" % (id,troe[0]))
    ##                if ntroe>1:
    ##                    self._write("user_data->troe_Tsss[%d] = %.17g;" % (id,troe[1]))
    ##                if ntroe>2:
    ##                    self._write("user_data->troe_Ts[%d]   = %.17g;" % (id,troe[2]))
    ##                if ntroe>3:
    ##                    self._write("user_data->troe_Tss[%d]  = %.17g;" % (id,troe[3]))
    ##                self._write("user_data->troe_len[%d]  = %d;" % (id,ntroe))
    ##            if reaction.sri:
    ##                sri = reaction.sri
    ##                nsri = len(sri)
    ##                is_sri = True
    ##                self._write("user_data->sri_a[%d]     = %.17g;" % (id,sri[0]))
    ##                if nsri>1:
    ##                    self._write("user_data->sri_b[%d]     = %.17g;" % (id,sri[1]))
    ##                if nsri>2:
    ##                    self._write("user_data->sri_c[%d]     = %.17g;" % (id,sri[2]))
    ##                if nsri>3:
    ##                    self._write("user_data->sri_d[%d]     = %.17g;" % (id,sri[3]))
    ##                if nsri>4:
    ##                    self._write("user_data->sri_e[%d]     = %.17g;" % (id,sri[4]))
    ##                self._write("user_data->sri_len[%d]   = %d;" % (id,nsri))

    ##        self._write("user_data->prefactor_units[%d]  = %.17g;" % (id,uc.value))
    ##        aeuc = self._activationEnergyUnits(reaction.units["activation"])
    ##        self._write("user_data->activation_units[%d] = %.17g;" % (id,aeuc / Rc / kelvin))
    ##        self._write("user_data->phase_units[%d]      = 1e-%d;" % (id,dim*6))

    ##        if low:
    ##            self._write("user_data->is_PD[%d] = 1;" % (id) )
    ##        else:
    ##            self._write("user_data->is_PD[%d] = 0;" % (id) )


    ##        if thirdBody:
    ##            efficiencies = reaction.efficiencies
    ##            self._write("user_data->nTB[%d] = %d;" % (id, len(efficiencies)))

    ##            if (len(efficiencies) > 0):
    ##                self._write("cudaMallocManaged(&user_data->TB[%d], %d * sizeof(double));"% (id, len(efficiencies)))
    ##                self._write("cudaMallocManaged(&user_data->TBid[%d], %d * sizeof(int));"% (id, len(efficiencies)))
    ##                for i, eff in enumerate(efficiencies):
    ##                    symbol, efficiency = eff
    ##                    self._write("user_data->TBid[%d][%d] = %.17g; user_data->TB[%d][%d] = %.17g; // %s"
    ##                                % (id, i, mechanism.species(symbol).id, id, i, efficiency, symbol ))
    ##        else:
    ##            self._write("user_data->nTB[%d] = 0;" % (id))

    ##        self._write()

    ##    self._write('return;')
    ##    self._outdent()
    ##    self._write('}')
    ##    return 

    def _productionRate_GPU(self, mechanism):

        nElement = len(mechanism.element())
        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        # main function
        self._write()
        self._write('#ifdef AMREX_USE_CUDA')
        self._write(self.line('GPU version of productionRate: no more use of thermo namespace vectors'))
        self._write(self.line('compute the production rate for each species'))
        self._write('AMREX_GPU_HOST_DEVICE inline void  productionRate(double * wdot, double * sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        if (nReactions == 0):
            self._write()
        else:
            self._write()
            self._write('double qdot, q_f[%d], q_r[%d];' % (nReactions,nReactions))
            self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

        self._write()
        self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
        self._indent()
        self._write('wdot[i] = 0.0;')
        self._outdent()
        self._write('}')

        for i in range(nReactions):
            self._write()
            self._write("qdot = q_f[%d]-q_r[%d];" % (i,i))
            reaction = mechanism.reaction(id=i)
            agents = list(set(reaction.reactants + reaction.products))
            agents = sorted(agents, key=lambda x: mechanism.species(x[0]).id)
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a:
                        if coefficient == 1.0:
                            self._write("wdot[%d] -= qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] -= %f * qdot;" % (mechanism.species(symbol).id, coefficient))
                for b in reaction.products: 
                    if b == a:
                        if coefficient == 1.0:
                            self._write("wdot[%d] += qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] += %f * qdot;" % (mechanism.species(symbol).id, coefficient))


        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        # k_f function
        ##self._write()
        ##self._write('AMREX_GPU_HOST_DEVICE void comp_k_f(double *  tc, double invT, double *  k_f)')
        ##self._write('{')
        ##self._indent()
        ##for j in range(nReactions):
        ##    reaction   = mechanism.reaction(id=j)
        ##    dim        = self._phaseSpaceUnits(reaction.reactants)
        ##    aeuc       = self._activationEnergyUnits(reaction.units["activation"])
        ##    A, beta, E = reaction.arrhenius
        ##    thirdBody  = reaction.thirdBody
        ##    low        = reaction.low
        ##    if not thirdBody:
        ##        uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB
        ##    elif not low:
        ##        uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
        ##    else:
        ##        uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
        ##    self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
        ##    self._write("k_f[%d] = %.17g * %.17g" % (reaction.id,uc.value,A))
        ##    self._write("            * exp(%.17g * tc[0] - %.17g * %.17g * invT);" % (beta,aeuc / Rc / kelvin,E))
        ##    self._write()
        ##self._write('return;')
        ##self._outdent()
        ##self._write('}')

        # Kc
        ##self._write()
        ##self._write('AMREX_GPU_HOST_DEVICE inline void comp_Kc(double *  tc, double invT, double *  Kc)')
        ##self._write('{')
        ##self._indent()

        ##self._write(self.line('compute the Gibbs free energy'))
        ##self._write('double g_RT[%d];' % (nSpecies))
        ##self._write('gibbs(g_RT, tc);')

        ##self._write()

        ##for reaction in mechanism.reaction():
        ##    KcExpArg = self._sortedKcExpArg(mechanism, reaction)
        ##    self._write("Kc[%d] = %s;" % (reaction.id-1,KcExpArg))

        ##self._write()
        ##
        ##self._outdent()
        ##self._write('#ifdef __INTEL_COMPILER')
        ##self._indent()
        ##self._write(' #pragma simd')
        ##self._outdent()
        ##self._write('#endif')
        ##self._indent()
        ##self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        ##self._indent()
        ##self._write("Kc[i] = exp(Kc[i]);")
        ##self._outdent()
        ##self._write("};")

        ##self._write()

        ##self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        ##self._write('double refC = %g / %g * invT;' % (atm.value, R.value))
        ##self._write('double refCinv = 1 / refC;')

        ##self._write()

        ##for reaction in mechanism.reaction():
        ##    KcConv = self._KcConv(mechanism, reaction)
        ##    if KcConv:
        ##        self._write("Kc[%d] *= %s;" % (reaction.id-1,KcConv))        
        ##
        ##self._write()

        ##self._write('return;')
        ##self._outdent()
        ##self._write('}')


        # qdot
        self._write()
        self._write('AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)')
        self._write('{')
        self._indent()

        if (nReactions > 0):

            nclassd = nReactions - nspecial
            #nCorr   = n3body + ntroe + nsri + nlindemann

            for i in range(nReactions):
                self._write()
                reaction = mechanism.reaction(id=i)
                self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))
                if (len(reaction.ford) > 0):
                    self._write("qf[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.ford)))
                else:
                    self._write("qf[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.reactants)))
                if reaction.reversible:
                    self._write("qr[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.products)))
                else:
                    self._write("qr[%d] = 0.0;" % (i))

            self._write()

            # Mixt concentration for PD & TB
            self._write(self.line('compute the mixture concentration'))
            self._write('double mixture = 0.0;')
            self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
            self._indent()
            self._write('mixture += sc[i];')
            self._outdent()
            self._write('}')
            self._write()

            # Kc stuff
            self._write(self.line('compute the Gibbs free energy'))
            self._write('double g_RT[%d];' % (nSpecies))
            self._write('gibbs(g_RT, tc);')

            self._write()

            self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
            self._write('double refC = %g / %g * invT;' % (atm.value, R.value))
            self._write('double refCinv = 1 / refC;')

            self._write()
            
            # kfs
            self._write("/* Evaluate the kfs */")
            #self._write("double k_f[%d];"% nclassd)
            #self._write("double Corr[%d];" % nclassd)
            self._write("double k_f, k_r, Corr;")
            if ntroe > 0:
                self._write("double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;")
            if nsri > 0:
                self._write("double redP, F, X, F_sri;")
            self._write()

            # build reverse reaction map
            rmap = {}
            for i, reaction in zip(range(nReactions), mechanism.reaction()):
                rmap[reaction.orig_id-1] = i

            for i in range(nReactions):
                reaction = mechanism.reaction()[rmap[i]]
                idx = reaction.id - 1

                KcExpArg = self._sortedKcExpArg(mechanism, reaction)
                KcConv = self._KcConv(mechanism, reaction)

                A, beta, E = reaction.arrhenius
                dim = self._phaseSpaceUnits(reaction.reactants)
                thirdBody = reaction.thirdBody
                low = reaction.low
                if not thirdBody:
                    uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB
                elif not low:
                    uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
                else:
                    uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
                    low_A, low_beta, low_E = low
                    if reaction.troe:
                        troe = reaction.troe
                        ntroe = len(troe)
                        is_troe = True
                    if reaction.sri:
                        sri = reaction.sri
                        nsri = len(sri)
                        is_sri = True
                aeuc = self._activationEnergyUnits(reaction.units["activation"])

                self._write("// (%d):  %s" % (reaction.orig_id - 1, reaction.equation()))
                self._write("k_f = %.17g * %.17g " % (uc.value,A)) 
                self._write("           * exp(%.17g * tc[0] - %.17g * (%.17g) * invT);" % (beta, aeuc / Rc / kelvin, E))

                if not thirdBody:
                    self._write("Corr  = 1.0;")
                    self._write("qf[%d] *= Corr * k_f;" % idx)
                elif not low:
                    alpha = self._enhancement_d(mechanism, reaction)
                    self._write("Corr  = %s;" %(alpha))
                    self._write("qf[%d] *= Corr * k_f;" % idx)
                else:
                    alpha = self._enhancement_d(mechanism, reaction)
                    self._write("Corr  = %s;" %(alpha))
                    self._write("redP = Corr / k_f * 1e-%d * %.17g " % (dim*6, low_A)) 
                    self._write("           * exp(%.17g  * tc[0] - %.17g  * (%.17g) *invT);" % (low_beta, aeuc / Rc / kelvin, low_E))
                    if reaction.troe:
                        self._write("F = redP / (1.0 + redP);")
                        self._write("logPred = log10(redP);")
                        self._write('logFcent = log10(')
                        if (abs(troe[1]) > 1.e-100):
                            if(troe[0] < 0):
                                self._write('    (1.+%.17g)*exp(-tc[1] / %.17g) ' % (-troe[0],troe[1]))
                            else:
                                self._write('    (1.-%.17g)*exp(-tc[1] / %.17g) ' % (troe[0],troe[1]))
                        else:
                            self._write('     0. ' )
                        if (abs(troe[2]) > 1.e-100):
                            self._write('    + %.17g * exp(-tc[1]/%.17g)  ' % (troe[0],troe[2]))
                        else:
                            self._write('    + 0. ')
                        if (ntroe == 4):
                            if(troe[3] < 0):
                                self._write('    + exp(%.17g * invT));' % -troe[3])
                            else:
                                self._write('    + exp(-%.17g * invT));' % troe[3])
                        else:
                            self._write('    + 0.);' )
                        self._write("troe_c = -.4 - .67 * logFcent;")
                        self._write("troe_n = .75 - 1.27 * logFcent;")
                        self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
                        self._write("F_troe = pow(10., logFcent / (1.0 + troe*troe));")
                        self._write("Corr = F * F_troe;")
                        self._write("qf[%d] *= Corr * k_f;" % idx)
                    elif reaction.sri:
                        self._write("F = redP / (1.0 + redP);")
                        self._write("logPred = log10(redP);")
                        self._write("X = 1.0 / (1.0 + logPred*logPred);")
                        if (sri[1] < 0):
                            self._write("F_sri = exp(X * log(%.17g * exp(%.17g*invT)" % (sri[0],-sri[1]))
                        else:
                            self._write("F_sri = exp(X * log(%.17g * exp(-%.17g*invT)" % (sri[0],sri[1]))
                        if (sri[2] > 1.e-100):
                            self._write("   +  exp(tc[0]/%.17g) " % sri[2])
                        else:
                            self._write("   +  0. ") 
                        self._write("   *  (%d > 3 ? %.17g*exp(%.17g*tc[0]) : 1.0);" % (nsri,sri[3],sri[4]))
                        self._write("Corr = F * F_sri;")
                        self._write("qf[%d] *= Corr * k_f;" % idx)
                    elif (nlindemann > 0):
                        self._write("Corr = redP / (1. + redP);")
                        self._write("qf[%d] *= Corr * k_f;" % idx)

                if reaction.rev:
                    Ar, betar, Er = reaction.rev
                    dim_rev       = self._phaseSpaceUnits(reaction.products)
                    if not thirdBody:
                        uc_rev = self._prefactorUnits(reaction.units["prefactor"], 1-dim_rev)
                    elif not low:
                        uc_rev = self._prefactorUnits(reaction.units["prefactor"], -dim_rev)
                    else:
                        print("REV reaction cannot be PD")
                        sys.exit(1)
                    self._write("k_r = %.17g * %.17g " % (uc_rev.value,Ar)) 
                    self._write("           * exp(%.17g * tc[0] - %.17g * %.17g * invT);" % (betar, aeuc / Rc / kelvin, Er))
                    self._write("qr[%d] *= Corr * k_r;" % idx)
                else:
                    if KcConv:
                        self._write("qr[%d] *= Corr * k_f / (exp(%s) * %s);" % (idx,KcExpArg,KcConv))        
                    else:
                        self._write("qr[%d] *= Corr * k_f / exp(%s);" % (idx,KcExpArg))

            self._write()


        #for reaction in mechanism.reaction():
        #    idx = reaction.id-1
        #    KcExpArg = self._sortedKcExpArg(mechanism, reaction)
        #    KcConv = self._KcConv(mechanism, reaction)
        #    if KcConv:
        #        self._write("qr[%d] *= qf[%d] / (exp(%s) * %s);" % (idx,idx,KcExpArg,KcConv))        
        #    else:
        #        self._write("qr[%d] *= qf[%d] / exp(%s);" % (idx,idx,KcExpArg))

        #self._write()
        
        self._write()
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('#endif')
        self._write()

        return

    def _productionRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        # OMP stuff
        self._write()
        self._write("#ifndef AMREX_USE_CUDA")
        self._write('static double T_save = -1;')
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(T_save)')
        self._write('#endif')
        self._write()
        self._write('static double k_f_save[%d];' % nReactions)
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(k_f_save)')
        self._write('#endif')
        self._write()
        self._write('static double Kc_save[%d];' % nReactions)
        self._write('#ifdef _OPENMP')
        self._write('#pragma omp threadprivate(Kc_save)')
        self._write('#endif')
        self._write()

        # main function
        self._write()
        self._write(self.line('compute the production rate for each species pointwise on CPU'))
        self._write('void productionRate(double *  wdot, double *  sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        self._write()
        self._write('if (T != T_save)')
        self._write('{')
        self._indent()
        self._write('T_save = T;')
        self._write('comp_k_f(tc,invT,k_f_save);');
        self._write('comp_Kc(tc,invT,Kc_save);');
        self._outdent()
        self._write("}")

        self._write()
        self._write('double qdot, q_f[%d], q_r[%d];' % (nReactions,nReactions))
        self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

        self._write()
        self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
        self._indent()
        self._write('wdot[i] = 0.0;')
        self._outdent()
        self._write('}')

        for i in range(nReactions):
            self._write()
            self._write("qdot = q_f[%d]-q_r[%d];" % (i,i))
            reaction = mechanism.reaction(id=i)
            agents = list(set(reaction.reactants + reaction.products))
            agents = sorted(agents, key=lambda x: mechanism.species(x[0]).id)
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d] -= qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] -= %f * qdot;" % (mechanism.species(symbol).id, coefficient))
                for b in reaction.products: 
                    if b == a:
                        if coefficient == 1:
                            self._write("wdot[%d] += qdot;" % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d] += %f * qdot;" % (mechanism.species(symbol).id, coefficient))


        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        # k_f function
        self._write()
        self._write('void comp_k_f(double *  tc, double invT, double *  k_f)')
        self._write('{')
        self._indent()
        self._outdent()
        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write('#pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()
        self._write("k_f[i] = prefactor_units[i] * fwd_A[i]")
        self._write("            * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);")
        self._outdent()
        self._write("};")
        self._write('return;')
        self._outdent()
        self._write('}')

        # Kc
        self._write()
        self._write('void comp_Kc(double *  tc, double invT, double *  Kc)')
        self._write('{')
        self._indent()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        for reaction in mechanism.reaction():
            KcExpArg = self._sortedKcExpArg(mechanism, reaction)
            self._write("Kc[%d] = %s;" % (reaction.id-1,KcExpArg))

        self._write()
        
        self._outdent()
        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write(' #pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<%d; ++i) {' % (nReactions))
        self._indent()
        self._write("Kc[i] = exp(Kc[i]);")
        self._outdent()
        self._write("};")

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g * invT;' % (atm.value, R.value))
        self._write('double refCinv = 1 / refC;')

        self._write()

        for reaction in mechanism.reaction():
            KcConv = self._KcConv(mechanism, reaction)
            if KcConv:
                self._write("Kc[%d] *= %s;" % (reaction.id-1,KcConv))        
        
        self._write()

        self._write('return;')
        self._outdent()
        self._write('}')

        # qdot
        self._write()
        self._write('void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)')
        self._write('{')
        self._indent()

        nclassd = nReactions - nspecial
        nCorr   = n3body + ntroe + nsri + nlindemann

        for i in range(nclassd):
            self._write()
            reaction = mechanism.reaction(id=i)
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))
            if (len(reaction.ford) > 0):
                self._write("qf[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.ford)))
            else:
                self._write("qf[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.reactants)))
            if reaction.reversible:
                self._write("qr[%d] = %s;" % (i, self._sortedPhaseSpace(mechanism, reaction.products)))
            else:
                self._write("qr[%d] = 0.0;" % (i))

        self._write()
        self._write('double T = tc[1];')
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int i = 0; i < %d; ++i) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write("double Corr[%d];" % nclassd)
        self._write('for (int i = 0; i < %d; ++i) {' % nclassd)
        self._indent()
        self._write('Corr[i] = 1.0;')
        self._outdent()
        self._write('}')

        if ntroe > 0:
            self._write()
            self._write(self.line(" troe"))
            self._write("{")
            self._indent()
            self._write("double alpha[%d];" % ntroe)
            alpha_d = {}
            for i in range(itroe[0],itroe[1]):
                ii = i - itroe[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha in alpha_d:
                        self._write("alpha[%d] = %s;" %(ii,alpha_d[alpha]))
                    else:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                        alpha_d[alpha] = "alpha[%d]" % ii

            if ntroe >= 4:
                self._outdent()
                self._outdent()
                # self._write('#ifdef __INTEL_COMPILER')
                # self._indent()
                # self._indent()
                # self._write(' #pragma simd')
                # self._outdent()
                # self._outdent()
                # self._write('#endif')
                self._indent()
                self._indent()
            self._write("for (int i=%d; i<%d; i++)" %(itroe[0],itroe[1]))
            self._write("{")
            self._indent()
            self._write("double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;")
            self._write("redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);" % itroe[0])
            self._write("F = redP / (1.0 + redP);")
            self._write("logPred = log10(redP);")
            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-T/troe_Tsss[i]) : 0.) ')
            self._write('    + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T/troe_Ts[i]) : 0.) ')
            self._write('    + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );')
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10., logFcent / (1.0 + troe*troe));")
            self._write("Corr[i] = F * F_troe;")
            self._outdent()
            self._write('}')

            self._outdent()
            self._write("}")

        if nsri > 0:
            self._write()
            self._write(self.line(" SRI"))
            self._write("{")
            self._indent()
            self._write("double alpha[%d];" % nsri)
            self._write("double redP, F, X, F_sri;")
            alpha_d = {}
            for i in range(isri[0],isri[1]):
                ii = i - isri[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha in alpha_d:
                        self._write("alpha[%d] = %s;" %(ii,alpha_d[alpha]))
                    else:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                        alpha_d[alpha] = "alpha[%d]" % ii

            if nsri >= 4:
                self._outdent()
                self._outdent()
                # self._write('#ifdef __INTEL_COMPILER')
                # self._indent()
                # self._indent()
                # self._write(' #pragma simd')
                # self._outdent()
                # self._outdent()
                # self._write('#endif')
                self._indent()
                self._indent()
            self._write("for (int i=%d; i<%d; i++)" %(isri[0],isri[1]))
            self._write("{")
            self._indent()
            self._write("redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);" % itroe[0])
            self._write("F = redP / (1.0 + redP);")
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[i] * exp(-sri_b[i]*invT)")
            self._write("   +  (sri_c[i] > 1.e-100 ? exp(T/sri_c[i]) : 0.0) )")
            self._write("   *  (sri_len[i] > 3 ? sri_d[i]*exp(sri_e[i]*tc[0]) : 1.0);")
            self._write("Corr[i] = F * F_sri;")
            self._outdent()
            self._write('}')

            self._outdent()
            self._write("}")

        if nlindemann > 0:
            self._write()
            self._write(self.line(" Lindemann"))
            self._write("{")
            self._indent()
            if nlindemann > 1:
                self._write("double alpha[%d];" % nlindemann)
            else:
                self._write("double alpha;")

            for i in range(ilindemann[0],ilindemann[1]):
                ii = i - ilindemann[0]
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if nlindemann > 1:
                        self._write("alpha[%d] = %s;" %(ii,alpha))
                    else:
                        self._write("alpha = %s;" %(alpha))

            if nlindemann == 1:
                self._write("double redP = alpha / k_f_save[%d] * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] * invT);" 
                            % (ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0],ilindemann[0]))
                self._write("Corr[%d] = redP / (1. + redP);" % ilindemann[0])
            else:
                if nlindemann >= 4:
                    self._outdent()
                    # self._write('#ifdef __INTEL_COMPILER')
                    # self._indent()
                    # self._write(' #pragma simd')
                    # self._outdent()
                    # self._write('#endif')
                    self._indent()
                self._write("for (int i=%d; i<%d; i++)" % (ilindemann[0], ilindemann[1]))
                self._write("{")
                self._indent()
                self._write("double redP = alpha[i-%d] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] * invT);"
                            % ilindemann[0])
                self._write("Corr[i] = redP / (1. + redP);")
                self._outdent()
                self._write('}')

            self._outdent()
            self._write("}")

        if n3body > 0:
            self._write()
            self._write(self.line(" simple three-body correction"))
            self._write("{")
            self._indent()
            self._write("double alpha;")
            alpha_save = ""
            for i in range(i3body[0],i3body[1]):
                reaction = mechanism.reaction(id=i)
                if reaction.thirdBody:
                    alpha = self._enhancement(mechanism, reaction)
                    if alpha != alpha_save:
                        alpha_save = alpha
                        self._write("alpha = %s;" % alpha)
                    self._write("Corr[%d] = alpha;" % i)
            self._outdent()
            self._write("}")

        self._write()
        self._write("for (int i=0; i<%d; i++)" % nclassd)
        self._write("{")
        self._indent()
        self._write("qf[i] *= Corr[i] * k_f_save[i];")
        self._write("qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];")
        self._outdent()
        self._write("}")
        
        if nspecial > 0:

            print "\n\n ***** WARNING: %d unclassified reactions\n" % nspecial

            self._write()
            self._write(self.line('unclassified reactions'))
            self._write('{')
            self._indent()

            self._write(self.line("reactions: %d to %d" % (ispecial[0]+1,ispecial[1])))

            #self._write('double Kc;                      ' + self.line('equilibrium constant'))
            self._write('double k_f;                     ' + self.line('forward reaction rate'))
            self._write('double k_r;                     ' + self.line('reverse reaction rate'))
            self._write('double q_f;                     ' + self.line('forward progress rate'))
            self._write('double q_r;                     ' + self.line('reverse progress rate'))
            self._write('double phi_f;                   '
                        + self.line('forward phase space factor'))
            self._write('double phi_r;                   ' + self.line('reverse phase space factor'))
            self._write('double alpha;                   ' + self.line('enhancement'))

            #self._write('double redP;                    ' + self.line('reduced pressure'))
            #self._write('double logPred;                 ' + self.line('log of above'))
            #self._write('double F;                       ' + self.line('fallof rate enhancement'))
            #self._write()
            #self._write('double F_troe;                  ' + self.line('TROE intermediate'))
            #self._write('double logFcent;                ' + self.line('TROE intermediate'))
            #self._write('double troe;                    ' + self.line('TROE intermediate'))
            #self._write('double troe_c;                  ' + self.line('TROE intermediate'))
            #self._write('double troe_n;                  ' + self.line('TROE intermediate'))

            for i in range(ispecial[0],ispecial[1]):
                self._write()
                reaction = mechanism.reaction(id=i)
                self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

                # compute the rates
                self._forwardRate(mechanism, reaction)
                self._reverseRate(mechanism, reaction)

                # store the progress rate
                self._write("qf[%d] = q_f;" % i)
                self._write("qr[%d] = q_r;" % i)

            self._outdent()
            self._write('}')

        self._write()
        self._write('return;')
        self._outdent()
        self._write('}')

        self._write("#endif")

        return

    def _DproductionRatePYJAC(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        self._write('#ifdef USE_PYJAC')
        self._write()
        self._write(self.line('compute the reaction Jacobian using PyJac'))
        self._write('void DWDOT_PYJAC(double *  J, double *  y, double *  Tp, double *  Press)')
        self._write('{')
        self._indent()

        self._write('double y_pyjac[%d];' % (nSpecies + 1))
        self._write('double J_reorg[%d];' % (nSpecies+1)**2)

        self._write()
        self._write(self.line(' INPUT Y'))
        self._write('y_pyjac[0] = *Tp;')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('y_pyjac[1+k] = y[k];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('double Press_MKS = *Press / 10.0;')

        self._write()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1)**2)
        self._indent()
        self._write('J[k] = 0.0;')
        self._write('J_reorg[k] = 0.0;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('eval_jacob(0, Press_MKS, y_pyjac, J);')

        self._write()
        self._write('/* Reorganization */')
        self._write('for (int k=0; k<%d; k++) {' %(nSpecies + 1))
        self._indent()
        self._write('J_reorg[k*%d + %d] = J[k*%d + 0];' % (nSpecies+1, nSpecies, nSpecies+1))
        self._write('for (int i=0; i<%d; i++) {' %(nSpecies))
        self._indent()
        self._write('J_reorg[k*%d + i] = J[k*%d + (i + 1)];' % (nSpecies+1, nSpecies+1))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write('for (int k=0; k<%d; k++) {' %(nSpecies + 1))
        self._indent()
        self._write('J[%d*%d + k] = J_reorg[k];' % (nSpecies, nSpecies+1))
        self._outdent()
        self._write('}')

        self._write('for (int k=0; k<%d; k++) {' %((nSpecies+1)*nSpecies))
        self._indent()
        self._write('J[k] = J_reorg[k + %d];' % (nSpecies + 1))
        self._outdent()
        self._write('}')

        #self._write()
        #self._write('/* Go back to CGS */')
        #self._write('for (int k=0; k<%d; k++) {' %(nSpecies))
        #self._indent()
        #self._write('/* dwdot[k]/dT */')
        #self._write('J[%d*%d + k] *= 1.e-3;' % (nSpecies, nSpecies+1))
        #self._write('/* dTdot/d[X] */')
        #self._write('J[k*%d + %d] *= 1.e3;' % (nSpecies+1, nSpecies))
        #self._outdent()
        #self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write('#endif')

        return

    def _DproductionRate(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        self._write()
        self._write(self.line('compute the reaction Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.e6 * sc[k];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, *Tp, *consP);')

        self._write()
        self._write('/* dwdot[k]/dT */')
        self._write('/* dTdot/d[X] */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('J[%d+k] *= 1.e-6;' % (nSpecies*(nSpecies+1)))
        self._write('J[k*%d+%d] *= 1.e6;' % (nSpecies+1, nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return

    def _DproductionRatePrecond(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        self._write()
        self._write(self.line('compute an approx to the reaction Jacobian (for preconditioning)'))
        self._write('AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(double *  J, double *  sc, double *  Tp, int * HP)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.e6 * sc[k];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian_precond(J, c, *Tp, *HP);')

        self._write()
        self._write('/* dwdot[k]/dT */')
        self._write('/* dTdot/d[X] */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('J[%d+k] *= 1.e-6;' % (nSpecies*(nSpecies+1)))
        self._write('J[k*%d+%d] *= 1.e6;' % (nSpecies+1, nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return

    def _DproductionRateSPSPrecond(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        self._write()
        self._write(self.line('compute an approx to the SPS Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(double *  Jsps, int * indx, int * len, double * sc, double * Tp, int * HP, double * gamma)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % ((nSpecies+1) * (nSpecies+1)))
        self._write('double mwt[%d];' % (nSpecies))
        self._write()
        self._write('get_mw(mwt);')
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.e6 * sc[k];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian_precond(J, c, *Tp, *HP);')

        self._write()
        self._write('/* Change of coord */')
        self._write('/* dwdot[k]/dT */')
        self._write('/* dTdot/d[X] */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('J[%d+k] = 1.e-6 * J[%d+k] * mwt[k];' % (nSpecies*(nSpecies+1), nSpecies*(nSpecies+1)))
        self._write('J[k*%d+%d] = 1.e6 * J[k*%d+%d] / mwt[k];' % (nSpecies+1, nSpecies, nSpecies+1, nSpecies))
        self._outdent()
        self._write('}')

        self._write('/* dTdot/dT */')
        self._write('/* dwdot[l]/[k] */')
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent() 
        self._write('for (int l=0; l<%d; l++) {' % nSpecies)
        self._indent()
        self._write('/* DIAG elem */')
        self._write('if (k == l){')
        self._indent()
        self._write('J[ %d * k + l] =  J[ %d * k + l] * mwt[l] / mwt[k];' % (nSpecies+1, nSpecies+1))
        self._outdent()
        self._write('/* NOT DIAG and not last column nor last row */')
        self._write('} else {')
        self._indent()
        self._write('J[ %d * k + l] =  J[ %d * k + l] * mwt[l] / mwt[k];' % (nSpecies+1, nSpecies+1))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write()

        self._write('for (int k=0; k<(*len); k++) {')
        self._indent()
        self._write('Jsps[k] = J[indx[k]];')
        self._outdent() 
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return

    def _sparsity(self, mechanism):

        species_list = [x.symbol for x in mechanism.species()]
        nSpecies = len(species_list)

        ####
        self._write()
        self._write(self.line('compute the sparsity pattern of the chemistry Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, 1500.0, *consP);')

        self._write()
        self._write('int nJdata_tmp = 0;')
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('if(J[ %d * k + l] != 0.0){' % (nSpecies+1))
        self._indent()
        self._write('nJdata_tmp = nJdata_tmp + 1;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*nJdata = NCELLS * nJdata_tmp;')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()
        self._write()

        ####
        self._write()
        self._write(self.line('compute the sparsity pattern of the system Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, 1500.0, *consP);')

        self._write()
        self._write('int nJdata_tmp = 0;')
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('if(k == l){')
        self._indent()
        self._write('nJdata_tmp = nJdata_tmp + 1;')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[ %d * k + l] != 0.0){' % (nSpecies+1))
        self._indent()
        self._write('nJdata_tmp = nJdata_tmp + 1;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('*nJdata = NCELLS * nJdata_tmp;')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()
        self._write()

        ####
        self._write()
        self._write(self.line('compute the sparsity pattern of the simplified (for preconditioning) system Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian_precond(J, c, 1500.0, *consP);')

        self._write()
        self._write('int nJdata_tmp = 0;')
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('if(k == l){')
        self._indent()
        self._write('nJdata_tmp = nJdata_tmp + 1;')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[ %d * k + l] != 0.0){' % (nSpecies+1))
        self._indent()
        self._write('nJdata_tmp = nJdata_tmp + 1;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('nJdata[0] = nJdata_tmp;')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()
        self._write()


        ####
        self._write(self.line('compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write('int offset_row;')
        self._write('int offset_col;')
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, 1500.0, *consP);')

        self._write()
        self._write('colPtrs[0] = 0;')
        self._write('int nJdata_tmp = 0;')
        self._write('for (int nc=0; nc<NCELLS; nc++) {')
        self._indent()
        self._write('offset_row = nc * %d;' % (nSpecies+1))
        self._write('offset_col = nc * %d;' % (nSpecies+1))
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('rowVals[nJdata_tmp] = l + offset_row; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('colPtrs[offset_col + (k + 1)] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        ####
        self._write(self.line('compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write('int offset;')
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, 1500.0, *consP);')

        self._write()
        self._write('if (base == 1) {')
        self._indent()
        self._write('rowPtrs[0] = 1;')
        self._write('int nJdata_tmp = 1;')
        self._write('for (int nc=0; nc<NCELLS; nc++) {')
        self._indent()
        self._write('offset = nc * %d;' % (nSpecies+1))
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp-1] = k+1 + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtrs[offset + (l + 1)] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
         
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('rowPtrs[0] = 0;')
        self._write('int nJdata_tmp = 0;')
        self._write('for (int nc=0; nc<NCELLS; nc++) {')
        self._indent()
        self._write('offset = nc * %d;' % (nSpecies+1))
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp] = k + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtrs[offset + (l + 1)] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        ####
        self._write(self.line('compute the sparsity pattern of the system Jacobian'))
        self._write(self.line('CSR format BASE is user choice'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write('int offset;')

        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian(J, c, 1500.0, *consP);')

        self._write()
        self._write('if (base == 1) {')
        self._indent()
        self._write('rowPtr[0] = 1;')
        self._write('int nJdata_tmp = 1;')
        self._write('for (int nc=0; nc<NCELLS; nc++) {')
        self._indent()
        self._write('offset = nc * %d;' % (nSpecies+1));
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if (k == l) {')
        self._indent()
        self._write('colVals[nJdata_tmp-1] = l+1 + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp-1] = k+1 + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtr[offset + (l + 1)] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('rowPtr[0] = 0;')
        self._write('int nJdata_tmp = 0;')
        self._write('for (int nc=0; nc<NCELLS; nc++) {')
        self._indent()
        self._write('offset = nc * %d;' % (nSpecies+1));
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if (k == l) {')
        self._indent()
        self._write('colVals[nJdata_tmp] = l + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp] = k + offset; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtr[offset + (l + 1)] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        ####
        self._write(self.line('compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU'))
        self._write(self.line('BASE 0'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian_precond(J, c, 1500.0, *consP);')

        self._write()
        self._write('colPtrs[0] = 0;')
        self._write('int nJdata_tmp = 0;')
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('if (k == l) {')
        self._indent()
        self._write('rowVals[nJdata_tmp] = l; ')
        self._write('indx[nJdata_tmp] = %d*k + l;' % (nSpecies+1))
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('rowVals[nJdata_tmp] = l; ')
        self._write('indx[nJdata_tmp] = %d*k + l;' % (nSpecies+1))
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('colPtrs[k+1] = nJdata_tmp;')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        ####
        self._write(self.line('compute the sparsity pattern of the simplified (for precond) system Jacobian'))
        self._write(self.line('CSR format BASE is under choice'))
        self._write('AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)')
        self._write('{')
        self._indent()

        self._write('double c[%d];' % (nSpecies))
        self._write('double J[%d];' % (nSpecies+1)**2)
        self._write()
        self._write('for (int k=0; k<%d; k++) {' % nSpecies)
        self._indent()
        self._write('c[k] = 1.0/ %f ;' % (nSpecies))
        self._outdent()
        self._write('}')

        self._write()
        self._write('aJacobian_precond(J, c, 1500.0, *consP);')

        self._write()
        self._write('if (base == 1) {')
        self._indent()
        self._write('rowPtr[0] = 1;')
        self._write('int nJdata_tmp = 1;')
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if (k == l) {')
        self._indent()
        self._write('colVals[nJdata_tmp-1] = l+1; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp-1] = k+1; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtr[l+1] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('rowPtr[0] = 0;')
        self._write('int nJdata_tmp = 0;')
        self._write('for (int l=0; l<%d; l++) {' % (nSpecies+1))
        self._indent()
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies+1))
        self._indent()
        self._write('if (k == l) {')
        self._indent()
        self._write('colVals[nJdata_tmp] = l; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('} else {')
        self._indent()
        self._write('if(J[%d*k + l] != 0.0) {' % (nSpecies+1))
        self._indent()
        self._write('colVals[nJdata_tmp] = k; ')
        self._write('nJdata_tmp = nJdata_tmp + 1; ')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')
        self._write('rowPtr[l+1] = nJdata_tmp;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')
        self._write()

        return


    def _ajacPrecond(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write(self.line('compute an approx to the reaction Jacobian'))
        self._write('AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)')
        self._write('{')
        self._indent()

        self._write('for (int i=0; i<%d; i++) {' % (nSpecies+1)**2)
        self._indent()
        self._write('J[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double wdot[%d];' % (nSpecies))
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies))
        self._indent()
        self._write('wdot[k] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        self._write('double invT2 = invT * invT;')

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1.0 / refC;')

        self._write()

        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[k];')
        self._outdent()
        self._write('}')

        self._write()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        self._write(self.line('compute the species enthalpy'))
        self._write('double h_RT[%d];' % (nSpecies))
        self._write('speciesEnthalpy(h_RT, tc);')

        self._write()

        self._write('double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;') 
        self._write('double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;')
        self._write('double dqdci, dcdc_fac, dqdc[%d];' % (nSpecies))
        self._write('double Pr, fPr, F, k_0, logPr;') 
        self._write('double logFcent, troe_c, troe_n, troePr_den, troePr, troe;')
        self._write('double Fcent1, Fcent2, Fcent3, Fcent;')
        self._write('double dlogFdc, dlogFdn, dlogFdcn_fac;')
        self._write('double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;')
        self._write('const double ln10 = log(10.0);')
        self._write('const double log10e = 1.0/log(10.0);')

        for i, reaction in zip(range(nReactions), mechanism.reaction()):

            lt = reaction.lt
            if lt:
                print "Landau-Teller reactions are not supported"
                sys.exit(1)

            self._write(self.line('reaction %d: %s' % (i+1, reaction.equation())))
            if reaction.low:  # case 1
                self._write(self.line('a pressure-fall-off reaction'))
                self._ajac_reaction_precond(mechanism, reaction, 1)
            elif reaction.thirdBody:  # case 2
                self._write(self.line('a third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction_precond(mechanism, reaction, 2)
            else:  # case 3
                self._write(self.line('a non-third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction_precond(mechanism, reaction, 3)
            self._write()

        self._write('double c_R[%d], dcRdT[%d], e_RT[%d];' % (nSpecies, nSpecies, nSpecies))
        self._write('double * eh_RT;')
        self._write('if (HP) {')
        self._indent()

        self._write('cp_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('eh_RT = &h_RT[0];');

        self._outdent()
        self._write('}')
        self._write('else {')
        self._indent()

        self._write('cv_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('speciesInternalEnergy(e_RT, tc);');
        self._write('eh_RT = &e_RT[0];');

        self._outdent()
        self._write('}')

        self._write()

        self._write('double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('cmix += c_R[k]*sc[k];')
        self._write('dcmixdT += dcRdT[k]*sc[k];')
        self._write('ehmix += eh_RT[k]*wdot[k];')
        self._write('dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];' % \
                        (nSpecies*(nSpecies+1)))
        self._outdent()
        self._write('}')

        self._write()
        self._write('double cmixinv = 1.0/cmix;')
        self._write('double tmp1 = ehmix*cmixinv;')
        self._write('double tmp3 = cmixinv*T;')
        self._write('double tmp2 = tmp1*tmp3;')
        self._write('double dehmixdc;')

        self._write('/* dTdot/d[X] */')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('dehmixdc = 0.0;')
        self._write('for (int m = 0; m < %d; ++m) {' % nSpecies)
        self._indent()
        self._write('dehmixdc += eh_RT[m]*J[k*%s+m];' % (nSpecies+1))
        self._outdent()
        self._write('}')        
        self._write('J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;' % (nSpecies+1,nSpecies))
        self._outdent()
        self._write('}')

        self._write('/* dTdot/dT */')
        self._write('J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;' % \
                        (nSpecies*(nSpecies+1)+nSpecies))

        self._outdent()
        self._write('}')
        return


    def _ajac_reaction_precond(self, mechanism, reaction, rcase):

        if rcase == 1: # pressure-dependent reaction
            isPD = True
            if reaction.thirdBody:
                has_alpha = True
                self._write('/* also 3-body */')
            else:
                has_alpha = False
                self._write('/* non 3-body */')
                print 'FIXME: pressure dependent non-3-body reaction in _ajac_reaction'
                sys.exit(1)
        elif rcase == 2: # third-body and non-pressure-dependent reaction
            isPD = False
            has_alpha = True
        elif rcase == 3: # simple non-third and non-pressure-dependent reaction
            isPD = False
            has_alpha = False
        else:
            print '_ajac_reaction: wrong case ', rcase
            exit(1)

        nSpecies = len(mechanism.species())
        rea_dict = {}
        pro_dict = {}
        all_dict = {}
        sumNuk = 0
        for symbol, coefficient in reaction.reactants:
            k = mechanism.species(symbol).id
            sumNuk -= coefficient
            if k in rea_dict:
                coe_old = rea_dict[k][1]
                rea_dict[k] = (symbol, coefficient+coe_old)
            else:
                rea_dict[k] = (symbol,  coefficient)
        for symbol, coefficient in reaction.products:
            k = mechanism.species(symbol).id
            sumNuk += coefficient
            if k in pro_dict:
                coe_old = pro_dict[k][1]
                pro_dict[k] = (symbol, coefficient+coe_old)
            else:
                pro_dict[k] = (symbol, coefficient)
        for k in range(nSpecies):
            if k in rea_dict and k in pro_dict:
                sr, nur = rea_dict[k]
                sp, nup = pro_dict[k]
                all_dict[k] = (sr, nup-nur)
            elif k in rea_dict:
                sr, nur = rea_dict[k]
                all_dict[k] = (sr, -nur)
            elif k in pro_dict:
                sp, nup = pro_dict[k]
                all_dict[k] = (sp, nup)

        sorted_reactants = sorted(rea_dict.values())
        sorted_products = sorted(pro_dict.values()) 

        if not reaction.reversible:
            if isPD or has_alpha:
                print 'FIXME: irreversible reaction in _ajac_reaction may not work'
                self._write('/* FIXME: irreversible reaction in _ajac_reaction may not work*/')
            for k in range(nSpecies):
                if k in sorted_reactants and k in sorted_products:
                    print 'FIXME: irreversible reaction in _ajac_reaction may not work'
                    self._write('/* FIXME: irreversible reaction in _ajac_reaction may not work*/')

        dim = self._phaseSpaceUnits(reaction.reactants)
        if isPD:
            Corr_s = 'Corr *'
            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
        elif has_alpha:
            Corr_s = 'alpha * '
            uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
        else:
            Corr_s = ''
            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB 
        aeuc = self._activationEnergyUnits(reaction.units["activation"])

        if has_alpha:
            self._write("/* 3-body correction factor */")
            self._write("alpha = %s;" % self._enhancement_d(mechanism, reaction))

        # forward
        A, beta, E = reaction.arrhenius
        self._write('/* forward */')
        self._write("phi_f = %s;" % self._sortedPhaseSpace(mechanism, sorted_reactants))
        #
        self._write("k_f = %.17g * %.17g" % (uc.value,A))
        self._write("            * exp(%.17g * tc[0] - %.17g * (%.17g) * invT);"
                    %(beta,aeuc / Rc / kelvin,E))
        #
        self._write("dlnkfdT = %.17g * invT + %.17g *  (%.17g)  * invT2;"
                    %(beta,aeuc / Rc / kelvin,E))

        if isPD:
            low_A, low_beta, low_E = reaction.low
            self._write('/* pressure-fall-off */')
            self._write("k_0 = %.17g * exp(%.17g * tc[0] - %.17g * (%.17g) * invT);"
                        %(low_A,low_beta,aeuc / Rc / kelvin,low_E))
            self._write('Pr = 1e-%d * alpha / k_f * k_0;' % (dim*6))
            self._write('fPr = Pr / (1.0+Pr);')
            self._write("dlnk0dT = %.17g * invT + %.17g * (%.17g) * invT2;"
                        %(low_beta,aeuc / Rc / kelvin,low_E))
            self._write('dlogPrdT = log10e*(dlnk0dT - dlnkfdT);')
            self._write('dlogfPrdT = dlogPrdT / (1.0+Pr);')
            #
            if reaction.sri:
                self._write('/* SRI form */')
                print "FIXME: sri not supported in _ajac_reaction yet"
                sys.exit(1)
            elif reaction.troe:
                self._write('/* Troe form */')
                troe = reaction.troe
                ntroe = len(troe)
                self._write("logPr = log10(Pr);")
                if (abs(troe[1]) > 1.e-100):
                    if (troe[0] < 0):
                        self._write('Fcent1 = (1.+%.17g)*exp(-T/%.17g);'
                                %(-troe[0],troe[1]))
                    else:
                        self._write('Fcent1 = (1.-%.17g)*exp(-T/%.17g);'
                                %(troe[0],troe[1]))
                else:
                    self._write('Fcent1 = 0.;')
                if (abs(troe[2]) > 1.e-100):
                    self._write('Fcent2 = %.17g * exp(-T/%.17g);'
                                %(troe[0],troe[2]))
                else:
                    self._write('Fcent2 = 0.;')
                if (ntroe == 4):
                    if (troe[3] < 0):
                        self._write('Fcent3 = exp(%.17g * invT);'
                                % -troe[3] )
                    else:
                        self._write('Fcent3 = exp(-%.17g * invT);'
                                % troe[3] )
                else:
                    self._write('Fcent3 = 0.;')
                self._write('Fcent = Fcent1 + Fcent2 + Fcent3;')
                self._write("logFcent = log10(Fcent);")
                self._write("troe_c = -.4 - .67 * logFcent;")
                self._write("troe_n = .75 - 1.27 * logFcent;")
                self._write("troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));")
                self._write("troePr = (troe_c + logPr) * troePr_den;")
                self._write("troe = 1.0 / (1.0 + troePr*troePr);")
                self._write("F = pow(10.0, logFcent * troe);")

                self._write("dlogFcentdT = log10e/Fcent*( ")
                if (abs(troe[1]) > 1.e-100):
                    self._write("    -Fcent1/%.17g"
                                % troe[1] )
                else:
                    self._write("    +0.")
                if (abs(troe[2]) > 1.e-100):
                    self._write("    -Fcent2/%.17g"
                                % troe[2] )
                else:
                    self._write("    +0.")
                if (ntroe == 4):
                    self._write("    + Fcent3*%.17g*invT2);"
                                % troe[3] )
                else:
                    self._write("    + 0.);")

                self._write("dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;")
                self._write('dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;')
                self._write('dlogFdn = dlogFdcn_fac * troePr;')
                self._write('dlogFdlogPr = dlogFdc;')
                self._write('dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;')
            else:
                self._write('/* Lindemann form */')
                self._write('F = 1.0;')
                self._write('dlogFdlogPr = 0.0;')
                self._write('dlogFdT = 0.0;')

        # reverse
        if not reaction.reversible:
            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f;')
            else:
                self._write('q_nocor = k_f*phi_f;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %sdlnkfdT*k_f*phi_f;' % Corr_s)
        else:
            self._write('/* reverse */')
            self._write("phi_r = %s;" % self._sortedPhaseSpace(mechanism, sorted_products))
            self._write('Kc = %s;' % self._sortedKc(mechanism, reaction))
            self._write('k_r = k_f / Kc;')

            dlnKcdT_s = 'invT * ('
            terms = []
            for symbol, coefficient in sorted(sorted_reactants,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += '-(' + ' + '.join(terms) + ')'
            terms = []
            for symbol, coefficient in sorted(sorted_products,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += ' + (' + ' + '.join(terms) + ')'
            if sumNuk > 0:
                dlnKcdT_s += ' - %f' % sumNuk
            elif sumNuk < 0:
                dlnKcdT_s += ' + %f' % (-sumNuk)
            dlnKcdT_s += ')'
            self._write('dlnKcdT = %s;' % dlnKcdT_s)

            self._write('dkrdT = (dlnkfdT - dlnKcdT)*k_r;')

            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f - k_r*phi_r;')
            else:
                self._write('q_nocor = k_f*phi_f - k_r*phi_r;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);' % Corr_s)

        self._write("/* update wdot */")
        for k in sorted(all_dict.keys()):
            s, nu = all_dict[k]
            if nu == 1:
                self._write('wdot[%d] += q; /* %s */' % (k, s))
            elif nu == -1:
                self._write('wdot[%d] -= q; /* %s */' % (k, s))
            elif nu > 0:
                self._write('wdot[%d] += %.17g * q; /* %s */' % (k, nu, s))
            elif nu < 0:
                self._write('wdot[%d] -= %.17g * q; /* %s */' % (k, -nu, s))

        if isPD:
            self._write('/* for convenience */')
            self._write('k_f *= Corr;')
            if reaction.reversible:
                self._write('k_r *= Corr;')
        elif has_alpha:
            self._write('/* for convenience */')
            self._write('k_f *= alpha;')
            if reaction.reversible:
                self._write('k_r *= alpha;')
            else:
                self._write('k_r = 0.0;')

        if isPD:
            self._write('dcdc_fac = 0.0;')
        #elif has_alpha:
        #    self._write('dcdc_fac = q_nocor;')

        def dqdc_simple_precond(dqdc_s, k):
            if dqdc_s ==  "0":
                dqdc_s = ''
            if k in sorted(rea_dict.keys()):
                dps = self._DphaseSpace(mechanism,sorted_reactants,rea_dict[k][0])
                if dps == "1.0":
                    dps_s = ''
                else:
                    dps_s = '*'+dps
                dqdc_s += ' + k_f%s' % dps_s
            if reaction.reversible:
                if k in sorted(pro_dict.keys()):
                    dps = self._DphaseSpace(mechanism,sorted_products,pro_dict[k][0])
                    if dps == "1.0":
                        dps_s = ''
                    else:
                        dps_s = '*'+dps
                    dqdc_s += ' - k_r%s' % dps_s
            return dqdc_s

        if has_alpha or isPD:

            #self._write('if (consP) {')
            #self._indent()

            #for k in range(nSpecies):
            #    dqdc_s = self._Denhancement(mechanism,reaction,k,True)
            #    if dqdc_s != "0":
            #        if isPD:
            #            if dqdc_s == "1":
            #                dqdc_s ='dcdc_fac'
            #            else:
            #                dqdc_s +='*dcdc_fac'
            #        elif has_alpha:
            #            if dqdc_s == "1":
            #                dqdc_s ='q_nocor'
            #            else:
            #                dqdc_s +='*q_nocor'

            #    dqdc_s = dqdc_simple_precond(dqdc_s,k)
            #    if dqdc_s:
            #        symb_k = self.species[k].symbol
            #        self._write('/* d()/d[%s] */' % symb_k)
            #        self._write('dqdci = %s;' % (dqdc_s))
            #        #
            #        for m in sorted(all_dict.keys()):
            #            if all_dict[m][1] != 0:
            #                s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
            #                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
            #                s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], symb_k)
            #                self._write(s1.ljust(30) + s2)

            #self._outdent()
            #self._write('}')
            #self._write('else {')
            #self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement_d(mechanism,reaction,k,False)
                if dqdc_s != '0':
                    if isPD:
                        if dqdc_s == '1':
                            dqdc_s ='dcdc_fac'
                        else:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == '1':
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple_precond(dqdc_s,k)
                if dqdc_s:
                    self._write('dqdc[%d] = %s;' % (k,dqdc_s))
                else:
                    self._write('dqdc[%d] = 0.0;' % k)

            self._write('for (int k=0; k<%d; k++) {' % nSpecies)
            self._indent()
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d*k+%d] += %.17g * dqdc[k];' % ((nSpecies+1), m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)
            #self._outdent()
            #self._write('}')

            self._outdent()
            self._write('}')

            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT; /* dwdot[%s]/dT */' % \
                        (nSpecies*(nSpecies+1)+m, all_dict[m][1], all_dict[m][0])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)

        else:

            for k in range(nSpecies):
                dqdc_s = dqdc_simple_precond('',k)
                if dqdc_s:
                    self._write('/* d()/d[%s] */' % all_dict[k][0])
                    self._write('dqdci = %s;' % (dqdc_s))
                    if reaction.reversible or k in rea_dict:
                        for m in sorted(all_dict.keys()):
                            if all_dict[m][1] != 0:
                                s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                                s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], all_dict[k][0])
                                self._write(s1.ljust(30) + s2)
            self._write('/* d()/dT */')
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT;' % (nSpecies*(nSpecies+1)+m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=').replace('+= -1 *', '-=')
                    s2 = '/* dwdot[%s]/dT */' % (all_dict[m][0])
                    self._write(s1.ljust(30) + s2)

        return




    def _ajac_GPU(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write('#ifdef AMREX_USE_CUDA')
        self._write(self.line('compute the reaction Jacobian on GPU'))
        self._write('AMREX_GPU_HOST_DEVICE')
        self._write('void aJacobian(double * J, double * sc, double T, int consP)')
        self._write('{')
        self._indent()

        self._write()

        self._write()

        self._write('for (int i=0; i<%d; i++) {' % (nSpecies+1)**2)
        self._indent()
        self._write('J[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double wdot[%d];' % (nSpecies))
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies))
        self._indent()
        self._write('wdot[k] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        self._write('double invT2 = invT * invT;')

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1.0 / refC;')

        self._write()

        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[k];')
        self._outdent()
        self._write('}')

        self._write()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        self._write(self.line('compute the species enthalpy'))
        self._write('double h_RT[%d];' % (nSpecies))
        self._write('speciesEnthalpy(h_RT, tc);')

        self._write()

        self._write('double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;') 
        self._write('double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;')
        self._write('double dqdci, dcdc_fac, dqdc[%d];' % (nSpecies))
        self._write('double Pr, fPr, F, k_0, logPr;') 
        self._write('double logFcent, troe_c, troe_n, troePr_den, troePr, troe;')
        self._write('double Fcent1, Fcent2, Fcent3, Fcent;')
        self._write('double dlogFdc, dlogFdn, dlogFdcn_fac;')
        self._write('double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;')
        self._write('const double ln10 = log(10.0);')
        self._write('const double log10e = 1.0/log(10.0);')

        for i, reaction in zip(range(nReactions), mechanism.reaction()):

            lt = reaction.lt
            if lt:
                print "Landau-Teller reactions are not supported"
                sys.exit(1)

            self._write(self.line('reaction %d: %s' % (i+1, reaction.equation())))
            if reaction.low:  # case 1
                self._write(self.line('a pressure-fall-off reaction'))
                self._ajac_reaction_d(mechanism, reaction, 1)
            elif reaction.thirdBody:  # case 2
                self._write(self.line('a third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction_d(mechanism, reaction, 2)
            else:  # case 3
                self._write(self.line('a non-third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction_d(mechanism, reaction, 3)
            self._write()

        self._write('double c_R[%d], dcRdT[%d], e_RT[%d];' % (nSpecies, nSpecies, nSpecies))
        self._write('double * eh_RT;')
        self._write('if (consP) {')
        self._indent()

        self._write('cp_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('eh_RT = &h_RT[0];');

        self._outdent()
        self._write('}')
        self._write('else {')
        self._indent()

        self._write('cv_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('speciesInternalEnergy(e_RT, tc);');
        self._write('eh_RT = &e_RT[0];');

        self._outdent()
        self._write('}')

        self._write()

        self._write('double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('cmix += c_R[k]*sc[k];')
        self._write('dcmixdT += dcRdT[k]*sc[k];')
        self._write('ehmix += eh_RT[k]*wdot[k];')
        self._write('dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];' % \
                        (nSpecies*(nSpecies+1)))
        self._outdent()
        self._write('}')

        self._write()
        self._write('double cmixinv = 1.0/cmix;')
        self._write('double tmp1 = ehmix*cmixinv;')
        self._write('double tmp3 = cmixinv*T;')
        self._write('double tmp2 = tmp1*tmp3;')
        self._write('double dehmixdc;')

        self._write('/* dTdot/d[X] */')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('dehmixdc = 0.0;')
        self._write('for (int m = 0; m < %d; ++m) {' % nSpecies)
        self._indent()
        self._write('dehmixdc += eh_RT[m]*J[k*%s+m];' % (nSpecies+1))
        self._outdent()
        self._write('}')        
        self._write('J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;' % (nSpecies+1,nSpecies))
        self._outdent()
        self._write('}')

        self._write('/* dTdot/dT */')
        self._write('J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;' % \
                        (nSpecies*(nSpecies+1)+nSpecies))

        self._outdent()
        self._write()
        self._write('return;')
        self._write('}')
        self._write('#endif')
        self._write()
        return


    def _ajac_reaction_d(self, mechanism, reaction, rcase):

        if rcase == 1: # pressure-dependent reaction
            isPD = True
            if reaction.thirdBody:
                has_alpha = True
                self._write('/* also 3-body */')
            else:
                has_alpha = False
                self._write('/* non 3-body */')
                print 'FIXME: pressure dependent non-3-body reaction in _ajac_reaction'
                sys.exit(1)
        elif rcase == 2: # third-body and non-pressure-dependent reaction
            isPD = False
            has_alpha = True
        elif rcase == 3: # simple non-third and non-pressure-dependent reaction
            isPD = False
            has_alpha = False
        else:
            print '_ajac_reaction: wrong case ', rcase
            exit(1)

        nSpecies = len(mechanism.species())
        rea_dict = {}
        pro_dict = {}
        all_dict = {}
        sumNuk = 0
        for symbol, coefficient in reaction.reactants:
            k = mechanism.species(symbol).id
            sumNuk -= coefficient
            if k in rea_dict:
                coe_old = rea_dict[k][1]
                rea_dict[k] = (symbol, coefficient+coe_old)
            else:
                rea_dict[k] = (symbol,  coefficient)
        for symbol, coefficient in reaction.products:
            k = mechanism.species(symbol).id
            sumNuk += coefficient
            if k in pro_dict:
                coe_old = pro_dict[k][1]
                pro_dict[k] = (symbol, coefficient+coe_old)
            else:
                pro_dict[k] = (symbol, coefficient)
        for k in range(nSpecies):
            if k in rea_dict and k in pro_dict:
                sr, nur = rea_dict[k]
                sp, nup = pro_dict[k]
                all_dict[k] = (sr, nup-nur)
            elif k in rea_dict:
                sr, nur = rea_dict[k]
                all_dict[k] = (sr, -nur)
            elif k in pro_dict:
                sp, nup = pro_dict[k]
                all_dict[k] = (sp, nup)

        sorted_reactants = sorted(rea_dict.values())
        sorted_products = sorted(pro_dict.values()) 

        if not reaction.reversible:
            if isPD or has_alpha:
                print 'FIXME: irreversible reaction in _ajac_reaction may not work'
                self._write('/* FIXME: irreversible reaction in _ajac_reaction may not work*/')
            for k in range(nSpecies):
                if k in sorted_reactants and k in sorted_products:
                    print 'FIXME: irreversible reaction in _ajac_reaction may not work'
                    self._write('/* FIXME: irreversible reaction in _ajac_reaction may not work*/')

        dim = self._phaseSpaceUnits(reaction.reactants)
        if isPD:
            Corr_s = 'Corr *'
            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 1 PD, TB
        elif has_alpha:
            Corr_s = 'alpha * '
            uc = self._prefactorUnits(reaction.units["prefactor"], -dim) # Case 2 !PD, TB
        else:
            Corr_s = ''
            uc = self._prefactorUnits(reaction.units["prefactor"], 1-dim) # Case 3 !PD, !TB 
        aeuc = self._activationEnergyUnits(reaction.units["activation"])

        if has_alpha:
            self._write("/* 3-body correction factor */")
            self._write("alpha = %s;" % self._enhancement_d(mechanism, reaction))

        # forward
        A, beta, E = reaction.arrhenius
        self._write('/* forward */')
        self._write("phi_f = %s;" % self._sortedPhaseSpace(mechanism, sorted_reactants))
        #
        self._write("k_f = %.17g * %.17g" % (uc.value,A))
        self._write("            * exp(%.17g * tc[0] - %.17g * (%.17g) * invT);"
                    %(beta,aeuc / Rc / kelvin,E))
        #
        self._write("dlnkfdT = %.17g * invT + %.17g *  %.17g  * invT2;"
                    %(beta,aeuc / Rc / kelvin,E))

        if isPD:
            low_A, low_beta, low_E = reaction.low
            self._write('/* pressure-fall-off */')
            self._write("k_0 = %.17g * exp(%.17g * tc[0] - %.17g * (%.17g) * invT);"
                        %(low_A,low_beta,aeuc / Rc / kelvin,low_E))
            self._write('Pr = 1e-%d * alpha / k_f * k_0;' % (dim*6))
            self._write('fPr = Pr / (1.0+Pr);')
            self._write("dlnk0dT = %.17g * invT + %.17g * (%.17g) * invT2;"
                        %(low_beta,aeuc / Rc / kelvin,low_E))
            self._write('dlogPrdT = log10e*(dlnk0dT - dlnkfdT);')
            self._write('dlogfPrdT = dlogPrdT / (1.0+Pr);')
            #
            if reaction.sri:
                self._write('/* SRI form */')
                print "FIXME: sri not supported in _ajac_reaction yet"
                sys.exit(1)
            elif reaction.troe:
                self._write('/* Troe form */')
                troe = reaction.troe
                ntroe = len(troe)
                self._write("logPr = log10(Pr);")
                if (abs(troe[1]) > 1.e-100):
                    if (troe[0] < 0):
                        self._write('Fcent1 = (1.+%.17g)*exp(-T/%.17g);'
                                %(-troe[0],troe[1]))
                    else:
                        self._write('Fcent1 = (1.-%.17g)*exp(-T/%.17g);'
                                %(troe[0],troe[1]))
                else:
                    self._write('Fcent1 = 0.;')
                if (abs(troe[2]) > 1.e-100):
                    self._write('Fcent2 = %.17g * exp(-T/%.17g);'
                                %(troe[0],troe[2]))
                else:
                    self._write('Fcent2 = 0.;')
                if (ntroe == 4):
                    if (troe[3] < 0):
                        self._write('Fcent3 = exp(%.17g * invT);'
                                % -troe[3] )
                    else:
                        self._write('Fcent3 = exp(-%.17g * invT);'
                                % troe[3] )
                else:
                    self._write('Fcent3 = 0.;')
                self._write('Fcent = Fcent1 + Fcent2 + Fcent3;')
                self._write("logFcent = log10(Fcent);")
                self._write("troe_c = -.4 - .67 * logFcent;")
                self._write("troe_n = .75 - 1.27 * logFcent;")
                self._write("troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));")
                self._write("troePr = (troe_c + logPr) * troePr_den;")
                self._write("troe = 1.0 / (1.0 + troePr*troePr);")
                self._write("F = pow(10.0, logFcent * troe);")

                self._write("dlogFcentdT = log10e/Fcent*( ")
                if (abs(troe[1]) > 1.e-100):
                    self._write("    -Fcent1/%.17g"
                                % troe[1] )
                else:
                    self._write("    0.")
                if (abs(troe[2]) > 1.e-100):
                    self._write("    -Fcent2/%.17g"
                                % troe[2] )
                else:
                    self._write("    0.")
                if (ntroe == 4):
                    self._write("    + Fcent3*%.17g*invT2);"
                                % troe[3] )
                else:
                    self._write("    + 0.);")

                self._write("dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;")
                self._write('dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;')
                self._write('dlogFdn = dlogFdcn_fac * troePr;')
                self._write('dlogFdlogPr = dlogFdc;')
                self._write('dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;')
            else:
                self._write('/* Lindemann form */')
                self._write('F = 1.0;')
                self._write('dlogFdlogPr = 0.0;')
                self._write('dlogFdT = 0.0;')

        # reverse
        if not reaction.reversible:
            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f;')
            else:
                self._write('q_nocor = k_f*phi_f;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %sdlnkfdT*k_f*phi_f;' % Corr_s)
        else:
            self._write('/* reverse */')
            self._write("phi_r = %s;" % self._sortedPhaseSpace(mechanism, sorted_products))
            self._write('Kc = %s;' % self._sortedKc(mechanism, reaction))
            self._write('k_r = k_f / Kc;')

            dlnKcdT_s = 'invT * ('
            terms = []
            for symbol, coefficient in sorted(sorted_reactants,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += '-(' + ' + '.join(terms) + ')'
            terms = []
            for symbol, coefficient in sorted(sorted_products,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += ' + (' + ' + '.join(terms) + ')'
            if sumNuk > 0:
                dlnKcdT_s += ' - %f' % sumNuk
            elif sumNuk < 0:
                dlnKcdT_s += ' + %f' % (-sumNuk)
            dlnKcdT_s += ')'
            self._write('dlnKcdT = %s;' % dlnKcdT_s)

            self._write('dkrdT = (dlnkfdT - dlnKcdT)*k_r;')

            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f - k_r*phi_r;')
            else:
                self._write('q_nocor = k_f*phi_f - k_r*phi_r;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);' % Corr_s)

        self._write("/* update wdot */")
        for k in sorted(all_dict.keys()):
            s, nu = all_dict[k]
            if nu == 1:
                self._write('wdot[%d] += q; /* %s */' % (k, s))
            elif nu == -1:
                self._write('wdot[%d] -= q; /* %s */' % (k, s))
            elif nu > 0:
                self._write('wdot[%d] += %.17g * q; /* %s */' % (k, nu, s))
            elif nu < 0:
                self._write('wdot[%d] -= %.17g * q; /* %s */' % (k, -nu, s))

        if isPD:
            self._write('/* for convenience */')
            self._write('k_f *= Corr;')
            if reaction.reversible:
                self._write('k_r *= Corr;')
        elif has_alpha:
            self._write('/* for convenience */')
            self._write('k_f *= alpha;')
            if reaction.reversible:
                self._write('k_r *= alpha;')
            else:
                self._write('k_r = 0.0;')

        if isPD:
            self._write('dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);')
        #elif has_alpha:
        #    self._write('dcdc_fac = q_nocor;')

        def dqdc_simple_d(dqdc_s, k):
            if dqdc_s ==  "0":
                dqdc_s = ''
            if k in sorted(rea_dict.keys()):
                dps = self._DphaseSpace(mechanism,sorted_reactants,rea_dict[k][0])
                if dps == "1.0":
                    dps_s = ''
                else:
                    dps_s = '*'+dps
                dqdc_s += ' + k_f%s' % dps_s
            if reaction.reversible:
                if k in sorted(pro_dict.keys()):
                    dps = self._DphaseSpace(mechanism,sorted_products,pro_dict[k][0])
                    if dps == "1.0":
                        dps_s = ''
                    else:
                        dps_s = '*'+dps
                    dqdc_s += ' - k_r%s' % dps_s
            return dqdc_s

        if has_alpha or isPD:

            self._write('if (consP) {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement_d(mechanism,reaction,k,True)
                if dqdc_s != "0":
                    if isPD:
                        if dqdc_s == "1":
                            dqdc_s ='dcdc_fac'
                        else:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == "1":
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple_d(dqdc_s,k)
                if dqdc_s:
                    symb_k = self.species[k].symbol
                    self._write('/* d()/d[%s] */' % symb_k)
                    self._write('dqdci = %s;' % (dqdc_s))
                    #
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                            s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                            s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], symb_k)
                            self._write(s1.ljust(30) + s2)

            self._outdent()
            self._write('}')
            self._write('else {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement_d(mechanism,reaction,k,False)
                if dqdc_s != '0':
                    if isPD:
                        if dqdc_s == '1':
                            dqdc_s ='dcdc_fac'
                        else:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == '1':
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple_d(dqdc_s,k)
                if dqdc_s:
                    self._write('dqdc[%d] = %s;' % (k,dqdc_s))

            self._write('for (int k=0; k<%d; k++) {' % nSpecies)
            self._indent()
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d*k+%d] += %.17g * dqdc[k];' % ((nSpecies+1), m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)
            self._outdent()
            self._write('}')

            self._outdent()
            self._write('}')

            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT; /* dwdot[%s]/dT */' % \
                        (nSpecies*(nSpecies+1)+m, all_dict[m][1], all_dict[m][0])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)

        else:

            for k in range(nSpecies):
                dqdc_s = dqdc_simple_d('',k)
                if dqdc_s:
                    self._write('/* d()/d[%s] */' % all_dict[k][0])
                    self._write('dqdci = %s;' % (dqdc_s))
                    if reaction.reversible or k in rea_dict:
                        for m in sorted(all_dict.keys()):
                            if all_dict[m][1] != 0:
                                s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                                s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], all_dict[k][0])
                                self._write(s1.ljust(30) + s2)
            self._write('/* d()/dT */')
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT;' % (nSpecies*(nSpecies+1)+m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=').replace('+= -1 *', '-=')
                    s2 = '/* dwdot[%s]/dT */' % (all_dict[m][0])
                    self._write(s1.ljust(30) + s2)

        return

    def _ajac(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line('compute the reaction Jacobian on CPU'))
        self._write('void aJacobian(double *  J, double *  sc, double T, int consP)')
        self._write('{')
        self._indent()

        self._write('for (int i=0; i<%d; i++) {' % (nSpecies+1)**2)
        self._indent()
        self._write('J[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double wdot[%d];' % (nSpecies))
        self._write('for (int k=0; k<%d; k++) {' % (nSpecies))
        self._indent()
        self._write('wdot[k] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        self._write('double invT2 = invT * invT;')

        self._write()

        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1.0 / refC;')

        self._write()

        self._write(self.line('compute the mixture concentration'))
        self._write('double mixture = 0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[k];')
        self._outdent()
        self._write('}')

        self._write()

        self._write(self.line('compute the Gibbs free energy'))
        self._write('double g_RT[%d];' % (nSpecies))
        self._write('gibbs(g_RT, tc);')

        self._write()

        self._write(self.line('compute the species enthalpy'))
        self._write('double h_RT[%d];' % (nSpecies))
        self._write('speciesEnthalpy(h_RT, tc);')

        self._write()

        self._write('double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;') 
        self._write('double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;')
        self._write('double dqdci, dcdc_fac, dqdc[%d];' % (nSpecies))
        self._write('double Pr, fPr, F, k_0, logPr;') 
        self._write('double logFcent, troe_c, troe_n, troePr_den, troePr, troe;')
        self._write('double Fcent1, Fcent2, Fcent3, Fcent;')
        self._write('double dlogFdc, dlogFdn, dlogFdcn_fac;')
        self._write('double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;')
        self._write('const double ln10 = log(10.0);')
        self._write('const double log10e = 1.0/log(10.0);')

        for i, reaction in zip(range(nReactions), mechanism.reaction()):

            lt = reaction.lt
            if lt:
                print "Landau-Teller reactions are not supported"
                sys.exit(1)

            self._write(self.line('reaction %d: %s' % (i+1, reaction.equation())))
            if reaction.low:  # case 1
                self._write(self.line('a pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 1)
            elif reaction.thirdBody:  # case 2
                self._write(self.line('a third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 2)
            else:  # case 3
                self._write(self.line('a non-third-body and non-pressure-fall-off reaction'))
                self._ajac_reaction(mechanism, reaction, 3)
            self._write()

        self._write('double c_R[%d], dcRdT[%d], e_RT[%d];' % (nSpecies, nSpecies, nSpecies))
        self._write('double * eh_RT;')
        self._write('if (consP) {')
        self._indent()

        self._write('cp_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('eh_RT = &h_RT[0];');

        self._outdent()
        self._write('}')
        self._write('else {')
        self._indent()

        self._write('cv_R(c_R, tc);')
        self._write('dcvpRdT(dcRdT, tc);')
        self._write('speciesInternalEnergy(e_RT, tc);');
        self._write('eh_RT = &e_RT[0];');

        self._outdent()
        self._write('}')

        self._write()

        self._write('double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('cmix += c_R[k]*sc[k];')
        self._write('dcmixdT += dcRdT[k]*sc[k];')
        self._write('ehmix += eh_RT[k]*wdot[k];')
        self._write('dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[%d+k];' % \
                        (nSpecies*(nSpecies+1)))
        self._outdent()
        self._write('}')

        self._write()
        self._write('double cmixinv = 1.0/cmix;')
        self._write('double tmp1 = ehmix*cmixinv;')
        self._write('double tmp3 = cmixinv*T;')
        self._write('double tmp2 = tmp1*tmp3;')
        self._write('double dehmixdc;')

        self._write('/* dTdot/d[X] */')
        self._write('for (int k = 0; k < %d; ++k) {' % nSpecies)
        self._indent()
        self._write('dehmixdc = 0.0;')
        self._write('for (int m = 0; m < %d; ++m) {' % nSpecies)
        self._indent()
        self._write('dehmixdc += eh_RT[m]*J[k*%s+m];' % (nSpecies+1))
        self._outdent()
        self._write('}')        
        self._write('J[k*%d+%d] = tmp2*c_R[k] - tmp3*dehmixdc;' % (nSpecies+1,nSpecies))
        self._outdent()
        self._write('}')

        self._write('/* dTdot/dT */')
        self._write('J[%d] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;' % \
                        (nSpecies*(nSpecies+1)+nSpecies))

        self._outdent()
        self._write('}')
        self._write('#endif')
        self._write()

        return

    def _ajac_reaction(self, mechanism, reaction, rcase):

        if rcase == 1: # pressure-dependent reaction
            isPD = True
            if reaction.thirdBody:
                has_alpha = True
                self._write('/* also 3-body */')
            else:
                has_alpha = False
                self._write('/* non 3-body */')
                print 'FIXME: pressure dependent non-3-body reaction in _ajac_reaction'
                sys.exit(1)
        elif rcase == 2: # third-body and non-pressure-dependent reaction
            isPD = False
            has_alpha = True
        elif rcase == 3: # simple non-third and non-pressure-dependent reaction
            isPD = False
            has_alpha = False
        else:
            print '_ajac_reaction: wrong case ', rcase
            exit(1)

        nSpecies = len(mechanism.species())
        rea_dict = {}
        pro_dict = {}
        all_dict = {}
        sumNuk = 0
        for symbol, coefficient in reaction.reactants:
            k = mechanism.species(symbol).id
            sumNuk -= coefficient
            if k in rea_dict:
                coe_old = rea_dict[k][1]
                rea_dict[k] = (symbol, coefficient+coe_old)
            else:
                rea_dict[k] = (symbol,  coefficient)
        for symbol, coefficient in reaction.products:
            k = mechanism.species(symbol).id
            sumNuk += coefficient
            if k in pro_dict:
                coe_old = pro_dict[k][1]
                pro_dict[k] = (symbol, coefficient+coe_old)
            else:
                pro_dict[k] = (symbol, coefficient)
        for k in range(nSpecies):
            if k in rea_dict and k in pro_dict:
                sr, nur = rea_dict[k]
                sp, nup = pro_dict[k]
                all_dict[k] = (sr, nup-nur)
            elif k in rea_dict:
                sr, nur = rea_dict[k]
                all_dict[k] = (sr, -nur)
            elif k in pro_dict:
                sp, nup = pro_dict[k]
                all_dict[k] = (sp, nup)

        sorted_reactants = sorted(rea_dict.values())
        sorted_products = sorted(pro_dict.values()) 

        if not reaction.reversible:
            if isPD or has_alpha:
                print 'FIXME: inreversible reaction in _ajac_reaction may not work'
                self._write('/* FIXME: inreversible reaction in _ajac_reaction may not work*/')
            for k in range(nSpecies):
                if k in sorted_reactants and k in sorted_products:
                    print 'FIXME: inreversible reaction in _ajac_reaction may not work'
                    self._write('/* FIXME: inreversible reaction in _ajac_reaction may not work*/')

        if isPD:
            Corr_s = 'Corr *'
        elif has_alpha:
            Corr_s = 'alpha * '
        else:
            Corr_s = ''

        if has_alpha:
            self._write("/* 3-body correction factor */")
            self._write("alpha = %s;" % self._enhancement(mechanism, reaction))

        # forward
        self._write('/* forward */')
        self._write("phi_f = %s;" % self._sortedPhaseSpace(mechanism, sorted_reactants))
        #
        self._write("k_f = prefactor_units[%d] * fwd_A[%d]" % (reaction.id-1,reaction.id-1))
        self._write("            * exp(fwd_beta[%d] * tc[0] - activation_units[%d] * fwd_Ea[%d] * invT);"
                    %(reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("dlnkfdT = fwd_beta[%d] * invT + activation_units[%d] * fwd_Ea[%d] * invT2;"
                    %(reaction.id-1,reaction.id-1,reaction.id-1))

        if isPD:
            self._write('/* pressure-fall-off */')
            self._write("k_0 = low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] * invT);"
                        %(reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('Pr = phase_units[%d] * alpha / k_f * k_0;' % (reaction.id-1))
            self._write('fPr = Pr / (1.0+Pr);')
            self._write("dlnk0dT = low_beta[%d] * invT + activation_units[%d] * low_Ea[%d] * invT2;"
                        %(reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('dlogPrdT = log10e*(dlnk0dT - dlnkfdT);')
            self._write('dlogfPrdT = dlogPrdT / (1.0+Pr);')
            if reaction.sri:
                self._write('/* SRI form */')
                print "FIXME: sri not supported in _ajac_reaction yet"
                sys.exit(1)
            elif reaction.troe:
                self._write('/* Troe form */')
                troe = reaction.troe
                self._write("logPr = log10(Pr);")
                self._write('Fcent1 = (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0.);'
                            %(reaction.id-1,reaction.id-1,reaction.id-1))
                self._write('Fcent2 = (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0.);'
                            %(reaction.id-1,reaction.id-1,reaction.id-1))
                self._write('Fcent3 = (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0.);'
                            %(reaction.id-1,reaction.id-1))
                self._write('Fcent = Fcent1 + Fcent2 + Fcent3;')
                self._write("logFcent = log10(Fcent);")
                self._write("troe_c = -.4 - .67 * logFcent;")
                self._write("troe_n = .75 - 1.27 * logFcent;")
                self._write("troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));")
                self._write("troePr = (troe_c + logPr) * troePr_den;")
                self._write("troe = 1.0 / (1.0 + troePr*troePr);")
                self._write("F = pow(10.0, logFcent * troe);")

                self._write("dlogFcentdT = log10e/Fcent*( ")
                self._write("    (fabs(troe_Tsss[%d]) > 1.e-100 ? -Fcent1/troe_Tsss[%d] : 0.)"
                            %(reaction.id-1,reaction.id-1))
                self._write("  + (fabs(troe_Ts[%d]) > 1.e-100 ? -Fcent2/troe_Ts[%d] : 0.)"
                            %(reaction.id-1,reaction.id-1))
                self._write("  + (troe_len[%d] == 4 ? Fcent3*troe_Tss[%d]*invT2 : 0.) );"
                            %(reaction.id-1,reaction.id-1))

                self._write("dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;")
                self._write('dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;')
                self._write('dlogFdn = dlogFdcn_fac * troePr;')
                self._write('dlogFdlogPr = dlogFdc;')
                self._write('dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;')
            else:
                self._write('/* Lindemann form */')
                self._write('F = 1.0;')
                self._write('dlogFdlogPr = 0.0;')
                self._write('dlogFdT = 0.0;')

        # reverse
        if not reaction.reversible:
            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f;')
            else:
                self._write('q_nocor = k_f*phi_f;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %sdlnkfdT*k_f*phi_f + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %sdlnkfdT*k_f*phi_f;' % Corr_s)
        else:
            self._write('/* reverse */')
            self._write("phi_r = %s;" % self._sortedPhaseSpace(mechanism, sorted_products))
            self._write('Kc = %s;' % self._sortedKc(mechanism, reaction))
            self._write('k_r = k_f / Kc;')

            dlnKcdT_s = 'invT * ('
            terms = []
            for symbol, coefficient in sorted(sorted_reactants,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += '-(' + ' + '.join(terms) + ')'
            terms = []
            for symbol, coefficient in sorted(sorted_products,
                                              key=lambda x: mechanism.species(x[0]).id):
                k = mechanism.species(symbol).id
                if coefficient == 1.0:
                    terms.append('h_RT[%d]' % (k))
                else:
                    terms.append('%f*h_RT[%d]' % (coefficient, k))
            dlnKcdT_s += ' + (' + ' + '.join(terms) + ')'
            if sumNuk > 0:
                dlnKcdT_s += ' - %f' % sumNuk
            elif sumNuk < 0:
                dlnKcdT_s += ' + %f' % (-sumNuk)
            dlnKcdT_s += ')'
            self._write('dlnKcdT = %s;' % dlnKcdT_s)

            self._write('dkrdT = (dlnkfdT - dlnKcdT)*k_r;')

            self._write('/* rate of progress */')
            if (not has_alpha) and (not isPD):
                self._write('q = k_f*phi_f - k_r*phi_r;')
            else:
                self._write('q_nocor = k_f*phi_f - k_r*phi_r;')
                if isPD:
                    self._write('Corr = fPr * F;')
                    self._write('q = Corr * q_nocor;')
                else: 
                    self._write('q = alpha * q_nocor;')

            if isPD:
                self._write('dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);')
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;' % Corr_s)
            else:
                self._write('dqdT = %s(dlnkfdT*k_f*phi_f - dkrdT*phi_r);' % Corr_s)

        self._write("/* update wdot */")
        for k in sorted(all_dict.keys()):
            s, nu = all_dict[k]
            if nu == 1:
                self._write('wdot[%d] += q; /* %s */' % (k, s))
            elif nu == -1:
                self._write('wdot[%d] -= q; /* %s */' % (k, s))
            elif nu > 0:
                self._write('wdot[%d] += %.17g * q; /* %s */' % (k, nu, s))
            elif nu < 0:
                self._write('wdot[%d] -= %.17g * q; /* %s */' % (k, -nu, s))

        if isPD:
            self._write('/* for convenience */')
            self._write('k_f *= Corr;')
            if reaction.reversible:
                self._write('k_r *= Corr;')
        elif has_alpha:
            self._write('/* for convenience */')
            self._write('k_f *= alpha;')
            if reaction.reversible:
                self._write('k_r *= alpha;')
            else:
                self._write('k_r = 0.0;')

        if isPD:
            self._write('dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);')
        #elif has_alpha:
        #    self._write('dcdc_fac = q_nocor;')

        def dqdc_simple(dqdc_s, k):
            if dqdc_s ==  "0":
                dqdc_s = ''
            if k in sorted(rea_dict.keys()):
                dps = self._DphaseSpace(mechanism,sorted_reactants,rea_dict[k][0])
                if dps == "1.0":
                    dps_s = ''
                else:
                    dps_s = '*'+dps
                dqdc_s += ' + k_f%s' % dps_s
            if reaction.reversible:
                if k in sorted(pro_dict.keys()):
                    dps = self._DphaseSpace(mechanism,sorted_products,pro_dict[k][0])
                    if dps == "1.0":
                        dps_s = ''
                    else:
                        dps_s = '*'+dps
                    dqdc_s += ' - k_r%s' % dps_s
            return dqdc_s

        if has_alpha or isPD:

            self._write('if (consP) {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement(mechanism,reaction,k,True)
                if dqdc_s != "0":
                    if isPD:
                        if dqdc_s == "1":
                            dqdc_s ='dcdc_fac'
                        else:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == "1":
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple(dqdc_s,k)
                if dqdc_s:
                    symb_k = self.species[k].symbol
                    self._write('/* d()/d[%s] */' % symb_k)
                    self._write('dqdci = %s;' % (dqdc_s))
                    #
                    for m in sorted(all_dict.keys()):
                        if all_dict[m][1] != 0:
                            s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                            s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                            s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], symb_k)
                            self._write(s1.ljust(30) + s2)

            self._outdent()
            self._write('}')
            self._write('else {')
            self._indent()

            for k in range(nSpecies):
                dqdc_s = self._Denhancement(mechanism,reaction,k,False)
                if dqdc_s != '0':
                    if isPD:
                        if dqdc_s == '1':
                            dqdc_s ='dcdc_fac'
                        elif isPD:
                            dqdc_s +='*dcdc_fac'
                    elif has_alpha:
                        if dqdc_s == '1':
                            dqdc_s ='q_nocor'
                        else:
                            dqdc_s +='*q_nocor'

                dqdc_s = dqdc_simple(dqdc_s,k)
                if dqdc_s:
                    self._write('dqdc[%d] = %s;' % (k,dqdc_s))

            self._write('for (int k=0; k<%d; k++) {' % nSpecies)
            self._indent()
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d*k+%d] += %.17g * dqdc[k];' % ((nSpecies+1), m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)
            self._outdent()
            self._write('}')

            self._outdent()
            self._write('}')

            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT; /* dwdot[%s]/dT */' % \
                        (nSpecies*(nSpecies+1)+m, all_dict[m][1], all_dict[m][0])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                    self._write(s1)

        else:

            for k in range(nSpecies):
                dqdc_s = dqdc_simple('',k)
                if dqdc_s:
                    self._write('/* d()/d[%s] */' % all_dict[k][0])
                    self._write('dqdci = %s;' % (dqdc_s))
                    if reaction.reversible or k in rea_dict:
                        for m in sorted(all_dict.keys()):
                            if all_dict[m][1] != 0:
                                s1 = 'J[%d] += %.17g * dqdci;' % (k*(nSpecies+1)+m, all_dict[m][1])
                                s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=')
                                s2 = '/* dwdot[%s]/d[%s] */' % (all_dict[m][0], all_dict[k][0])
                                self._write(s1.ljust(30) + s2)
            self._write('/* d()/dT */')
            for m in sorted(all_dict.keys()):
                if all_dict[m][1] != 0:
                    s1 = 'J[%d] += %.17g * dqdT;' % (nSpecies*(nSpecies+1)+m, all_dict[m][1])
                    s1 = s1.replace('+= 1 *', '+=').replace('+= -1 *', '-=').replace('+= -1 *', '-=')
                    s2 = '/* dwdot[%s]/dT */' % (all_dict[m][0])
                    self._write(s1.ljust(30) + s2)

        return


    def _vproductionRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write('#ifndef AMREX_USE_CUDA')
        self._write(self.line('compute the production rate for each species'))
        self._write('void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)')
        self._write('{')
        self._indent()

        self._write('double k_f_s[%d*npt], Kc_s[%d*npt], mixture[npt], g_RT[%d*npt];'
                    % (nReactions, nReactions, nSpecies))
        self._write('double tc[5*npt], invT[npt];')

        self._write()

        self._outdent()
        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write(' #pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('tc[0*npt+i] = log(T[i]);')
        self._write('tc[1*npt+i] = T[i];')
        self._write('tc[2*npt+i] = T[i]*T[i];')
        self._write('tc[3*npt+i] = T[i]*T[i]*T[i];')
        self._write('tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];')
        self._write('invT[i] = 1.0 / T[i];')
        self._outdent()
        self._write('}')

        self._write()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('mixture[i] = 0.0;')
        self._outdent()
        self._write('}')
        
        self._write()
        self._write('for (int n=0; n<%d; n++) {' % nSpecies)
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('mixture[i] += sc[n*npt+i];')
        self._write('wdot[n*npt+i] = 0.0;')
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()
        self._write('vcomp_k_f(npt, k_f_s, tc, invT);')
        self._write()
        self._write('vcomp_gibbs(npt, g_RT, tc);')
        self._write()
        self._write('vcomp_Kc(npt, Kc_s, g_RT, invT);')
        self._write()
        if nReactions <= 50:
            self._write('vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);')
        else:
            for i in range(0,nReactions,50):
                self._write('vcomp_wdot_%d_%d(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);' % (i+1,min(i+50,nReactions)))

        self._outdent()
        self._write('}')

        self._write()

        self._write('void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)')
        self._write('{')
        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write('#pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        for reaction in mechanism.reaction():
            self._write("k_f_s[%d*npt+i] = prefactor_units[%d] * fwd_A[%d] * exp(fwd_beta[%d] * tc[i] - activation_units[%d] * fwd_Ea[%d] * invT[i]);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1,
                           reaction.id-1,reaction.id-1,reaction.id-1))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')        

        self._write()

        self._write('void vcomp_gibbs(int npt, double *  g_RT, double *  tc)')
        self._write('{')
        self._indent()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write('double tg[5], g[%d];' % nSpecies)
        self._write('tg[0] = tc[0*npt+i];')
        self._write('tg[1] = tc[1*npt+i];')
        self._write('tg[2] = tc[2*npt+i];')
        self._write('tg[3] = tc[3*npt+i];')
        self._write('tg[4] = tc[4*npt+i];')
        self._write()
        self._write('gibbs(g, tg);')
        self._write()
        for ispec in range(nSpecies):
            self._write('g_RT[%d*npt+i] = g[%d];' % (ispec, ispec))
        self._outdent()
        self._write('}')
        self._outdent()
        self._write('}')

        self._write()

        self._write('void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)')
        self._write('{')
        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write('#pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = (101325. / 8.31451) * invT[i];');
        self._write('double refCinv = 1.0 / refC;');
        self._write()
        for reaction in mechanism.reaction():
            K_c = self._vKc(mechanism, reaction)
            self._write("Kc_s[%d*npt+i] = %s;" % (reaction.id-1,K_c))
        self._outdent()
        self._write('}')        
        self._outdent()
        self._write('}')        

        self._write()
        if nReactions <= 50:
            self._write('void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,')
            self._write('		double *  k_f_s, double *  Kc_s,')
            self._write('		double *  tc, double *  invT, double *  T)')
            self._write('{')
            self._vcomp_wdot(mechanism,0,nReactions)
            self._write('}')
        else:
            for i in range(0,nReactions,50):
                nr = min(50, nReactions-i)
                self._write('void vcomp_wdot_%d_%d(int npt, double *  wdot, double *  mixture, double *  sc,' % (i+1,i+nr))
                self._write('		double *  k_f_s, double *  Kc_s,')
                self._write('		double *  tc, double *  invT, double *  T)')
                self._write('{')
                self._vcomp_wdot(mechanism,i,nr)
                self._write('}')
                self._write()

        self._write('#endif')

        return

    def _vcomp_wdot(self, mechanism, istart, nr):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        # self._write('#ifdef __INTEL_COMPILER')
        # self._indent()
        # self._write('#pragma simd')
        # self._outdent()
        # self._write('#endif')
        self._indent()
        self._write('for (int i=0; i<npt; i++) {')
        self._indent()

        self._write('double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;')
        self._write('double alpha;')
        #if istart < isimple[0]:
        #    self._write('double alpha;')
        if istart < i3body[0]:
            self._write('double redP, F;') 
        if istart < ilindemann[0]:
            self._write('double logPred;')
            if ntroe>0:
                self._write('double logFcent, troe_c, troe_n, troe, F_troe;')
            if nsri>0:
                self._write('double X, F_sri;')

        first_id = istart + 1
        last_id  = istart + nr

        for reaction in mechanism.reaction():

            if reaction.id < first_id or reaction.id > last_id:
                continue

            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            # compute the rates
            self._vforwardRate(mechanism, reaction)
            self._vreverseRate(mechanism, reaction)

            # store the progress rate
            self._write("qdot = q_f - q_r;")

            agents = list(set(reaction.reactants + reaction.products))
            agents = sorted(agents, key=lambda x: mechanism.species(x[0]).id)
            # note that a species might appear as both reactant and product
            # a species might alos appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a:
                        if coefficient == 1.0:
                            self._write("wdot[%d*npt+i] -= qdot;" 
                                        % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d*npt+i] -= %f * qdot;" 
                                        % (mechanism.species(symbol).id, coefficient))
                for b in reaction.products: 
                    if b == a:
                        if coefficient == 1.0:
                            self._write("wdot[%d*npt+i] += qdot;" 
                                        % (mechanism.species(symbol).id))
                        else:
                            self._write("wdot[%d*npt+i] += %f * qdot;" 
                                        % (mechanism.species(symbol).id, coefficient))

        self._outdent()
        self._write('}')
        self._outdent()        

        return


    def _progressRate(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write(self.line('compute the progress rate for each reaction'))
        self._write('AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  sc, double T)')
        self._write('{')
        self._indent()

        self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
        self._write('double invT = 1.0 / tc[1];')
        
        self._write()
        self._outdent()
        self._write('#ifndef AMREX_USE_CUDA')
        self._indent()
        self._write('if (T != T_save)')
        self._write('{')
        self._indent()
        self._write('T_save = T;')
        self._write('comp_k_f(tc,invT,k_f_save);');
        self._write('comp_Kc(tc,invT,Kc_save);');
        self._outdent()
        self._write("}")
        self._outdent()
        self._write('#endif')
        self._indent()

        if (nReactions == 0):
            self._write()
        else:
            self._write()
            self._write('double q_f[%d], q_r[%d];' % (nReactions,nReactions))
            self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');
            self._write()
            self._write('for (int i = 0; i < %d; ++i) {' % nReactions)
            self._indent()
            self._write('qdot[i] = q_f[i] - q_r[i];')
            self._outdent()
            self._write('}')

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _initializeRateCalculation(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # declarations
        self._write()
        self._write('int id; ' + self.line('loop counter'))

        self._write('double mixture;                 '
                    + self.line('mixture concentration'))
        self._write('double g_RT[%d];                ' % nSpecies
                    + self.line('Gibbs free energy'))

        self._write('double Kc;                      ' + self.line('equilibrium constant'))
        self._write('double k_f;                     ' + self.line('forward reaction rate'))
        self._write('double k_r;                     ' + self.line('reverse reaction rate'))
        self._write('double q_f;                     ' + self.line('forward progress rate'))
        self._write('double q_r;                     ' + self.line('reverse progress rate'))
        self._write('double phi_f;                   '
                    + self.line('forward phase space factor'))
        self._write('double phi_r;                   '
                    + self.line('reverse phase space factor'))
        self._write('double alpha;                   ' + self.line('enhancement'))


        self._write('double redP;                    ' + self.line('reduced pressure'))
        self._write('double logPred;                 ' + self.line('log of above'))
        self._write('double F;                       '
                    + self.line('fallof rate enhancement'))
        self._write()
        self._write('double F_troe;                  ' + self.line('TROE intermediate'))
        self._write('double logFcent;                ' + self.line('TROE intermediate'))
        self._write('double troe;                    ' + self.line('TROE intermediate'))
        self._write('double troe_c;                  ' + self.line('TROE intermediate'))
        self._write('double troe_n;                  ' + self.line('TROE intermediate'))
        self._write()

        self._write(
            'double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; '
            + self.line('temperature cache'))

        self._write()
        self._write('double invT = 1.0 / tc[1];')

        # compute the reference concentration
        self._write()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))
        self._write('double refCinv = 1 / refC;')

        # compute the mixture concentration
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('mixture = 0.0;')
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[id];')
        self._outdent()
        self._write('}')

        # compute the Gibbs free energies
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(g_RT, tc);')
        
        return


    def _forwardRate(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._phaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)

        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_save[%d];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return
            
        alpha = self._enhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_save[%d];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return

        self._write("k_f = k_f_save[%d];" % (reaction.id-1))

        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] *invT);"
                    %(reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T)" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[0]) : 1);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0) );' 
                        % (reaction.id-1,reaction.id-1))
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10, logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f = phi_f * k_f;")
        return
        

    def _vforwardRate(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._vphaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)
                
        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_s[%d*npt+i];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return
            
        alpha = self._venhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_s[%d*npt+i];" % (reaction.id-1))
            self._write("q_f = phi_f * k_f;")
            return

        self._write("k_f = k_f_s[%d*npt+i];" % (reaction.id-1))
        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[i] - activation_units[%d] * low_Ea[%d] * invT[i]);" 
                    % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T[i])" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T[i]/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[i]) : 1.);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T[i]/troe_Tsss[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T[i]/troe_Ts[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT[i]) : 0.) );' 
                        % (reaction.id-1,reaction.id-1))
            
            d = .14
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10., logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f = phi_f * k_f;")
        return
        

    def _reverseRate(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r = 0.0;")
            return

        phi_r = self._phaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:

            self._write("k_r = prefactor_units_rev[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[0] - activation_units_rev[%d] * rev_Ea[%d] * invT);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_r = phi_r * k_r;")
            return
        
        self._write("Kc = Kc_save[%d];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")
        self._write("q_r = phi_r * k_r;")

        return


    def _vreverseRate(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r = 0.0;")
            return

        phi_r = self._vphaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:

            self._write("k_r = prefactor_units[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[i] - activation_units[%d] * rev_Ea[%d] * invT[i]);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_r = phi_r * k_r;")
            return
        
        self._write("Kc = Kc_s[%d*npt+i];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")
        self._write("q_r = phi_r * k_r;")

        return


    def _progressRateFR(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        itroe      = self.reactionIndex[0:2]
        isri       = self.reactionIndex[1:3]
        ilindemann = self.reactionIndex[2:4]
        i3body     = self.reactionIndex[3:5] 
        isimple    = self.reactionIndex[4:6]
        ispecial   = self.reactionIndex[5:7]

        if len(self.reactionIndex) != 7:
            print '\n\nCheck this!!!\n'
            sys.exit(1)
        
        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        self._write()
        self._write()
        self._write(self.line('compute the progress rate for each reaction'))
        self._write('AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)')
        self._write('{')
        self._indent()

        if (nReactions > 0):

            self._write('double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */')
            self._write('double invT = 1.0 / tc[1];')
            
            self._outdent()
            self._write('#ifndef AMREX_USE_CUDA')
            self._indent()
            self._write()
            self._write('if (T != T_save)')
            self._write('{')
            self._indent()
            self._write('T_save = T;')
            self._write('comp_k_f(tc,invT,k_f_save);');
            self._write('comp_Kc(tc,invT,Kc_save);');
            self._outdent()
            self._write("}")
            self._outdent()
            self._write('#endif')
            self._indent()

            self._write()
            self._write('comp_qfqr(q_f, q_r, sc, tc, invT);');

            self._write()

        self._write('return;')
        self._outdent()

        self._write('}')

        return

    def _getCriticalParameters(self, mechanism):
      
        
        TabulatedCriticalParams = {
          "H2":{'Tci':33.145,"Pci":12.964,"wt":2.01588,"acentric_factor":-0.219},
            "O2":{'Tci':154.581,"Pci":50.4304658,"wt":31.9988,"acentric_factor":0.0222},
            "H2O":{'Tci':647.096,"Pci":220.640,"wt":18.015340,"acentric_factor":0.3443},
            "N2":{'Tci':126.192,"Pci":33.958,"wt":28.013400,"acentric_factor":0.0372},
            "CH4":{'Tci':190.56,"Pci":45.99,"wt":16.043030,"acentric_factor":0.011},
            "C2H6":{'Tci':305.32,"Pci":48.72,"wt":30.070120,"acentric_factor":0.099},
            "C3H8":{'Tci':369.83,"Pci":42.48,"wt":44.097210,"acentric_factor":0.152},
            "CO2":{'Tci':304.12,"Pci":73.74,"wt":44.009950,"acentric_factor":0.225},
            "He":{'Tci':5.1953,"Pci":2.2746,"wt":4.002602,"acentric_factor":-0.382},
            "CO":{'Tci':132.85,"Pci":34.94,"wt":28.010,"acentric_factor":0.045},
            "AR":{'Tci':150.86,"Pci":48.98,"wt":39.948,"acentric_factor":-0.002},
            "NO":{'Tci':180.0,"Pci":64.80,"wt":30.006,"acentric_factor":0.582},
            "CH3OH":{'Tci':512.64,"Pci":80.97,"wt":32.042,"acentric_factor":0.565},
            "C2H2":{'Tci':308.30,"Pci":61.14,"wt":26.038,"acentric_factor":0.189},
            "C2H4":{'Tci':282.34,"Pci":50.41,"wt":28.054,"acentric_factor":0.087},
            "N2O":{'Tci':309.60,"Pci":72.55,"wt":44.013,"acentric_factor":0.162}
                    }
      
        nSpecies = len(mechanism.species()) 
        self._write()
        self._write()
        self._write(self.line('compute the critical parameters for each species'))
        self._write('void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)')
        self._write('{')
        self._write()
        self._indent()
        
        
        self._write('double   EPS[%d];'%nSpecies)
        self._write('double   SIG[%d];' %nSpecies)
        self._write('double    wt[%d];' %nSpecies)
        self._write('double avogadro = 6.02214199e23;')
        self._write('double boltzmann = 1.3806503e-16; //we work in CGS')
        self._write('double Rcst = 83.144598; //in bar [CGS] !')
        
        self._write()

        self._write('egtransetEPS(EPS);')
        self._write('egtransetSIG(SIG);')
        self._write('get_mw(wt);')


        for species in mechanism.species():

          if species.symbol in TabulatedCriticalParams:
          
            self._write()
            self._write(self.line('species %d: %s' % (species.id, species.symbol)))
            self._write(self.line('Imported from NIST'))
            self._write('Tci[%d] = %f ; ' % (
                species.id,TabulatedCriticalParams[species.symbol]["Tci"]))
            self._write('ai[%d] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[%d],2.0) / (pow(%f,2.0) * %f); ' % (
                species.id,species.id,TabulatedCriticalParams[species.symbol]["wt"],TabulatedCriticalParams[species.symbol]["Pci"]))
            self._write('bi[%d] = 0.08664 * Rcst * Tci[%d] / (%f * %f); ' % (
                species.id,species.id,TabulatedCriticalParams[species.symbol]["wt"],TabulatedCriticalParams[species.symbol]["Pci"]))
            self._write('acentric_i[%d] = %f ;'
                %(species.id,TabulatedCriticalParams[species.symbol]["acentric_factor"]))
          else:
                      
            self._write()
            self._write(self.line('species %d: %s' % (species.id, species.symbol)))
            self._write('Tci[%d] = 1.316 * EPS[%d] ; ' % (
                species.id,species.id))
            self._write('ai[%d] = (5.55 * pow(avogadro,2.0) * EPS[%d]*boltzmann * pow(1e-8*SIG[%d],3.0) ) / (pow(wt[%d],2.0)); ' % (
                species.id,species.id,species.id,species.id))
            self._write('bi[%d] = 0.855 * avogadro * pow(1e-8*SIG[%d],3.0) / (wt[%d]); ' % (
                species.id,species.id,species.id))
            self._write('acentric_i[%d] = 0.0 ;'
                %(species.id))
        
        self._write()
        self._write('return;')
        self._outdent() 
        self._write('}')

        return

    def _initializeRateCalculationFR(self, mechanism):

        nSpecies = len(mechanism.species())
        nReactions = len(mechanism.reaction())

        # declarations
        self._write()
        self._write('int id; ' + self.line('loop counter'))

        self._write('double mixture;                 '
                    + self.line('mixture concentration'))
        self._write('double g_RT[%d];                ' % nSpecies
                    + self.line('Gibbs free energy'))

        self._write('double Kc;                      ' + self.line('equilibrium constant'))
        self._write('double k_f;                     ' + self.line('forward reaction rate'))
        self._write('double k_r;                     ' + self.line('reverse reaction rate'))
        self._write('double phi_f;                   '
                    + self.line('forward phase space factor'))
        self._write('double phi_r;                   '
                    + self.line('reverse phase space factor'))
        self._write('double alpha;                   ' + self.line('enhancement'))


        self._write('double redP;                    ' + self.line('reduced pressure'))
        self._write('double logPred;                 ' + self.line('log of above'))
        self._write('double F;                       '
                    + self.line('fallof rate enhancement'))
        self._write()
        self._write('double F_troe;                  ' + self.line('TROE intermediate'))
        self._write('double logFcent;                ' + self.line('TROE intermediate'))
        self._write('double troe;                    ' + self.line('TROE intermediate'))
        self._write('double troe_c;                  ' + self.line('TROE intermediate'))
        self._write('double troe_n;                  ' + self.line('TROE intermediate'))
        self._write()

        self._write(
            'double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; '
            + self.line('temperature cache'))

        self._write()
        self._write('double invT = 1.0 / tc[1];')

        # compute the reference concentration
        self._write()
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))

        # compute the mixture concentration
        self._write()
        self._write(self.line('compute the mixture concentration'))
        self._write('mixture = 0.0;')
        self._write('for (id = 0; id < %d; ++id) {' % nSpecies)
        self._indent()
        self._write('mixture += sc[id];')
        self._outdent()
        self._write('}')

        # compute the Gibbs free energies
        self._write()
        self._write(self.line('compute the Gibbs free energy'))
        self._write('gibbs(g_RT, tc);')
        
        return


    def _forwardRateFR(self, mechanism, reaction):

        lt = reaction.lt
        if lt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return self._landau(reaction)

        dim = self._phaseSpaceUnits(reaction.reactants)

        phi_f = self._phaseSpace(mechanism, reaction.reactants)
        self._write("phi_f = %s;" % phi_f)

        thirdBody = reaction.thirdBody
        if not thirdBody:
            self._write("k_f = k_f_save[%d];" % (reaction.id-1))
            self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
            return
            
        alpha = self._enhancement(mechanism, reaction)
        self._write("alpha = %s;" % alpha)

        sri = reaction.sri
        low = reaction.low
        troe = reaction.troe

        if not low:
            self._write("k_f = alpha * k_f_save[%d];" % (reaction.id-1))
            self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
            return

        self._write("k_f = k_f_save[%d];" % (reaction.id-1))

        self._write("redP = alpha / k_f * phase_units[%d] * low_A[%d] * exp(low_beta[%d] * tc[0] - activation_units[%d] * low_Ea[%d] *invT);"
                    % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))
        self._write("F = redP / (1 + redP);")

        if sri:
            self._write("logPred = log10(redP);")
            self._write("X = 1.0 / (1.0 + logPred*logPred);")
            self._write("F_sri = exp(X * log(sri_a[%d] * exp(-sri_b[%d]/T)" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   +  (sri_c[%d] > 1.e-100 ? exp(T/sri_c[%d]) : 0.) )" 
                        % (reaction.id-1,reaction.id-1))
            self._write("   *  (sri_len[%d] > 3 ? sri_d[%d]*exp(sri_e[%d]*tc[0]) : 1.);" 
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write("F *= F_sri;")

        elif troe:
            self._write("logPred = log10(redP);")

            self._write('logFcent = log10(')
            self._write('    (fabs(troe_Tsss[%d]) > 1.e-100 ? (1.-troe_a[%d])*exp(-T/troe_Tsss[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (fabs(troe_Ts[%d]) > 1.e-100 ? troe_a[%d] * exp(-T/troe_Ts[%d]) : 0.) '
                        % (reaction.id-1,reaction.id-1,reaction.id-1))
            self._write('    + (troe_len[%d] == 4 ? exp(-troe_Tss[%d] * invT) : 0) );'
                        % (reaction.id-1,reaction.id-1))

            d = .14
            self._write("troe_c = -.4 - .67 * logFcent;")
            self._write("troe_n = .75 - 1.27 * logFcent;")
            self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
            self._write("F_troe = pow(10, logFcent / (1.0 + troe*troe));")
            self._write("F *= F_troe;")

        self._write("k_f *= F;")
        self._write("q_f[%d] = phi_f * k_f;" % (reaction.id - 1))
        return
        

    def _reverseRateFR(self, mechanism, reaction):
        if not reaction.reversible:
            self._write("q_r[%d] = 0.0;" % (reaction.id - 1))
            return

        phi_r = self._phaseSpace(mechanism, reaction.products)
        self._write("phi_r = %s;" % phi_r)

        if reaction.rlt:
            import pyre
            pyre.debug.Firewall.hit("Landau-Teller reactions are not supported yet")
            return

        if reaction.rev:
            self._write("k_r = prefactor_units[%d] * rev_A[%d] * exp(rev_beta[%d] * tc[0] - activation_units[%d] * rev_Ea[%d] * invT);"
                        % (reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1,reaction.id-1))

            thirdBody = reaction.thirdBody
            if thirdBody:
                self._write("k_r *= alpha;")

            self._write("q_r[%d] = phi_r * k_r;" % (reaction.id - 1))
            return
        
        self._write("Kc = Kc_save[%d];" % (reaction.id-1))

        self._write("k_r = k_f / Kc;")

        self._write("q_r[%d] = phi_r * k_r;" % (reaction.id - 1))

        return

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

        #print "UNITS/SECOND/EXP ", units.value, second.value, exponent
        #print " units ** exponent / second (value) " , units.value ** exponent / second.value
        #print " units ** exponent / second (returned) " , units ** exponent / second
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


    def _equilibriumConstants(self, mechanism):
        self._write()
        self._write()
        self._write(self.line('compute the equilibrium constants for each reaction'))
        self._write('void equilibriumConstants(double *  kc, double *  g_RT, double T)')
        self._write('{')
        self._indent()

        # compute the reference concentration
        self._write(self.line('reference concentration: P_atm / (RT) in inverse mol/m^3'))
        self._write('double refC = %g / %g / T;' % (atm.value, R.value))

        # compute the equilibrium constants
        for reaction in mechanism.reaction():
            self._write()
            self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))

            K_c = self._Kc(mechanism, reaction)
            self._write("kc[%d] = %s;" % (reaction.id - 1, K_c))

        self._write()
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _phaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in reagents:
            if (coefficient == "1.0"):
                conc = "sc[%d]" % mechanism.species(symbol).id
            else:
                conc = "pow(sc[%d], %f)" % (mechanism.species(symbol).id, coefficient)
            phi += [conc]

        return "*".join(phi)

    def _sortedPhaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            if (float(coefficient) == 1.0):
                conc = "sc[%d]" % mechanism.species(symbol).id
            else:
                conc = "pow(sc[%d], %f)" % (mechanism.species(symbol).id, float(coefficient))
            phi += [conc]

        return "*".join(phi)


    def _DphaseSpace(self, mechanism, reagents, r):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            if symbol == r:
                if coefficient > 1:
                    phi += ["%f" % coefficient]
                    if ((coefficient-1) == 1.0):
                        conc = "sc[%d]" % (mechanism.species(symbol).id)
                    else:
                        conc = "pow(sc[%d],%f)" % (mechanism.species(symbol).id, (coefficient-1))
                    phi += [conc]
            else:
                if (coefficient == 1.0):
                    conc = "sc[%d]" % mechanism.species(symbol).id
                else:
                    conc = "pow(sc[%d], %f)" % (mechanism.species(symbol).id, coefficient)
                phi += [conc]

        if phi:
            return "*".join(phi)
        else:
            return "1.0"


    def _vphaseSpace(self, mechanism, reagents):

        phi = []

        for symbol, coefficient in sorted(reagents,key=lambda x:mechanism.species(x[0]).id):
            if (coefficient == 1.0):
                conc = "sc[%d*npt+i]" % (mechanism.species(symbol).id)
            else:
                conc = "pow(sc[%d*npt+i], %f)" % (mechanism.species(symbol).id, coefficient)
            phi += [conc]

        return "*".join(phi)


    def _phaseSpaceUnits(self, reagents):
        dim = 0.0
        for symbol, coefficient in reagents:
            dim += float(coefficient)

        return dim


    def _enhancement_d(self, mechanism, reaction):
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
            #return "sc[%d]" % mechanism.species(species).id
            return "mixture"

        alpha = ["mixture"]
        for i, eff in enumerate(efficiencies):
            symbol, efficiency = eff
            factor = "( %.17g - 1)" % (efficiency)
            conc = "sc[%d]" % mechanism.species(symbol).id
            alpha.append("%s*%s" % (factor, conc))

        return " + ".join(alpha).replace('+ -','- ')


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
        for i, eff in enumerate(efficiencies):
            symbol, efficiency = eff
            factor = "(TB[%d][%d] - 1)" % (reaction.id-1, i)
            conc = "sc[%d]" % mechanism.species(symbol).id
            alpha.append("%s*%s" % (factor, conc))

        return " + ".join(alpha).replace('+ -','- ')

    def _Denhancement(self, mechanism, reaction, kid, consP):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                if consP:
                    return "0"
                else:
                    return "1"
            elif mechanism.species(species).id == kid:
                return "1"
            else:
                return "0"
        else:
            if consP:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "(TB[%d][%d] - 1)" % (reaction.id-1, i)
                return "0"
            else:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "TB[%d][%d]" % (reaction.id-1,i)
                return "1"

    def _Denhancement_d(self, mechanism, reaction, kid, consP):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                if consP:
                    return "0"
                else:
                    return "1"
            elif mechanism.species(species).id == kid:
                return "1"
            else:
                return "0"
        else:
            if consP:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "(%.17g - 1)" % (efficiency)
                return "0"
            else:
                for i, eff in enumerate(efficiencies):
                    symbol, efficiency = eff
                    if mechanism.species(symbol).id == kid:
                        return "%.17g" % (efficiency)
                return "1"

    def _venhancement(self, mechanism, reaction):
        thirdBody = reaction.thirdBody
        if not thirdBody:
            import pyre
            pyre.debug.Firewall.hit("_enhancement called for a reaction without a third body")
            return

        species, coefficient = thirdBody
        efficiencies = reaction.efficiencies

        if not efficiencies:
            if species == "<mixture>":
                return "mixture[i]"
            return "sc[%d*npt+i]" % mechanism.species(species).id

        alpha = ["mixture[i]"]
        for i, eff in enumerate(efficiencies):
            symbol, efficiency = eff
            factor = "(TB[%d][%d] - 1)" % (reaction.id-1,i)
            conc = "sc[%d*npt+i]" % mechanism.species(symbol).id
            alpha.append("%s*%s" % (factor, conc))

        return " + ".join(alpha)


    def _cv(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cv/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("cv_R", self._cvNASA, speciesInfo)

        return

    def _cv_GPU_H(self):

        self._write(self.line('compute Cv/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("cv_R_d")

        return

    def _cv_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cv/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("cv_R", self._cvNASA, speciesInfo)

        return
    
    def _cp(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cp/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("cp_R", self._cpNASA, speciesInfo)

        return

    def _cp_GPU_H(self):

        self._write(self.line('compute Cp/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("cp_R")

        return

    def _cp_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute Cp/R at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("cp_R", self._cpNASA, speciesInfo)

        return

    def _dcvpdT(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("dcvpRdT", self._dcpdTNASA, speciesInfo)

        return

    def _dcvpRdT_GPU_H(self):

        self._write(self.line('compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("dcvpRdT")

        return

    def _dcvpRdT_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("dcvpRdT", self._dcpdTNASA, speciesInfo)

        return


    def _gibbs(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the g/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("gibbs", self._gibbsNASA, speciesInfo, 1)

        return

    def _gibbs_GPU_H(self):

        self._write(self.line('compute the g/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("gibbs")

        return

    def _gibbs_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the g/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("gibbs", self._gibbsNASA, speciesInfo, 1)

        return

    def _helmholtz(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the a/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("helmholtz", self._helmholtzNASA, speciesInfo, 1)

        return

    def _helmholtz_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the a/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("helmholtz", self._helmholtzNASA, speciesInfo, 1)

        return

    def _speciesEntropy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the S/R at the given temperature (Eq 21)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesEntropy", self._entropyNASA, speciesInfo)

        return

    def _speciesEntropy_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the S/R at the given temperature (Eq 21)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("speciesEntropy", self._entropyNASA, speciesInfo)

        return

    def _speciesInternalEnergy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the e/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesInternalEnergy", self._internalEnergy, speciesInfo, 1)

        return

    def _speciesInternalEnergy_GPU_H(self):

        self._write(self.line('compute the e/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("speciesInternalEnergy_d")

        return

    def _speciesInternalEnergy_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the e/(RT) at the given temperature'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("speciesInternalEnergy", self._internalEnergy, speciesInfo, 1)

        return

    def _speciesEnthalpy(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the h/(RT) at the given temperature (Eq 20)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine("speciesEnthalpy", self._enthalpyNASA, speciesInfo, 1)

        return

    def _speciesEnthalpy_GPU_H(self):

        self._write(self.line('compute the h/(RT) at the given temperature (Eq 20)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU_H("speciesEnthalpy")

        return

    def _speciesEnthalpy_GPU(self, speciesInfo):

        self._write()
        self._write()
        self._write(self.line('compute the h/(RT) at the given temperature (Eq 20)'))
        self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
        self._generateThermoRoutine_GPU("speciesEnthalpy", self._enthalpyNASA, speciesInfo, 1)

        return

    def _generateThermoRoutine_GPU_H(self, name):

        self._write('AMREX_GPU_HOST_DEVICE void %s(double * species, double *  tc);' % name)

        return

    def _generateThermoRoutine_GPU(self, name, expressionGenerator, speciesInfo, needsInvT=0):

        lowT, highT, midpoints = speciesInfo
        
        self._write('AMREX_GPU_HOST_DEVICE void %s(double * species, double *  tc)' % name)
        self._write('{')

        self._indent()

        # declarations
        self._write()
        self._write(self.line('temperature'))
        self._write('double T = tc[1];')
        if needsInvT != 0:
           self._write('double invT = 1 / T;')
        if needsInvT == 2:
           self._write('double invT2 = invT*invT;')

        # temperature check
        # self._write()
        # self._write(self.line('check the temperature value'))
        # self._write('if (T < %g || T > %g) {' % (lowT, highT))
        # self._indent()
        # self._write(
        #     'fprintf(stderr, "temperature %%g is outside the range (%g, %g)", T);'
        #     % (lowT, highT))
        # self._write('return;')
        # self._outdent()
        # self._write('}')
                    
        for midT, speciesList in midpoints.items():

            self._write('')
            self._write(self.line('species with midpoint at T=%g kelvin' % midT))
            self._write('if (T < %g) {' % midT)
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(lowRange.parameters)
                self._outdent()

            self._outdent()
            self._write('} else {')
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(highRange.parameters)
                self._outdent()

            self._outdent()
            self._write('}')
            
        self._write('return;')
        self._outdent()

        self._write('}')

        return

    def _generateThermoRoutine(self, name, expressionGenerator, speciesInfo, needsInvT=0):

        lowT, highT, midpoints = speciesInfo
        
        self._write('void %s(double *  species, double *  tc)' % name)
        self._write('{')

        self._indent()

        # declarations
        self._write()
        self._write(self.line('temperature'))
        self._write('double T = tc[1];')
        if needsInvT != 0:
           self._write('double invT = 1 / T;')
        if needsInvT == 2:
           self._write('double invT2 = invT*invT;')

        # temperature check
        # self._write()
        # self._write(self.line('check the temperature value'))
        # self._write('if (T < %g || T > %g) {' % (lowT, highT))
        # self._indent()
        # self._write(
        #     'fprintf(stderr, "temperature %%g is outside the range (%g, %g)", T);'
        #     % (lowT, highT))
        # self._write('return;')
        # self._outdent()
        # self._write('}')
                    
        for midT, speciesList in midpoints.items():

            self._write('')
            self._write(self.line('species with midpoint at T=%g kelvin' % midT))
            self._write('if (T < %g) {' % midT)
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(lowRange.parameters)
                self._outdent()

            self._outdent()
            self._write('} else {')
            self._indent()

            for species, lowRange, highRange in speciesList:
                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                self._write('species[%d] =' % species.id)
                self._indent()
                expressionGenerator(highRange.parameters)
                self._outdent()

            self._outdent()
            self._write('}')
            
        self._write('return;')
        self._outdent()

        self._write('}')

        return


    def _miscTransInfo(self, KK, NLITE, do_declarations, NO=4):

        self._write()
        self._write()
        LENIMC = 4*KK+NLITE
        self._generateTransRoutineInteger(["egtransetLENIMC", "EGTRANSETLENIMC", "egtransetlenimc", "egtransetlenimc_", "LENIMC"], LENIMC, do_declarations)

        self._write()
        self._write()
        LENRMC = (19+2*NO+NO*NLITE)*KK+(15+NO)*KK**2
        self._generateTransRoutineInteger(["egtransetLENRMC", "EGTRANSETLENRMC", "egtransetlenrmc", "egtransetlenrmc_", "LENRMC"], LENRMC, do_declarations)

        self._write()
        self._write()
        self._generateTransRoutineInteger(["egtransetNO", "EGTRANSETNO", "egtransetno", "egtransetno_", "NO"], NO, do_declarations)

        self._write()
        self._write()
        self._generateTransRoutineInteger(["egtransetKK", "EGTRANSETKK", "egtransetkk", "egtransetkk_", "KK"], KK, do_declarations)

        self._write()
        self._write()
        self._generateTransRoutineInteger(["egtransetNLITE", "EGTRANSETNLITE", "egtransetnlite", "egtransetnlite_", "NLITE"], NLITE, do_declarations)

        self._write()
        self._write()
        self._write(self.line('Patm in ergs/cm3'))

        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetPATM EGTRANSETPATM')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetPATM egtransetpatm')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetPATM egtransetpatm_')
            self._write('#endif')

        self._write('void egtransetPATM(double* PATM) {')
        self._indent()
        self._write('*PATM =   0.1013250000000000E+07;}')
        self._outdent()

        return


    def _wt(self, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the molecular weights in g/mol'))

        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetWT EGTRANSETWT')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetWT egtransetwt')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetWT egtransetwt_')
            self._write('#endif')

        self._write('void %s(double* %s ) {' % ("egtransetWT", "WT"))
        self._indent()

        for species in self.species:
            self._write('%s[%d] = %.8E;' % ("WT", species.id, float(species.weight)))

        self._outdent()
        self._write('}')

        return


    def _eps(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the lennard-jones potential well depth eps/kb in K'))

        #i=0
        #expression=[]
        #for species in mechanism.species():
        #    expression[i] = float(species.trans[0].eps)
        #    i++
        self._generateTransRoutineSimple(mechanism, ["egtransetEPS", "EGTRANSETEPS", "egtranseteps", "egtranseteps_", "EPS"], 1, speciesTransport, do_declarations)

        return
    

    def _sig(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the lennard-jones collision diameter in Angstroms'))
        self._generateTransRoutineSimple(mechanism, ["egtransetSIG", "EGTRANSETSIG", "egtransetsig", "egtransetsig_", "SIG"], 2, speciesTransport, do_declarations)

        return


    def _dip(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the dipole moment in Debye'))
        self._generateTransRoutineSimple(mechanism, ["egtransetDIP", "EGTRANSETDIP", "egtransetdip", "egtransetdip_", "DIP"], 3, speciesTransport, do_declarations)

        return


    def _pol(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the polarizability in cubic Angstroms'))
        self._generateTransRoutineSimple(mechanism, ["egtransetPOL", "EGTRANSETPOL", "egtransetpol", "egtransetpol_", "POL"], 4, speciesTransport, do_declarations)

        return


    def _zrot(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('the rotational relaxation collision number at 298 K'))
        self._generateTransRoutineSimple(mechanism, ["egtransetZROT", "EGTRANSETZROT", "egtransetzrot", "egtransetzrot_", "ZROT"], 5, speciesTransport, do_declarations)

        return


    def _nlin(self, mechanism, speciesTransport, do_declarations):

        self._write()
        self._write()
        self._write(self.line('0: monoatomic, 1: linear, 2: nonlinear'))
        #self._generateTransRoutineSimple(["egtransetNLIN", "EGTRANSETNLIN", "egtransetNLIN", "egtransetNLIN_", "NLIN"], 0, speciesTransport)

        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetNLIN EGTRANSETNLIN')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetNLIN egtransetnlin')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetNLIN egtransetnlin_')
            self._write('#endif')

        self._write('void egtransetNLIN(int* NLIN) {')
        self._indent()

        for species in mechanism.species():
            self._write('%s[%d] = %d;' % ("NLIN", species.id, int(speciesTransport[species][0])))

        self._outdent()
        self._write('}')

        return


    def _viscosity(self, speciesTransport, do_declarations, NTFit):

        #compute single constants in g/cm/s
        kb = 1.3806503e-16
        Na = 6.02214199e23 
        RU = 8.31447e7
        #conversion coefs
        AtoCM = 1.0e-8
        DEBYEtoCGS = 1.0e-18
        #temperature increment  
        dt = (self.highT-self.lowT) / (NTFit-1)
        #factor dependent upon the molecule
        m_crot = np.zeros(self.nSpecies)
        m_cvib = np.zeros(self.nSpecies)
        isatm = np.zeros(self.nSpecies)
        for spec in speciesTransport:
            if int(speciesTransport[spec][0]) == 0:
                m_crot[spec.id] = 0.0
                m_cvib[spec.id] = 0.0
                isatm[spec.id] = 0.0
            elif int(speciesTransport[spec][0]) == 1:
                m_crot[spec.id] = 1.0
                m_cvib[spec.id] = 5.0 / 2.0
                isatm[spec.id] = 1.0
            else:
                m_crot[spec.id] = 1.5
                m_cvib[spec.id] = 3.0
                isatm[spec.id] = 1.0
        #viscosities coefs (4 per spec)
        cofeta = {}
        #conductivities coefs (4 per spec)
        coflam = {}
        for spec in speciesTransport:
            spvisc = []
            spcond = []
            tlog = []
            for n in range(NTFit):
                t = self.lowT + dt*n
                #variables
                #eq. (2)
                tr = t/ float(speciesTransport[spec][1])
                conversion = DEBYEtoCGS * DEBYEtoCGS / AtoCM / AtoCM / AtoCM / kb 
                dst = 0.5 * conversion * float(speciesTransport[spec][3])**2 / (float(speciesTransport[spec][1]) \
                        * float(speciesTransport[spec][2])**3)
                #viscosity of spec at t
                #eq. (1)
                conversion = AtoCM * AtoCM
                visc = (5.0 / 16.0) * np.sqrt(np.pi * self.species[spec.id].weight * kb * t / Na) / \
                    (self.om22_CHEMKIN(tr,dst) * np.pi * \
                    float(speciesTransport[spec][2]) * float(speciesTransport[spec][2]) * conversion)
                #conductivity of spec at t
                #eq. (30)
                conversion = AtoCM * AtoCM
                m_red = self.species[spec.id].weight / (2.0 * Na)
                diffcoef = (3.0 / 16.0) * np.sqrt(2.0 * np.pi * kb**3 * t**3 / m_red) /  \
                        (10.0 * np.pi * self.om11_CHEMKIN(tr,dst) * float(speciesTransport[spec][2]) * \
                        float(speciesTransport[spec][2]) * conversion)
                #eq. (19)
                cv_vib_R = (self._getCVdRspecies(t, spec) - m_cvib[spec.id]) * isatm[spec.id]
                rho_atm = 10.0 * self.species[spec.id].weight /(RU * t)
                f_vib = rho_atm * diffcoef / visc
                #eq. (20)
                A = 2.5 - f_vib
                #eqs. (21) + (32-33)
                cv_rot_R = m_crot[spec.id]
                #note: the T corr is not applied in CANTERA
                B = (float(speciesTransport[spec][5]) \
                        * self.Fcorr(298.0, float(speciesTransport[spec][1])) / self.Fcorr(t, float(speciesTransport[spec][1])) \
                        + (2.0 / np.pi) * ((5.0 / 3.0 ) * cv_rot_R  + f_vib))
                #eq. (18)
                f_rot = f_vib * (1.0 + 2.0 / np.pi * A / B )
                #eq. (17) 
                cv_trans_R = 3.0 / 2.0 
                f_trans = 5.0 / 2.0 * (1.0 - 2.0 / np.pi * A / B * cv_rot_R / cv_trans_R )
                if (int(speciesTransport[spec][0]) == 0):
                    cond = (visc * RU / self.species[spec.id].weight) * \
                            (5.0 / 2.0) * cv_trans_R
                else:
                    cond = (visc * RU / self.species[spec.id].weight) * \
                        (f_trans * cv_trans_R + f_rot * cv_rot_R + f_vib * cv_vib_R)

                #log transformation for polyfit
                tlog.append(np.log(t))
                spvisc.append(np.log(visc))
                spcond.append(np.log(cond))

            cofeta[spec.id] = np.polyfit(tlog, spvisc, 3)
            coflam[spec.id] = np.polyfit(tlog, spcond, 3)

        #header for visco
        self._write()
        self._write()
        self._write(self.line('Poly fits for the viscosities, dim NO*KK'))
        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetCOFETA EGTRANSETCOFETA')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetCOFETA egtransetcofeta')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetCOFETA egtransetcofeta_')
            self._write('#endif')

        #visco coefs
        self._write('void egtransetCOFETA(double* COFETA) {')
        self._indent()

        for spec in self.species:
            for i in range(4):
                self._write('%s[%d] = %.8E;' % ('COFETA', spec.id*4+i, cofeta[spec.id][3-i]))

        self._outdent()
        self._write('}')

        #header for cond
        self._write()
        self._write()
        self._write(self.line('Poly fits for the conductivities, dim NO*KK'))
        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetCOFLAM EGTRANSETCOFLAM')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetCOFLAM egtransetcoflam')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetCOFLAM egtransetcoflam_')
            self._write('#endif')

        #visco coefs
        self._write('void egtransetCOFLAM(double* COFLAM) {')

        self._indent()

        for spec in self.species:
            for i in range(4):
                self._write('%s[%d] = %.8E;' % ('COFLAM', spec.id*4+i, coflam[spec.id][3-i]))

        self._outdent()

        self._write('}')

        return


    def _thermaldiffratios(self, speciesTransport, lightSpecList, do_declarations, NTFit):

        # This is an overhaul of CHEMKIN version III
        #REORDERING OF SPECS
        specOrdered = []
        for i in range(self.nSpecies):
            for spec in speciesTransport:
                if spec.id == i:
                    specOrdered.append(spec)
                    break

        #compute single constants in g/cm/s
        kb = 1.3806503e-16
        #conversion coefs
        DEBYEtoCGS = 1.0e-18
        AtoCM = 1.0e-8
        #temperature increment  
        dt = (self.highT-self.lowT) / (NTFit-1)
        #diff ratios (4 per spec pair involving light species) 
        coftd = []
        k = -1
        for i,spec1 in enumerate(specOrdered):
            if (i != spec1.id):
                print "Problem in _thermaldiffratios computation"
                stop
            if spec1.id in lightSpecList:
                k = k + 1
                if (lightSpecList[k] != spec1.id):
                    print "Problem in  _thermaldiffratios computation"
                    stop
                coftd.append([])
                epsi = float(speciesTransport[spec1][1]) * kb
                sigi = float(speciesTransport[spec1][2]) * AtoCM
                poli = float(speciesTransport[spec1][4]) * AtoCM * AtoCM * AtoCM
                #eq. (12)
                poliRed = poli / sigi**3
                for j,spec2 in enumerate(specOrdered):
                    if (j != spec2.id):
                        print "Problem in _thermaldiffratios computation"
                        stop
                    #eq. (53)
                    Wji = (self.species[spec2.id].weight - self.species[spec1.id].weight) / \
                            (self.species[spec1.id].weight + self.species[spec2.id].weight) 
                    epsj = float(speciesTransport[spec2][1]) * kb
                    sigj = float(speciesTransport[spec2][2]) * AtoCM
                    dipj = float(speciesTransport[spec2][3]) * DEBYEtoCGS
                    #eq. (13)
                    dipjRed = dipj / np.sqrt(epsj*sigj**3)
                    epsRatio = epsj / epsi
                    tse = 1.0 + 0.25*poliRed*dipjRed**2*np.sqrt(epsRatio)
                    eok = tse**2 * np.sqrt(float(speciesTransport[spec1][1]) * float(speciesTransport[spec2][1]))
                    #enter the loop on temperature
                    spthdiffcoef = []
                    tTab = []
                    for n in range(NTFit):
                       t = self.lowT + dt*n
                       tslog = np.log(t) - np.log(eok)
                       #eq. (53)
                       thdifcoeff = 15.0 / 2.0 * Wji * (2.0 * self.astar(tslog) + 5.0) * (6.0 * self.cstar(tslog) - 5.0) / \
                               (self.astar(tslog) * (16.0 * self.astar(tslog) - 12.0 * self.bstar(tslog) + 55.0))

                       #log transformation for polyfit
                       tTab.append(t)
                       spthdiffcoef.append(thdifcoeff)

                    coftd[k].append(np.polyfit(tTab, spthdiffcoef, 3))

        #header for thermal diff ratios
        self._write()
        self._write()
        self._write(self.line('Poly fits for thermal diff ratios, dim NO*NLITE*KK'))
        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetCOFTD EGTRANSETCOFTD')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetCOFTD egtransetcoftd')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetCOFTD egtransetcoftd_')
            self._write('#endif')

        #visco coefs
        self._write('void egtransetCOFTD(double* COFTD) {')
        self._indent()

        for i in range(len(coftd)):
            for j in range(self.nSpecies):
                for k in range(4):
                    self._write('%s[%d] = %.8E;' % ('COFTD', i*4*self.nSpecies+j*4+k, coftd[i][j][3-k]))

        self._outdent()
        self._write('}')

        return


    def _diffcoefs(self, speciesTransport, do_declarations, NTFit) :

        #REORDERING OF SPECS
        specOrdered = []
        for i in range(self.nSpecies):
            for spec in speciesTransport:
                if spec.id == i:
                    specOrdered.append(spec)
                    break
        #checks
        #for spec in speciesTransport:
        #    print spec.symbol, spec.id
        #for i in range(self.nSpecies):
        #    print i, specOrdered[i].id, specOrdered[i].symbol
        #stop

        #compute single constants in g/cm/s
        kb = 1.3806503e-16
        Na = 6.02214199e23 
        #conversion coefs
        AtoCM = 1.0e-8
        DEBYEtoCGS = 1.0e-18
        PATM = 0.1013250000000000E+07
        #temperature increment  
        dt = (self.highT-self.lowT) / (NTFit-1)
        #diff coefs (4 per spec pair) 
        cofd = []
        for i,spec1 in enumerate(specOrdered):
            cofd.append([])
            if (i != spec1.id):
                print "Problem in _diffcoefs computation"
                stop
            for j,spec2 in enumerate(specOrdered[0:i+1]):
                if (j != spec2.id):
                    print "Problem in _diffcoefs computation"
                    stop
                #eq. (9)
                sigm = (0.5 * (float(speciesTransport[spec1][2]) + float(speciesTransport[spec2][2])) * AtoCM)\
                        * self.Xi(spec1, spec2, speciesTransport)**(1.0/6.0)
                #eq. (4)
                m_red = self.species[spec1.id].weight * self.species[spec2.id].weight / \
                        (self.species[spec1.id].weight + self.species[spec2.id].weight) / Na
                #eq. (8) & (14)
                epsm_k = np.sqrt(float(speciesTransport[spec1][1]) * float(speciesTransport[spec2][1])) \
                        * self.Xi(spec1, spec2, speciesTransport)**2.0

                #eq. (15)
                conversion = DEBYEtoCGS * DEBYEtoCGS / kb  
                dst = 0.5 * conversion * float(speciesTransport[spec1][3]) * float(speciesTransport[spec2][3]) / \
                    (epsm_k * sigm**3)
                if self.Xi_bool(spec1, spec2, speciesTransport)==False:
                    dst = 0.0
                #enter the loop on temperature
                spdiffcoef = []
                tlog = []
                for n in range(NTFit):
                   t = self.lowT + dt*n
                   tr = t/ epsm_k
                   #eq. (3)
                   #note: these are "corrected" in CHEMKIN not in CANTERA... we chose not to
                   difcoeff = 3.0 / 16.0 * 1 / PATM * (np.sqrt(2.0 * np.pi * t**3 * kb**3 / m_red) / \
                           ( np.pi * sigm * sigm * self.om11_CHEMKIN(tr,dst)))

                   #log transformation for polyfit
                   tlog.append(np.log(t))
                   spdiffcoef.append(np.log(difcoeff))

                cofd[i].append(np.polyfit(tlog, spdiffcoef, 3))

        #use the symmetry for upper triangular terms 
        #note: starting with this would be preferable (only one bigger loop)
        #note2: or write stuff differently !
        #for i,spec1 in enumerate(specOrdered):
        #    for j,spec2 in enumerate(specOrdered[i+1:]):
        #        cofd[i].append(cofd[spec2.id][spec1.id])

        #header for diffusion coefs
        self._write()
        self._write()
        self._write(self.line('Poly fits for the diffusion coefficients, dim NO*KK*KK'))
        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetCOFD EGTRANSETCOFD')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetCOFD egtransetcofd')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetCOFD egtransetcofd_')
            self._write('#endif')

        #coefs
        self._write('void egtransetCOFD(double* COFD) {')

        self._indent()

        for i,spec1 in enumerate(specOrdered):
            #for j,spec2 in enumerate(specOrdered):
            for j,spec2 in enumerate(specOrdered[0:i+1]):
                for k in range(4):
                    #self._write('%s[%d] = %.8E;' % ('COFD', i*self.nSpecies*4+j*4+k, cofd[j][i][3-k]))
                    self._write('%s[%d] = %.8E;' % ('COFD', i*self.nSpecies*4+j*4+k, cofd[i][j][3-k]))
            for j,spec2 in enumerate(specOrdered[i+1:]):
                for k in range(4):
                    self._write('%s[%d] = %.8E;' % ('COFD', i*self.nSpecies*4+(j+i+1)*4+k, cofd[j+i+1][i][3-k]))

        self._outdent()

        self._write('}')
        
        return


    def _lightSpecs(self, speclist, do_declarations):
        
        #header 
        self._write()
        self._write()
        self._write(self.line('List of specs with small weight, dim NLITE'))
        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define egtransetKTDIF EGTRANSETKTDIF')
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define egtransetKTDIF egtransetktdif')
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define egtransetKTDIF egtransetktdif_')
            self._write('#endif')

        #coefs
        self._write('void egtransetKTDIF(int* KTDIF) {')
        self._indent()

        for i in range(len(speclist)):
            self._write('%s[%d] = %d;' % ('KTDIF', i, speclist[i]+1))

        self._outdent()
        self._write('}')
        
        return


    def astar(self, tslog):

        aTab = [.1106910525E+01, -.7065517161E-02,-.1671975393E-01,
                .1188708609E-01,  .7569367323E-03,-.1313998345E-02,
                .1720853282E-03]

        B = aTab[6]
        for i in range(6):
            B = aTab[5-i] + B*tslog

        return B


    def bstar(self, tslog):

        bTab = [.1199673577E+01, -.1140928763E+00,-.2147636665E-02,
                .2512965407E-01, -.3030372973E-02,-.1445009039E-02,
                .2492954809E-03]

        B = bTab[6]
        for i in range(6):
            B = bTab[5-i] + B*tslog

        return B


    def cstar(self, tslog):

        cTab = [ .8386993788E+00,  .4748325276E-01, .3250097527E-01,
                -.1625859588E-01, -.2260153363E-02, .1844922811E-02,
                -.2115417788E-03]

        B = cTab[6]
        for i in range(6):
            B = cTab[5-i] + B*tslog

        return B


    def Xi(self, spec1, spec2, speciesTransport):

        dipmin = 1e-20
        #1 is polar, 2 is nonpolar
        #err in eq. (11) ?
        if (float(speciesTransport[spec2][3]) < dipmin) and (float(speciesTransport[spec1][3]) > dipmin):
            xi = 1.0 + 1.0/4.0 * self.redPol(spec2, speciesTransport)*self.redDip(spec1, speciesTransport) *\
                    self.redDip(spec1, speciesTransport) *\
                    np.sqrt(float(speciesTransport[spec1][1]) / float(speciesTransport[spec2][1]))
        #1 is nonpolar, 2 is polar
        elif (float(speciesTransport[spec2][3]) > dipmin) and (float(speciesTransport[spec1][3]) < dipmin):
            xi = 1.0 + 1.0/4.0 * self.redPol(spec1, speciesTransport)*self.redDip(spec2, speciesTransport) *\
                    self.redDip(spec2, speciesTransport) *\
                    np.sqrt(float(speciesTransport[spec2][1]) / float(speciesTransport[spec1][1]))
        #normal case, either both polar or both nonpolar
        else:
            xi = 1.0

        return xi


    def Xi_bool(self, spec1, spec2, speciesTransport):

        dipmin = 1e-20
        #1 is polar, 2 is nonpolar
        #err in eq. (11) ?
        if (float(speciesTransport[spec2][3]) < dipmin) and (float(speciesTransport[spec1][3]) > dipmin):
            xi_b = False
        #1 is nonpolar, 2 is polar
        elif (float(speciesTransport[spec2][3]) > dipmin) and (float(speciesTransport[spec1][3]) < dipmin):
            xi_b = False
        #normal case, either both polar or both nonpolar
        else:
            xi_b = True

        return xi_b


    def redPol(self, spec, speciesTransport): 

        return (float(speciesTransport[spec][4]) / float(speciesTransport[spec][2])**3.0)


    def redDip(self, spec, speciesTransport): 

        #compute single constants in g/cm/s
        kb = 1.3806503e-16
        #conversion coefs
        AtoCM = 1.0e-8
        DEBYEtoCGS = 1.0e-18
        convert = DEBYEtoCGS / np.sqrt( kb * AtoCM**3.0 )
        return convert * float(speciesTransport[spec][3]) / \
                np.sqrt(float(speciesTransport[spec][1]) * float(speciesTransport[spec][2])**3.0)


    def Fcorr(self, t, eps_k):

        thtwo = 3.0 / 2.0
        return 1 + np.pi**(thtwo) / 2.0 * np.sqrt(eps_k / t) + \
                (np.pi**2 / 4.0 + 2.0) * (eps_k / t) + \
                (np.pi * eps_k / t)**(thtwo)


    def om11(self, tr, dst):

        # This is an overhaul of CANTERA version 2.3
        #range of dst
        dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

        #range of tr
        trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

        #tab of astar corresp. to (tr, dst)
        #CANTERA
        astarTab = [1.0065, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840,
                1.0231, 1.0660, 1.0380, 1.0400, 1.0430, 1.0500, 1.0520, 1.0510,
                1.0424, 1.0450, 1.0480, 1.0520, 1.0560, 1.0650, 1.0660, 1.0640,
                1.0719, 1.0670, 1.0600, 1.0550, 1.0580, 1.0680, 1.0710, 1.0710,
                1.0936, 1.0870, 1.0770, 1.0690, 1.0680, 1.0750, 1.0780, 1.0780,
                1.1053, 1.0980, 1.0880, 1.0800, 1.0780, 1.0820, 1.0840, 1.0840,
                1.1104, 1.1040, 1.0960, 1.0890, 1.0860, 1.0890, 1.0900, 1.0900,
                1.1114, 1.1070, 1.1000, 1.0950, 1.0930, 1.0950, 1.0960, 1.0950,
                1.1104, 1.1070, 1.1020, 1.0990, 1.0980, 1.1000, 1.1000, 1.0990,
                1.1086, 1.1060, 1.1020, 1.1010, 1.1010, 1.1050, 1.1050, 1.1040,
                1.1063, 1.1040, 1.1030, 1.1030, 1.1040, 1.1080, 1.1090, 1.1080,
                1.1020, 1.1020, 1.1030, 1.1050, 1.1070, 1.1120, 1.1150, 1.1150,
                1.0985, 1.0990, 1.1010, 1.1040, 1.1080, 1.1150, 1.1190, 1.1200,
                1.0960, 1.0960, 1.0990, 1.1030, 1.1080, 1.1160, 1.1210, 1.1240,
                1.0943, 1.0950, 1.0990, 1.1020, 1.1080, 1.1170, 1.1230, 1.1260,
                1.0934, 1.0940, 1.0970, 1.1020, 1.1070, 1.1160, 1.1230, 1.1280,
                1.0926, 1.0940, 1.0970, 1.0990, 1.1050, 1.1150, 1.1230, 1.1300,
                1.0934, 1.0950, 1.0970, 1.0990, 1.1040, 1.1130, 1.1220, 1.1290,
                1.0948, 1.0960, 1.0980, 1.1000, 1.1030, 1.1120, 1.1190, 1.1270,
                1.0965, 1.0970, 1.0990, 1.1010, 1.1040, 1.1100, 1.1180, 1.1260,
                1.0997, 1.1000, 1.1010, 1.1020, 1.1050, 1.1100, 1.1160, 1.1230,
                1.1025, 1.1030, 1.1040, 1.1050, 1.1060, 1.1100, 1.1150, 1.1210,
                1.1050, 1.1050, 1.1060, 1.1070, 1.1080, 1.1110, 1.1150, 1.1200,
                1.1072, 1.1070, 1.1080, 1.1080, 1.1090, 1.1120, 1.1150, 1.1190,
                1.1091, 1.1090, 1.1090, 1.1100, 1.1110, 1.1130, 1.1150, 1.1190,
                1.1107, 1.1110, 1.1110, 1.1110, 1.1120, 1.1140, 1.1160, 1.1190,
                1.1133, 1.1140, 1.1130, 1.1140, 1.1140, 1.1150, 1.1170, 1.1190,
                1.1154, 1.1150, 1.1160, 1.1160, 1.1160, 1.1170, 1.1180, 1.1200,
                1.1172, 1.1170, 1.1170, 1.1180, 1.1180, 1.1180, 1.1190, 1.1200,
                1.1186, 1.1190, 1.1190, 1.1190, 1.1190, 1.1190, 1.1200, 1.1210,
                1.1199, 1.1200, 1.1200, 1.1200, 1.1200, 1.1210, 1.1210, 1.1220,
                1.1223, 1.1220, 1.1220, 1.1220, 1.1220, 1.1230, 1.1230, 1.1240,
                1.1243, 1.1240, 1.1240, 1.1240, 1.1240, 1.1240, 1.1250, 1.1250,
                1.1259, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260,
                1.1273, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1280,
                1.1297, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1290,
                1.1339, 1.1340, 1.1340, 1.1350, 1.1350, 1.1340, 1.1340, 1.1320,
                1.1364, 1.1370, 1.1370, 1.1380, 1.1390, 1.1380, 1.1370, 1.1350,
                1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187,
                1.14187]


        #Find for each fixed tr the poly of deg 6 in dst approx astar values
        #store the poly coefs in m_apoly
        m_apoly = []
        for i in range(37):
            dstDeg = 6
            #Polynomial coefficients, highest power first 
            polycoefs = np.polyfit(dstTab,astarTab[8*(i+1):8*(i+2)],dstDeg)
            m_apoly.append(polycoefs)

        #Find 3 referenced temp points around tr
        for i in range(37):
            if tr<trTab[i]:
                break
        i1 = max(i-1, 0)
        i2 = i1+3
        if (i2 > 36):
            i2 = 36
            i1 = i2 - 3
        #compute astar value for these 3 points
        values = []
        for j in range(i1,i2):
            if (dst == 0.0):
                values.append(astarTab[8*(j+1)])
            else:
                poly6 = np.poly1d(m_apoly[j])
                values.append(poly6(dst))

        #interpolate to find real tr value
        trTab_log = []
        for j in range(len(trTab)):
            trTab_log.append(np.log(trTab[j]))

        astar_interp = self.quadInterp(np.log(tr), trTab_log[i1:i2], values)

        return self.om22(tr,dst)/astar_interp


    def om11_CHEMKIN(self, tr, dst):

        # This is an overhaul of CANTERA version 2.3
        #range of dst
        dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

        #range of tr
        trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

        #tab of omega11 corresp. to (tr, dst)
        #CANTERA
        omegaTab = [4.008, 4.002, 4.655, 5.52, 6.454, 8.214, 9.824, 11.31,
                3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160,
                2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483 ,
                2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539 ,
                2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918 ,
                1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747 ,
                1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137 ,
                1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872 ,
                1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657 ,
                1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780,
                1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199,
                1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99,
                1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827,
                1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7,
                1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592,
                1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397,
                 .95  ,  .9538,  .9656,  .9852, 1.012 , 1.082 , 1.168 , 1.265,
                 .9131,  .9162,  .9256,  .9413,  .9626, 1.019 , 1.09  , 1.17,
                 .8845,  .8871,  .8948,  .9076,  .9252,  .972 , 1.03  , 1.098,
                 .8428,  .8446,  .850 ,  .859 ,  .8716,  .9053,  .9483,  .9984,
                 .813 ,  .8142,  .8183,  .825 ,  .8344,  .8598,  .8927,  .9316,
                 .7898,  .791 ,  .794 ,  .7993,  .8066,  .8265,  .8526,  .8836,
                 .7711,  .772 ,  .7745,  .7788,  .7846,  .8007,  .822 ,  .8474,
                 .7555,  .7562,  .7584,  .7619,  .7667,  .78  ,  .7976,  .8189,
                 .7422,  .743 ,  .7446,  .7475,  .7515,  .7627,  .7776,  .796 ,
                 .72022, .7206,  .722 ,  .7241,  .7271,  .7354,  .7464,  .76  ,
                 .7025,  .703 ,  .704 ,  .7055,  .7078,  .7142,  .7228,  .7334,
                 .68776, .688,   .6888,  .6901,  .6919,  .697 ,  .704 ,  .7125,
                 .6751,  .6753,  .676 ,  .677 ,  .6785,  .6827,  .6884,  .6955,
                 .664 ,  .6642,  .6648,  .6657,  .6669,  .6704,  .6752,  .681,
                 .6414,  .6415,  .6418,  .6425,  .6433,  .6457,  .649 ,  .653,
                 .6235,  .6236,  .6239,  .6243,  .6249,  .6267,  .629 ,  .632,
                 .60882, .6089,  .6091,  .6094,  .61  ,  .6112,  .613 ,  .6154,
                 .5964,  .5964,  .5966,  .597 ,  .5972,  .5983,  .600 ,  .6017,
                 .5763,  .5763,  .5764,  .5766,  .5768,  .5775,  .5785,  .58,
                 .5415,  .5415,  .5416,  .5416,  .5418,  .542 ,  .5424,  .543,
                 .518 ,  .518 ,  .5182,  .5184,  .5184,  .5185,  .5186,  .5187]

        #First test on tr
        if (tr > 75.0):
            omeg12 = 0.623 - 0.136e-2*tr + 0.346e-5*tr*tr - 0.343e-8*tr*tr*tr
        else:
            #Find tr idx in trTab
            if (tr <= 0.2):
                ii = 1
            else:
                ii = 36
            for i in range(1,37):
                if (tr > trTab[i-1]) and (tr <= trTab[i]):
                    ii = i
                    break
            #Find dst idx in dstTab 
            if (abs(dst) >= 1.0e-5):
                if (dst <= 0.25):
                    kk = 1
                else:
                    kk = 6
                for i in range(1,7):
                    if (dstTab[i-1] < dst) and (dstTab[i] >= dst):
                        kk = i
                        break
                #Find surrounding values and interpolate
                #First on dst
                vert = np.zeros(3) 
                for i in range(3):
                    arg = np.zeros(3) 
                    val = np.zeros(3) 
                    for k in range(3):
                      arg[k] = dstTab[kk-1+k]
                      val[k] = omegaTab[8*(ii-1+i) + (kk-1+k)]
                    vert[i] = self.qinterp(dst, arg, val)
                #Second on tr
                arg = np.zeros(3) 
                for i in range(3):
                   arg[i] = trTab[ii-1+i]
                omeg12 = self.qinterp(tr, arg, vert)
            else:
                arg = np.zeros(3) 
                val = np.zeros(3) 
                for i in range(3):
                   arg[i] = trTab[ii-1+i]
                   val[i] = omegaTab[8*(ii-1+i)]
                omeg12 =self. qinterp(tr, arg, val)

        return omeg12


    def om22_CHEMKIN(self, tr, dst):

        # This is an overhaul of CANTERA version 2.3
        #range of dst
        dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

        #range of tr
        trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

        #tab of omega22 corresp. to (tr, dst)
        #CANTERA
        omegaTab = [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89,
                3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618,
                2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874,
                2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895,
                2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249,
                2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786,
                1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435,
                1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156,
                1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933,
                1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746,
                1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451,
                1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228,
                1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053,
                1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912,
                1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795,
                1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578,
                1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428,
                0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319,
                0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236,
                0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121,
                0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044,
                0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
                0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
                0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
                0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
                0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
                0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
                0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
                0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
                0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642,
                0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339,
                0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
                0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
                0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
                0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
                0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
                0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885]

        #First test on tr
        if (tr > 75.0):
            omeg12 = 0.703 - 0.146e-2*tr + 0.357e-5*tr*tr - 0.343e-8*tr*tr*tr
        else:
            #Find tr idx in trTab
            if (tr <= 0.2):
                ii = 1
            else:
                ii = 36
            for i in range(1,37):
                if (tr > trTab[i-1]) and (tr <= trTab[i]):
                    ii = i
                    break
            #Find dst idx in dstTab 
            if (abs(dst) >= 1.0e-5):
                if (dst <= 0.25):
                    kk = 1
                else:
                    kk = 6
                for i in range(1,7):
                    if (dstTab[i-1] < dst) and (dstTab[i] >= dst):
                        kk = i
                        break
                #Find surrounding values and interpolate
                #First on dst
                vert = np.zeros(3) 
                for i in range(3):
                    arg = np.zeros(3) 
                    val = np.zeros(3) 
                    for k in range(3):
                      arg[k] = dstTab[kk-1+k]
                      val[k] = omegaTab[8*(ii-1+i) + (kk-1+k)]
                    vert[i] = self.qinterp(dst, arg, val)
                #Second on tr
                arg = np.zeros(3) 
                for i in range(3):
                   arg[i] = trTab[ii-1+i]
                omeg12 = self.qinterp(tr, arg, vert)
            else:
                arg = np.zeros(3) 
                val = np.zeros(3) 
                for i in range(3):
                   arg[i] = trTab[ii-1+i]
                   val[i] = omegaTab[8*(ii-1+i)]
                omeg12 =self. qinterp(tr, arg, val)

        return omeg12


    def om22(self, tr, dst):

        # This is an overhaul of CANTERA version 2.3
        #range of dst
        dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

        #range of tr
        trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

        #tab of omega22 corresp. to (tr, dst)
        #CANTERA
        omegaTab = [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89,
                3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618,
                2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874,
                2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895,
                2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249,
                2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786,
                1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435,
                1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156,
                1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933,
                1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746,
                1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451,
                1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228,
                1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053,
                1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912,
                1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795,
                1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578,
                1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428,
                0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319,
                0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236,
                0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121,
                0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044,
                0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
                0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
                0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
                0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
                0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
                0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
                0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
                0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
                0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642,
                0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339,
                0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
                0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
                0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
                0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
                0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
                0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885]


        #Find for each fixed tr the poly of deg 6 in dst approx omega22 values
        #store the poly coefs in m_o22poly
        m_o22poly = []
        for i in range(37):
            dstDeg = 6
            #Polynomial coefficients, highest power first 
            polycoefs = np.polyfit(dstTab,omegaTab[8*i:8*(i+1)],dstDeg)
            m_o22poly.append(polycoefs)

        #Find 3 referenced temp points around tr
        for i in range(37):
            if tr<trTab[i]:
                break
        i1 = max(i-1, 0)
        i2 = i1+3
        if (i2 > 36):
            i2 = 36
            i1 = i2 - 3
        #compute omega22 value for these 3 points
        values = []
        for j in range(i1,i2):
            if (dst == 0.0):
                values.append(omegaTab[8*j])
            else:
                poly6 = np.poly1d(m_o22poly[j])
                values.append(poly6(dst))

        #interpolate to find real tr value
        trTab_log = []
        for j in range(len(trTab)):
            trTab_log.append(np.log(trTab[j]))
        #print trTab_log[i1:i2], values
        om22_interp = self.quadInterp(np.log(tr), trTab_log[i1:i2], values)

        return om22_interp


    def qinterp(self, x0, x, y):

        val1 = y[0] + (x0-x[0])*(y[1]-y[0]) / (x[1]-x[0])
        val2 = y[1] + (x0-x[1])*(y[2]-y[1]) / (x[2]-x[1])
        fac1 = (x0-x[0]) / (x[1]-x[0]) / 2.0
        fac2 = (x[2]-x0) / (x[2]-x[1]) / 2.0
        if (x0 >= x[1]):
           val = (val1*fac2+val2) / (1.0+fac2)
        else:
           val = (val1+val2*fac1) / (1.0+fac1)
        return val


    def quadInterp(self, x0, x, y):

        dx21 = x[1] - x[0] 
        dx32 = x[2] - x[1]
        dx31 = dx21 + dx32
        dy32 = y[2] - y[1]
        dy21 = y[1] - y[0]
        a = (dx21*dy32 - dy21*dx32)/(dx21*dx31*dx32)
        return a*(x0 - x[0])*(x0 - x[1]) + (dy21/dx21)*(x0 - x[1]) + y[1]


    def _generateTransRoutineSimple(self, mechanism, nametab, id, speciesTransport, do_declarations):

        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define %s %s' % (nametab[0], nametab[1]))
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define %s %s' % (nametab[0], nametab[2]))
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define %s %s' % (nametab[0], nametab[3]))
            self._write('#endif')

        self._write('void %s(double* %s ) {' % (nametab[0], nametab[4]))
        self._indent()

        for species in mechanism.species():
            self._write('%s[%d] = %.8E;' % (nametab[4], species.id, float(speciesTransport[species][id])))

        self._outdent()
        self._write('}')

        return

    def _generateTransRoutineInteger(self, nametab, expression, do_declarations):

        if (do_declarations):
            self._write('#if defined(BL_FORT_USE_UPPERCASE)')
            self._write('#define %s %s' % (nametab[0], nametab[1]))
            self._write('#elif defined(BL_FORT_USE_LOWERCASE)')
            self._write('#define %s %s' % (nametab[0], nametab[2]))
            self._write('#elif defined(BL_FORT_USE_UNDERSCORE)')
            self._write('#define %s %s' % (nametab[0], nametab[3]))
            self._write('#endif')

        self._write('void %s(int* %s ) {' % (nametab[0], nametab[4]))
        self._indent()

        self._write('*%s = %d;}' % (nametab[4], expression ))
        self._outdent()

        return


    def _getCVdRspecies(self, t, species):

        models = species.thermo
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

        if t < mid:
            parameters = lowRange.parameters
        else:
            parameters = highRange.parameters

        return ((parameters[0] - 1.0) + parameters[1] * t + parameters[2] * t * t \
                + parameters[3] * t * t * t + parameters[4] * t * t * t * t)


    def _analyzeThermodynamics(self, mechanism):
        lowT = 0.0
        highT = 1000000.0

        midpoints = {}

        for species in mechanism.species():

            models = species.thermo
            if len(models) > 2:
                print 'species: ', species
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
        
        self.lowT = lowT
        self.highT = highT
        return lowT, highT, midpoints
    

    def _analyzeTransport(self, mechanism):

        transdata = {}

        for species in mechanism.species():

            models = species.trans
            if len(models) > 2:
                print 'species: ', species
                import pyre
                pyre.debug.Firewall.hit("unsupported configuration in species.trans")
                return
            
            m1 = models[0]

            lin = m1.parameters[0]
            eps = m1.eps
            sig = m1.sig
            dip = m1.dip
            pol = m1.pol
            zrot = m1.zrot

            transdata[species] = [lin, eps, sig, dip, pol, zrot]

        return transdata


    def _Kc(self, mechanism, reaction):

        dim = 0
        dG = ""

        terms = []
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
                    
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim -= coefficient
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in reaction.products:
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
            dim += coefficient
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        if dim == 0:
            conversion = ""
        elif dim > 0:
            if (dim == 1.0):
                conversion = "*".join(["refC"]) + ' * '
            else:
                conversion = "*".join(["pow(refC,%f)" % dim]) + ' * '
        else:
            if (dim == -1.0):
                conversion = "1.0 / (" + "*".join(["refC"]) + ') * '
            else:
                conversion = "1.0 / (" + "*".join(["pow(refC,%f)" % abs(dim)]) + ') * '

        K_c = conversion + K_p

        return K_c

    def _KcConv(self, mechanism, reaction):
        dim = 0
        for symbol, coefficient in reaction.reactants:
            dim -= coefficient
        # flip the signs
        for symbol, coefficient in reaction.products:
            dim += coefficient

        if dim == 0:
            conversion = ""
        elif dim > 0:
            if (dim == 1.0):
                conversion = "*".join(["refC"])
            else:
                conversion = "*".join(["pow(refC,%f)" % dim])
        else:
            if (dim == -1.0):
                conversion = "*".join(["refCinv"])
            else:
                conversion = "*".join(["pow(refCinv,%f)" % abs(dim)])

        return conversion

    def _sortedKcExpArg(self, mechanism, reaction):

        nSpecies = len(mechanism.species())

        terms = []
        for i in range(nSpecies):
            terms.append('')
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1.0:
                factor = " + "
            else:
                factor = " + %f*" % coefficient
            i = mechanism.species(symbol).id
            terms[i] += "%sg_RT[%d]"%(factor,i)

        for symbol, coefficient in reaction.products:
            if coefficient == 1.0:
                factor = " - "    # flip the signs
            else:
                factor = " - %f*" % coefficient
            i = mechanism.species(symbol).id
            terms[i] += "%sg_RT[%d]"%(factor,i)

        dG = ""
        for i in range(nSpecies):
            if terms[i]:
                dG += terms[i]
        if dG[0:3] == " + ":
            return dG[3:]
        else:
            return "-"+dG[3:]


    def _sortedKc(self, mechanism, reaction):
        conv = self._KcConv(mechanism,reaction)
        exparg = self._sortedKcExpArg(mechanism,reaction)
        if conv:
            return conv + ' * exp('+exparg+')'
        else:
            return 'exp('+exparg+')'


    def _vKc(self, mechanism, reaction):

        dim = 0
        dG = ""

        terms = []
        for symbol, coefficient in sorted(reaction.reactants, 
                                          key=lambda x: mechanism.species(x[0]).id):
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
                    
            terms.append("%sg_RT[%d*npt+i]" % (factor, mechanism.species(symbol).id))
            dim -= coefficient
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in sorted(reaction.products,
                                          key=lambda x: mechanism.species(x[0]).id):
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
            terms.append("%sg_RT[%d*npt+i]" % (factor, mechanism.species(symbol).id))
            dim += coefficient
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        if dim == 0:
            conversion = ""
        elif dim > 0:
            if (dim == 1.0):
                conversion = "*".join(["refC"]) + ' * '
            else:
                conversion = "*".join(["pow(refC,%f)" % dim]) + ' * '
        else:
            if (dim == -1.0):
                conversion = "*".join(["refCinv"]) + ' * '
            else:
                conversion = "*".join(["pow(refCinv,%f)" % abs(dim)]) + ' * '

        K_c = conversion + K_p

        return K_c


    def _Kc_exparg(self, mechanism, reaction):

        dG = ""

        terms = []
        for symbol, coefficient in reaction.reactants:
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
                    
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
        dG += '(' + ' + '.join(terms) + ')'

        # flip the signs
        terms = []
        for symbol, coefficient in reaction.products:
            if coefficient == 1.0:
                factor = ""
            else:
                factor = "%f * " % coefficient
            terms.append("%sg_RT[%d]" % (factor, mechanism.species(symbol).id))
        dG += ' - (' + ' + '.join(terms) + ')'

        K_p = 'exp(' + dG + ')'

        return dG

    def _cpNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4];' % parameters[4])
        return

    def _dcpdTNASA(self, parameters):
        self._write('%+15.8e' % parameters[1])
        self._write('%+15.8e * tc[1]' % (parameters[2]*2.))
        self._write('%+15.8e * tc[2]' % (parameters[3]*3.))
        self._write('%+15.8e * tc[3];' % (parameters[4]*4.))
        return

    def _cvNASA(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % parameters[1])
        self._write('%+15.8e * tc[2]' % parameters[2])
        self._write('%+15.8e * tc[3]' % parameters[3])
        self._write('%+15.8e * tc[4];' % parameters[4])
        return


    def _enthalpyNASA(self, parameters):
        self._write('%+15.8e' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e * invT;' % (parameters[5]))
        return


    def _internalEnergy(self, parameters):
        self._write('%+15.8e' % (parameters[0] - 1.0))
        self._write('%+15.8e * tc[1]' % (parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (parameters[2]/3))
        self._write('%+15.8e * tc[3]' % (parameters[3]/4))
        self._write('%+15.8e * tc[4]' % (parameters[4]/5))
        self._write('%+15.8e * invT;' % (parameters[5]))
        return

    
    def _gibbsNASA(self, parameters):
        self._write('%+20.15e * invT' % parameters[5])
        self._write('%+20.15e' % (parameters[0] - parameters[6]))
        self._write('%+20.15e * tc[0]' % (-parameters[0]))
        self._write('%+20.15e * tc[1]' % (-parameters[1]/2))
        self._write('%+20.15e * tc[2]' % (-parameters[2]/6))
        self._write('%+20.15e * tc[3]' % (-parameters[3]/12))
        self._write('%+20.15e * tc[4];' % (-parameters[4]/20))
        return
    
    def _helmholtzNASA(self, parameters):
        self._write('%+15.8e * invT' % parameters[5])
        self._write('%+15.8e' % (parameters[0] - parameters[6] - 1.0))
        self._write('%+15.8e * tc[0]' % (-parameters[0]))
        self._write('%+15.8e * tc[1]' % (-parameters[1]/2))
        self._write('%+15.8e * tc[2]' % (-parameters[2]/6))
        self._write('%+15.8e * tc[3]' % (-parameters[3]/12))
        self._write('%+15.8e * tc[4];' % (-parameters[4]/20))
        return

    def _entropyNASA(self, parameters):
        self._write('%+15.8e * tc[0]' % parameters[0])
        self._write('%+15.8e * tc[1]' % (parameters[1]))
        self._write('%+15.8e * tc[2]' % (parameters[2]/2))
        self._write('%+15.8e * tc[3]' % (parameters[3]/3))
        self._write('%+15.8e * tc[4]' % (parameters[4]/4))
        self._write('%+15.8e ;' % (parameters[6]))
        return

    def _T_given_ey(self, mechanism):
        self._write()
        self._write(self.line(' get temperature given internal energy in mass units and mass fracs'))
        self._write('AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)')
        self._write('{')
        self._write('#ifdef CONVERGENCE')
        self._indent()
        self._write('const int maxiter = 5000;')
        self._write('const double tol  = 1.e-12;')
        self._outdent()
        self._write('#else')
        self._indent()
        self._write('const int maxiter = 200;')
        self._write('const double tol  = 1.e-6;')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('double ein  = *e;')
        self._write('double tmin = 90;'+self.line('max lower bound for thermo def'))
        self._write('double tmax = 4000;'+self.line('min upper bound for thermo def'))
        self._write('double e1,emin,emax,cv,t1,dt;')
        self._write('int i;'+self.line(' loop counter'))
        self._write('CKUBMS(&tmin, y, &emin);')
        self._write('CKUBMS(&tmax, y, &emax);')
        self._write('if (ein < emin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCVBS(&tmin, y, &cv);')
        self._write('*t = tmin - (emin-ein)/cv;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('if (ein > emax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCVBS(&tmax, y, &cv);')
        self._write('*t = tmax - (emax-ein)/cv;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('t1 = *t;')
        self._write('if (t1 < tmin || t1 > tmax) {')
        self._indent()
        self._write('t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);')
        self._outdent()
        self._write('}')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKUBMS(&t1,y,&e1);')
        self._write('CKCVBS(&t1,y,&cv);')
        self._write('dt = (ein - e1) / cv;')
        self._write('if (dt > 100.) { dt = 100.; }')
        self._write('else if (dt < -100.) { dt = -100.; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('else if (t1+dt == t1) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        self._write('*t = t1;')
        self._write('*ierr = 0;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write()

    def _T_given_hy(self, mechanism):
        self._write(self.line(' get temperature given enthalpy in mass units and mass fracs'))
        self._write('AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)')
        self._write('{')
        self._write('#ifdef CONVERGENCE')
        self._indent()
        self._write('const int maxiter = 5000;')
        self._write('const double tol  = 1.e-12;')
        self._outdent()
        self._write('#else')
        self._indent()
        self._write('const int maxiter = 200;')
        self._write('const double tol  = 1.e-6;')
        self._outdent()
        self._write('#endif')
        self._indent()
        self._write('double hin  = *h;')
        self._write('double tmin = 90;'+self.line('max lower bound for thermo def'))
        self._write('double tmax = 4000;'+self.line('min upper bound for thermo def'))
        self._write('double h1,hmin,hmax,cp,t1,dt;')
        self._write('int i;'+self.line(' loop counter'))
        self._write('CKHBMS(&tmin, y, &hmin);')
        self._write('CKHBMS(&tmax, y, &hmax);')
        self._write('if (hin < hmin) {')
        self._indent()
        self._write(self.line('Linear Extrapolation below tmin'))
        self._write('CKCPBS(&tmin, y, &cp);')
        self._write('*t = tmin - (hmin-hin)/cp;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('if (hin > hmax) {')
        self._indent()
        self._write(self.line('Linear Extrapolation above tmax'))
        self._write('CKCPBS(&tmax, y, &cp);')
        self._write('*t = tmax - (hmax-hin)/cp;')
        self._write('*ierr = 1;')
        self._write('return;')
        self._outdent()
        self._write('}')
        self._write('t1 = *t;')
        self._write('if (t1 < tmin || t1 > tmax) {')
        self._indent()
        self._write('t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);')
        self._outdent()
        self._write('}')
        self._write('for (i = 0; i < maxiter; ++i) {')
        self._indent()
        self._write('CKHBMS(&t1,y,&h1);')
        self._write('CKCPBS(&t1,y,&cp);')
        self._write('dt = (hin - h1) / cp;')
        self._write('if (dt > 100.) { dt = 100.; }')
        self._write('else if (dt < -100.) { dt = -100.; }')
        self._write('else if (fabs(dt) < tol) break;')
        self._write('else if (t1+dt == t1) break;')
        self._write('t1 += dt;')
        self._outdent()
        self._write('}')
        self._write('*t = t1;')
        self._write('*ierr = 0;')
        self._write('return;')
        self._outdent()
        self._write('}')

    def _emptygjs(self, mechanism):
        self._write()
        self._write(self.line(' Replace this routine with the one generated by the Gauss Jordan solver of DW'))
        self._write('AMREX_GPU_HOST_DEVICE void sgjsolve(double* A, double* x, double* b) {')
        self._indent()
        self._write('amrex::Abort("sgjsolve not implemented, choose a different solver ");')
        self._outdent()
        self._write('}')

        self._write()
        self._write(self.line(' Replace this routine with the one generated by the Gauss Jordan solver of DW'))
        self._write('AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(double* A, double* x, double* b) {')
        self._indent()
        self._write('amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");')
        self._outdent()
        self._write('}')


    ####################
    #unused
    ####################
    def _end(self):
        self._timestamp()
        self._rep += self.footer()
        return


# version
__id__ = "$Id$"

#  End of file 
