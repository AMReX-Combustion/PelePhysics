"""Symbolic math for symbolic differentiation."""
from collections import OrderedDict
import sympy as smp
import re

def convertToCPP(symStr):
    cppStr = symStr
    # Fix indices
    symbReplace = re.findall(r'_smp\d+', cppStr)
    for symb in symbReplace:
        ind = symb[4:]
        cppStr = cppStr.replace(symb,"["+ind+"]",1)
    # Fix superfluous _smp
    symbReplace = re.findall(r'_smp', cppStr)
    for symb in symbReplace:
        cppStr = cppStr.replace(symb,"",1)
    # Fix power
    #symbReplace = re.findall(r'\*\*\d+', cppStr)
    return cppStr

class SymbolicMath:
    """Symbols to carry throughout operations."""

    def __init__(self,species_info,reaction_info):
        self.T_smp = smp.symbols('T_smp')
        self.tc_smp = [ smp.log(self.T_smp), self.T_smp, self.T_smp**2, self.T_smp**3, self.T_smp**4 ]
        self.invT_smp = 1.0 / self.tc_smp[1]
        
        n_species = species_info.n_species
        n_qssa_species = species_info.n_qssa_species 
        
        self.sc_smp = smp.symbols('sc_smp:'+str(n_species))
        self.g_RT_smp = [ smp.symbols('g_RT_smp'+str(i)) for i in range(n_species)]
        self.h_RT_smp = [ smp.symbols('h_RT_smp'+str(i)) for i in range(n_species)]
        self.wdot_smp = [ smp.symbols('wdot_smp'+str(i)) for i in range(n_species)]
 
        if n_qssa_species > 0:

            n_qssa_reactions = reaction_info.n_qssa_reactions
            self.sc_qss_smp = [ smp.symbols('sc_qss_smp'+str(i)) for i in range(n_qssa_species)]
            self.kf_qss_smp = [ smp.symbols('kf_qss_smp'+str(i)) for i in range(n_qssa_reactions)]
            self.qf_qss_smp = [ smp.symbols('qf_qss_smp'+str(i)) for i in range(n_qssa_reactions)]
            self.qr_qss_smp = [ smp.symbols('qr_qss_smp'+str(i)) for i in range(n_qssa_reactions)]
            self.g_RT_qss_smp = [ smp.symbols('g_RT_qss_smp'+str(i)) for i in range(n_qssa_species)] 
            self.h_RT_qss_smp = [ smp.symbols('h_RT_qss_smp'+str(i)) for i in range(n_qssa_species)]
    

