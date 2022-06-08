"""Test function for the SymJac functionality."""
import sympy as smp

T, A, Ea, R = smp.symbols("T A Ea R")
kf = A * smp.exp(-Ea / (R * T))
dkf_dt = smp.diff(kf, T)
print(dkf_dt)


Aarr = smp.symbols("Aarr:198")
Eaarr = smp.symbols("Eaarr:198")

ind1 = 87
ind2 = 90

kffromarr = Aarr[ind1] * smp.exp(-Eaarr[ind2] / (R * T))
dkffromarrdt = smp.diff(kffromarr, T)
print(dkffromarrdt)

len_arr = 198

A = smp.symbols("A:" + str(len_arr))
Ea = smp.symbols("Ea:" + str(len_arr))
kffromarr = 1.0

kffromarr = smp.symbols("kffromarr")
kffromarr = A[87] * smp.exp(-Ea[90] / (R * T))
dkffromarrdt = smp.diff(kffromarr, T)
print(dkffromarrdt)
