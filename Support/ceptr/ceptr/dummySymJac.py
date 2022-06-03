import sympy as smp


T, A, Ea, R = smp.symbols('T A Ea R')
kf = A * smp.exp (-Ea/(R*T))
dkfdT = smp.diff (kf, T)
print(dkfdT)


Aarr = smp.symbols('Aarr:198')
Eaarr = smp.symbols('Eaarr:198')

kffromarr = Aarr[87] * smp.exp(-Eaarr[90]/(R*T))
dkffromarrdT = smp.diff (kffromarr, T)
print(dkffromarrdT)


A = smp.symbols('A:198')
Ea = smp.symbols('Ea:198')

kffromarr = A[87] * smp.exp(-Ea[90]/(R*T))
dkffromarrdT = smp.diff (kffromarr, T)
print(dkffromarrdT)
