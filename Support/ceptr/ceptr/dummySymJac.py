import sympy as smp


T, A, Ea, R = smp.symbols('T A Ea R')
kf = A * smp.exp (-Ea/(R*T))
dkfdT = smp.diff (kf, T)
print(dkfdT)


Aarr = smp.symbols('Aarr:198')
Eaarr = smp.symbols('Eaarr:198')

ind1 = 87
ind2 = 90

kffromarr = Aarr[ind1] * smp.exp(-Eaarr[ind2]/(R*T))
dkffromarrdT = smp.diff (kffromarr, T)
print(dkffromarrdT)

lenArr = 198

A = smp.symbols('A:'+str(lenArr))
Ea = smp.symbols('Ea:'+str(lenArr))

kffromarr = smp.symbols('kffromarr')
kffromarr = A[87] * smp.exp(-Ea[90]/(R*T))
dkffromarrdT = smp.diff (kffromarr, T)
print(dkffromarrdT)
