C       CVS Revision:$Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
C	this is CHEMKIN-III file ckstrt.h V.3.0 January 1997;
C	it contains pointers for the gas-phase kinetics
C       subroutines' data storage arrays
C
C
C     include file for CHEMKIN-III cklib.f, dated: Oct. 2, 1966
C
      INTEGER
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD, 
     4   MXORD, KEL,  NKKI,
     5   IcMM, IcKK,
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL,
     7   IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, IcEI, 
     8   IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, IcTD, 
     9   IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,
     2   NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4


      COMMON /CKSTRT/ 
C
C     Integer constants
C
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD, 
     4   MXORD, KEL,  NKKI,
C
C     Integer pointers to character arrays in CCKWRK
C
     5   IcMM, IcKK,
C
C     Integer pointers to integer arrays in ICKWRK
C
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL,
     7   IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, IcEI, 
     8   IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, IcTD, 
     9   IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
C
C     Integer pointers to real variables and arrays in RCKWRK
C
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,
     2   NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4

C Logical global variables
C
C      NPERT = TRUE if any of the perturbation factors for the rate
C              constants have been changed from the value of one.
C
       LOGICAL         LPERT
       COMMON /CKGLBL/ LPERT
C
C     END include file for cklib.f
C
