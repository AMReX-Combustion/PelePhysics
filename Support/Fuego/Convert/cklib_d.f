C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $

      SUBROUTINE CKABS
C
C  START PROLOGUE
C
C  SUBROUTINE CKABS is the comment and upgrades information for the
C  CHEMKIN-III subroutine library.
C
C  START PROLOGUE
C
C///////////////////////////////////////////////////////////////////////
C
C  V 5.15 98/07/26 E. Meeks
C  1. Action#0185:  Corrected indexing in CKIORD; In IF statement when
C     KSPEC < 0, 'FORD(KSPEC,NFORD)' should be ' FORD(-KSPEC,NFORD)'
C  2. Action#0185: Changed call list, routine description and 
C     function of CKIORD so that it just returns available information 
C     about the NORD reactions: reaction indices and species orders, 
C     rather than information about forward vs. reverse order 
C     specifications, which is not available on the linking file.
C  V.5.14 98/04/08 E. Meeks
C  1. Action#0166:  Added subroutine CKIREV, to return explicitly
C     defined reverse reaction-rate coefficients A,B, and E
C  V.5.13 98/04/01 E. Meeks
C  1. Fix bug#0163: Correct pointer definition I_NK in soubroutine
C     CKIRNU:  I_NK = IcNK + (I-1)*MXSP.
C  V.5.12 97/10/28 E. Meeks
C  1. Fix bug#054b: improve overflow handling for FIT1 option in
C     CKRATT.
C  V.5.11 97/08/28 F. Rupley
C  1. Fix bug#005a: in CKRATX ensure non-negative concentrations to
c     calculate rates for real stoichiometric coefficients (DO 160),
C     and for changed orders (DO 200)
C  1. Fix bug#104: remove unused CKFRCH from EXTERNAL and INTEGER
C     statements in CKPARR.
C  V.5.10 97/07/23 (F. Rupley)
C  1. Fix bug#012: removed spacing from logical operators for f77 stnd.
C     Spaces found in CKEQ, line 3658 and CKRATT, line 9742.
C  V.5.9 97/04/15
C  1. move SAVE statements before DATA statements
C  V.5.8 97/03/01
C  1. fix NTR=1 case;
C     "final" version for Reaction Design
C  V.5.7, 96/12/24
C  1. fix indexing error, CKPY
C  V.5.6, 96/12/11
C  1. several bugfixes
C  V.5.5, 96/11/11
C  1. new time-saving logic in most thermodynamic and rate routines
C  V.5.4, 96/09/04
C  1. add INTEGER FUNCTION CKLKUP, returns index of an integer if
C     present in an integer array
C  V.5.3, 96/08/19
C  1. in CKRATT, set default EQK(I)=1.0 to cover SUMSMH=0.0.
C  VERSION 5.2, 96/05/24
C  1. initial sccs version
C  VERSION 5.1 (F. Rupley, April 26, 1996)
C  CAUTION...THIS WORK IS STILL IN PROGRESS:
C  1. CKRATX - IF (KFAL(N) .EQ. 0) should be IF (KFAL(N) .GT. 0)
C  2. delete CKRATN and CKQNCK (per E. Meeks)
C  3. reworked CKSYMR
C  CHANGES FOR VERSION 5.0 (CHEMKIN-III initial version,
C                           F. Rupley, March 1 1996)
C  1. add, modify comments sections throughout in preparation for
C     updated documentation.
C  2. add ascii linkfile options.
C  3. VERS=5.0 is linkfile version number, not to be changed unless
C     linkfile changes.
C  4. additional IFLO (new pointer IiFLO) for pressure-dependent
C     reactions, and additional logic for high- vs. low-pressure in
C     CKRATX.
C  5. eliminate some excess arguments from CKRATT list?
C  6. since many more linkfile records, in case of READ error,
C     set IFLAG=NREC, the index number of the record
C  CHANGES VERSION 4.8 (F. Rupley 11/27/95)
C  1. added or modified comments for CKSTRT, CKABS, and other subroutine
C     prolog sections (per E. Meeks) for upgrading documentation
C  2. reversed chronological or of VERSION change comments.
C  CHANGES FOR VERSION 4.7 (8/14.95 E. Meeks)
C  1.  Change maximum number of species in a reaction to 12
C      with corrections per R. Larson
C  2.  Added subroutine CKKTFL to initialize species temperature
C      array flag for nonthermal systems.  CKINIT now initializes
C      this array to ones for thermal systems.  The array only gets
C      set otherwise if CKKTFL is called by the application code.
C  3.  Generalized all temperature-dependent subroutines calculations
C      to allow for multiple species temperatures.  The default
C      temperature is T(1) for background gas temperature in plasma
C      cases or the temperature in thermal systems.
C  4.  Added subroutines for determining electron and ion transport
C      properties from reaction rate specifications for momentum-
C      transfer collisions.
C     CHANGES FOR VERSION 4.6 (2/27/95 F. Rupley)
C     1.  Change character index "(:" to "(1:"
C     CHANGES for VERSION 4.5 (1/19/95 per M. Coltrin)
C     1.  Add integer flag IFLAG to CKLEN and CKINIT and allow
C         RETURN instead of STOP for IFLAG error conditions.
C     CHANGES for VERSION 4.4 (11/17/94 per R. Steeper)
C     1.  Simplify CKRHOC,PKRHOC (superflous calculation for RHO)
C     2.  Simplify CKPC   (superflous calculations)
C     CHANGES for VERSION 4.3 (10/3/94 F. Rupley per E. Meeks)
C     1.  Correct calculation of RHO in PKRHOX.
C     CHANGES for VERSION 4.23 (8/26/94)
C     1.  Correct value of RUC (RCKWRK(NcRC)) in CKINIT.
C     CHANGES for VERSION 4.22 (8/15/94)
C     1.  Remove NSPEC(*) from CKRATX call list (ICKWRK(IcNS))
C     CHANGES for VERSION 4.21 (8/10/94)
C     1.  Accepts version 3.9 linkfile
C     CHANGES for VERSION 4.2 (6/13/94 F. Rupley, per E. Meeks)
C     1.  Modify CKRATT for plasma options.
C     2.  Add SUBROUTINES CKHRX, CKIEIM, CKIEXC
C     CHANGES for VERSION 4.10c (6/3/94 F. Rupley)
C     1.  Accept linkfile 3.6c (bugfixes per H. Moffat)
C     CHANGES for VERSION 4.10b (5/20/94 F. Rupley per E. Meeks)
C     1.  Incorporate plasma options (linkfile 3.6b)
C     CHANGES for VERSION 4.10 (4/28/94 F. Rupley, per M. Coltrin)
C     1.  New subroutines CKINU, CKIRNU, and CKIORD for real
C         stoichiometric coefficients and change of order reactions.
C     2.  Recognize linkfile
C     CHANGES FOR VERSION 4.9 (4/20/94 F. Rupley)
C     1.  Accept binary file V.3.6 (correction to CKUNIT)
C     CHANGES FOR VERSION 4.8 (4/19/94 F. Rupley)
C     1.  Accept binary file V.3.5 (correction to CKBAL, CKRBAL)
C     CHANGES FOR VERSION 4.7 (4/14/94 F. Rupley, suggested by E. Meeks)
C     1.  use INCLUDE 'ckstrt.h' instead of having the CKSTRT common
C         block in every subroutine.
C     CHANGES FOR VERSION 4.6 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements (updating interpreter), removing
C         unusued but possibly initialized variables.
C     CHANGES FOR V.4.5 (1/26/94 F. Rupley per R. Kee)
C     1. Implement real stoichometric coefficients; binary file V.3.3.
C     CHANGES FOR V.4.4 (11/10/93 F. Rupley)
C     1. Accept binary file V.3.2 (correction to CKUNIT)
C     CHANGES FOR V.4.3 (11/9/93 F. Rupley per E. Meeks)
C     1. Min/max single-precision exponent in CKR2CH should be
C        be 30 instead of 38.
C     CHANGES FOR V.4.2 (9/14/93 F. Rupley)
C     1. Move perturbation factoring from CKRATT to CKRATX
C     CHANGES FOR V.4.1 (2/24/93 F. Rupley)
C     1. Accept binary file V.3.1 (correction to CKREAC)
C     CHANGES FOR V.4.0 (10/1/92 F. Rupley per M. Coltrin)
C     1. COMMON /CKCONS/ VERS, PREC, KERR, LENI, LENR, LENC
C        eliminates need for LINKCK in argument list of CKSAVE
C     CHANGES FOR V.3.9 (4/17/92 F. Rupley)
C     1. Bugfix in CKSAVE (did not write new pointers NcKF,NcKR)
C     CHANGES FOR V.3.8 (4/15/92 F. Rupley)
C     1. Accept binary file V.3.0 (correction to CKDUP)
C     CHANGES FOR VERSION 3.7 (3/10/92 F. Rupley per Kee/Grcar)
C     1. Calls to CKRAT replaced by calls to CKRATT and CKRATX.
C     2. New subroutine CKKFRT returns the forward and reverse
C        rates (RKFT, RKRT) calculated by CKRATT (does not consider
C        pressure dependencies).
C     3. New subroutine CKWYPK returns the rates of production
C        given the RKFT and RKRT from (2).
C     CHANGES FOR VERSION 3.6 (2/24/92 F. Rupley per E. Meeks)
C     1. Accept binary file V.2.9 (additional error checking for
C        reverse T-L reactions, 2*II additional real work space)
C     2. Correct calculation for reverse T-L reaction rates
C     3. New subroutines CKRATT, CKRATX (subsets of CKRAT)
C     4. New pointers NcKF,NcKR to store intermediate temperature-
C        dependent rates.
CC     CHANGES FOR VERSION 3.4 (2/19/92 F. Rupley)
C     1. Correct error in CKITR (IcNR should be IcNS)
C     CHANGES FOR VERSION 3.3 (6/27/91 F. Rupley)
C     1. Accept binary file V.2.8 (modified interpreter output to
C        print all 16 characters of species names)
C     CHANGES FOR VERSION 3.2 (6/10/91 H. Moffat)
C     1. Added Subroutine CKFAL, which returns the fall-off parameters
C        for the mechanism.
C     2. Added Subroutine CKNUF, which returns the reactant
C        stoichiometric coefficients.
C     3. Fixed an error in CKSNUM, which caused an error condition when
C        the input string did not have any blanks inbetween words.
C     4. Fixed two errors in CKTHB. The default third body efficiency
C        should be equal to 1.
C     CHANGES FOR VERSION 3.1 (5/9/91 F. Rupley)
C     1. Add Subroutine CKMXTP to return number of temperatures used
C        in thermodynamic fits.
C     CHANGES FOR VERSION 3.0 (4/1/91 F. Rupley)
C     1. Accept binary file V.2.7 (modification of CKDUP)
C        coefficient a6.
C    CHANGES TO VERSION 2.9 (2/15/91, F. Rupley per R. Kee)
C     2. Subroutine CKRHEX allows perturbation of thermodynamic
C     1. Add a fourth parameter to the array of Arhennius coefficients
C        for the II reactions;
C        increase the value of NPAR in COMMON /CKSTRT/ by one (this
C        also increases the length of the array of reverse Arhennius
C        parameters);
C        initialize the value of the fourth parameter to 1.0 in
C        CKINIT;
C        use this value as a "perturbation factor" for the forward
C        rates in CKRAT;
C        add SUBROUTINE CKRDEX to allow applications codes to change
C        the perturbation factor RD(I) in sensitivity calculations.
C     2. Accept binary file V.2.6 (LENRCK was increased by II+NREV to
C        reflect above changes in RCKWRK array.
C    CHANGES TO VERSION 2.8 (1/18/91, F. Rupley)
C     1. Accept binary file V.2.5
C    CHANGES TO VERSION 2.7 (12/20/90, F. Rupley)
C     1. Accept binary file V.2.4
C    CHANGES TO VERSION 2.6 (12/15/90, F. Rupley)
C     1. Accept binary file V.2.3
C    CHANGES TO VERSION 2.5 (11/15/90, F. Rupley)
C     1. Accept binary file V.2.2
C    CHANGES TO VERSION 2.4
C     1. Accept binary  file V.2.1
C    CHANGES TO VERSION 2.3
C     1. Accept binary file V.2.0
C    CHANGES TO VERSION 2.2
C     1. Bugfix in CKABML
C     2. In CKXNUM (and CKSNUM), if NEXP is negative, it is not an
C        error to find fewer values.
C    CHANGES TO VERSION 2.1
C     1. New binary file has an additional record to indicate its
C        version, machine precision, and error status
C     2. SUBROUTINE CKPNT reads a binary file to get COMMON /CKSTRT/
C        pointers.
C     3. SUBROUTINE CKSAVE writes pointers and work arrays to a
C        binary file.
C     4. Add COMMON /MACH/ and initialization of BIG,SMALL,EXPARG to
C        SUBROUTINE CKPNT
C     5. Change LOG(*) to LOG(MAX(*,SMALL)) in several subroutines.
C    CHANGES TO VERSION 2.0
C     1. Subroutine CKLEN to provide necessary lengths of work arrays.
C     2. Subroutine CKKFKR provides arrays of forward and reverse
C        reaction rates.
C    CHANGES FROM VERSION 1.8
C     1. vax/cray change blocks for machine constants changed to
C        smallexp, bigexp change blocks
C     2. add ERR= to first read statement in CKINIT
C    CHANGES FROM VERSION 1.7
C     1. Get rid of non-standard comment statements.
C    CHANGES FROM VERSION 1.6
C     1. Added error checking and additional arguments to character
C        manipulation subroutines.
C     2. Fixed an error with IFIRCH(LINE) in IPPLEN
C    CHANGES FROM VERSION 1.5
C     1. Implement Landau-Teller rate expression
C     2. Add utility routine CKR2CH
C    CHANGES FROM VERSION 1.4
C     1. New versions of CKABML/CKABMS, CKGBML/CKGBMS, CKSBML/CKSBMS
C        for mixture-averaging require additional argument P
C     2. Replace SDOT's and DDOT's by loops
C     3. Reaction strings now have "=>" in irreversible reactions,
C                                 "<=>" in reversible reactions.
C    CHANGES FROM VERSION 1.3
C     1. added PROLOGUES
C     2. binary file now includes minimum length of arrays
C    CHANGES FROM VERSION 1.2
C     1. change SCOPY for integer arrays to DO loops
C    CHANGES FROM VERSION 1.1
C     1. add change block to CKCHRG
C     2. correct change block in CKHORT
C    CHANGES FROM VERSION 1.0
C     1. REPLACE "REAL*8" WITH "DOUBLE PRECISION"
C
C//////////////////////////////////////////////////////////////////////C
C
C  START PROLOGUE
C
C  SUBROUTINE CKABS
C
C  Work arrays ICKWRK, RCKWRK, and CCKWRK contain information about the
C  elements, species and reactions in the mechanism; they also contain
C  some work space needed for internal manipulations.  A user wishing
C  to modify a subroutine or to write new routines will probably want
C  to use the work arrays directly.  The pointers described below are
C  starting addresses for information stored in the work arrays, and
C  are found in the labeled common block COMMON /CKSTRT/, declared
C  by the use of the include file ckstrt.h.
C
C  COMMON /CKSTRT/
C
C  Integer constants
C
C 1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
C 2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
C 3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
C 4   MXORD,KEL,  NKKI,
C
C  Integer pointers to character string arrays in CCKWRK
C
C 5   IcMM, IcKK,
C
C  Integer pointers to integer arrays in ICKWRK
C
C 6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL,
C 7   IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, IcEI,
C 8   IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, IcTD,
C 9   IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
C
C  Integer pointers to real variables and arrays in RCKWRK
C
C *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
C 1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,
C 2   NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
C
C  INTEGER CONSTANTS:
C
C  NMM,   Total count, elements in problem.
C  NKK,   Total count, species in problem.
C  NII,   Total count, reactions in problem.
C  MXSP,  Maximum number of species (reactants plus products) allowed
C         in any reaction; unless changed in the interpreter, MXSP=12.
C  MXTB,  Maximum number of enhanced third-bodies allowed in any
C         reaction; unless changed in the interpreter, MXTB=10.
C  MXTP,  Maximum number of temperatures allowed in fits of
C         thermodynamic properties for any species;
C         unless changed in the interpreter and the thermodynamic
C         database, MXTP=3.
C  NCP,   Number of polynomial coefficients to fits of CP/R for a
C         species; unless changed in the interpreter and the
C         thermodynamic database, NCP=5.
C  NCP1,  NCP + 1
C  NCP2,  NCP + 2
C  NCP2T, Total number of thermodynamic fit coefficients for species;
C         unless changed, NCP2T = (MXTP-1)*NCP2 = 14.
C  NPAR,  Number of parameters required in the rate expression for
C         reactions;  in the current formulation NPAR=3, however,
C         a 4th parameter is used for purposes of scaling.
C  NLAR,  Number of parameters required for Landau-Teller reactions;
C         NLAR=4.
C  NFAR,  Number of parameters allowed for pressure-dependent
C         reactions; NFAR=8.
C  NLAN,  Total count, Landau-Teller reactions.
C  NFAL,  Total count, pressure-dependent reactions.
C  NREV,  Total count, reactions with reverse parameters.
C  NTHB,  Total count, reactions with third-bodies.
C  NRLT,  Total count, Landau-Teller reactions with reverse parameters.
C  NWL,   Total count, reactions with radiation wavelength enhancement.
C  NEIM,  Total count, electron-impact reactions.
C  NJAN,  Total count, Janev-Langer,Evans,Post reactions.
C  NJAR,  Number of parameters required for an NJAN reaction.
C  NFT1,  Total count, reactions using fit#1.
C  NF1R,  Number of parameters required for an NFT1 reaction.
C  NEXC,  Total count, excitation-only reactions.
C  NMOM,  Total count, electron momentum-transfer reactions.
C  NXSM,  Total count, ion momentum-transfer reactions.
C  NTDE,  Total count, non-thermal-equilibrium reactions.
C  NRNU,  Total count, real stoichiometry reactions.
C  NORD,  Total count, changed-order reactions.
C  MXORD, Maximum number of order changes allowed for above.
C  KEL,   Species index of the electron species if present.
C  NKKI,  Total count, ion species in the mechanism.
C
C  STARTING ADDRESSES FOR THE CHARACTER WORK SPACE, CCKWRK.
C
C  IcMM,  CCKWRK(I = IcMM) starts an array of element names;
C         CCKWRK(I + M - 1) is the name of element M.
C  IcKK,  CCKWRK(I = IcKK) starts an array of species names;
C         CCKWRK(I + K - 1) is the name of species K.
C
C  STARTING ADDRESSES FOR THE INTEGER WORK SPACE, ICKWRK.
C
C  IcNC, ICKWRK(I = IcNC) starts a matrix of elemental composition
C        for the species;
C        ICKWRK(I + (K-1)*NMM + M - 1) is the quantity of element M
C        in species K.
C  IcPH, ICKWRK(I = IcPH) starts an array of physical phases for the
C        species;
C        ICKWRK(I + K - 1) = -1, species K is solid
C                          =  0, species K is gaseous
C                          = +1, species K is liquid
C  IcCH, ICKWRK(I = IcCH) starts an array of electronic charges for
C        the species;
C        ICKWRK(I + K - 1) = -2, species K has two excess electrons.
C  IcNT, ICKWRK(I = IcNT) starts an array of the total number of
C        temperatures dividing the ranges of thermodynamic fits of
C        the species;
C        ICKWRK(I + K - 1) is the number of dividing temperatures
C        for thermodynamic fits for species K.
C  IcNU, ICKWRK(I = IcNU) starts a matrix of stoichiometric coefficients
C        for the MXSP species in the reactions;
C        ICKWRK(I + (N-1)*MXSP + L - 1) is the coefficient of the Lth
C        participant species in the Nth reaction.
C  IcNK, ICKWRK(I = IcNK) starts a matrix of indices for the MXSP
C        species in the reactions;
C        ICKWRK(I + (N-1)*MXSP + L -1) is the species index for the
C        Lth participant species in the Nth reaction.
C  IcNS, ICKWRK(I = IcNS) starts an array of the total counts of
C        participant species for the reactions,
C        and indicates the reversibility of the reactions;
C        ICKWRK(I + N - 1) = +L, reaction N is reversible and has
C                             L participant species (reactants+products)
C                          = -L, reaction N is irreversible and has
C                             L participant species (reactants+products)
C  IcNR, ICKWRK(I = IcNR) starts an array of the total count of
C        reactants only for the reactions;
C        ICKWRK(I + N - 1) is the total reactant count for reaction N.
C  IcLT, ICKWRK(I = IcLT) starts an array of reaction indices for
C        Landau-Teller reactions;
C        ICKWRK(I + N - 1) is the reaction index of the Nth LT reaction.
C  IcRL, ICKWRK(I = IcRL) starts an array of reaction indices for
C        Landau-Teller reactions with explicit reverse parameters;
C        ICKWRK(I + N - 1) is the reaction index of the Nth reaction
C        with reverse Landau-Teller parameters.
C  IcRV, ICKWRK(I = IcRV) starts an array of reaction indices for those
C        with explicit reverse Arrhenius coefficients;
C        ICKWRK(I + N - 1) is the reaction index of the Nth reaction
C        with reverse coefficients.
C  IcWL, ICKWRK(I = IcWL) starts an array of reaction indices for those
C        with radiation wavelength enhancement;
C        ICKWRK(I + N - 1) is the reaction index of the Nth reaction
C        with wavelength enhancement.
C  IcFL, ICKWRK(I = IcFL) starts an array of reaction indices for those
C        with pressure-dependent formulations;
C        ICKWRK(I + N - 1) is the reaction index of the Nth pressure-
C        dependent reaction.
C  IcFO, ICKWRK(I = IcFO) starts an array of formulation types for
C        pressure-dependent reactions;
C        ICKWRK(I + N - 1) is the type of the Nth pressure-dependent
C        reaction,  1 for 3-parameter Lindemann Form
C                   2 for 6- or 8-parameter SRI Form
C                   3 for 6-parameter Troe Form
C                   4 for 7-parameter Troe form
C  IcFT, ICKWRK(I = IcFT) starts an array of option types for pressure-
C        dependent reactions;
C        ICKWRK(I + N - 1) is an option for the Nth pressure-dependent
C        reaction, 0 for unimolecular fall-off,
C                  1 for chemically activated.
C  IcKF, ICKWRK(I = IcKF) starts an array of third-body species flags
C        for pressure-dependent reactions;
C        ICKWRK(I + N - 1) is the third-body species flag for the Nth
C        pressure-dependent reaction,
C        0, the concentration of the third-body is the sum of the
C           concentrations of all species in the problem
C        K, the concentration of the third-body is the concentration
C           of species K.
C  IcTB, ICKWRK(I = IcTB) starts an array of reaction indices for those
C        with enhanced third-bodies;
C        ICKWRK(I + N - 1) is the reaction index of the Nth third-body
C        reaction.
C  IcKN, ICKWRK(I = IcKN) starts an array of enhanced species counts
C        for 3rd-body reactions;
C        ICKWRK(I + N - 1) is the total enhanced species count for the
C        Nth third-body reaction.
C  IcKT, ICKWRK(I = IcTB) starts a matrix of species indices for
C        enhanced 3rd bodies in third-body reactions;
C        ICKWRK(I + (N-1)*MXTB + L - 1) is the species index of the
C        Lth enhanced species in the Nth third-body reaction.
C  IcEI, ICKWRK(I = IcEI) starts an array of reaction indices for
C        electron-impact reactions;
C        ICKWRK(I + N - 1) is the reaction index of the Nth electron-
C        impact reaction.
C  IcET, ICKWRK(I = IcET) starts an array of temperature-dependence
C        flags for electron-impact reactions;
C        ICKWRK(I + N - 1) is a pointer to the temperature in an array
C        which is used to compute the reaction's rate.
C  IcJN, ICKWRK(I = IcJN) starts an array of reaction indices for
C        Janev-Langer reactions;
C        ICKWRK(I + N - 1) is the reaction index of the Nth Janev-
C        Langer reaction.
C  IcF1, ICKWRK(I = IcF1) starts an array of reaction indices for
C        fit-type reactions;
C        ICKWRK(I + N - 1) is the reaction index of the Nth fit-type
C        reaction.
C  IcEX, ICKWRK(I = IcEX) starts an array of reaction indices for
C        excitation-only reactions;
C        ICKWRK(I + N - 1) is the reaction index of the Nth excitation-
C        only reaction.
C  IcMO, ICKWRK(I = IcMO) starts an array of reaction indices for those
C        with electron momentum-transfer;
C        ICKWRK(I + N - 1) is the reaction index of the Nth electron
C        momentum-transfer reaction.
C  IcMK, ICKWRK(I = IcMK) starts an array of species indices for an
C        electron's collision partner in the electron momentum-transfer
C        reactions;
C        ICKWRK(I + N - 1) is the species index of the collision
C        partner in the Nth electron momentum-transfer reaction.
C  IcXS, ICKWRK(I = IcXS) starts an array of reaction indices for those
C        with ion momentum-transfer cross-section;
C        ICKWRK(I + N - 1) is the reaction index of the Nth ion
C        momentum-transfer cross-section reaction.
C  IcXI, ICKWRK(I = IcXI) starts an array of species indices for the
C        ion collision partner in ion momentum-transfer reactions;
C        ICKWRK(I + N - 1) is the species index of the ion collision
C        partner in the Nth ion momentum-transfer reaction.
C  IcXK, ICKWRK(I = IcXK) starts an array of species indices for the
C        non-ion collision partner in ion momentum-transfer reactions;
C        ICKWRK(I + N - 1) is the species index of the non-ion
C        collision partner for the Nth.
C  IcTD, ICKWRK(I = IcTD) starts an array of reaction indices for those
C        with non-thermal-equilibrium temperature-dependence;
C        ICKWRK(I + N - 1) is the reaction index of the Nth non-
C        thermal-equilibrium reaction.
C  IcTK, ICKWRK(I = IcTK) starts an array of temperature-dependent
C        species indices for the non-thermal-equilibrium reactions;
C        ICKWRK(I + N - 1) is the index of the species which determines
C        the temperature for the Nth non-thermal-equilibrium reaction.
C  IcRNU,ICKWRK(I = IcRNU) starts an array of reaction indices for those
C        with real stoichiometry;
C        ICKWRK(I + N - 1) is the reaction index of the Nth reaction
C        with real stoichiometry.
C  IcORD,ICKWRK(I = IcORD) starts an array of reaction indices for those
C        with changed-order species;
C        ICKWRK(I + N - 1) is the reaction index of the Nth reaction
C        with changes of order.
C  IcKOR,ICKWRK(I = IcKOR) starts a matrix of species indices for the
C        changed-order reactions;
C        K = ICKWRK(I + (N-1)*MXORD + L - 1) is the species number of
C            the Lth change-of-order for the Nth changed-order
C            reaction;
C          > 0, K participates in the reverse direction (product),
C          < 0, K participates in the forward direction (reactant).
C  IcKI, ICKWRK(I = IcKI) starts an array of species indices for those
C        which are ions;
C        ICKWRK(I + N - 1) is the species index for the Nth ion.
C  IcKTF,ICKWRK(I = IcKTF) starts an array of temperature array pointers
C        for the species;
C        ICKWRK(I + K - 1) is the pointer to the temperature array
C        for species K.
C  IcK1, ICKWRK(I = IcK1) starts scratch storage space of length NKK.
C  IcK2, ditto
C
C  STARTING ADDRESSES FOR THE REAL WORK SPACE, RCKWRK.
C
C  NcAW, RCKWRK(I = NcAW) starts an array of atomic weights (gm/mole);
C        RCKWRK(I + M - 1) is the atomic weight of element M.
C  NcWT, RCKWRK(I = NcWT) starts an array of molecule weights (gm/mole);
C        RCKWRK(I + K - 1) is the molecular weight of species K.
C  NcTT, RCKWRK(I = NcTT) starts an array of temperatures (Kelvin)
C        used to fit thermodynamic properties for the species;
C        RCKWRK(I + (K-1)*MXTP + N - 1) is the Nth temperature for
C        species K.
C  NcAA, RCKWRK(I = NcAA) starts a three-dimensional array of polynomial
C        coefficients for thermodynamic properties of the species;
C        RCKWRK(I + (L-1)*NCP2 + (K-1)*NCP2T + N - 1) is the Nth
C        polynomial coefficient A(N,L,K) for species K, in the Lth
C        temperature range.
C  NcCO, RCKWRK(I = NcCO) starts a matrix of Arrhenius parameters for
C        reactions;
C        RCKWRK(I + (N-1)*(NPAR+1) + L -1) is the Lth parameter of
C        reaction N, where
C           L=1 is the pre-exponential factor (mole-cm-sec-K),
C           L=2 is the temperature exponent,
C           L=3 is the activation energy (Kelvins), and
C           L=4 is used as a scalar for sensitivity analysis.
C  NcRV, RCKWRK(I = NcRV) starts a matrix of reverse Arrhenius
C        parameters for reactions which give them explicitly;
C        RCKWRK(I + (N-1)*(NPAR+1) + L - 1) is the Lth reverse parameter
C        for the Nth reaction with reverse parameters declared, where
C           L=1 is the pre-exponential factor (mole-cm-sec-K),
C           L=2 is the temperature exponent,
C           L=3 is the activation energy (Kelvins),
C           L=4 is used as a scalar for sensitivity analysis.
C           The reaction index is ICKWRK(IcRV + N - 1).
C  NcLT, RCKWRK(I = NcLT) starts a matrix of parameters for the Landau-
C        Teller reactions;
C        RCKWRK(I + (N-1)*NLAR + L - 1) is the Lth Landau-Teller
C        parameter for the Nth Landau-Teller reaction, where
C           L=1 is B(I) (Eq. 73) (Kelvins**1/3), and
C           L=2 is C(I) (Eq. 73) (Kelvins**2/3).
C           The reaction index is ICKWRK(IcLT + N - 1).
C  NcRL, RCKWRK(I = NcRL) starts a matrix of explicitly-given reverse
C        Landau-Teller parameters;
C        RCKWRK(I + (N-1)*NLAR + L - 1) is the Lth reverse parameter
C        for the Nth reaction with reverse Landau-Teller parameters,
C        where
C           L=1 is B(I) (Eq. 73) (Kelvins**1/3), and
C           L=2 is C(I) (Eq. 73) (Kelvins**2/3).
C           The reaction index is ICKWRK(IcRL + N - 1).
C  NcFL, RCKWRK(I = NcFL) starts a matrix of parameters for the
C        pressure-dependent reactions;
C        RCKWRK(I + (N-1)*NFAR + L - 1) is the Lth parameter for the
C        Nth pressure-dependent reaction, where the low pressure limits
C        limits are defined by
C           L=1 is the pre-exponential factor (mole-cm-sec-K),
C           L=2 is the temperature exponent, and
C           L=3 is the activation energy (Kelvins).
C        Additional parameters define the centering, depending on
C        the type of formulation -
C           Troe: L=4 is the Eq. 68 parameter a,
C                 L=5 is the Eq. 68 parameter T*** (Kelvins),
C                 L=6 is the Eq. 68 parameter T*   (Kelvins), and
C                 L=7 is the Eq. 68 parameter T**  (Kelvins).
C           SRI:  L=4 is the Eq. 69 parameter a,
C                 L=5 is the Eq. 69 parameter b (Kelvins),
C                 L=6 is the Eq. 69 parameter c (kelvins),
C                 L=7 is the Eq. 69 parameter d, and
C                 L=8 is the Eq. 69 parameter e.
C        The reaction index is ICKWRK(IcFL+N-1) and the type of
C        formulation is ICKWRK(IcFO+N-1).
C  NcKT, RCKWRK(I = NcKT) starts a matrix of enhancement factors for
C        third-body reactions;
C        RCKWRK(I + (N-1)*MXTB + L - 1) is an enhancement factor for
C        the Lth enhanced species in the Nth third-body reaction;
C        the reaction index is ICKWRK(IcTB+N-1), and the Lth
C        enhanced species index is ICKWRK(IcKT+(N-1)*MXTB+L-1).
C  NcWL, RCKWRK(I = NcWL) starts an array of wavelengths for wavelength-
C        enhanced reactions;
C        RCKWRK(I + N - 1) is the wavelength enhancement (angstrom)
C        for the Nth wavelength-enhanced reaction;
C        the reaction index is ICKWRK(IcWL+N-1).
C  NcJN, RCKWRK(I = NcJN) starts a matrix of parameters for the Janev-
C        Langer reactions;
C        RCKWRK(I + (N-1)*NJAR + L - 1) is Lth parameter for the Nth
C        Janev-Langer reaction, where
C           L=1 is
C           L=2 is
C        The reaction index is ICKWRK(IcJN+N-1).
C  NcF1, RCKWRK(I = NcF1) starts an array of parameters for fit#1
C        reactions;
C        RCKWRK(I + (N-1)*NF1R + L - 1) is the Lth parameter for the
C        Nth fit-type reaction, where
C           L=1 is
C           L=2 is
C        The reaction index is ICKWRK(IcF1+N-1).
C  NcEX, RCKWRK(I = NcEX) starts an array of enervy losses for
C        excitation-only reactions;
C        RCKWRK(I + N - 1) is the excitation energy loss per event
C        for the Nth excitation only reaction.
C        The reaction index is ICKWRK(IcEX+N-1).
C  NcRU, RCKWRK(I = NcRU) is the universal gas constant (ergs/mole-K).
C  NcRC, RCKWRK(I = NcRC) is the universal gas constant (cal/mole-K).
C  NcPA, RCKWRK(I = NcPA) is the pressure of one standard atmosphere
C          (dynes/cm**2).
C  NcKF, RCKWRK(I = NcKF) starts an array of the temperature-dependent
C        forward rate components for reactions.
C  NcKR, RCKWRK(I = NcKR) starts an array of the temperature-dependent
C        reverse rate components for reactions.
C  NcRNU,RCKWRK(I = NcRNU) starts a matrix of stoichiometric
C        coefficients for reactions with real stoichiometry;
C        RCKWRK(I + (N-1)*MXSP + L - 1) is the coefficient for
C        the Lth species in the Nth real stoichiometry reaction.
C        The reaction index is ICKWRK(IcRNU+N-1).
C        The species index is ICKWRK(IcNUNK+(N-1)*MXSP+L-1).
C  NcKOR,RCKWRK(I = NcKOR) starts a matrix of order values for changed-
C        order species reactions;
C        RCKWRK(I + (N-1)*MXORD + L - 1) is the order for the Lth
C        changed-order species in the Nth change-order reaction.
C        The reaction index is ICKWRK(IcKOR+N-1).
C        The change-order species index is
C        ICKWRK(IcKOR+(N-1)*MXORD+L-1).
C  NcK1, RCKWRK(I = NcK1) starts species scratch workspace.
C  NcK2, RCKWRK(I = NcK2) starts species scratch workspace.
C  NcK3, RCKWRK(I = NcK3) starts species scratch workspace.
C  NcK4, RCKWRK(I = NcK4) starts species scratch workspace.
C  NcI1, RCKWRK(I = NcI1) starts reaction scratch workspace.
C  NcI2, RCKWRK(I = NcI2) starts reaction scratch workspace.
C  NcI3, RCKWRK(I = NcI3) starts reaction scratch workspace.
C  NcI4, RCKWRK(I = ncI4) starts reaction scratch workspace.
C
C  STORING DATA INTO THE ARRAYS is usally accomplished by a
C  CALL SKINIT, which reads a linkfile generated by the gas-phase
C  mechanism interpreter;
C  the linkfile consists of the following records:
C
C  Linkfile information:
C  1. FILVER
C  2. PRVERS
C  3. PREC
C  4. KERR
C  5. LENI, LENR, LENC
C
C  Parameters and constants:
C  6. MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,
C     NIFAR, NJA, MXORD, NF1
C  7. MM, KK, II, NREV, NFAL, NTHB, NLAN, NRLT, NWL,
C     NCHRG, NEIM, NJAN, NFT1, NEXC, NMOM, NXSM,
C     NTDE, NRNU, NORD, KEL, KKION
C  8. CKMIN
C
C  Element data:
C  9. ( CCKWRK(IcMM + M - 1), M = 1, MM)      element names
C 10. ( RCKWRK(NcAW + M - 1), M = 1, MM)      atomic weights
C
C  Species data:
C 11. ( CCKWRK(IcKK + K - 1), K = 1, KK)      species names
C 12. ( RCKWRK(NcWT + K - 1), K = 1, KK)      molecular weights
C 13. (( ICKWRK(IcNC + (K-1)*MM + M - 1),     composition
C        M = 1, MM), K = 1, KK)
C 14. ( ICKWRK(IcCH + K - 1), K = 1, KK)      electronic charge
C 15. ( ICKWRK(IcNT + K - 1), K = 1, KK)      #fit temperatures
C 16. ( ICKWRK(IcPH + K - 1), K = 1, KK)      physical phase
C 17. (( RCKWRK(NcTT + (K-1)*MAXTP + L - 1),  fit temperatures
C        L = 1, MAXTP), K = 1, KK)
C 18. ((( RCKWRK(NcAA + (L-1)*NCP2 + (K-1)*NCP2T + N - 1),
C         N = 1, NCP2), L = 1, (MAXTP-1)), K = 1, KK)  thermodynamics
C
C  Ion data (if NKKI > 0):
C 19. NKKI
C 20. ( ICKWRK(IcKI + K - 1), K = 1, NKKI)    ion species indices
C
C  Reaction data (if II > 0:
C 21. ( ICKWRK(IcNS + I - 1), I = 1, II)      species counts
C 22. ( ICKWRK(IcNR + I - 1), I = 1, II)      reactant counts
C 23. (( ICKWRK(IcNU + (I-1)*MAXSP + N - 1),  stoichiometric coeff'nts
C        ICKWRK(IcNK + (I-1)*MAXSP + N - 1),  species indices
C        N = 1, MAXSP), I = 1, II)
C 24. (( RCKWRK(NcCO + (I-1)*(NPAR+1) + N - 1),
C        N = 1, NPAR), I = 1, II)             Arrhenius coefficients
C
C  Explicit reverse parameter reaction data (if NREV > 0):
C 25. NREV
C 26. ( ICKWRK(IcRV + N - 1), N = 1, NREV)    reaction indices
C 27. (( RCKWRK(NcRV + (N-1)*(NPAR+1) + L - 1),
C        L = 1, NPAR), N = 1, NREV)           reverse coefficients
C
C  Pressure-dependent reaction data (if NFAL > 0):
C 28. NFAL, NFAR
C 29. ( ICKWRK(IcFL+N-1),                     reaction indices
C       ICKWRK(IcFO+N-1),                     option type
C       ICKWRK(IcFT+N-1),                     flow type
C       ICKWRK(IcKF+N-1), N = 1, NFAL)        3rd-body species index
C 30. (( RCKWRK(NcFL + (N-1)*NFAR + L - 1),   option parameters
C        L = 1, NFAR), N = 1, NFAL)
C
C  Third-body reaction data (if NTHB > 0):
C 31. NTHB
C 32. ( ICKWRK(IcTB + N - 1),                 reaction indices
C       ICKWRK(IcKN + N - 1), N = 1, NTHB)    3rd body count
C 33. (( ICKWRK(IcKT + (N-1)*MAXTB + L - 1),  3rd body species indices
C        L = 1, MAXTB), N = 1, NTHB)
C 34. (( RCKWRK(NcKT + (N-1)*MAXTB + L - 1),  enhancement factors
C        L = 1, MAXTB), N = 1, NTHB)
C
C  Landau-Teller reaction data (if NLAN > 0):
C 35. NLAN, NLAR
C 36. ( ICKWRK(IcLT + N - 1), N = 1, NLAN)    reaction indices
C 37. (( RCKWRK(NcLT + (N-1)*NLAR + L - 1),   L-T parameters
C        L = 1, NLAR), N = 1, NLAN)
C
C  Landau-Teller reverse reaction data (if NRLT > 0):
C 38. NRLT
C 39. ( ICKWRK(IcRL + N - 1), N = 1, NRL)     reaction indices
C 40. ((RCKWRK(NcRL + (N-1)*NLAR + L - 1),    reverse L-T parameters
C       L = 1, NLAR), N = 1, NRLT)
C
C  Wavelength enhancement reaction data (if NWL > 0):
C 41. NWL
C 42. ( ICKWRK(IcWL + N - 1), N = 1, NWL)     reaction indices
C 43. ( RCKWRK(NcWL + N - 1), N = 1, NWL)     enhancement factors
C
C  Electron-impact reaction data (if NEIM > 0):
C 44. NEIM
C 45. ( ICKWRK(IcEI + N - 1), N = 1, NEIM)    reaction indices
C 46. ( ICKWRK(IcET + N - 1), N = 1, NEIM)    electron energy
C
C  Janev-Langer reaction data (if NJAN > 0):
C 47. NJAN, NJAR
C 48. ( ICKWRK(IcJN + N - 1), N = 1, NJAN)    reaction indices
C 49. ((RCKWRK(NcJN + (N-1)*NJAR + L - 1),    J-L parameters
C       L = 1, NJAR), N = 1, NJAN)
C
C  Fit #1 reaction data (if NFT1 > 0):
C 50. NFT1, NF1R
C 51. ( ICKWRK(IcF1 + N - 1), N = 1, NFT1)    reaction indices
C 52. ((RCKWRK(NcF1 + (N-1)*NF1R + L - 1),    fit#1 parameters
C       N = 1, NF1R), N = 1, NFT1)
C
C  Excitation-only reaction data (if NEXC > 0):
C 53. NEXC
C 54. ( ICKWRK(IcEX + N - 1), N = 1, NEXC)    reaction indices
C 55. ( RCKWRK(NcEX + N - 1), N = 1, NEXC)
C
C  Electron momentum-transfer collision reaction data (if NMOM > 0):
C 56. NMOM
C 57. ( ICKWRK(IcMO + N - 1), N = 1, NMOM)    reaction indices
C 58. ( ICKWRK(IcMK + N - 1), N = 1, NMOM)    partner species indices
C
C  Ion momentum-transfer cross-section reaction data (if NXSM > 0):
C 59. NXSM
C 60. ( ICKWRK(IcXS + N - 1), N = 1, NXSM)    reaction indices
C 61. ( ICKWRK(IcXI + N - 1), N = 1, NXSM)    ion species indices
C 62. ( ICKWRK(IcXK + N - 1), N = 1, NXSM)    partner species indices
C
C  Non-thermal-equilibium reaction data (if NTDE > 0):
C 63. NTDE
C 64. ( ICKWRK(IcTD + N - 1), N = 1, NTDE)    reaction indices
C 65. ( ICKWRK(IcTK + N - 1), N = 1, NTDE)    species indices
C
C  Real stoichiometry reaction data (if NRNU > 0):
C 66. NRNU
C 67. ( ICKWRK(IcRNU + N - 1), N = 1, NRNU)   reaction indices
C 68. ((RCKWRK(NCRNU + (N-1)*MAXSP + L - 1),  stoichiometric coeff'nts
C       L = 1, MAXSP), N = 1, NRNU)
C
C  Changed-order reaction data (if NORD > 0):
C 69. NORD
C 70. ( ICKWRK(IcORD + N - 1), N = 1, NORD)   reaction indices
C 71. ((ICKWRK(IcKOR + (N-1)*MXORD + L - 1),  change-order species
C       L = 1, MXORD), N = 1, NORD)
C 72. ((RCKWRK(NcKOR + (N-1)*MXORD + L - 1),  change-order values
C       L = 1, MXORD), N = 1, NORD)
C
C  END PROLOGUE
C
C     end of SUBROUTINE CKABS
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKABE  (ICKWRK, RCKWRK, RA, RB, RE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKABE  (ICKWRK, RCKWRK, RA, RB, RE)
C  Returns the Arrhenius coefficients of the reactions; see Eq. (52).
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RA(*)     - Real array, pre-exponential constants for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, mole-cm-sec-K
C  RB(*)     - Real array, temperature dependence exponents for
C              reactions;
C              dimension at least II, total reaction count.
C                 cgs units none
C  RE(*)     - Real array, activation energies for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, K
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), RA(*), RB(*), RE(*)
C
      IND = NcCO
      DO 100 I = 1, NII
         RA(I) = RCKWRK(IND)
         RB(I) = RCKWRK(IND+1)
         RE(I) = RCKWRK(IND+2)
         IND = IND + NPAR + 1
  100 CONTINUE
C
C     end of SUBROUTINE CKABE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKABML (P, T, X, ICKWRK, RCKWRK, ABML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKABML (P, T, X, ICKWRK, RCKWRK, ABML)*
C  Returns the Helmholtz free energy of the mixture in molar units
C  given pressure, temperature(s), and mole fractions; see Eq. (46).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ABML      - Real scalar, mean Helmholtz free energy.
C                 cgs units, ergs/mole
C
C   END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      RU = RCKWRK(NcRU)
      RLNP = RU * LOG(P / RCKWRK(NcPA))
C
      ABML = 0.0
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         ABML = ABML + X(K+1)
     1        * ( RCKWRK(NcK2 + K) - T(ICKWRK(IcKTF + K)) *
     2        (RCKWRK(NcK1 + K) - RU * LOG(MAX(X(K+1),SMALL)) - RLNP) )
  100 CONTINUE
C
C     end of SUBROUTINE CKABML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKABMS (P, T, Y, ICKWRK, RCKWRK, ABMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKABMS (P, T, Y, ICKWRK, RCKWRK, ABMS)*
C  Returns the mean Helmholtz free energy of the mixture in mass units
C  given pressure, temperature(s) and mass fractions; see Eq. (47).
C
C  INPUT
C  P         - Real scalar, pressure.
C                cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ABMS      - Real scalar, mean Helmholtz free energy.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK3))
      CALL CKMMWY(Y, ICKWRK, RCKWRK, WTM)
      RU = RCKWRK(NcRU)
      RLNP = RU * LOG (P / RCKWRK(NcPA))
C
      SUM = 0.0
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         SUM = SUM + RCKWRK(NcK3 + K) *
     1             ( RCKWRK(NcK2 + K) - T(ICKWRK(IcKTF + K)) *
     2             ( RCKWRK(NcK1 + K) - RU *
     4               LOG(MAX(RCKWRK(NcK3 + K),SMALL)) - RLNP))
  100 CONTINUE
      ABMS = SUM / WTM
C
C     end of SUBROUTINE CKABMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKAML  (T, ICKWRK, RCKWRK, AML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKAML  (T, ICKWRK, RCKWRK, AML)
C  Returns the standard state Helmholtz free energies in molar units;
C  see Eq. (25).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  AML(*)    - Real array, standard state Helmholtz free energies
C              for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), AML(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
C
      RU = RCKWRK(NcRU)
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         AML(K) = RCKWRK(NcK2 + K) -
     1            T(ICKWRK(IcKTF + K)) * (RU + RCKWRK(NcK1 + K))
150   CONTINUE
C
C     end of SUBROUTINE CKAML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKAMS  (T, ICKWRK, RCKWRK, AMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKAMS  (T, ICKWRK, RCKWRK, AMS)
C  Returns the standard state Helmholtz free energies in mass units;
C  see Eq. (32).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  AMS(*)    - Real array, standard state Helmholtz free energies
C              for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), AMS(*)
C
      CALL CKSMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
C
      RU = RCKWRK(NcRU)
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TEMP = T(ICKWRK(IcKTF + K))
         AMS(K+1) = RCKWRK(NcK2 + K)
     1          - TEMP * (RU/RCKWRK(NcWT + K) + RCKWRK(NcK1 + K))
150   CONTINUE
C
C     end of SUBROUTINE CKAMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C*****precision > double
      DOUBLE PRECISION FUNCTION CKATOM(ENAME)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION CKATOM(ENAME)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C
C  FUNCTION CKATOM(ENAME)
C  returns atomic weight, given character-string element name.
C
C  Arguments:
C  ENAME   - Character string, element name.
C  AWT     - Real scalar, element atomic weight.
C
C  END PROLOGUE
C
      PARAMETER (NATOM = 102)
      DIMENSION ATOM(NATOM)
      CHARACTER ENAME*(*), IATOM(NATOM)*2, CKCHUP*2
      EXTERNAL CKCHUP
C
      DATA (IATOM(I),ATOM(I),I=1,40) /
     *'H ',  1.00797, 'HE',  4.00260, 'LI',  6.93900, 'BE',  9.01220,
     *'B ', 10.81100, 'C ', 12.01115, 'N ', 14.00670, 'O ', 15.99940,
     *'F ', 18.99840, 'NE', 20.18300, 'NA', 22.98980, 'MG', 24.31200,
     *'AL', 26.98150, 'SI', 28.08600, 'P ', 30.97380, 'S ', 32.06400,
     *'CL', 35.45300, 'AR', 39.94800, 'K ', 39.10200, 'CA', 40.08000,
     *'SC', 44.95600, 'TI', 47.90000, 'V ', 50.94200, 'CR', 51.99600,
     *'MN', 54.93800, 'FE', 55.84700, 'CO', 58.93320, 'NI', 58.71000,
     *'CU', 63.54000, 'ZN', 65.37000, 'GA', 69.72000, 'GE', 72.59000,
     *'AS', 74.92160, 'SE', 78.96000, 'BR', 79.90090, 'KR', 83.80000,
     *'RB', 85.47000, 'SR', 87.62000, 'Y ', 88.90500, 'ZR', 91.22000/
C
      DATA (IATOM(I),ATOM(I),I=41,80) /
     *'NB', 92.90600, 'MO', 95.94000, 'TC', 99.00000, 'RU',101.07000,
     *'RH',102.90500, 'PD',106.40000, 'AG',107.87000, 'CD',112.40000,
     *'IN',114.82000, 'SN',118.69000, 'SB',121.75000, 'TE',127.60000,
     *'I ',126.90440, 'XE',131.30000, 'CS',132.90500, 'BA',137.34000,
     *'LA',138.91000, 'CE',140.12000, 'PR',140.90700, 'ND',144.24000,
     *'PM',145.00000, 'SM',150.35000, 'EU',151.96000, 'GD',157.25000,
     *'TB',158.92400, 'DY',162.50000, 'HO',164.93000, 'ER',167.26000,
     *'TM',168.93400, 'YB',173.04000, 'LU',174.99700, 'HF',178.49000,
     *'TA',180.94800, 'W ',183.85000, 'RE',186.20000, 'OS',190.20000,
     *'IR',192.20000, 'PT',195.09000, 'AU',196.96700, 'HG',200.59000/
C
      DATA (IATOM(I),ATOM(I),I=81,NATOM) /
     *'TL',204.37000, 'PB',207.19000, 'BI',208.98000, 'PO',210.00000,
     *'AT',210.00000, 'RN',222.00000, 'FR',223.00000, 'RA',226.00000,
     *'AC',227.00000, 'TH',232.03800, 'PA',231.00000, 'U ',238.03000,
     *'NP',237.00000, 'PU',242.00000, 'AM',243.00000, 'CM',247.00000,
     *'BK',249.00000, 'CF',251.00000, 'ES',254.00000, 'FM',253.00000,
     *'D ',002.01410, 'E',5.48578E-4/
C
      CKATOM = 0.0
      CALL CKCOMP ( CKCHUP(ENAME, 2), IATOM, NATOM, L)
      IF (L .GT. 0) CKATOM = ATOM(L)
C
C     end of FUNCTION CKATOM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKATHM (NDIM1, NDIM2, ICKWRK, RCKWRK, MAXTP, NT, TMP,
     1                   A)
C
C  START PROLOGUE
C
C  SUBROUTINE CKATHM (NDIM1, NDIM2, ICKWRK, RCKWRK, MAXTP, NT, TMP,
C                     A)
C  Returns the coefficients of the fits for thermodynamic properties
C  of species; see Eqns. (19)-(21).
C
C  INPUT
C  NDIM1     - Integer scalar, first dimension of A, the three-
C              dimensional array of thermodynamic fit coefficients;
C              NDIM1 must be at least NPCP2, the total number of
C              coefficients for one temperature range.
C  NDIM2     - Integer scalar, second dimension of A; NDIM2 must be
C              at least MXTP-1, the total number of temperature ranges.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  MAXTP     - Integer scalar, number of temperatures used to divide
C              the temperature ranges of thermodynamic fits.
C
C  OUTPUT
C  NT(*)     - Integer array, total number of temperatures used in
C              fitting coefficients of thermodynamic properties for
C              the species;
C              dimension at least KK, the total species count.
C  TMP(*,*)  - Real matrix, temperatures for dividing the
C              thermodynamic fits for species; dimension at least
C              MAXTP for the first, and at least KK for the second,
C              the total species count.
C                 cgs units, K
C  A(*,*,*)  - Real three-dimensioned array of fit coefficients to the
C              thermodynamic data for species;
C              dimension exactly NPCP2 for the first, exactly MAXTP-1
C              for the second, and at least KKTOT for the third, the
C              total species count.
C              The indicies in  A(N,L,K) mean-
C              N = 1,NN represent polynomial coefficients in CP/R
C                CP/R(K)=A(1,L,K) + A(2,L,K)*T + A(3,L,K)*T**2 + ...
C              N = NN+1 is for the formation enthalpies, i.e.,
C                HO/R = A(NN+1,L,K)
C              N = NN+2 is for the formation entropies, i.e.,
C                SO/R = A(NN+2,L,K)
C              L = 1 is for temperature <= TMP(2,K)
C              L = 2 is for TMP(2,K) < temperature <= TMP(3)
C                :
C              L = (NTMP-1) is for TMP(NTMP-1) <= temperature;
C              K  is  the  species index
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NT(*), TMP(MAXTP,*),
     1          A(NDIM1,NDIM2,*)
C
      DO 150 K = 1, NKK
         NT(K) = ICKWRK(IcNT + K - 1)
         DO 120 L = 1, MXTP
            TMP(L,K) = RCKWRK(NcTT + (K-1)*MXTP + L - 1)
  120    CONTINUE
         DO 140 L = 1, MXTP-1
            NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
            DO 130 M = 1, NCP2
               A(M, L, K) = RCKWRK(NA1 + M - 1)
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
C
C     end of SUBROUTINE CKATHM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKAWT  (ICKWRK, RCKWRK, AWT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKAWT  (ICKWRK, RCKWRK, AWT)
C  Returns the atomic weights of the elements
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  AWT(*)    - Real array, atomic weights of the elements;
C              dimension at least MM, the total element count.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), AWT(*)
C
      DO 100 M = 1, NMM
         AWT(M) = RCKWRK(NcAW + M - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKAWT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C*****precision > double
      DOUBLE PRECISION FUNCTION CKBSEC (NPTS, X, XX, TT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION CKBSEC (NPTS, X, XX, TT)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C  START PROLOGUE
C  CKBSEC uses bisection to interpolate f(X), given X and other pairs
C  of X and f(X).
C
C  INPUT
C  NPTS	  - Integer scalar, total pairs of data.
C  X      - Real scalar, location for which f(X) is required.
C  XX(*)  - Real array, locations for which data is given.
C  TT(*)  - Real array, function values for locations given.
C
C  END PROLOGUE
C
      DIMENSION XX(NPTS), TT(NPTS)
C
      ZERO = 0.0
C
C     X is outside (1..NPTS)?
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
         NLO = 1
         NHI = NPTS
         S   = ZERO
C
C        bisect interval
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               N = N-1
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               N = N + 1
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
  100 CONTINUE
      CKBSEC      = TT(N) + S * (X - XX(N))
C
C     end of FUNCTION CKBSEC
      RETURN
      END
C                                                                      C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKAVG (NN, S1, S2, SAVG)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER(I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER(I-N)
C*****END precision > single
C
C     START PROLOGUE
C     For arrays of length nn,
C     SAVG(n) is the average value of S1(n) and S2(n).
C     END PROLOGUE
C
      DIMENSION S1(NN), S2(NN), SAVG(NN)
C
      DO 10 N = 1, NN
         SAVG(N) = 0.5 * (S1(N) + S2(N))
   10 CONTINUE
C
C     end of SUBROUTINE CKAVG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCDC  (T, C, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCDC  (T, C, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the nolar creation and destruction rates of the species
C  given temperature(s) and molar concentrations;  see Eq. (76).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), CDOT(*), DDOT(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKRATX (T, C, ICKWRK(IcNU), ICKWRK(IcNK), RCKWRK(NcCO),
     2             ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcFT),
     3             ICKWRK(IcKF), RCKWRK(NcFL), ICKWRK(IcTB),
     3             ICKWRK(IcKN), RCKWRK(NcKT), ICKWRK(IcKT),
     4             RCKWRK(NcKF), RCKWRK(NcKR), RCKWRK(NcI1),
     5             RCKWRK(NcI2), RCKWRK(NcI3), ICKWRK(IcRNU),
     6             RCKWRK(NcRNU), ICKWRK(IcORD), ICKWRK(IcKOR),
     7             RCKWRK(NcKOR), ICKWRK(IcMO), ICKWRK(IcXS))
      CALL CKDOT (RCKWRK(NcI1), RCKWRK(NcI2), ICKWRK, RCKWRK,
     1            CDOT, DDOT)
C
C     end of SUBROUTINE CKCDC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCOPY (NN, X1, X2)
C
C  START PROLOGUE
C  SUBROUTINE CKCOPY (NN, X1, X2)
C  Copy X1(*) array members into X2(*) array.
C
C  INPUT
C  NN        - Integer scalar; number of elements to copy.
C  X1(*)     - Real array; dimension at least NN.
C
C  OUTPUT
C  X2(*)     - Real array; dimension at least NN.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X1(NN), X2(NN)
      DO 10 N = 1, NN
         X2(N) = X1(N)
   10 CONTINUE
C
C     end of CKCOPY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKDOT (RKF, RKR, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDOT (RKF, RKR, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the molar creation and destruction rates of the species
C  given reactions' rates of progress.
C
C  INPUT
C  RKF(*)    - Real array, reactions' forward rates of progress;
C              dimension at least II, the total reaction count.
C  RKR(*)    - Real array, reactions' reverse rates of progress;
C              dimension at least II, the total reaction count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), RKF(*), RKR(*), CDOT(*), DDOT(*)
      INCLUDE 'ckstrt.h'
C
      DO 100 K = 1, NKK
         CDOT(K) = 0.0
         DDOT(K) = 0.0
  100 CONTINUE
C
      I_NK = IcNK
      I_NU = IcNU
C
      DO 200 I = 1, NII
C
C        first and all integer coefficients zero
C        if this reaction has real coefficients
         IF (ICKWRK(I_NU) .EQ. 0) GO TO 200
C
C        there is at least one reactant and one product
         K = ICKWRK(I_NK)
         NU = IABS(ICKWRK(I_NU))
         CDOT(K) = CDOT(K) + RKR(I) * NU
         DDOT(K) = DDOT(K) + RKF(I) * NU
         K = ICKWRK(I_NK + 1)
         IF (K .GT. 0) THEN
            NU = IABS(ICKWRK(I_NU + 1))
            CDOT(K) = CDOT(K) + RKR(I) * NU
            DDOT(K) = DDOT(K) + RKF(I) * NU
            K = ICKWRK(I_NK + 2)
            IF (K .GT. 0) THEN
               NU = IABS(ICKWRK(I_NU + 2))
               CDOT(K) = CDOT(K) + RKR(I) * NU
               DDOT(K) = DDOT(K) + RKF(I) * NU
               K = ICKWRK(I_NK + 3)
               IF (K .GT. 0) THEN
                  NU = IABS(ICKWRK(I_NU + 3))
                  CDOT(K) = CDOT(K) + RKR(I) * NU
                  DDOT(K) = DDOT(K) + RKF(I) * NU
                  K = ICKWRK(I_NK + 4)
                  IF (K .GT. 0) THEN
                     NU = IABS(ICKWRK(I_NU + 4))
                     CDOT(K) = CDOT(K) + RKR(I) * NU
                     DDOT(K) = DDOT(K) + RKF(I) * NU
                     K = ICKWRK(I_NK + 5)
                     IF (K .GT. 0) THEN
                        NU = IABS(ICKWRK(I_NU + 5))
                        CDOT(K) = CDOT(K) + RKR(I) * NU
                        DDOT(K) = DDOT(K) + RKF(I) * NU
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         K = ICKWRK(I_NK + 6)
         NU = ICKWRK(I_NU + 6)
         CDOT(K) = CDOT(K) + RKF(I) * NU
         DDOT(K) = DDOT(K) + RKR(I) * NU
         K = ICKWRK(I_NK + 7)
         IF (K .GT. 0) THEN
            NU = ICKWRK(I_NU + 7)
            CDOT(K) = CDOT(K) + RKF(I) * NU
            DDOT(K) = DDOT(K) + RKR(I) * NU
            K = ICKWRK(I_NK + 8)
            IF (K .GT. 0) THEN
               NU = ICKWRK(I_NU + 8)
               CDOT(K) = CDOT(K) + RKF(I) * NU
               DDOT(K) = DDOT(K) + RKR(I) * NU
               K = ICKWRK(I_NK + 9)
               IF (K .GT. 0) THEN
                  NU = ICKWRK(I_NU + 9)
                  CDOT(K) = CDOT(K) + RKF(I) * NU
                  DDOT(K) = DDOT(K) + RKR(I) * NU
                  K = ICKWRK(I_NK + 10)
                  IF (K .GT. 0) THEN
                     NU = ICKWRK(I_NU + 10)
                     CDOT(K) = CDOT(K) + RKF(I) * NU
                     DDOT(K) = DDOT(K) + RKR(I) * NU
                     K = ICKWRK(I_NK + 11)
                     IF (K .GT. 0) THEN
                        NU = ICKWRK(I_NU + 11)
                        CDOT(K) = CDOT(K) + RKF(I) * NU
                        DDOT(K) = DDOT(K) + RKR(I) * NU
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
  200 CONTINUE
C
      IF (NRNU .LE. 0) RETURN
      I_NU = NcRNU
      DO 300 L = 1, NRNU
         I   = ICKWRK(IcRNU + L - 1)
         I_NK = IcNK + (I-1)*MXSP
         K = ICKWRK(I_NK)
         RNU = ABS(RCKWRK(I_NU))
         CDOT(K) = CDOT(K) + RKR(I) * RNU
         DDOT(K) = DDOT(K) + RKF(I) * RNU
         K = ICKWRK(I_NK + 1)
         IF (K .GT. 0) THEN
            RNU = ABS(RCKWRK(I_NU + 1))
            CDOT(K) = CDOT(K) + RKR(I) * RNU
            DDOT(K) = DDOT(K) + RKF(I) * RNU
            K = ICKWRK(I_NK + 2)
            IF (K .GT. 0) THEN
               RNU = ABS(RCKWRK(I_NU + 2))
               CDOT(K) = CDOT(K) + RKR(I) * RNU
               DDOT(K) = DDOT(K) + RKF(I) * RNU
               K = ICKWRK(I_NK + 3)
               IF (K .GT. 0) THEN
                  RNU = ABS(RCKWRK(I_NU + 3))
                  CDOT(K) = CDOT(K) + RKR(I) * RNU
                  DDOT(K) = DDOT(K) + RKF(I) * RNU
                  K = ICKWRK(I_NK + 4)
                  IF (K .GT. 0) THEN
                     RNU = ABS(RCKWRK(I_NU + 4))
                     CDOT(K) = CDOT(K) + RKR(I) * RNU
                     DDOT(K) = DDOT(K) + RKF(I) * RNU
                     K = ICKWRK(I_NK + 5)
                     IF (K .GT. 0) THEN
                        RNU = ABS(RCKWRK(I_NU + 5))
                        CDOT(K) = CDOT(K) + RKR(I) * RNU
                        DDOT(K) = DDOT(K) + RKF(I) * RNU
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         K = ICKWRK(I_NK + 6)
         RNU = RCKWRK(I_NU + 6)
         CDOT(K) = CDOT(K) + RKF(I) * RNU
         DDOT(K) = DDOT(K) + RKR(I) * RNU
         K = ICKWRK(I_NK + 7)
         IF (K .GT. 0) THEN
            RNU = RCKWRK(I_NU + 7)
            CDOT(K) = CDOT(K) + RKR(I) * RNU
            DDOT(K) = DDOT(K) + RKF(I) * RNU
            K = ICKWRK(I_NK + 8)
            IF (K .GT. 0) THEN
               RNU = RCKWRK(I_NU + 8)
               CDOT(K) = CDOT(K) + RKR(I) * RNU
               DDOT(K) = DDOT(K) + RKF(I) * RNU
               K = ICKWRK(I_NK + 9)
               IF (K .GT. 0) THEN
                  RNU = RCKWRK(I_NU + 9)
                  CDOT(K) = CDOT(K) + RKR(I) * RNU
                  DDOT(K) = DDOT(K) + RKF(I) * RNU
                  K = ICKWRK(I_NK + 10)
                  IF (K .GT. 0) THEN
                     RNU = RCKWRK(I_NU + 10)
                     CDOT(K) = CDOT(K) + RKR(I) * RNU
                     DDOT(K) = DDOT(K) + RKF(I) * RNU
                     K = ICKWRK(I_NK + 11)
                     IF (K .GT. 0) THEN
                        RNU = RCKWRK(I_NU + 11)
                        CDOT(K) = CDOT(K) + RKR(I) * RNU
                        DDOT(K) = DDOT(K) + RKF(I) * RNU
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         I_NU = I_NU + MXSP
  300 CONTINUE
C
C     end of SUBROUTINE CKDOT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCDXP (P, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCDXP (P, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the molar creation and destruction rates of the species
C  given pressure, temperature(s) and mole fractions;  see Eq. (76).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), CDOT(*), DDOT(*)
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKCDC  (T, RCKWRK(NcK4), ICKWRK, RCKWRK, CDOT, DDOT)
C
C     end of SUBROUTINE CKCDXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCDXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCDXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the molar creation and destruction rates of the species
C  given mass density, temperature(s) and mole fractions; see Eq. (76).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), CDOT(*), DDOT(*)
C
      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKCDC  (T, RCKWRK(NcK4), ICKWRK, RCKWRK, CDOT, DDOT)
C
C     end of SUBROUTINE CKCDXR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCDYP (P, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCDYP (P, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the molar creation and destruction rates of the species
C  given mass density, temperature(s) and mass fractions; see Eq. (76).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), CDOT(*), DDOT(*)
C
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKCDC  (T, RCKWRK(NcK4), ICKWRK, RCKWRK, CDOT, DDOT)
C
C     end of SUBROUTINE CKCDYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCDYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCDYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, DDOT)
C  Returns the molar creation and destruction rates of the species
C  given mass density, temperature(s) and mass fractions; see Eq. (76).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  DDOT(*)   - Real array, chemical destruction rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), CDOT(*), DDOT(*)
C
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKCDC  (T, RCKWRK(NcK4), ICKWRK, RCKWRK, CDOT, DDOT)
C
C     end of SUBROUTINE CKCDYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCHRG (ICKWRK, RCKWRK, KCHARG)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCHRG (ICKWRK, RCKWRK, KCHARG)
C  Returns the electronic charges of the species.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  KCHARG(*) - Integer array, electronic charges of the species;
C              dimension at least KK, the total species count.
C              KCHARG(K)=-2 indicates that species K has two
C              excess electrons.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), KCHARG(*)
C
      DO 100 K = 1, NKK
         KCHARG(K) = ICKWRK(IcCH + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKCHRG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCOMP (IST, IRAY, II, I)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCOMP (IST, IRAY, II, I)*
C  Returns the index of an element of a reference character string
C  array which corresponds to a character string;
C  leading and trailing blanks are ignored.
C
C
C  INPUT
C  IST      - Character string; length determined by application
C             program.
C  IRAY(*)  - Character string array; dimension at least II, the total
C             number of character strings for be searched.
C  II       - Integer scalar, the length of IRAY to be searched.
C
C  OUTPUT
C  I        - Integer scalar, the first array index in IRAY of a
C             character string IST, or 0 if IST is not found.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) IST, IRAY(*)
      INTEGER CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
C
      I = 0
      IS1 = CKFRCH(IST)
      IS2 = CKLSCH(IST)
      ISLEN = IS2 - IS1 + 1
      ILEN = LEN(IRAY(1))
      IF (ILEN .LT. ISLEN) RETURN
      DO 10 N = 1, II
         IR1 = CKFRCH(IRAY(N)(1:ILEN))
         IR2 = CKLSCH(IRAY(N)(1:ILEN))
         IF (IR2-IR1+1 .NE. ISLEN) GO TO 10
         IF (IST(IS1:IS2) .EQ. IRAY(N)(IR1:IR2)) THEN
            I = N
            RETURN
         ENDIF
   10 CONTINUE
C
C     end of SUBROUTINE CKCOMP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNCMP (STR, IRAY, II, I, NF)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNCMP (STR, IRAY, II, I, NF)
C  Returns the first index of the character string STR if it occurs
C  in the character string IRAY, and returns the total number of
C  times STR occurs in IRAY.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) STR, IRAY(II)
      INTEGER CKFRCH, CKLSCH
      EXTERNAL CKFRCH, CKLSCH
C
      I = 0
      NF = 0
      IS1 = CKFRCH(STR)
      IS2 = CKLSCH(STR)
      DO 10 N = II, 1, -1
         IR1 = CKFRCH(IRAY(N))
         IR2 = CKLSCH(IRAY(N))
         IF (IS2.GE.IS1 .AND. IS2.GT.0 .AND.
     1       IR2.GE.IR1 .AND. IR2.GT.0 .AND.
     2       STR(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) THEN
             I = N
             NF = NF + 1
         ENDIF
   10 CONTINUE
C
C     end of SUBROUTINE CKNCMP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCONT (K, Q, ICKWRK, RCKWRK, CIK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCONT (K, Q, ICKWRK, RCKWRK, CIK)
C  Returns the contributions of the reactions to the molar production
C  rate of a species;  see Eqs. (49) and (51).
C
C  INPUT
C  K         - Integer scalar; species index number.
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CIK(*)    - Real array, contributions of the reactions to the
C              production rate of species K;
C              dimension least II, the total reaction count.
C                 cgs units, mole/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION Q(*), ICKWRK(*), RCKWRK(*), CIK(*)
C
      DO 100 I = 1, NII
         CIK(I) = 0.0
  100 CONTINUE
C
      I_NK = IcNK - 1
      I_NU = IcNU - 1
      DO 200 I = 1, NII
         DO 180 N = 1, MXSP
            IF (ICKWRK(I_NK+N) .EQ. K)
     1            CIK(I) = CIK(I) + ICKWRK(I_NU+N) * Q(I)
  180    CONTINUE
         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
  200 CONTINUE

      IF (NRNU .GT. 0) THEN
         I_NU = NcRNU - 1
         DO 300 L = 1, NRNU
            I = ICKWRK(IcRNU + L - 1)
            I_NK = IcNK + MXSP*(I-1) - 1
            DO 280 N = 1, MXSP
               IF (ICKWRK(I_NK+N) .EQ. K)
     1            CIK(I) = CIK(I) + RCKWRK(I_NU+N) * Q(I)
  280       CONTINUE
            I_NU = I_NU + MXSP
  300    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKCONT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCPBL (T, X, ICKWRK, RCKWRK, CPBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPBL (T, X, ICKWRK, RCKWRK, CPBML)
C  Returns the mean specific heat at constant pressure in molar units;
C  see Eq. (33).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CPBML     - Real scalar, mean specific heat at constant pressure.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKCPML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      CPBML = 0.0
      DO 100 K = 1, NKK
         CPBML = CPBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKCPBL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCPBS (T, Y, ICKWRK, RCKWRK, CPBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPBS (T, Y, ICKWRK, RCKWRK, CPBMS)
C  Returns the mean specific heat at constant pressure; see Eq. (34).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CPBMS     - Real scalar, mean specific heat at constant pressure.
C                 cgs units - ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKCPMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      CPBMS = 0.0
      DO 100 K = 1, NKK
         CPBMS = CPBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKCPBS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
C  Returns the specific heats at constant pressure in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CPML(*)   - Real array, specific heats at constant pressure for
C              the species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), CPML(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPML(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1)
     1             + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2             + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4)
250   CONTINUE
C
C     end of SUBROUTINE CKCPML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCPMS (T, ICKWRK, RCKWRK, CPMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPMS (T, ICKWRK, RCKWRK, CPMS)
C  Returns the specific heats at constant pressure in mass units;
C  see Eq. (26).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CPMS(*)   - Real array, specific heats at constant pressure for
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), CPMS(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPMS(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1)
     1             + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2             + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4)
     3             / RCKWRK(NcWT + K)
250   CONTINUE
C
C     end of SUBROUTINE CKCPMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCPOR (T, ICKWRK, RCKWRK, CPOR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPOR (T, ICKWRK, RCKWRK, CPOR)
C  Returns the nondimensional specific heats at constant pressure;
C  see Eq. (19).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CPOR(*)   - Real array, nondimensional specific heats at constant
C              pressure for species;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), CPOR(*)
      SAVE TN1, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         CPOR(K+1) = RCKWRK(NA1)
     1             + RCKWRK(NA1+1)*TK1 + RCKWRK(NA1+2)*TK2
     2             + RCKWRK(NA1+3)*TK3 + RCKWRK(NA1+4)*TK4
250   CONTINUE
C
C     end of SUBROUTINE CKCPOR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCRAY (LINE, NN, KRAY, LOUT, NDIM, NRAY, NF, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCRAY (LINE, NN, KRAY, LOUT, NDIM, NRAY, NF, KERR)
C  Searches a character string, LINE, and compares the space-delimited
C  substrings in LINE, to an array of character strings, KRAY;
C  if a substring in LINE is located in KRAY, the index of its location
C  in KRAY is stored in the integer array NRAY.  For example, the
C  subroutine might be called to assign Chemkin species indices to a
C  given list of species names.  This application is illusgrated in the
C  following example:
C
C     input:  LINE    = "OH  N2  NO"
C             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
C             NN      = 9, the number of entries in KRAY(*)
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C             NDIM    = 10, the dimension of array NRAY(*)
C     output: NRAY(*) = 7, 3, 9, the index numbers of the entries
C                       in KRAY(*) corresponding to the substrings
C                       in LINE
C             NF      = 3, the number of correspondences found.
C             KERR    = .FALSE.
C
C  INPUT
C  LINE    - Character string.
C  KRAY(*) - Character string array; dimension at least NN.
C  NN      - Integer scalar, total character string count of KRAY.
C  LOUT    - Integer scalar, formatted output file unit.
C  NDIM    - Integer scalar, dimension of the integer array NRAY.
C
C  OUTPUT
C  NRAY(*) - Integer array, indices of the elements of KRAY
C            which correspond to the substrings in LINE;
C            dimension at least NDIM.
C  NF      - Integer scalar, count of correspondences found.
C  KERR    - Logical, syntax or dimensioning Error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), KRAY(*)*(*), SUB(80)*80
      DIMENSION NRAY(*)
      LOGICAL KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      NF = 0
C
      IDIM = 80
      CALL CKSUBS (LINE, LOUT, IDIM, SUB, NFOUND, IERR)
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT,*) ' Error in CKCRAY...'
         RETURN
      ENDIF
C
      DO 50 N = 1, NFOUND
         CALL CKCOMP (SUB(N), KRAY, NN, K)
         IF (K .LE. 0) THEN
            LT = MAX (CKLSCH(SUB(N)), 1)
            WRITE (LOUT,'(A)')
     1      ' Error in CKCRAY...'//SUB(N)(1:LT)//' not found...'
            KERR = .TRUE.
         ELSE
            IF (NF+1 .GT. NDIM) THEN
               WRITE (LOUT,'(A)')
     1       ' Error in CKCRAY...dimension of NRAY too small...'
               KERR = .TRUE.
            ELSE
               NF = NF + 1
               NRAY(NF) = K
            ENDIF
         ENDIF
   50 CONTINUE
C
C     end of SUBROUTINE CKCRAY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTC  (T, C, ICKWRK, RCKWRK, CDOT, TAU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTC  (T, C, ICKWRK, RCKWRK, CDOT, TAU)
C  Returns the molar creation rates and characteristic destruction
C  times of the species given temperature(s) and molar concentrations;
C  see Eqs. (79) and (81).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  TAU(*)    - Real array, characteristic destruction times of species;
C              dimension at least KK, the total species count.
C                 cgs units, sec
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), CDOT(*), TAU(*)
C
      CALL CKCDC (T, C, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      DO 150 K = 1, NKK
         TAU(K) = C(K) / (RCKWRK(NcK1 + K - 1)+SMALL)
150   CONTINUE
C
C     end of SUBROUTINE CKCTC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTX  (C, ICKWRK, RCKWRK, X)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTX  (C, ICKWRK, RCKWRK, X)
C  Returns the mole fractions given molar concentrations; see Eq. (13).
C
C  INPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                   cgs units - mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  X(*)      - Real array, mole fraction of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION C(*), ICKWRK(*), RCKWRK(*), X(*)
C
      CTOT = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
  100 CONTINUE
      DO 200 K = 1, NKK
         X(K) = C(K)/CTOT
200   CONTINUE
C
C     end of SUBROUTINE CKCTX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTXP (P, T, X, ICKWRK, RCKWRK, CDOT, TAU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTXP (P, T, X, ICKWRK, RCKWRK, CDOT, TAU)
C  Returns the molar creation rates and characteristic destruction
C  times of the species given pressure, temperature(s) and mole
C  fractions;  see Eqs. (79) and (81).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  TAU(*)    - Real array, characteristic destruction times of the
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, sec
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), CDOT(*), TAU(*)
C
      CALL CKCDXP (P, T, X, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TAU(K+1) = RCKWRK(NcK2 + K) / (RCKWRK(NcK1 + K)+SMALL)
150   CONTINUE
C
C     end of SUBROUTINE CKTCXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, TAU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, TAU)
C  Returns the molar creation rates and characteristic destruction
C  times of the species given mass density, temperature(s) and mole
C  fractions;  see Eqs. (79) and (81).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  TAU(*)    - Real array, characteristic destruction times of species;
C              dimension at least KK, the total species count.
C                 cgs units, sec
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), CDOT(*), TAU(*)
C
      CALL CKCDXR (RHO, T, X, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TAU(K+1) = RCKWRK(NcK2 + K) / (RCKWRK(NcK1 + K)+SMALL)
150   CONTINUE
C
C     end of SUBROUTINE CKCTXR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTY  (C, ICKWRK, RCKWRK, Y)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTY  (C, ICKWRK, RCKWRK, Y)
C  Returns the mass fractions given molar concentrations; see Eq. (12).
C
C  INPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION C(*), ICKWRK(*), RCKWRK(*), Y(*)
C
      RHO = 0.0
      DO 100 K = 1, NKK
         RHO = RHO + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      DO 200 K = 1, NKK
         Y(K) = C(K) * RCKWRK(NcWT + K - 1)/RHO
200   CONTINUE
C
C     end of SUBROUTINE CKCTY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTYP (P, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTYP (P, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
C  Returns the molar creation rates and characteristic destruction
C  times of the species given mass density, temperature(s) and mass
C  fractions;  see Eqs. (79) and (81).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  TAU(*)    - Real array, characteristic destruction times of the
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, sec
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), CDOT(*), TAU(*)
C
      CALL CKCDYP (P, T, Y, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TAU(K+1) = RCKWRK(NcK2 + K) / (RCKWRK(NcK1 + K)+SMALL)
150   CONTINUE
C
C     end of SUBROUTINE CKCTYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCTYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCTYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, TAU)
C  Returns the molar creation rates and characteristic destruction
C  times of the species given mass density, temperature(s) and mass
C  fractions;  see Eqs. (79) and (81).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CDOT(*)   - Real array, chemical creation rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/(cm**3*sec)
C  TAU(*)    - Real array, characteristic destruction times of the
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, sec
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), CDOT(*), TAU(*)
C
      CALL CKCDYR (RHO, T, Y, ICKWRK, RCKWRK, CDOT, RCKWRK(NcK1))
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TAU(K+1) = RCKWRK(NcK2 + K) / (RCKWRK(NcK1 + K)+SMALL)
150   CONTINUE
C
C     end of SUBROUTINE CKCTYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCVBL (T, X, ICKWRK, RCKWRK, CVBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCVBL (T, X, ICKWRK, RCKWRK, CVBML)
C  Returns the mean specific heat at constant volume in molar units;
C  see Eq. (35).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CVBML     - Real scalar, mean specific heat at constant volume.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKCVML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      CVBML = 0.0
      DO 100 K = 1, NKK
         CVBML = CVBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKCVBL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCVBS (T, Y, ICKWRK, RCKWRK, CVBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCVBS (T, Y, ICKWRK, RCKWRK, CVBMS)
C  Returns the mean specific heat at constant volume in mass units;
C  see Eq. (36).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CVBMS     - Real scalar, mean specific heat at constant volume.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKCVMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      CVBMS = 0.0
      DO 100 K = 1, NKK
         CVBMS = CVBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKCVBS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCVML (T, ICKWRK, RCKWRK, CVML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCVML (T, ICKWRK, RCKWRK, CVML)
C  Returns the specific heats in constant volume in molar units;
C  see Eq. (22).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CVML(*)   - Real array, specific heats at constant volume for
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), CVML(*)
C
      CALL CKCPML (T, ICKWRK, RCKWRK, CVML)
C
      DO 150 K = 1, NKK
         CVML(K) = CVML(K) - RCKWRK(NcRU)
150   CONTINUE
C
C     end of SUBROUTINE CKCVML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKCVMS (T, ICKWRK, RCKWRK, CVMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCVMS (T, ICKWRK, RCKWRK, CVMS)
C  Returns the specific heats at constant volume in mass units;
C  see Eq. (29).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  CVMS(*)   - Real array, specific heats at constant volume for
C              species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), CVMS(*)
C
      CALL CKCPMS (T, ICKWRK, RCKWRK, CVMS)
C
      RU = RCKWRK(NcRU)
      DO 150 K = 1, NKK
         CVMS(K) = CVMS(K) - RU / RCKWRK(NcWT + K - 1)
150   CONTINUE
C
C     end of SUBROUTINE CKCVMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKECON (P, T, X, XNUEH, ICKWRK, RCKWRK, CONE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKECON (P, T, X, XNUEH, ICKWRK, RCKWRK, CONE)
C  Returns the electron species thermal conductivity given collision
C  frequency.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUEH     - Real scalar, total of the momentum-transfer
C              collision frequencies for electrons.
C
C  OUTPUT
C  CONE      - Real scalar, electron thermal conductivity
C                 cgs units, ERG/CM*K*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      EDEN = TDEN * MAX(X(KEL),SMALL)
      EMASS = RCKWRK(NcWT+KEL-1)/AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      CONE = 2.4*(BOLTZ**2)*EDEN*TE/(EMASS*XNUEH)
C
C     end of SUBROUTINE CKECON
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEDIF (T, XNUES, ICKWRK, RCKWRK, DEK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEDIF (T, XNUES, ICKWRK, RCKWRK, DEK)
C  Returns the electron species binary diffusion coefficients given
C  collision frequency.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  XNUES(*)  - Real array, momentum-transfer collision frequency
C              for electrons with each species;
C              dimension at least KK, the total species count.
C                 cgs units, /S
C  OUTPUT
C  DEK(*)    - Real array, electron binary diffusion coefficients;
C              dimension at least KK, the total species count.
C                 cgs units, ERG/CM*K*S
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), XNUES(*), ICKWRK(*), RCKWRK(*), DEK(*)
C
      EMASS = RCKWRK(NcWT+KEL-1)/AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      DO 150 K = 1, NKK
         IF (EMASS*XNUES(K) .GT. SMALL) THEN
            DEK(K) = (BOLTZ*TE)/(EMASS*XNUES(K))
         ELSE
            DEK(K) = BIG
         ENDIF
150   CONTINUE
C
C     end of SUBROUTINE CKEDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQ (RCKWRK, ICKWRK, T, NU, NUNK, PAR, IREV, RPAR,
     1                 ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU, IEIM,
     2                 IEIMT, IJAN, PJAN, IFT1, PF1, ITDE, ITDK,
     3                 KTFL, RKFT, RKRT, EQK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQ (RCKWRK, ICKWRK, T, NU, NUNK, PAR, IREV, RPAR,
C                   ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU, IEIM,
C                   IEIMT, IJAN, PJAN, IFT1, PF1, ITDE, ITDK,
C                   KTFL, RKFT, RKRT, EQK)
C
C  INPUT
C
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C              IREV(N) is the reaction index for the Nth reaction
C              with explicit reverse parameters.
C  RPAR(*,*) - Real matrix,  reverse Arrhenius rate coefficients for
C              the NREV reactions; dimension at least NPAR for the
C              first, and at least NREV for the second.
C              RPAR(L,N) is the Lth coefficient for the Nth reaction
C              reaction with explicit reverse rate coefficients.
C                 cgs units same as PAR(*,*)
C  ILAN(*)   - Integer array, reaction indices for the NLAN reactions;
C              dimension at least NLAN.
C              ILAN(N) is the reaction index for the Nth Landau-
C              Teller reaction.
C  PLT(*,*)  - Real matrix, the additional parameters for the NLAN
C              reactions; dimension at least NLAR for the first and at
C              least NLAN for the second.
C              PLAN(L,N) is the Lth parameter for the Nth
C              Landau-Teller reaction.
C  IRLT(*)   - Integer array, reaction indices for the NRLT reactions;
C              dimension at least NRLT.
C              IRLT(N) is the reaction index for the Nth reaction
C              with Landau-Teller reverse parameters.
C  RPLT(*,*) - Real matrix, the additional rate parameters for the
C              NRLT reactions; dimension at least NLAR for the first
C              and at least NRLT for the second.
C              RPLT(L,N) is the Lth reverse parameter for the Nth
C              Landau-Teller reaction with reverse parameters.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  IEIM(*)   - Integer array, reaction indices for the NEIM reactions;
C              dimension at least NEIM.
C              IEIM(N) is the reaction index for the Nth electron-
C              impact reaction.
C  IEIMT(*)  - Integer array, temperature-dependence indices for the
C              NEIM reactions; dimension at least NEIM.
C              IEIMT(N) is a pointer into a temperature array for
C              the Nth electron-impact reaction.
C  IJAN(*)   - Integer array, reaction indices of the NJAN reactions;
C              dimension at least NJAN.
C              IJAN(N) is the reaction index for the Nth Janev et al.
C              reaction.
C  PJAN(*,*) - Real matrix, rate parameters for the NJAN reactions;
C              dimension at least NJAR for the first and NJAN for the
C              second.
C              PJAN(L,N) is the Lth parameter for the Nth Janev et al.
C              reaction.
C  IFT1(*)   - Integer array, reaction indices for the NFT1 reactions;
C              dimension at least NFT1.
C              IFT1(N) is the reaction index for the Nth fit-type
C              reaction.
C  PF1(*,*)  - Real matrix, the additional rate parameters for the
C              NFT1 reactions; dimension at least NF1R for the first
C              and at least NFT1 for the second.
C              PF1(L,N) is the Lth fit parameter for the Nth fit-type
C              reaction.
C  ITDE(*)   - Integer array, reaction indices for the NTDE reactions;
C              dimension at least NTDE.
C              ITDE(N) is the reaction index for the Nth non-
C              thermal-equilibrium reaction.
C  ITDK(*)   - Integer array, special species for the NTDE reactions;
C              dimension at least NTDE.
C              ITDK(N) is the indentifying species, K, whose associated
C              temperature (KTFL(K)) should be used in calculating
C              quantities associated with the reaction.
C  KTFL(*)   - Integer array, indices into the temperature array,
C              for the KK species, as required by the NTDE reactions;
C              dimension at least KK.
C              KTFL(K) is the temperature array index for species K.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT; dimension at least KK.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units depend on the reaction
C                 NOTE: The values in RKFT(I) returned by this routine
C                       aren't correct
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units depend on the reaction
C                 NOTE: The values in RKRT(I) returned by this routine
C                       aren't correct
C  EQK(*)    - Real array, equilibrium constants in concentration
C              units for reactions; dimension at least II.
C              EQKC(I) is the equilibrium constant for reaction I.
C                 cgs units (mole/cm**3)**some power,
C                 depends on reaction
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
C     Integer arrays
      DIMENSION ICKWRK(*), NU(MXSP,NII), NUNK(MXSP,NII), IREV(NREV),
     1          ILAN(NLAN), IRLT(NRLT), IRNU(NRNU), IEIM(NEIM),
     2          IEIMT(NEIM), IJAN(NJAN), IFT1(NFT1), ITDE(NTDE),
     3          ITDK(NTDE), KTFL(*),
C     Real arrays
     4          RCKWRK(*), PAR(NPAR+1,NII), RPAR(NPAR+1,NREV),
     2          PLT(NLAR,NLAN), RPLT(NLAR,NRLT), SMH(NKK),
     3          RKFT(NII), RKRT(NII), EQK(NII), RNU(MXSP,NRNU),
     4          T(*), PJAN(NJAR,NJAN), PF1(NF1R,NFT1)
      COMMON /MACH/ SMALL, BIG, EXPARG
      INTEGER CKLKUP
      EXTERNAL CKSMH, CKLKUP
C
C     Find Gibbs/ RT for all species in the mechanism
C     Note: CKSMH takes a vector of temperatures, and the G/RT
C           value for each species will be evaluated at the
C           temperature associated with that species.
C
      CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
      TAV = T(1)
      ALOGT = LOG(TAV)
      RU = RCKWRK(NcRU)
      PATM = RCKWRK(NcPA)
      PFAC = PATM / (RU * TAV)
C
C.....Default way to calculate the equilibrium constant.................
C
C     Loop over NII -> largest expense in the routine, so optimize)
C
      DO 50 I = 1, NII
C
         IF (NUNK(1,I) .EQ. 0) GO TO 50
C         Initialize the net mole change number and the
C         the net DeltaG value with the contributions from
C         the first reactant and product
C
        NUSUMK = NU(1,I) + NU(7,I)
        SUMSMH = NU(1,I)*SMH(NUNK(1,I)) + NU(7,I)*SMH(NUNK(7,I))
        IF (NUNK(2,I) .NE. 0) THEN
          NUSUMK = NUSUMK + NU(2,I)
          SUMSMH = SUMSMH + NU(2,I)*SMH(NUNK(2,I))
          IF (NUNK(3,I) .NE. 0) THEN
            NUSUMK = NUSUMK + NU(3,I)
            SUMSMH = SUMSMH + NU(3,I)*SMH(NUNK(3,I))
            IF (NUNK(4,I) .NE. 0) THEN
              NUSUMK = NUSUMK + NU(4,I)
              SUMSMH = SUMSMH + NU(4,I)*SMH(NUNK(4,I))
              IF (NUNK(5,I) .NE. 0) THEN
                NUSUMK = NUSUMK + NU(5,I)
                SUMSMH = SUMSMH + NU(5,I)*SMH(NUNK(5,I))
                IF (NUNK(6,I) .NE. 0) THEN
                  NUSUMK = NUSUMK + NU(6,I)
                  SUMSMH = SUMSMH + NU(6,I)*SMH(NUNK(6,I))
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF (NUNK(8,I) .NE. 0) THEN
          NUSUMK = NUSUMK + NU(8,I)
          SUMSMH = SUMSMH + NU(8,I)*SMH(NUNK(8,I))
          IF (NUNK(9,I) .NE. 0) THEN
            NUSUMK = NUSUMK + NU(9,I)
            SUMSMH = SUMSMH + NU(9,I)*SMH(NUNK(9,I))
            IF (NUNK(10,I) .NE. 0) THEN
              NUSUMK = NUSUMK + NU(10,I)
              SUMSMH = SUMSMH + NU(10,I)*SMH(NUNK(10,I))
              IF (NUNK(11,I) .NE. 0) THEN
                NUSUMK = NUSUMK + NU(11,I)
                SUMSMH = SUMSMH + NU(11,I)*SMH(NUNK(11,I))
                IF (NUNK(12,I) .NE. 0) THEN
                  NUSUMK = NUSUMK + NU(12,I)
                  SUMSMH = SUMSMH + NU(12,I)*SMH(NUNK(12,I))
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C
C        Calculate the concentration equilibrium constant,
C        Protecting against overflow in the exponential
C
         IF (NUSUMK .NE. 0) THEN
           EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**NUSUMK)
         ELSE
           EQK(I) = EXP(MIN(SUMSMH,EXPARG))
         ENDIF
   50 CONTINUE
C
C.....Fix-up look for rxn's with real-stoichiometry.....................
C
C     Must completely redo the calculation of the Gibbs free energy
C     of reaction because the stoichiometric coefficients have been
C     overridden.
C
      DO 70 N = 1, NRNU
         I = IRNU(N)
         RNUSUM = RNU(1,N) + RNU(7,N)
         SUMSMH = RNU(1,N)*SMH(NUNK(1,I)) + RNU(7,N)*SMH(NUNK(7,I))
         DO 60 L = 2, 6
            IF (NUNK(L,I) .EQ. 0) GO TO 61
            SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
            RNUSUM = RNUSUM + RNU(L,N)
   60    CONTINUE
   61    CONTINUE
         DO 62 L = 8, 12
            IF (NUNK(L,I) .EQ. 0) GO TO 63
            SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
            RNUSUM = RNUSUM + RNU(L,N)
   62    CONTINUE
   63    CONTINUE
         IF (RNUSUM .NE. 0.0) THEN
            EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**RNUSUM)
         ELSE
            EQK(I) = EXP(MIN(SUMSMH,EXPARG))
         ENDIF
   70 CONTINUE
C
C.....Rxns with explicitly defined reverse rate constants...............
C
C     Fix-up loop for reactions which have a defined reverse
C     rate constant. We will return the ratio of the forward and
C     reverse rate constants in this case
C
      DO 90 N = 1, NREV
         I = IREV(N)
C        Default K_f (forward rate) is Arrhenius expression
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/TAV)
C        K_r for reaction I, with explicit parameters
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)/TAV)
C        new K_eq = K_f / new K_r
         EQK(I)  = RKFT(I) / MAX(RKRT(I), SMALL)
   90 CONTINUE
C
C.....Landau-Teller reactions...........................................
C
C     Need to fix up equilibirum constant for REV keyword case only
C
      DO 200 N = 1, NLAN
         I = ILAN(N)
C
C        Lookup whether this is an Landau-Teller Rxn with an
C        explicitly given reverse reaction rate
C        - Note: If this is such a rxn, The forward and reverse rate
C                constants have already been calculated in the loop
C                above. However, additional modifications to RKFT and
C                RKRT must be made, and then the equilibrium constant
C                must be recalculated,
C
         ISRLT = CKLKUP(I, IRLT, NRLT)
         IF (ISRLT .NE. 0) THEN
C
C           Calculate the Modified forward rxn rate for Landau
C           Teller Rxn's
C
            TFAC = PLT(1,N)/TAV**(1.0/3.0) + PLT(2,N)/TAV**(2.0/3.0)
            RKFT(I) = RKFT(I) * EXP(TFAC)
C
C           Calculate the Modified, Explicitly Given Reverse rxn
C           rate for Landau Teller Rxn's
C
            TFAC =   RPLT(1,ISRLT)/TAV**(1.0/3.0)
     1             + RPLT(2,ISRLT)/TAV**(2.0/3.0)
            RKRT(I) = RKRT(I) * EXP(TFAC)
C           new K_eq = new K_f / new K_r
            EQK(I) = RKFT(I) / MAX(RKRT(I), SMALL)
         ENDIF
  200 CONTINUE
C
C.....Electron-impact reactions.........................................
C
C     Need to fix up the equilibrium constant for REV and non-REV case
C     The non-REV case may have used concentrations based on the wrong
C     temperature. The REV case may have calculated the forward and
C     reverse rate constants based on the wrong temperature.
C
      DO 300 N = 1, NEIM
C
C       No change in equilibrium constant unless temperature other
C       than first temperature is associated with this reaction.
C
         IF (IEIMT(N) .EQ. 1) GO TO 300
C
C        Look up the temperature and global rxn number for this
C        reaction
C
         TEMP = T(IEIMT(N))
         IF (TEMP .EQ. TAV) GO TO 300
         I = IEIM(N)
C
C        Lookup whether this rxn has explicit Reverse Parameter
C        parameters and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C           new K_f based on the correct temperature
            RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C           new K_r with explicit parameters and new temperature
            RKRT(I) = RPAR(1,ISREV)
     1               * EXP(RPAR(2,ISREV)*LOG(TEMP) - RPAR(3,ISREV)/TEMP)
C           new K_eq = new K_f / new K_r
            EQK(I) = RKFT(I) / MAX(RKRT(I),SMALL)
         ELSE
C            modify concentration based on new temperature
           PFAC2 = PATM / (RU*TEMP)
C
C          Fix concentration term in equilibrium constant, using the
C          previously calculated Gibbs Free energy part, which is
C          still good.
C
           ISREAL = CKLKUP (I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
              L = ISREAL
              RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1             + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2             + RNU(11,L)+RNU(12,L)
              IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
              NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1             + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
              IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF

         ENDIF
C
  300 CONTINUE
C
C.....Non-thermal-equilibrium, species-temperature-dependent reactions..
C
C       Need to fix up the equilibrium constant for REV and non-REV case
C       The non-REV case may have used concentrations based on the wrong
C       temperature. The REV case may have calculated the forward and
C       reverse rate constants based on the wrong temperature.
C
      DO 400 N = 1, NTDE
C
C          Find the temperature index for the special species in the
C          reaction
C
         KTEMP = KTFL(ITDK(N))
C
C          No change in equilibrium constant unless temperature other
C          than first temperature is associated with this reaction.
C
         IF (KTEMP .EQ. 1) GO TO 400
         TEMP = T(KTEMP)
         IF (TEMP .EQ. TAV) GO TO 400
         I = ITDE(N)
C
C          Lookup whether the reaction has explicit reverse rate
C          constants and branch accordingly.
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C            new K_f for reaction I, Arrhenius expression with new
C            temperature
           RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C             new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
C          new K_eq = new K_f / new K_r
           EQK(I) = RKFT(I) / MAX(RKRT(I),SMALL)
         ELSE
C
C          modify concentration based on new temperature
C
           PFAC2 = PATM / (RU*TEMP)
C
C            Look up whether the reaction has real coefficients
C            and branch accordingly
C
           ISREAL = CKLKUP(I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
              L = ISREAL
              RNUSUM=RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1            +RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2            +RNU(11,L) +RNU(12,L)
              IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
              NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1             + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
              IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF
         ENDIF
  400 CONTINUE
C
C.....Reactions using fit#1:  k = A * T^B * exp(v1/T+v2/T^2+v3/T^3...)..
C
C       Only need to handle the REV keyword case in this block. The
C       possibility that the concentration part of the equilibrium
C       constant may have a different temperature associated with it
C       has already been handled by the 300 and 400 do loops.
C
      DO 500 N = 1, NFT1
         IF (N .EQ. 1) BLOG = LOG(BIG)
         I = IFT1(N)
C
C          Lookup rxn in IREV to see if it has explicit reverse rate
C          constant parameters -> return the position in IREV or 0
C          If this isn't a REV rxn, bail out
C
         ISREV = CKLKUP (I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C
C            Species temperature-dependent reaction?
C
           N2 = CKLKUP(I, ITDE, NTDE)
           IF (N2 .EQ. 0) THEN
             TEMP = TAV
           ELSE
             TEMP = T(KTFL(ITDK(N2)))
           ENDIF
C
C            Electron-impact reaction?
C
           N2 = CKLKUP(I, IEIM, NEIM)
           IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C            new K_f based on fit#1, possibly different temperature
C
           RKFT(I) = PAR(1,I) * TEMP**PAR(2,I)
           SUMJ = 0.0
           DO 470 J = 1, NF1R
            ROOTJ = 1./J
            IF (TEMP.GE.BIG**ROOTJ .OR. SUMJ .GE. BLOG) THEN
               SUMJ = BLOG
            ELSE
               SUMJ = SUMJ + PF1(J,N)/TEMP**J
            ENDIF
 470       CONTINUE
           SUMJ = MIN (SUMJ, BLOG)
           RKFT(I) = MIN(BIG, RKFT(I) * EXP(SUMJ))
C
           IF (TEMP .NE. TAV) THEN
C              new K_r with explicit parameters and new temperature
             RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1               - RPAR(3,ISREV)/TEMP)
           ENDIF
C          new K_eq = new K_f / new K_r
           EQK(I) = RKFT(I) / MAX(RKRT(I),SMALL)
         ENDIF
  500 CONTINUE
C
C.....jannev, langer, evans & post - type reactions.....................
C
C       Only need to handle the REV keyword case in this block. The
C       possibility that the concentration part of the equilibrium
C       constant may have a different temperature associated with it
C       has already been handled by the 300 and 400 do loops.
C
      DO 600 N = 1, NJAN
         I = IJAN(N)
C
C        Lookup rxn in IREV to see if it has explicit reverse rate
C        constant parameters -> return the position in IREV or 0
C        Bail out if this isn't a REV reaction.
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C
C            Specify the temperature to be used below by looking up
C            whether this rxn has a special temperature variable
C            associated with it. If it doesn't, use the default first
C            temperature variable.
C
           N2 = CKLKUP(I, ITDE, NTDE)
           IF (N2 .EQ. 0) THEN
             TEMP = TAV
           ELSE
             TEMP = T(KTFL(ITDK(N2)))
           ENDIF
C
C            Lookup up whether this rxn is an electron impact rxn
C            If it is, override the temperature variable determined
C            in the previous loop.
C
           N2 = CKLKUP (I, IEIM, NEIM)
           IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C            Re-evaluate Arrhenius K_f, possibly different temperature,
C            then modify for jannev expression
C
           RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C            convert E- temperature to eV's
           TEV = TEMP / 11595.0
           SUMJ = 0.0
           DO 530 J = 1, NJAR
              SUMJ = SUMJ + PJAN(J,N) * (LOG(TEV))**(J-1)
  530      CONTINUE
           RKFT(I) =  MIN(BIG, RKFT(I) * EXP(SUMJ))

           IF (TEMP .NE. TAV) THEN
C              new K_r with explicit parameters and new temperature
             RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1               - RPAR(3,ISREV)/TEMP)
           ENDIF
C            new K_eq = new K_f / new K_r
           EQK(I) = RKFT(I) / MAX(RKRT(I),SMALL)
         ENDIF
 600  CONTINUE
C
C     end of SUBROUTINE CKEQ
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQC  (T, C, ICKWRK, RCKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQC  (T, C, ICKWRK, RCKWRK, EQKC)
C  Returns the equilibrium constants of the reactions given
C  temperature(s) and molar concentrations;  see Eq. (54).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants in concentration units
C              for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, (mole/cm**3)**some power, depending on
C                            the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), C(*), EQKC(*), T(*)
      INCLUDE 'ckstrt.h'
C
      CALL CKEQ (RCKWRK, ICKWRK, T, ICKWRK(IcNU), ICKWRK(IcNK),
     1           RCKWRK(NcCO), ICKWRK(IcRV), RCKWRK(NcRV), ICKWRK(IcLT),
     2           RCKWRK(NcLT), ICKWRK(IcRL),   RCKWRK(NcRL),
     3           RCKWRK(NcK1), ICKWRK(IcRNU),  RCKWRK(NcRNU),
     4           ICKWRK(IcEI), ICKWRK(IcET),   ICKWRK(IcJN),
     5           RCKWRK(NcJN), ICKWRK(IcF1),   RCKWRK(NcF1),
     6           ICKWRK(IcTD), ICKWRK(IcTK),   ICKWRK(IcKTF),
     7           RCKWRK(NcKF), RCKWRK(NcKR),   EQKC)
C
C     end of SUBROUTINE CKEQC
      RETURN
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQXP (P, T, X, ICKWRK, RCKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQXP (P, T, X, ICKWRK, RCKWRK, EQKC)
C  Returns the equilibrium constants for reactions given pressure,
C  temperature(s) and mole fractions;  see Eq. (54).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, (mole/cm**3)**some power, depending on
C                            the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), EQKC(*)
C
      CALL CKEQ (RCKWRK, ICKWRK, T, ICKWRK(IcNU), ICKWRK(IcNK),
     1           RCKWRK(NcCO), ICKWRK(IcRV), RCKWRK(NcRV), ICKWRK(IcLT),
     2           RCKWRK(NcLT), ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     3           ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     4           ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     5           ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     5           ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     6           RCKWRK(NcKR), EQKC)
C
C     end of SUBROUTINE CKEQXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQXR (RHO, T, X, ICKWRK, RCKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQXR (RHO, T, X, ICKWRK, RCKWRK, EQKC)
C  Returns the equilibrium constants of the reactions given mass
C  density, temperature(s) and mole fractions;  see Eq. (54).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C                   species.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants in concentration units
C              for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, (mole/cm**3)**some power, depending on
C                             the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), EQKC(*)
C
      CALL CKEQ (RCKWRK, ICKWRK, T, ICKWRK(IcNU), ICKWRK(IcNK),
     1           RCKWRK(NcCO), ICKWRK(IcRV), RCKWRK(NcRV), ICKWRK(IcLT),
     2           RCKWRK(NcLT), ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     3           ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     4           ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     5           ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     5           ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     6           RCKWRK(NcKR), EQKC)
C
C     end of SUBROUTINE CKEQXR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQYP (P, T, Y, ICKWRK, RCKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQYP (P, T, Y, ICKWRK, RCKWRK, EQKC)
C  Returns the equilibrium constants for reactions given pressure
C  temperature(s) and mass fractions;  see Eq. (54).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants in concentration units
C              for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, (mole/cm**3)**some power, depending on
C                              the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), EQKC(*)
C
      CALL CKEQ (RCKWRK, ICKWRK, T, ICKWRK(IcNU), ICKWRK(IcNK),
     1           RCKWRK(NcCO), ICKWRK(IcRV), RCKWRK(NcRV), ICKWRK(IcLT),
     2           RCKWRK(NcLT), ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     3           ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     4           ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     5           ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     5           ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     6           RCKWRK(NcKR), EQKC)
C
C     end of SUBROUTINE CKEQYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEQYR (RHO, T, Y, ICKWRK, RCKWRK, EQKC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEQYR (RHO, T, Y, ICKWRK, RCKWRK, EQKC)
C  Returns the equilibrium constants of the reactions given mass
C  density, temperature(s) and mass fractions;  see Eq. (54).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units; gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  EQKC(*)   - Real array, equilibrium constants in concentration units
C              for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units; (mole/cm**3)**some power, depending on
C                             the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), EQKC(*)
C
      CALL CKEQ (RCKWRK, ICKWRK, T, ICKWRK(IcNU), ICKWRK(IcNK),
     1           RCKWRK(NcCO), ICKWRK(IcRV), RCKWRK(NcRV), ICKWRK(IcLT),
     2           RCKWRK(NcLT), ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     3           ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     4           ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     5           ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     5           ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     6           RCKWRK(NcKR), EQKC)
C
C     end of SUBROUTINE CKEQYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKESIG (P, T, X, XNUEH, ICKWRK, RCKWRK, SIGE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKESIG (P, T, X, XNUEH, ICKWRK, RCKWRK, SIGE)
C  Returns the electron species electrical conductivity given
C  collision frequency.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUEH     - Real scalar, total of the momentum-transfer collision
C              frequencies of the electronix.
C
C  OUTPUT
C  SIGE      - Real scalar, electron electrical conductivity (DC)
C                   cgs units, ERG/CM*V*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23, ECHRG=1.6022D-12, SMALLX=1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23, ECHRG=1.6022E-12, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      EDEN = TDEN * MAX(X(KEL),SMALLX)
      EMASS = RCKWRK(NcWT+KEL-1) / AVAG
      SIGE = EDEN * (ECHRG**2) / (EMASS*XNUEH)
C
C     end of SUBROUTINE CKESIG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKEVIS (P, T, X, XNUEH, ICKWRK, RCKWRK, VISE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKEVIS (P, T, X, XNUEH, ICKWRK, RCKWRK, VISE)
C  Returns the electron species viscosity given collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUEH     - Real scalar, total of the momentum-transfer collision
C              frequencies of the elctronis
C  OUTPUT
C  VISE      - Real scalar, electron viscosity
C                   cgs units, GM/CM*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      EDEN = TDEN * MAX(X(KEL),SMALL)
      BOLTZ = RCKWRK(NcRU)/AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      VISE = (4.*BOLTZ*TE*EDEN)/(PI*MAX(XNUEH,SMALL))
C
C     end of SUBROUTINE CKEVIS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKFAL  (NDIM, ICKWRK, RCKWRK, IFOP, IFLO, KFAL, FPAR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKFAL  (NDIM, ICKWRK, RCKWRK, IFOP, IFLO, KFAL, FPAR)
C  Returns a set of flags indicating whether a reaction has pressure-
C  dependent behavior and an array of parameters.
C
C  INPUT
C  NDIM      - Integer scalar, first dimension of the matrix FPAR;
C              NDIM must be greater than or equal to NFAR, the
C              maximum number of supplemental rate  parameters, which
C              is currently equal to 8.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  IFOP(*)   - Integer array,  flags indicating pressure-dependent
C              behavior;
C              dimension at least II, the total reaction count.
C              IFOP(I) indicates the pressure-dependent behavior
C              of reaction I:
C               0 - No pressure dependency
C               1 - Lindeman form (3 parameters)
C               2 - SRI form      (8 parameters)
C               3 - Troe form     (6 parameters)
C               4 - Troe form     (7 parameters)
C  IFLO(*)   - Integer array, flags indication pressure-depencency;
C              dimension at least II, the total reaction count.
C              IFLO(I) indicates
C              0 - unimolecular fall-off,
C              1 - chemically activated bi-molecular.
C  KFAL(*)   - Integer array, flags indicating type of bath-gas
C              concentration to be used in expressions
C              (see footnote on page 27);
C              dimension at least II, the total reaction count.
C              KFAL(I) indicates the type of reaction I:
C               0 - Use total concentration of gas mixture
C                    (with the added capability of using enhanced
C                     third body coefficients) (default)
C               K - Use the concentration of species K
C  FPAR(*,*) - Real matrix, pressure dependency parameters;
C              dimension at least NFAR for the first, the maximum
C              number of parameters (currently 8), and
C              at least II for the second, the total reaction
C              count.
C              The number of parameters depends on the
C              particular functional form indicated by the IFOP array:
C              FPAR(1,I), FPAR(2,I), FPAR(3,I) are always the
C              parameters entered on the LOW auxiliary keyword line
C              in the CHEMKIN interpretor input file.
C                 FPAR(1,I) = Pre-exponential for low pressure
C                             limiting rate constant
C                             cgs units, mole-cm-sec-K
C                 FPAR(2,I) = Temperature dependence exponents
C                             for the low pressure limiting rate
C                             constants.
C                 FPAR(3,I) = Activation energy for the low
C                             pressure limiting rate constant.
C                             cgs units, K
C              Additional FPAR values depend on IFOP:
C              IFOP(I) = 2:
C                 FPAR(4,I) = a           (See Eqn. (69))
C                 FPAR(5,I) = b (Kelvin)  (See Eqn. (69))
C                 FPAR(6,I) = c (Kelvin)  (See Eqn. (69))
C                 FPAR(7,I) = d           (See Eqn. (69))
C                 FPAR(8,I) = e           (See Eqn. (69))
C              IFOP(I) = 3:
C                 FPAR(4,I) = a             (See Eqn. (68))
C                 FPAR(5,I) = T*** (Kelvin) (See Eqn. (68))
C                 FPAR(6,I) = T*   (Kelvin) (See Eqn. (68))
C              IFOP(I) = 4:
C                 FPAR(4,I) = a             (See Eqn. (68))
C                 FPAR(5,I) = T*** (Kelvin) (See Eqn. (68))
C                 FPAR(6,I) = T*   (Kelvin) (See Eqn. (68))
C                 FPAR(7,I) = T**  (Kelvin) (See Eqn. (68))
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), IFOP(NII), IFLO(NII),
     1          KFAL(NII), FPAR(NDIM,NII)
C
      DO 100 I = 1, NII
        IFOP(I) = 0
        IFLO(I) = 0
        KFAL(I) = 0
        DO 50 N = 1, NFAR
           FPAR(N,I) = 0.0
   50   CONTINUE
  100 CONTINUE
C
      I_FPAR = NcFL - NFAR
      DO 250 N = 0, NFAL - 1
        I       = ICKWRK(IcFL + N)
        IFOP(I) = ICKWRK(IcFO + N)
        IFLO(I) = ICKWRK(IcFT + N)
        KFAL(I) = ICKWRK(IcKF + N)
        IF (IFOP(I) .EQ. 1) THEN
           NF = 3
        ELSEIF (IFOP(I) .EQ. 2) THEN
           NF = 8
        ELSEIF (IFOP(I) .EQ. 3) THEN
           NF = 6
        ELSEIF (IFOP(I) .EQ. 4) THEN
           NF = 7
        ENDIF
        I_FPAR = I_FPAR + NFAR
        DO 200 L = 1, NF
          FPAR(L,I) = RCKWRK(I_FPAR + L - 1)
  200   CONTINUE
  250 CONTINUE
C
C     end of SUBROUTINE CKFAL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKFALP  (P, T, X, ICKWRK, RCKWRK, I,
     $                    RKLOW, CTB, PR, FC, PCOR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKFALP  (P, T, X, ICKWRK, RCKWRK, I,
C                      RKLOW, CTB, PR, FC, PCOR)
C
C     This subroutine returns details concerning the reaction rate
C     constant for fall-off reactions.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T(*)   - Temperature array.
C                   cgs units - K
C                   Data type - real vector
C     X      - Mole fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension X(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C
C     RKLOW    - Low Pressure forward reaction rate for fall-off
C                reactions. It is defined to be zero for
C                non-fall-off reactions.
C                   cgs units - 1/(sec) *
C                           (cm**3/mole)**(sum of forward stoich. coeff)
C                   Data type - real
C     CTB      - Effective concentration for reaction, I_SAVE.
C                This takes into account the effectiveness factors
C                for the reaction, applicable to third body
C                and fall-off reactions.  It is defined to be equal
C                to the total concentration for other fall-off
C                or third body reactions, and to be equal to one
C                for reactions which don't use it
C                Units are moles/cm**3.
C                   cgs units - mole/(cm**3)
C                   Data type - real
C
C     PR       - Reduced Pressure for fall-off reactions.  This is
C                defined to be equal to CTB*RKLOW_SAVE/RCF_INF.
C                where RCF_INF is the high pressure rate constant.
C                This is a dimensionless quantity.  For non-fall-off
C                reactions, this quantity is defined to be 0.
C                   cgs units - unitless
C                   Data type - real
C     FC       - Correction to L-H rate constant for fall-off
C                reactions.  It is defined to be 0 for non-fall-off
C                reactions.
C                   cgs units - unitless
C                   Data type - real
C     PCOR     - This is equal to the pressure correction ratio for
C                fall-off reactions, i.e., RC(T,P) / RC(T)_inf
C                the ratio of the actual reaction rate to the high
C                pressure reaction rate constant.
C                   cgs units - unitless
C                   Data type - real array
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), X(*), T(*)
C
*
*       Call the normal temperature dependent part of the rate constant
*       using the internal work space,  RCKWRK(NcKF), RCKWRK(NcKR),
*       to store the complete forward and reverse rate constant
*       vectors.
*
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF),  RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
C
C     Calculate the concentrations given the mole fractions
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
C       Extract various fall-off parameters
      CALL CKFAL1(RCKWRK, ICKWRK, NII, NKK, MXTB, RCKWRK(NcRU),
     1            RCKWRK(NcPA), T, RCKWRK(NcK1), NPAR+1, RCKWRK(NcCO),
     2            NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF),
     3            NFAR, RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),
     4            RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),
     5            RCKWRK(NcKR), I, RKLOW, CTB, PR, FC, PCOR)
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKFAL1 (RCKWRK, ICKWRK, II, KK, MAXTB, RU, PATM, T, C,
     1                   NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR,
     5                   NTHB, ITHB, NTBS, AIK, NKTB, RCFT, RCRT,
     8                   I, RKLOW, CTB, PR, FC, PCOR)
C
C  START PROLOGUE
C
C   This subroutine modifies the forward and reverse rate constants
C obtained from CKRATT to account for those parts of CKRATX that
C don't involve multiplications of the concentrations of reactants
C or products. This specifically includes third body effects and
C and fall-off effects.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION  (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL  (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MACH/ SMALL, BIG, EXPARG
      DIMENSION ICKWRK(*), RCKWRK(*), IFAL(*), IFOP(*), KFAL(*),
     1          ITHB(*), NTBS(*), NKTB(MAXTB,*), C(*), PAR(NPAR,*),
     2          AIK(MAXTB,*), RCFT(*), RCRT(*), FPAR(NFAR,*), T(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
C         PRELIMINARIES
C
      TEMP = T(1)
      ALOGT = LOG(TEMP)
      AINVT = 1.0/TEMP
      PFAC  = PATM/(RU*TEMP)
      RKLOW = 0.0
      CTB   = 0.0
      PR    = 0.0
      FC    = 1.0
      PCOR  = 1.0
C
C  Find the total concentration
C
      CTOT = 0.0
      DO 10 K = 1, KK
         CTOT = CTOT + C(K)
   10 CONTINUE
*
*       Third-body reactions - find the correct one for the current
*       reaction
*
      ISTHB = CKLKUP (I, ITHB, NTHB)
      IF (ISTHB .GT. 0) THEN
         CTB = CTOT
         DO 75 N = 1, NTBS(ISTHB)
            CTB = CTB + (AIK(N, ISTHB)-1.0) * C(NKTB(N, ISTHB))
   75    CONTINUE
      ENDIF
C
C  Corrections for fall-off reactions
C         Only process the correct reaction
C
      ISFAL = CKLKUP (I, IFAL, NFAL)
      IF (ISFAL .GT. 0) THEN
         FP1 = FPAR(1,ISFAL)
         FP2 = FPAR(2,ISFAL)
         FP3 = FPAR(3,ISFAL)
         RKLOW = FP1 * EXP(FP2 * ALOGT - FP3 * AINVT)
C
C         Store the special species, if there is one
         K = KFAL(ISFAL)
         IF (K .EQ. 0) THEN
            PR = RKLOW * CTB  / RCFT(I)
         ELSE
            PR = RKLOW * C(K) / RCFT(I)
         ENDIF
C
C        This is the Lindemann form , i.e., IFOP(N) = 1
         PCOR = PR / (1.0 + PR)
C
         IF (IFOP(ISFAL) .GT. 1) THEN
            PRLOG = LOG10(MAX(PR,SMALL))
C
            FP4 = FPAR(4,ISFAL)
            FP5 = FPAR(5,ISFAL)
            FP6 = FPAR(6,ISFAL)
            FP7 = FPAR(7,ISFAL)
            FP8 = FPAR(8,ISFAL)
            IF (IFOP(ISFAL) .EQ. 2) THEN
C              SRI form
               XP = 1.0/(1.0 + PRLOG**2)
               FC = ( FP4 * EXP(-FP5*AINVT) + EXP(-TEMP/FP6) )**XP
     1              * FP7 * TEMP**FP8
C
            ELSE
C              6-parameter TROE form
               FCENT = (1.0-FP4) * EXP(-TEMP/FP5) + FP4 * EXP(-TEMP/FP6)
C
C              7-parameter TROE form
               IF (IFOP(ISFAL) .EQ. 4) FCENT = FCENT + EXP(-FP7*AINVT)
C
               FCLOG = LOG10(MAX(FCENT,SMALL))
               XN    = 0.75 - 1.27*FCLOG
               CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
               FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
               FC = 10.0**FLOG
            ENDIF
            PCOR = FC * PR/(1.0+PR)
         ENDIF
C
C           Correct both the forward and reverse rate constant
C
         RCFT(I) = RCFT(I) * PCOR
         RCRT(I) = RCRT(I) * PCOR
      ENDIF
C
C     Multiply the rate constant by the third body factor and
C     PAR(4,I), perturbation factor.
C
      RCFT(I) = RCFT(I)*CTB*PAR(4,I)
      RCRT(I) = RCRT(I)*CTB*PAR(4,I)
C
      RETURN
      END
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKGBML (P, T, X, ICKWRK, RCKWRK, GBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKGBML (P, T, X, ICKWRK, RCKWRK, GBML)*
C  Returns the mean Gibbs free energy of the mixture in molar units
C  given pressure, temperature(s) and mole fractions; see Eq. (44).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  GBML      - Real scalar, mean Gibbs free energy.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
C
      RU = RCKWRK(NcRU)
      RLNP = RU * LOG(P / RCKWRK(NcPA))
      GBML = 0.0
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         TK = T(ICKWRK(IcKTF + K))
         GBML = GBML + X(K+1) * ( RCKWRK(NcK2 + K) - TK *
     1          (RCKWRK(NcK1 + K) - RU *
     2           LOG(MAX(X(K+1),SMALL)) - RLNP))
  100 CONTINUE
C
C     end of SUBROUTINE CKGBML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKGBMS (P, T, Y, ICKWRK, RCKWRK, GBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKGBMS (P, T, Y, ICKWRK, RCKWRK, GBMS)*
C  Returns the mean Gibbs free energy of the mixture in mass units
C  given pressure, temperature(s), and mass fractions; see Eq. (45).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  GBMS      - Real scalar, mean Gibbs free energy.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK3))
      CALL CKMMWY(Y, ICKWRK, RCKWRK, WTM)
      RU = RCKWRK(NcRU)
      RLNP = RU * LOG(P / RCKWRK(NcPA))
C
      SUM = 0.0
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         SUM = SUM + RCKWRK(NcK3 + K) *
     1             ( RCKWRK(NcK2 + K) - T(ICKWRK(IcKTF + K)) *
     2             ( RCKWRK(NcK1 + K) - RU *
     3               LOG(MAX(RCKWRK(NcK3 + K),SMALL)) - RLNP))
  100 CONTINUE
      GBMS = SUM / WTM
C
C     end of SUBROUTINE CKGBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKGML  (T, ICKWRK, RCKWRK, GML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKGML  (T, ICKWRK, RCKWRK, GML)
C  Returns the standard state Gibbs free energies in molar units;
C  see Eq. (24).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  GML(*)    - Real array, standard state gibbs free energies for
C              the species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), GML(*)
C
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TK = T(ICKWRK(IcKTF + K))
         GML(K+1) = RCKWRK(NcK1 + K) - TK*RCKWRK(NcK2 + K)
150   CONTINUE
C
C     end of SUBROUTINE CKGML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKGMS  (T, ICKWRK, RCKWRK, GMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKGMS  (T, ICKWRK, RCKWRK, GMS)
C  Returns the standard state Gibbs free energies in mass units;
C  see Eq. (31).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  GMS(*)    - Real array, standard state Gibbs free energies for
C              the species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), GMS(*)
C
      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKSMS (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         TK = T(ICKWRK(IcKTF + K))
         GMS(K+1) = RCKWRK(NcK1 + K) - TK*RCKWRK(NcK2 + K)
150   CONTINUE
C
C     end of SUBROUTINE CKGMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHBML (T, X, ICKWRK, RCKWRK, HBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHBML (T, X, ICKWRK, RCKWRK, HBML)
C  Returns the mean enthalpy of the mixture in molar units;
C  see Eq. (37).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HBML      - Real scalar, mean enthalpy.
C                   cgs units - ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKHML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      HBML = 0.0
      DO 100 K = 1, NKK
         HBML = HBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKHBM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHBMS (T, Y, ICKWRK, RCKWRK, HBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHBMS (T, Y, ICKWRK, RCKWRK, HBMS)
C  Returns the mean enthalpy of the mixture in mass units; see Eq. (38).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HBMS      - Real scalar, mean enthalpy.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKHMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      HBMS = 0.0
      DO 100 K = 1, NKK
         HBMS = HBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKHBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHML  (T, ICKWRK, RCKWRK, HML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHML  (T, ICKWRK, RCKWRK, HML)
C  Returns the enthalpies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HML(*)    - Real array, enthalpies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), HML(*)
      SAVE TN1, TN2, TN3, TN4, TN5
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN5 = TN1*TN4
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
         TN5 = TN5 / 5
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK5 = TK1*TK4
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
            TK5 = TK5 / 5
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
            TK5 = TN5
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HML(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1)*TK1
     1            + RCKWRK(NA1+1) * TK2 + RCKWRK(NA1+2) * TK3
     2            + RCKWRK(NA1+3) * TK4 + RCKWRK(NA1+4) * TK5
     3            + RCKWRK(NA1 + NCP1 - 1))
250   CONTINUE
C
C     end of SUBROUTINE CKHML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHMS  (T, ICKWRK, RCKWRK, HMS)
C  Returns the enthalpies in mass units;  see Eq. (27).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HMS(*)    - Real array, enthalpies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), HMS(*)
      SAVE TN1, TN2, TN3, TN4, TN5
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN5 = TN1*TN4
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
         TN5 = TN5 / 5
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK5 = TK1*TK4
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
            TK5 = TK5 / 5
         ELSE
            TK1 = TN1
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
            TK5 = TN5
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HMS(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1) * TK1
     1            + RCKWRK(NA1+1) * TK2 + RCKWRK(NA1+2) * TK3
     2            + RCKWRK(NA1+3) * TK4 + RCKWRK(NA1+4) * TK5
     3            + RCKWRK(NA1 + NCP1 - 1)) / RCKWRK(NcWT + K)
250   CONTINUE
C
C     end of SUBROUTINE CKHMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHORT (T, ICKWRK, RCKWRK, HORT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHORT (T, ICKWRK, RCKWRK, HORT)
C  Returns the nondimensional enthalpies;  see Eq. (20).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HORT(*)   - Real array, nondimensional enthalpies for species;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), HORT(*)
      SAVE TN1, TNHALF, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNHALF = TN1 / 2
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 3
         TN3 = TN3 / 4
         TN4 = TN4 / 5
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKHALF = TK1 / 2
            TK2 = TK1 * TK1
            TK3 = TK1 * TK2
            TK4 = TK1 * TK3
            TK2 = TK2 / 3
            TK3 = TK3 / 4
            TK4 = TK4 / 5
         ELSE
            TK1 = TN1
            TKHALF = TNHALF
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         HORT(K+1) = RCKWRK(NA1)
     1             + RCKWRK(NA1+1) * TKHALF + RCKWRK(NA1+2) * TK2
     2             + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3             + RCKWRK(NA1 + NCP1 - 1) / TK1
250   CONTINUE
C
C     end of SUBROUTINE CKHORT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKHRX  (I, HML, ICKWRK, RCKWRK, HRXI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKHRX  (I, HML, ICKWRK, RCKWRK, HRXI)
C  Returns the molar heat of reaction I.
C
C  INPUT
C  I          - Integer scalar, reaction index.
C  HML(*)     - Real array, molar enthalpies for species;
C               dimension at lest KK, the total species count.
C                  cgs units, ergs/mole
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  HRXI      - Real scalar, molar heat of reaction I.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      DIMENSION ICKWRK(*), RCKWRK(*), HML(*)
      INTEGER CKLKUP
      EXTERNAL CKLKUP
C
      HRXI = 0.0
      ISRNU = CKLKUP (I, ICKWRK(IcRNU), NRNU)
      I_NK = IcNK + (I-1)*MXSP - 1
C
      IF (ISRNU .LE. 0) THEN
         I_NU = IcNU + (I-1)*MXSP - 1
         DO 30 N = 1, MXSP
            K = ICKWRK(I_NK + N)
            NUKI = ICKWRK(I_NU + N)
            IF (K.NE.0 .AND. NUKI.NE.0) HRXI = HRXI + NUKI*HML(K)
30       CONTINUE
      ELSE
         I_NU = NcRNU + (ISRNU-1)*MXSP - 1
         DO 35 N = 1, MXSP
            K = ICKWRK(I_NK + N)
            RNUKI = RCKWRK(I_NU + N)
            IF (K.NE.0 .AND. RNUKI.GT.CKMIN) HRXI=HRXI + RNUKI*HML(K)
   35    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKHRX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKI2CH (NUM, STR, I, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKI2CH (NUM, STR, I, KERR)
C  Returns a character string representation of an integer and the
C  character count of the string.
C
C  INPUT
C  NUM   - Integer scalar, to be converted to a character string;
C          the maximum magnitude of NUM is machine-dependent.
C
C  OUTPUT
C  STR   - Character string, left-justified character representation
C          of NUM.
C  I     - Integer scalar, the non-blank character count of STR.
C  KERR  - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STR*(*), IST(10)*(1)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
      DATA IST/'0','1','2','3','4','5','6','7','8','9'/
      BIGI = 2147483647.
C
      I = 0
      STR = ' '
      ILEN = LEN(STR)
      KERR = .FALSE.
      IF (ILEN.LT.1 .OR. IABS(NUM).GT.BIGI) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IF (NUM .EQ. 0) THEN
         STR = '0'
         I = 1
         RETURN
      ELSEIF (NUM .LT. 0) THEN
         STR(1:) = '-'
      ENDIF
C
      INUM = IABS(NUM)
      NCOL = NINT(LOG10(REAL(INUM))) + 1
C
      DO 10 J = NCOL, 1, -1
         IDIV = INUM / 10.0**(J-1)
         IF (J.EQ.NCOL .AND. IDIV.EQ.0) GO TO 10
         LT = CKLSCH(STR)
         IF (LT .EQ. ILEN) THEN
            STR = ' '
            KERR = .TRUE.
            RETURN
         ENDIF
         STR(LT+1:) = IST(IDIV+1)
         INUM = INUM - IDIV*10.0**(J-1)
   10 CONTINUE
      I = CKLSCH(STR)
C
C     end of SUBROUTINE CKI2CH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKICON (P, T, X, XNUIM, K, ICKWRK, RCKWRK, CONI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKICON (P, T, X, XNUIM, K, ICKWRK, RCKWRK, CONI)
C  Returns the ion species conductivities given collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUIM     - Real scalar, total momentum-transfer collision
C              frequency for an ion
C  K         - Integer scalar, species index of an ion
C
C  OUTPUT
C  CONI      - Real scalar, ion thermal conductivity
C                 cgs units, ERG/CM*K*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      CALL CKCPMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKCVMS (T, ICKWRK, RCKWRK, RCKWRK(NcK2))
      TDEN = (RHO/WTM)*AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      TI = T(ICKWRK(IcKTF + K - 1))
      XDEN = TDEN * MAX(X(K),SMALL)
      GAMI = RCKWRK(NcK1 + K - 1)/RCKWRK(NcK2 + K - 1)
      XMSS = RCKWRK(NcWT + K - 1)/AVAG
      CONI = ((9.*GAMI-5.)/(GAMI-1.)) *
     1             (BOLTZ**2) * XDEN * TI /(XMSS*MAX(XNUIM,SMALL))
C
C     end of SUBROUTINE CKICON
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIDIF (T, XNUIK, KIK, KDIM, ICKWRK, RCKWRK, DIK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIDIF (T, XNUIK, KIK, KDIM, ICKWRK, RCKWRK, DIK)
C  Returns the ion species binary diffusion coefficients given
C  collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  XNUIK(*,*)- Real matrix, momentum-transfer collision frequency
C              for the ion with each species;
C              dimension at least KK for the first, the total species
C              count, and at least NKKI for the second, the total
C              ion count.
C                 cgs units - /S
C  KIK       - Integer scalar, species index of the ion.
C  KDIM      - Integer scalar, first dimension of the matrix XNUIK;
C              KDIM must be greater than or equal to KK, the total
C              species count.
C  OUTPUT
C  DIK(*)    - Real array, ion viscosity;
C              dimension at least KK, the total species count.
C                 cgs units - GM/CM*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), XNUIK(KDIM,*), ICKWRK(*), RCKWRK(*), DIK(*)
C
C Find the Ion index
C
      DO 50 K1 = 1, NKKI
         IF (KIK .EQ. ICKWRK(IcKI+K1-1)) KI = K1
50    CONTINUE
      BOLTZ = RCKWRK(NcRU)/AVAG
      TI = T(ICKWRK(IcKTF + KIK - 1))
      XMSS = RCKWRK(NcWT + KIK - 1)/AVAG
C
      DO 100 K = 1, NKK
         IF (XMSS*XNUIK(K,KI) .GT. SMALL) THEN
            DIK(K)  = (BOLTZ*TI)/(XMSS*XNUIK(K,KI))
         ELSE
            DIK(K) = BIG
         ENDIF
 100  CONTINUE
C
C     end of SUBROUTINE CKIDIF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIEIM  (ICKWRK, RCKWRK, IEIM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIEIM (ICKWRK, RCKWRK, IEIM)
C  Returns a set of flags indicating whether the reactions are
C  electron-impact, and if so, the temperature dependence
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  IEIM(*)   - Integer array, electron-impact flags for reactions;
C              dimension at least II, the total reaction count.
C              IEIM(I)= -1  reaction I is not a third-body reactions
C              IEIM(I)=  N  reaction I is a third-body reaction with
C                        temperature dependence N
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), IEIM(*)
C
      DO 100 I = 1, NII
         IEIM(I) = -1
  100 CONTINUE
      DO 150 N = 0, NEIM - 1
         IEIM(ICKWRK(IcEI + N)) = ICKWRK(IcET + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKIEIM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIEXC  (ICKWRK, RCKWRK, IEXC, EEXC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIEXC (ICKWRK, RCKWRK, IEXC, EEXC)
C  Returns a set of flags indicating whether the reactions are
C  excitation reactions and, if so, the energy loss per event in eV.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  IEXC(*)   - Integer array, excitation-only reaction flag;
C              dimension at least II, the total reaction count.
C              IEXC(I)= -1  reaction I is not an excitation-only reax
C              IEXC(I)=  1  reaction I is an excitation reaction
C  EEXC(*)   - Real array, excitation energy loss per event in forward
C              direction for reactions;
C              dimension at least II, the total reaction count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), IEXC(*), EEXC(*)
C
      DO 100 I = 1, NII
         IEXC(I) = -1
         EEXC(I) = 0.0
  100 CONTINUE
      DO 150 N = 0, NEXC - 1
         IEXC(ICKWRK(IcEX + N)) = 1
         EEXC(ICKWRK(IcEX + N)) = RCKWRK(NcEX + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKIEXC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIMOM  (ICKWRK, IMOM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIMOM (ICKWRK, IMOM)
C  Returns a set of flags indicating whether the reactions are
C  electron momentum-transfer collision frequencies and, if so,
C  the index of the species with which the electron collides.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  IMOM(*)   - Integer array, electron momentum-transfer collision
C              frequency flags for reactions;
C              dimension at least II, the total reaction count.
C              IMOM(I)= -1  reaction I is not a mom-transfer coll freq
C              IMOM(I)=  K  reaction I is a mom-transfer coll frequency
C                        and K is species index of the electron's
C                        collision partner
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), IMOM(*)
C
      DO 100 I = 1, NII
         IMOM(I) = -1
  100 CONTINUE
      DO 150 N = 0, NMOM - 1
         IMOM(ICKWRK(IcMO + N)) = ICKWRK(IcMK + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKIMOM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)*
C  Returns a group of indices defining the size of the particular
C  reaction mechanism
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  MM        - Integer scalar, mechanism total element count.
C  KK        - Integer scalar, mechanism total species count.
C  II        - Integer scalar, mechanism total reaction count.
C  NFIT      - Integer scalar, number of coefficients in fits to
C              thermodynamic data for a temperature range;
C              NFIT=number of coefficients in polynomial fits to CP/R
C              plus 2.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      MM = NMM
      KK = NKK
      II = NII
      NFIT = NCP2
C
C     end of SUBROUTINE CKINDX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKINIT (LENICK, LENRCK, LENCCK, LINC, LOUT, ICKWRK,
     1                   RCKWRK, CCKWRK, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE CKINIT (LENICK, LENRCK, LENCCK, LINC, LOUT, ICKWRK,
C                     RCKWRK, CCKWRK, IFLAG)**
C  Reads the linkfile and creates the internal work arrays ICKWRK,
C  RCKWRK and CCKWRK.  CKINIT must be called before any other CHEMKIN
C  subroutine can be used, as the work arrays must be available as
C  their input.
C
C  INPUT
C  LENICK - Integer scalar, length of the integer work array, ICKWRK.
C  LENRCK - Integer scalar, length of the real work array, RCKWRK.
C  LENCCK - Integer scalar, length of the character work array, CCKWRK.
C  LINC  -  Integer scalar, linkfile input file unit number.
C  LOUT  -  Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C  IFLAG     - Integer scalar to indicate successful reading of
C              linkfile; IFLAG>0 is an error type.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
C     Data about the machine dependent constants is carried in:
      COMMON /MACH/ SMALL,BIG,EXPARG
C
C     Set the gas constant to the 1986 CODATA recommended
C     value - Gas constant in cals/mole is determined by
C     division of the later value by 4.184.
C     The number for conversion to one atmosphere is exact.
C*****precision > double
      PARAMETER (RU=8.314510D7, RUC=RU/4.184D7, PA=1.01325D6)
C*****END precision > double
C*****precision > single
C      PARAMETER (RU=8.314510E7, RUC=RU/4.184E7, PA=1.01325E6)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER*(*) CCKWRK(*)
      CHARACTER*16 FILVER, PRVERS, PREC, IFMT, RFMT, CFMT, LFMT,
     1             CKVERS, CKPREC, CKDATE
      PARAMETER (IFMT='(10I12)', CFMT='(8A16)', RFMT='(1P,5E24.16)',
     1           LFMT='(L8)')
      LOGICAL IOK, ROK, COK, KERR, LBIN
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
      EXPARG = LOG(BIG)
C
      CKVERS = '5.15'
      CKDATE = '98/07/30'
C*****precision > double
      CKPREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      CKPREC = 'SINGLE'
C*****END precision > single
C
      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1   ' CKLIB: CHEMKIN-III GAS-PHASE CHEMICAL KINETICS LIBRARY,',
     2   CKPREC(1:CKLSCH(CKPREC)), ' PRECISION Vers. ',
     3   CKVERS(1:CKLSCH(CKVERS)+1), CKDATE,
     4   ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
      CALL CKLEN (LINC, LOUT, LI, LR, LC, IFLAG)
C
      IF (IFLAG .GT. 0) RETURN
C
      IOK = (LENICK .GE. LI)
      ROK = (LENRCK .GE. LR)
      COK = (LENCCK .GE. LC)
      IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
         IF (.NOT. IOK) WRITE (LOUT, '(10X,A,I5)')
     1      'ICKWRK MUST BE DIMENSIONED AT LEAST ', LI
         IF (.NOT. ROK) WRITE (LOUT, '(10X,A,I5)')
     1      'RCKWRK MUST BE DIMENSIONED AT LEAST ', LR
         IF (.NOT. COK) WRITE (LOUT, '(10X,A,I5)')
     1      'CCKWRK MUST BE DIMENSIONED AT LEAST ', LC
         IFLAG = 20
         RETURN
      ENDIF
      IF (LEN(CCKWRK(1)) .LT. 16) THEN
         WRITE (LOUT, '(10X,A)')
     1      'CHARACTER LENGTH OF CCKWRK MUST BE AT LEAST 16'
         IFLAG = 21
         RETURN
      ENDIF
C
      REWIND LINC
C*****linkfile (gas) > binary
C      LBIN = .TRUE.
C*****END linkfile (gas) > binary
C*****linkfile (gas) > ascii
      LBIN = .FALSE.
C*****END linkfile (gas) > ascii
C
      IF (LBIN) THEN
         NREC = 1
         READ (LINC, ERR=100) FILVER
         NREC = 2
         READ (LINC, ERR=100) PRVERS
         NREC = 3
         READ (LINC, ERR=100) PREC
         NREC = 4
         READ (LINC, ERR=100) KERR
         NREC = 5
         READ (LINC, ERR=100) LENI, LENR, LENC
         NREC = 6
         READ (LINC, ERR=100) MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,
     1                        NIFAR, NJA, MAXORD, NF1
         NREC = 7
         READ (LINC, ERR=100) MM, KK, II, NRV, NFL, NTB, NLT, NRL,NW,
     1                        NCHRG, NEI, NIJAN, NIF1, NEX, NMO, NXS,
     2                        NTD, NSTO, NOR, KELECT, KKI
         NREC = 8
         READ (LINC, ERR=100) CKMN
      ELSE
         NREC = 1
         READ (LINC, CFMT, ERR=100) FILVER
         NREC = 2
         READ (LINC, CFMT, ERR=100) PRVERS
         NREC = 3
         READ (LINC, CFMT, ERR=100) PREC
         NREC = 4
         READ (LINC, LFMT,  ERR=100) KERR
         NREC = 5
         READ (LINC, IFMT, ERR=100) LENI, LENR, LENC
         NREC = 6
         READ (LINC, IFMT, ERR=100) MAXSP, MAXTB, MAXTP, NTHCF, NIPAR,
     1                              NITAR, NIFAR, NJA, MAXORD, NF1
         NREC = 7
         READ (LINC, IFMT, ERR=100) MM, KK, II, NRV, NFL, NTB, NLT,
     1                              NRL, NW, NCHRG, NEI, NIJAN,
     1                              NIF1, NEX, NMO, NXS, NTD, NSTO,
     2                              NOR, KELECT, KKI
         NREC = 8
         READ (LINC, RFMT, ERR=100) CKMN
      ENDIF
C
      NMM = MM
      NKK = KK
      NII = II
      MXSP = MAXSP
      MXTB = MAXTB
      MXTP = MAXTP
      NCP  = NTHCF
      NCP1 = NTHCF+1
      NCP2 = NTHCF+2
      NCP2T = NCP2*(MAXTP-1)
      NPAR = NIPAR
      NLAR = NITAR
      NFAR = NIFAR
      NTHB = NTB
      NLAN = NLT
      NFAL = NFL
      NREV = NRV
      NRLT = NRL
      NWL  = NW
      NEIM = NEI
      NJAR = NJA
      NJAN = NIJAN
      NFT1 = NIF1
      NF1R = NF1
      NEXC = NEX
      NMOM = NMO
      NXSM = NXS
      NTDE = NTD
      NRNU= NSTO
      NORD = NOR
      MXORD= MAXORD
      KEL  = KELECT
      NKKI = KKI
      CKMIN= CKMN
C
C             APPORTION work arrays
C
C             SET  ICKWRK(*)=1  TO FLAG THAT CKINIT HAS BEEN CALLED
C
      ICKWRK(1) = 1
C
C             STARTING LOCATIONS OF INTEGER SPACE
C
C! elemental composition of species
      IcNC = 2
C! species phase array
      IcPH = IcNC + KK*MM
C! species charge array
      IcCH = IcPH + KK
C! ion species index array
      IcKI = IcCH + KK
C! # of temperatures for fit
      IcNT = IcKI + KKI
C! stoichiometric coefficients
      IcNU = IcNT + KK
C! species indices for coefficients
      IcNK = IcNU + MAXSP*II
C! # of non-zero coefficients  (<0=reversible, >0=irreversible)
      IcNS = IcNK + MAXSP*II
C! # of reactants
      IcNR = IcNS + II
C! Landau-Teller reaction indices
      IcLT = IcNR + II
C! Reverse Landau-Teller reactions
      IcRL = IcLT + NLAN
C! Fall-off reaction indices
      IcFL = IcRL + NRLT
C! Fall-off formulation option numbers
      IcFO = IcFL + NFAL
C! Fall-off option (uni-molecular vs. chemically-activated)
      IcFT = IcFO + NFAL
C! Fall-off enhanced species
      IcKF = IcFT + NFAL
C! Third-body reaction indices
      IcTB = IcKF + NFAL
C! count of 3rd bodies for above
      IcKN = IcTB + NTHB
C! array of species #'s for above
      IcKT = IcKN + NTHB
C! Reverse parameter reaction indices
      IcRV = IcKT + MAXTB*NTHB
C! Radiation wavelength reactions
      IcWL = IcRV + NREV
C! Electon-impact reaction indices
      IcEI = IcWL + NWL
C! Electron-impact temperature dependence flags
      IcET = IcEI + NEIM
C! Non thermal-equilibrium reaction indices
      IcTD = IcET + NEIM
C! Non thermal-equilibrium temperature species indices
      IcTK = IcTD + NTDE
C! Janev-Langer_Evans&Post type reaction indices
      IcJN = IcTK + NTDE
C! Reaction indices using fit#1
      IcF1 = IcJN + NJAN
C! Reaction indices for excitation-only reactions
      IcEX = IcF1 + NFT1
C! Reaction indices for momentum transfer collision frequencies
      IcMO = IcEX + NEXC
C! Heavy collision partner for electron momentum-transfer freq
      IcMK = IcMO + NMOM
C! Reaction indices for ion momentum-transfer cross sections
      IcXS = IcMK + NMOM
C! Ion species indices for ion momentum-transfer cross sections
      IcXI = IcXS + NXSM
C! Collision partner species index for ion mom-transf x-sections
      IcXK = IcXI + NXSM
C! Real stoichometry reactions
      IcRNU= IcXK + NXSM
C! Change of order reactions
      IcORD= IcRNU + NRNU
C! Species for which there is a change of order
      IcKOR= IcORD + NORD
C! Array indicating which temperature (in the temperature array)
C  corresponds to each species
      IcKTF= IcKOR + NORD*MXORD
C! Internal workspace of lengh kk
      IcK1 = IcKTF + KK
C!      ditto
      IcK2 = IcK1 + KK
C
      ITOT = IcK2 + KK - 1
C
C             STARTING LOCATIONS OF CHARACTER SPACE
C
C! start of element names
      IcMM = 1
C! start of species names
      IcKK = IcMM + MM
      ITOC = IcKK + KK - 1
C
C             STARTING LOCATIONS OF REAL SPACE
C
C! atomic weights
      NcAW = 1
C! molecular weights
      NcWT = NcAW + MM
C! temperature fit array for species
      NcTT = NcWT + KK
C! thermodynamic coefficients
      NcAA = NcTT + MAXTP*KK
C! Arrhenius coefficients (3)
      NcCO = NcAA + (MAXTP-1)*NCP2*KK
C! Reverse coefficients
      NcRV = NcCO + (NPAR+1)*II
C! Landau-Teller #'s for NLT reactions
      NcLT = NcRV + (NPAR+1)*NREV
C! Reverse Landau-Teller #'s
      NcRL = NcLT + NLAR*NLAN
C! Fall-off parameters for NFL reactions
      NcFL = NcRL + NLAR*NRLT
C! 3rd body coef'nts for NTHB reactions
      NcKT = NcFL + NFAR*NFAL
C! wavelength
      NcWL = NcKT + MAXTB*NTHB
C! Janev-type coefficients
      NcJN = NcWL + NWL
C! Fit#1 parameters
      NcF1 = NcJN + NJAR*NJAN
C! Excitation-only reaction energy loss
      NcEX = NcF1 + NF1R*NFT1
C! real stoichometric coefficients
      NcRNU= NcEX + NEXC
C! change of order for species/reactions
      NcKOR= NcRNU + NRNU*MXSP
C! universal gas constant
      NcRU = NcKOR + NORD*MXORD
C! universal gas constant in units
      NcRC = NcRU + 1
C! pressure of one atmosphere
      NcPA = NcRC + 1
C! intermediate temperature-dependent forward rates
      NcKF = NcPA + 1
C! intermediate temperature-dependent reverse rates
      NcKR = NcKF + II
C! internal work space of length kk
      NcK1 = NcKR + II
C!          'ditto'
      NcK2 = NcK1 + KK
C!          'ditto'
      NcK3 = NcK2 + KK
C!          'ditto'
      NcK4 = NcK3 + KK
      NcI1 = NcK4 + KK
      NcI2 = NcI1 + II
      NcI3 = NcI2 + II
      NcI4 = NcI3 + II
      NTOT = NcI4 + II - 1
C
C        SET UNIVERSAL CONSTANTS IN CGS UNITS
C
      RCKWRK(NcRU) = RU
      RCKWRK(NcRC) = RUC
      RCKWRK(NcPA) = PA
C
C     element and species records
      IF (LBIN) THEN
         NREC = 9
         READ (LINC, ERR=100) (CCKWRK(IcMM+M-1), M = 1, MM)
         NREC = 10
         READ (LINC, ERR=100) (RCKWRK(NcAW+M-1), M = 1, MM)
         NREC = 11
         READ (LINC, ERR=100) (CCKWRK(IcKK+K-1), K = 1, KK)
         NREC = 12
         READ (LINC, ERR=100) (RCKWRK(NcWT+K-1), K = 1, KK)
         NREC = 13
         READ (LINC, ERR=100)
     1   ((ICKWRK(IcNC + (K-1)*MM + M - 1),  M=1,MM), K=1,KK)
         NREC = 14
         READ (LINC, ERR=100) (ICKWRK(IcCH+K-1), K = 1, KK)
         NREC = 15
         READ (LINC, ERR=100) (ICKWRK(IcNT+K-1), K = 1, KK)
         NREC = 16
         READ (LINC, ERR=100) (ICKWRK(IcPH+K-1), K = 1, KK)
         NREC = 17
         READ (LINC, ERR=100)
     1   ((RCKWRK(NcTT + (K-1)*MAXTP + L - 1), L=1,MAXTP),K=1,KK)
         NREC = 18
         READ (LINC, ERR=100)
     1   ( ( (RCKWRK(NcAA + (L-1)*NCP2 + (K-1)*NCP2T+N-1),
     2        N=1,NCP2), L=1,(MAXTP-1)), K = 1, KK)
      ELSE
         NREC = 9
         READ (LINC, CFMT, ERR=100) (CCKWRK(IcMM+M-1), M=1,MM)
         NREC = 10
         READ (LINC, RFMT, ERR=100) (RCKWRK(NcAW+M-1), M=1,MM)
         NREC = 11
         READ (LINC, CFMT, ERR=100) (CCKWRK(IcKK+K-1), K=1,KK)
         NREC = 12
         READ (LINC, RFMT, ERR=100) (RCKWRK(NcWT+K-1), K=1,KK)
         NREC = 13
         READ (LINC, IFMT, ERR=100)
     1   ((ICKWRK(IcNC + (K-1)*MM + M - 1), M=1,MM), K=1,KK)
         NREC = 14
         READ (LINC, IFMT, ERR=100) (ICKWRK(IcCH+K-1), K=1,KK)
         NREC = 15
         READ (LINC, IFMT, ERR=100) (ICKWRK(IcNT+K-1), K=1,KK)
         NREC = 16
         READ (LINC, IFMT, ERR=100) (ICKWRK(IcPH+K-1), K=1,KK)
         NREC = 17
         READ (LINC, RFMT, ERR=100)
     1   ((RCKWRK(NcTT + (K-1)*MAXTP + L - 1), L=1,MAXTP),K=1,KK)
         NREC = 18
         READ (LINC, RFMT, ERR=100)
     1   (((RCKWRK(NcAA + (L-1)*NCP2 + (K-1)*NCP2T+N-1),
     2      N=1,NCP2), L=1,(MAXTP-1)), K=1,KK)
      ENDIF
C
      DO 14 K = 1, NKK
         ICKWRK(IcKTF + K - 1) = 1
   14 CONTINUE
C     ion records
      IF (KKI .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 19
            READ (LINC, ERR=100)
            NREC = 20
            READ (LINC, ERR=100) (ICKWRK(IcKI + K - 1), K = 1, KKI)
         ELSE
            NREC = 19
            READ (LINC, IFMT, ERR=100)
            NREC = 20
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcKI+K-1), K=1,KKI)
         ENDIF
      ENDIF
C
      IF (II .EQ. 0) THEN
         REWIND LINC
         RETURN
      ENDIF
C     reaction records
      IF (LBIN) THEN
         NREC = 21
         READ (LINC, ERR=100) (ICKWRK(IcNS + I - 1), I = 1, II)
         NREC = 22
         READ (LINC, ERR=100) (ICKWRK(IcNR + I - 1), I = 1, II)
         NREC = 23
         READ (LINC, ERR=100)
     1   ((ICKWRK(IcNU+(I-1)*MAXSP+N-1),
     2     ICKWRK(IcNK+(I-1)*MAXSP+N-1), N=1,MAXSP), I=1,II)
         NREC = 24
         READ (LINC, ERR=100)
     1     ((RCKWRK(NcCO+(I-1)*(NPAR+1)+N-1), N=1,NPAR), I=1,II)
      ELSE
         NREC = 21
         READ (LINC, IFMT, ERR=100) (ICKWRK(IcNS+I-1), I=1,II)
         NREC = 22
         READ (LINC, IFMT, ERR=100) (ICKWRK(IcNR+I-1), I=1,II)
         NREC = 23
         READ (LINC, IFMT, ERR=100)
     1   ((ICKWRK(IcNU+(I-1)*MAXSP+N-1),
     2     ICKWRK(IcNK+(I-1)*MAXSP+N-1), N=1,MAXSP), I=1,II)
         NREC = 24
         READ (LINC, RFMT, ERR=100)
     1     ((RCKWRK(NcCO+(I-1)*(NPAR+1)+N-1), N=1,NPAR), I=1,II)
      ENDIF
C
      DO 10 I = 1, II
         RCKWRK(NcCO + (I-1)*(NPAR+1) + NPAR) = 1.0
   10 CONTINUE
C
      IF (NREV .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 25
            READ (LINC, ERR=100)
            NREC = 26
            READ (LINC, ERR=100) (ICKWRK(IcRV+N-1), N = 1, NREV)
            NREC = 27
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcRV+(N-1)*(NPAR+1)+L-1),L=1,NPAR), N=1,NREV)
         ELSE
            NREC = 25
            READ (LINC, IFMT, ERR=100)
            NREC = 26
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcRV+N-1), N=1,NREV)
            NREC = 27
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcRV+(N-1)*(NPAR+1)+L-1),L=1,NPAR), N=1,NREV)
         ENDIF
      ENDIF
      IF (NFAL .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 28
            READ (LINC, ERR=100)
            NREC = 29
            READ (LINC, ERR=100) (ICKWRK(IcFL+N-1), ICKWRK(IcFO+N-1),
     1                ICKWRK(IcFT+N-1), ICKWRK(IcKF+N-1), N = 1, NFAL)
            NREC = 30
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcFL+(N-1)*NFAR + L-1), L=1,NFAR), N=1,NFAL)
         ELSE
            NREC = 28
            READ (LINC, IFMT, ERR=100)
            NREC = 28
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcFL+N-1),
     1      ICKWRK(IcFO+N-1), ICKWRK(IcFT+N-1), ICKWRK(IcKF+N-1),
     2      N = 1, NFAL)
            NREC = 30
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NCFL+(N-1)*NFAR + L-1), L=1,NFAR), N=1,NFAL)
         ENDIF
      ENDIF
      IF (NTHB .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 31
            READ (LINC, ERR=100)
            NREC = 32
            READ (LINC, ERR=100)
     1      (ICKWRK(IcTB+N-1), ICKWRK(IcKN+N-1), N=1,NTHB)
            NREC = 33
            READ (LINC, ERR=100)
     1      ((ICKWRK(IcKT+(N-1)*MAXTB+L-1), L=1,MAXTB), N=1,NTHB)
            NREC = 34
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcKT+(N-1)*MAXTB+L-1), L=1,MAXTB), N=1,NTHB)
         ELSE
            NREC = 31
            READ (LINC, IFMT, ERR=100)
            NREC = 32
            READ (LINC, IFMT, ERR=100)
     1      (ICKWRK(IcTB+N-1), ICKWRK(IcKN+N-1), N=1,NTHB)
            NREC = 33
            READ (LINC, IFMT, ERR=100)
     1      ((ICKWRK(IcKT+(N-1)*MAXTB+L-1), L=1,MAXTB), N=1,NTHB)
            NREC = 34
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcKT+(N-1)*MAXTB+L-1), L=1,MAXTB), N=1,NTHB)
         ENDIF
      ENDIF
      IF (NLAN .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 35
            READ (LINC, ERR=100)
            NREC = 36
            READ (LINC, ERR=100) (ICKWRK(IcLT+N-1), N = 1, NLAN)
            NREC = 37
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcLT+(N-1)*NLAR+L-1), L=1,NLAR), N=1,NLAN)
         ELSE
            NREC = 35
            READ (LINC, IFMT, ERR=100)
            NREC = 36
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcLT+N-1),N=1,NLAN)
            NREC = 37
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcLT+(N-1)*NLAR+L-1), L=1,NLAR), N=1,NLAN)
         ENDIF
      ENDIF
      IF (NRLT .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 38
            READ (LINC, ERR=100)
            NREC = 39
            READ (LINC, ERR=100) (ICKWRK(IcRL+N-1), N = 1, NRLT)
            NREC = 40
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcRL+(N-1)*NLAR+L-1), L=1,NLAR), N=1,NRLT)
         ELSE
            NREC = 38
            READ (LINC, IFMT, ERR=100)
            NREC = 39
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcRL+N-1), N=1,NRLT)
            NREC = 40
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcRL+(N-1)*NLAR+L-1), L=1,NLAR), N=1,NRLT)
         ENDIF
      ENDIF
      IF (NWL .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 41
            READ (LINC, ERR=100)
            NREC = 42
            READ (LINC, ERR=100) (ICKWRK(IcWL+N-1), N = 1, NWL)
            NREC = 43
            READ (LINC, ERR=100) (RCKWRK(NcWL+N-1), N = 1, NWL)
         ELSE
            NREC = 41
            READ (LINC, IFMT, ERR=100)
            NREC = 42
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcWL+N-1), N=1,NWL)
            NREC = 43
            READ (LINC, RFMT, ERR=100) (RCKWRK(NcWL+N-1), N=1,NWL)
         ENDIF
      ENDIF
      IF (NEIM .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 44
            READ (LINC, ERR=100)
            NREC = 45
            READ (LINC, ERR=100) (ICKWRK(IcEI+N-1), N = 1, NEIM)
            NREC = 46
            READ (LINC, ERR=100) (ICKWRK(IcET+N-1), N = 1, NEIM)
         ELSE
            NREC = 44
            READ (LINC, IFMT, ERR=100)
            NREC = 45
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcEI+N-1), N=1,NEIM)
            NREC = 46
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcET+N-1), N=1,NEIM)
         ENDIF
      ENDIF
      IF (NJAN .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 47
            READ (LINC, ERR=100)
            NREC = 48
            READ (LINC, ERR=100) (ICKWRK(IcJN+N-1), N = 1, NJAN)
            NREC = 49
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcJN+(N-1)*NJAR+L-1),L=1,NJAR), N=1,NJAN)
         ELSE
            NREC = 47
            READ (LINC, IFMT, ERR=100)
            NREC = 48
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcJN+N-1), N = 1, NJAN)
            NREC = 49
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcJN+(N-1)*NJAR+L-1),L=1,NJAR), N=1,NJAN)
         ENDIF
      ENDIF
      IF (NFT1 .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 50
            READ (LINC, ERR=100)
            NREC = 51
            READ (LINC, ERR=100) (ICKWRK(IcF1+N-1), N = 1, NFT1)
            NREC = 52
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcF1+(N-1)*NF1R+L-1),L=1,NF1R), N=1,NFT1)
         ELSE
            NREC = 50
            READ (LINC, IFMT, ERR=100)
            NREC = 51
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcF1+N-1), N=1,NFT1)
            NREC = 52
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcF1+(N-1)*NF1R+L-1),L=1,NF1R), N=1,NFT1)
         ENDIF
      ENDIF
      IF (NEXC .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 53
            READ (LINC, ERR=100)
            NREC = 54
            READ (LINC, ERR=100) (ICKWRK(IcEX+N-1), N = 1, NEXC)
            NREC = 55
            READ (LINC, ERR=100) (RCKWRK(NcEX+N-1), N = 1, NEXC)
         ELSE
            NREC = 53
            READ (LINC, IFMT, ERR=100)
            NREC = 54
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcEX+N-1), N = 1, NEXC)
            NREC = 55
            READ (LINC, RFMT, ERR=100) (RCKWRK(NcEX+N-1), N = 1, NEXC)
         ENDIF
      ENDIF
      IF (NMOM .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 56
            READ (LINC, ERR=100)
            NREC = 57
            READ (LINC, ERR=100) (ICKWRK(IcMO+N-1), N = 1, NMOM)
            NREC = 58
            READ (LINC, ERR=100) (ICKWRK(IcMK+N-1), N = 1, NMOM)
         ELSE
            NREC = 56
            READ (LINC, IFMT, ERR=100)
            NREC = 57
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcMO+N-1), N = 1, NMOM)
            NREC = 58
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcMK+N-1), N = 1, NMOM)
         ENDIF
      ENDIF
      IF (NXSM .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 59
            READ (LINC, ERR=100)
            NREC = 60
            READ (LINC, ERR=100) (ICKWRK(IcXS+N-1), N = 1, NXSM)
            NREC = 61
            READ (LINC, ERR=100) (ICKWRK(IcXI+N-1), N = 1, NXSM)
            NREC = 62
            READ (LINC, ERR=100) (ICKWRK(IcXK+N-1), N = 1, NXSM)
         ELSE
            NREC = 59
            READ (LINC, IFMT, ERR=100)
            NREC = 60
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcXS+N-1), N = 1, NXSM)
            NREC = 61
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcXI+N-1), N = 1, NXSM)
            NREC = 62
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcXK+N-1), N = 1, NXSM)
         ENDIF
      ENDIF
      IF (NTDE .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 63
            READ (LINC, ERR=100)
            NREC = 64
            READ (LINC, ERR=100) (ICKWRK(IcTD+N-1), N = 1, NTDE)
            NREC = 65
            READ (LINC, ERR=100) (ICKWRK(IcTK+N-1), N = 1, NTDE)
         ELSE
            NREC = 63
            READ (LINC, IFMT, ERR=100)
            NREC = 64
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcTD+N-1), N = 1, NTDE)
            NREC = 65
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcTK+N-1), N = 1, NTDE)
         ENDIF
      ENDIF
      IF (NRNU .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 66
            READ (LINC, ERR=100)
            NREC = 67
            READ (LINC, ERR=100) (ICKWRK(IcRNU+N-1), N = 1, NRNU)
            NREC = 68
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcRNU+(N-1)*MAXSP+L-1),L=1,MAXSP), N=1,NRNU)
         ELSE
            NREC = 66
            READ (LINC, IFMT, ERR=100)
            NREC = 67
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcRNU+N-1), N = 1, NRNU)
            NREC = 68
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcRNU+(N-1)*MAXSP+L-1),L=1,MAXSP), N=1,NRNU)

         ENDIF
      ENDIF
      IF (NORD .GT. 0) THEN
         IF (LBIN) THEN
            NREC = 69
            READ (LINC, ERR=100)
            NREC = 70
            READ (LINC, ERR=100) (ICKWRK(IcORD+N-1), N = 1, NORD)
            NREC = 71
            READ (LINC, ERR=100)
     1      ((ICKWRK(IcKOR+(N-1)*MXORD+L-1),L=1,MXORD), N=1,NORD)
            NREC = 72
            READ (LINC, ERR=100)
     1      ((RCKWRK(NcKOR+(N-1)*MXORD+L-1),L=1,MXORD), N=1,NORD)
         ELSE
            NREC = 69
            READ (LINC, IFMT, ERR=100)
            NREC = 70
            READ (LINC, IFMT, ERR=100) (ICKWRK(IcORD+N-1), N = 1, NORD)
            NREC = 71
            READ (LINC, IFMT, ERR=100)
     1      ((ICKWRK(IcKOR+(N-1)*MXORD+L-1),L=1,MXORD), N=1,NORD)
            NREC = 72
            READ (LINC, RFMT, ERR=100)
     1      ((RCKWRK(NcKOR+(N-1)*MXORD+L-1),L=1,MXORD), N=1,NORD)
         ENDIF
      ENDIF
C
      REWIND LINC
      RETURN
C
  100 CONTINUE
      IFLAG = NREC
      WRITE (LOUT, '(/A,/A,I5)')
     1   ' Error reading gas-phase linkfile,',
     2   ' SUBROUTINE CKINIT record index #', IFLAG
      REWIND LINC
C
C     Generic Formats - limit lines to 132 characters
C
8001  FORMAT (10I12)
8002  FORMAT (1P,5E24.16)
8003  FORMAT (8A16)
C
C     end of SUBROUTINE CKINIT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKINU   (I, NDIM, ICKWRK, RCKWRK, NSPEC, KI, NU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKINU   (I, NDIM, ICKWRK, RCKWRK, NSPEC, KI, NU)
C  Returns a count of species in a reaction, and their indices
C  and stoichiometric coefficients; see Eq. (50).
C
C  INPUT
C  I         - Integer scalar, index of a reaction;
C              I must be positive, and less than or equal to NII,
C              the total reaction count.
C  NDIM      - Integer scalar, dimension of the arrays KI and NU;
C              NDIM must be at least MAXSP, the maximum number of
C              species allowed in a reaction.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NSPEC     - Integer scalar, the total species count for reaction I.
C  KI(*)     - Integer array, species indices for those in
C              reaction I; dimension at least MAXSP, the maximum
C              number of species allowed in a reaction.
C              KI(N) is the index of the Nth species in reaction I.
C  NU(*)     - Integer array, stoichiometric coefficients for those
C              in reaction I;
C              dimension at least MAXSP, the maximum number of
C              species allowed in a reaction.
C              NU(N) is the stoichiometric coefficient of the Nth
C              Nth species in reaction I, and
C              NU < 0 if the Nth species is a reactant;
C              NU > 0 if the Nth species is a product.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), KI(*), NU(*)
C
      NSPEC = 0
      DO 50 N = 1, NDIM
         KI(N) = 0
         NU(N) = 0
   50 CONTINUE
C
      IF (NII.LE.0 .OR. NDIM.LT.MXSP .OR. I.LE.0 .OR.
     1    I.GT.NII) RETURN
C
      I_NK = IcNK + (I-1)*MXSP - 1
      I_NU = IcNU + (I-1)*MXSP - 1
      DO 200 N = 1, MXSP
         K = ICKWRK(I_NK + N)
         IF (K .NE. 0) THEN
            NSPEC = NSPEC + 1
            KI(NSPEC) = K
            NU(NSPEC) = ICKWRK(I_NU + N)
         ENDIF
200   CONTINUE
C
C     end of SUBROUTINE CKINU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKION (ICKWRK, KION)
C
C  START PROLOGUE
C
C  SUBROUTINE CKION (ICKWRK, KION)
C  Returns the ion species indices
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  KION(*)   - Integer array, ion species indices;
C              dimension at least NKKI, the total ion count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), KION(*)
C
      DO 100 KI = 1, NKKI
         KION(KI) = ICKWRK(IcKI+KI-1)
 100  CONTINUE
C
C     end of SUBROUTINE CKION
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIORD (IDIM, KDIM, ICKWRK, RCKWRK, NIORD, IORD, FORD,
     1                   RORD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIORD (IDIM, KDIM, ICKWRK, RCKWRK, NIORD, IORD, FORD,
C                     RORD)
C  Returns the count and indices of reactions with modified species
C  order and the order values for the species.
C
C  INPUT
C  IDIM      - Integer scalar, dimension of arrays IFORD and IRORD;
C              IDIM must be at least NIORD, the total number of
C              reactions with modified species orders.
C  KDIM      - Integer scalar, first dimension of the arrays FORD and
C              RORD;
C              KDIM must be at least NKK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NIORD      - Integer scalar, total number of reactions with modified
C              species orders.
C  IORD(*)   - Integer array, indices of reactions with modified 
C              species orders; dimension at least NIORD.
C  FORD(*,*) - Real matrix, the modified forward species orders for the
C              NIORD reactions;
C              dimension at least NKK for the first, the total species
C              count, and at least NIORD for the second.
C              FORD(K,N) is the forward order of species K for the Nth
C              change-order reaction.
C  RORD(*,*) - Real matrix, the modified reverse species orders for the
C              NIORD reactions;
C              dimension at least NKK for the first, the total species
C              count, and at least NRORD for the second.
C              RORD(K,N) is the reverse order of species K for the Nth
C              change-order reaction.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), IORD(*), FORD(KDIM,*),
     1          RORD(KDIM,*)
C
      DO 100 N = 1, IDIM
         IORD(N) = 0
         DO 50 K = 1, KDIM
            FORD(K,N) = 0.0
            RORD(K,N) = 0.0
   50    CONTINUE
  100 CONTINUE
C
      IF (IDIM.LE.0 .OR. IDIM.LT.NORD .OR. KDIM.LT.NKK)
     1   RETURN
C
      NIORD = NORD
      I_NK = IcKOR - 1
      I_ORD = NcKOR - 1
C
      DO 200 N = 1, NORD
         IORD(N)  = ICKWRK(IcORD + N - 1)
         DO 150 K = 1, MXORD
            KSPEC = ICKWRK(I_NK + K)
            IF (KSPEC .EQ. 0) THEN
               I_NK = I_NK + MXORD
               I_ORD = I_ORD + MXORD
               GO TO 200
            ELSEIF (KSPEC .LT. 0) THEN
               FORD(-KSPEC,N) = RCKWRK(I_ORD + K)
            ELSE
               RORD(KSPEC,N) = RCKWRK(I_ORD + K)
            ENDIF
  150    CONTINUE
  200 CONTINUE
C
C     end of SUBROUTINE CKIORD
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIREV  (IR, ICKWRK, RCKWRK, IREV, RAR, RBR, RER)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIREV  (IR, ICKWRK, RCKWRK, IREV, RAR, RBR, RER)
C  Returns an integer flag to indicate whether reaction IR has an
C  explicitly assigned reverse rate constant.  It also returns the
C  reverse Arrhenius expression values for reaction IR,
C  if it was explicitly assigned in the Chemkin interpreter.
C  If reverse Arrhenius values were not explicitly assigned,
C  RAR, RBR and RER will be zero.
C
C  INPUT
C  IR        - Integer scalar, reaction index.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  IREV      - Integer scalar, 
C              1, reaction IR has explicit reverse rate parameters
C              0, no.
C  RAR       - Real scalar, explicit pre-exponential constants
C              for reaction IR.
C                 cgs units, mole-cm-sec-K
C  RBR       - Real scalar, explicit temperature dependence exponents
C              for reaction IR.
C  RER       - Real scalar, explicit activation energy for reaction IR.
C                 cgs units, Kelvins
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      IF (NREV.GT.0) THEN
         DO 100 N = 1, NREV
            I = ICKWRK(IcRV+N-1)
            IF (I .EQ. IR) THEN
               IREV = 1
               RAR = RCKWRK(NcRV+(N-1)*(NPAR+1))
               RBR = RCKWRK(NcRV+(N-1)*(NPAR+1)+1)
               RER = RCKWRK(NcRV+(N-1)*(NPAR+1)+2)
            ENDIF
 100     CONTINUE
      ELSE
         IREV = 0
         RAR = 0.0
         RBR = 0.0
         RER = 0.0
      ENDIF
C
C     end of SUBROUTINE CKIREV
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIRNU (IDIM, NDIM, ICKWRK, RCKWRK, NIRNU, IRNU, NSPEC,
     1                   KI, RNU)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIRNU (IDIM, NDIM, ICKWRK, RCKWRK, NIRNU, IRNU, NSPEC,
C                     KI, RNU)
C  Returns the count and indices of reactions with real stoichiometric
C  coefficients, counts of species in the reactions, and the species
C  indices and coefficients; see Eq. (50).
C
C  INPUT
C  IDIM      - Integer scalar, dimension of the arrays IRNU and NSPEC,
C              and the second dimension of matrices KI and RNU;
C              IDIM must be at least NIRNU, the number of reactions
C              with real stoichiometric coefficients.
C  NDIM      - Integer scalar, first dimension of matrices KI and RNU;
C              NDIM must be at least MAXSP, the maximum number of
C              species allowed in a reaction.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NIRNU     - Integer scalar, total number of reactions with real
C              stoichiometric coefficients.
C  IRNU(*)   - Integer array, indices of reactions with real
C              stoichiometric coefficients; dimension at least NIRNU.
C  NSPEC(*)  - Integer array, total number of species in a reaction;
C              dimension at least NIRNU.
C  KI(*,*)   - Integer matrix, species indices for species in the
C              NIRNU reactions; dimension at least MAXSP for the first,
C              and at least NIRNU for the second.
C              KI(M,N) is the species index of the Mth species in the
C              Nth real coefficient reaction.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for species
C              in the NIRNU reactions; dimension at least MAXSP for
C              the first, and at least NIRNU for the second.
C              RNU(M,N) is the stoichiometric coefficient of the Mth
C              species in the Nth real coefficient reaction, and
C              RNU(M,*) < 0 if the Mth species is a reactant;
C              RNU(M,*) > 0 if the Mth species is a product.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), IRNU(*), NSPEC(*),
     1          KI(NDIM,*), RNU(NDIM,*)
C
      NIRNU = NRNU
      DO 100 N = 1, IDIM
         NSPEC(N) = 0
         IRNU(N) = 0
         DO 50 M = 1, NDIM
            KI(M,N) = 0
            RNU(M,N) = 0.0
   50    CONTINUE
  100 CONTINUE
C
      IF (NRNU.LE.0 .OR. IDIM.LT.NRNU .OR. NDIM.LT.MXSP)
     1   RETURN
C
      I_NU = NcRNU - 1
      DO 200 N = 1, NRNU
         I = ICKWRK(IcRNU + N - 1)
         IRNU(N) = I
         NSPEC(N) = 0
C
         I_NK = IcNK + (I-1)*MXSP - 1
         DO 150 M = 1, MXSP
            K = ICKWRK(I_NK + M)
            IF (K .NE. 0) THEN
               NSPEC(N) = NSPEC(N) + 1
               KI(NSPEC(N),N) = K
               RNU(NSPEC(N),N) = RCKWRK(I_NU + M)
            ENDIF
  150    CONTINUE
         I_NU = I_NU + MXSP
  200 CONTINUE
C
C     end of SUBROUTINE CKIRNU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKISIG (P, T, X, XNUIK, KK, ICKWRK, RCKWRK, SIGI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKISIG (P, T, X, XNUIK, KK, ICKWRK, RCKWRK, SIGI)
C  Returns the ion species electrical conductivities given
C  collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUIK(*,*)- Real matrix, momentum-transfer collision frequencies
C              for the ions with the species;
C              dimension at least KK for the first, the total species
C              count, and at least NKKI for the second, the ion count.
C  KK        - Integer scalar, first dimension of XNUIK.
C
C  OUTPUT
C  SIGI(*)   - Real array, ion electrical conductivities (DC);
C              dimension at least NKKI, the total ion count.
C                 cgs units, GM/CM*S
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23, ECHRG=1.6022D-12, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23, ECHRG=1.6022E-12, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL, BIG, EXPARG
      DIMENSION T(*), X(*), XNUIK(KK,*), ICKWRK(*), RCKWRK(*),
     1          SIGI(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      DO 100 KI = 1, NKKI
         K = ICKWRK(IcKI + KI - 1)
         XDEN = TDEN * MAX(X(K),SMALLX)
         SIGI(KI) = 0.0
         DO 50 J = 1, NKK
            WTI = RCKWRK(NcWT + K - 1)
            WTJ = RCKWRK(NcWT + J - 1)
            RMASS = (WTI*WTJ)/(WTI+WTJ) / AVAG
            SIGI(KI) = SIGI(KI) + RMASS*XNUIK(J,KI)/(XDEN*ECHRG**2)
50       CONTINUE
         SIGI(KI) = 1.0/SIGI(KI)
100   CONTINUE
C
C     end of SUBROUTINE CKISIG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKITDE  (ICKWRK, RCKWRK, ITDE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKITDE (ICKWRK, RCKWRK, ITDE)
C  Returns a set of flags indicating whether the reactions are
C  non-thermal, and if so, returns the index of the species on
C  which the reaction depends.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ITDE(*)   - Integer array, electron-impact flags for reactions;
C              dimension at least II, the total reaction count.
C              ITDE(I)= -1  reaction I is not a third-body reactions
C              ITDE(I)=  K  reaction I is a third-body reaction with
C                        temperature dependence on species # K
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), ITDE(*)
C
      DO 100 I = 1, NII
         ITDE(I) = -1
  100 CONTINUE
      DO 150 N = 0, NTDE - 1
         ITDE(ICKWRK(IcTD + N)) = ICKWRK(IcTK + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKITDE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKITR  (ICKWRK, RCKWRK, ITHB, IREV)
C
C  START PROLOGUE
C
C  SUBROUTINE CKITR  (ICKWRK, RCKWRK, ITHB, IREV)
C  Returns a set of flags indicating whether the reactions are
C  reversible or whether they contain arbitrary third bodies
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ITHB(*)   - Integer array, third-body indices for reactions;
C              dimension at least II, the total reaction count.
C              ITHB(I)= -1  reaction I is not a third-body reactions
C              ITHB(I)=  0  reaction I is is a third-body reaction with
C                           no enhanced third body efficiencies
C              ITHB(I)=  N  reaction I is a third-body reaction with
C                        N species enhanced third-body efficiencies.
C
C  IREV(*)   - Integer array, reversibility indices and species
C              count (reactants plus products) for reactions;
C              dimension at least II, the total reaction count.
C              IREV(I)=+N, reversible reaction I has N species
C              IREV(I)=-N, irreversible reaction I has N species
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), ITHB(*), IREV(*)
C
      DO 100 I = 1, NII
         IREV(I) = ICKWRK(IcNS + I - 1)
         ITHB(I) = -1
  100 CONTINUE
      DO 150 N = 0, NTHB - 1
         ITHB(ICKWRK(IcTB + N)) = ICKWRK(IcKN + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKITR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIVIS (P, T, X, XNUIM, K, ICKWRK, RCKWRK, VISI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIVIS (P, T, X, XNUIM, K, ICKWRK, RCKWRK, VISI)
C  Returns the ion species viscosities given collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array of mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUIM     - Real scalar, momentum-transfer collision frequency
C              for an ion
C  K         - Integer scalar, species index of the ion
C  OUTPUT
C  VISI      - Real scalar, ion viscosity
C                 cgs units, GM/CM*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
C
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      TI = T(ICKWRK(IcKTF + K - 1))
      XDEN = TDEN * MAX(X(K),SMALL)
      VISI = (4.*BOLTZ*TI*XDEN)/(PI*MAX(XNUIM,SMALL))
C
C     end of SUBROUTINE CKIVIS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIXSM  (ICKWRK, IXSM, IXSK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIXSM (ICKWRK, IXSM, IXSK)
C  Returns a set of flags indicating whether the reactions are ion
C  momentum-transfer cross sections.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  IXSM(*)   - Integer array, ion momentum-transfer cross-section flag;
C              dimension at least II, the total reaction count.
C              IXSM(I)= -1  reaction I is not a ion mom-transfer x-sec
C              IXSM(I)=  KI reaction I is a ion mom-trans cross-section
C                        and KI is the ion species index
C  IXSK(*)   - Integer array, species indices for the collision partner
C              of the ion momentum-transfer cross-section reactions;
C              dimension at least II, the total reaction count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), IXSM(*), IXSK(*)
C
      DO 100 I = 1, NII
         IXSM(I) = -1
         IXSK(I) = 0
  100 CONTINUE
      DO 150 N = 0, NXSM - 1
         IXSM(ICKWRK(IcXS + N)) = ICKWRK(IcXI + N)
         IXSK(ICKWRK(IcXS + N)) = ICKWRK(IcXK + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKIXSM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKFKR (P, T, X, ICKWRK, RCKWRK, FWDK, REVK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKFKR (P, T, X, ICKWRK, RCKWRK, FWDK, REVK)
C  Returns the forward and reverse reaction rates for reactions
C  given pressure, temperature(s) and mole fractions.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  FWDK(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  REVK(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), FWDK(*), REVK(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKRATX (T, RCKWRK(NcK1), ICKWRK(IcNU), ICKWRK(IcNK),
     1             RCKWRK(NcCO), ICKWRK(IcFL), ICKWRK(IcFO),
     2             ICKWRK(IcFT), ICKWRK(IcKF), RCKWRK(NcFL),
     3             ICKWRK(IcTB), ICKWRK(IcKN), RCKWRK(NcKT),
     4             ICKWRK(IcKT), RCKWRK(NcKF), RCKWRK(NcKR), FWDK,
     5             REVK, RCKWRK(NcI3), ICKWRK(IcRNU), RCKWRK(NcRNU),
     6             ICKWRK(IcORD), ICKWRK(IcKOR), RCKWRK(NcKOR),
     7             ICKWRK(IcMO), ICKWRK(IcXS))
C
C     end of SUBROUTINE CKKFKR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKFRT (P, T, ICKWRK, RCKWRK, RKFT, RKRT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKFRT (P, T, ICKWRK, RCKWRK, RKFT, RKRT)
C  Returns the forward and reverse reaction rates for reactions
C  given pressure and temperature(s).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), RKFT(*), RKRT(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RKFT, RKRT,
     9             RCKWRK(NcI1))
C
C     end of SUBROUTINE CKKFRT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKTFL (ICKWRK, KTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKTFL (ICKWRK, KTFL)
C  Allows the user to assign a location in the temperature array
C  to use for each gas-phase species.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  KTFL(*)   - Integer array, indices into the temperature(s) for
C              species;
C              dimension at least KK, the total species count.
C              Default value stored in ICKWRK is set to 1 in CKINIT.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), KTFL(*)
C
      DO 100 K = 1, NKK
         ICKWRK(IcKTF + K - 1) = KTFL(K)
  100 CONTINUE
C
C     end of SUBROUTINE CKKTFL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKLEN (LINC, LOUT, LI, LR, LC, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE CKLEN (LINC, LOUT, LENI, LENR, LENC, IFLAG)
C   Returns the lengths required for work arrays.
C
C  INPUT
C  LINC     - Integer scalar, input file unit for the linkfile.
C  LOUT     - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  LENI     - Integer scalar, minimum length required for the
C             integer work array.
C  LENR     - Integer scalar, minimum length required for the
C             real work array.
C  LENC     - Integer scalar, minimum length required for the
C             character work array.
C  IFLAG    - Integer scalar, indicates successful reading of
C             linkfile; IFLAG>0 indicates error type.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NLIST = 1)
      LOGICAL KERR, VOK, POK, LBIN
      CHARACTER*16 LIST(NLIST), FILVER, PREC, PRVERS, IFMT, CFMT,
     1             RFMT, LFMT
      PARAMETER (IFMT='(10I12)', CFMT='(8A16)', RFMT='(1P,5E24.16)',
     1           LFMT='(L8)')
C
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
      DATA LIST(1) /'1.0'/
C
      FILVER = ' '
      PRVERS   = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0
C
      KERR = .FALSE.
      IFLAG = 0
      REWIND LINC
      NREC = 1
C*****linkfile (gas) > binary
C      LBIN = .TRUE.
C*****END linkfile (gas) > binary
C*****linkfile (gas) > ascii
      LBIN = .FALSE.
C*****END linkfile (gas) > ascii
C
      IF (LBIN) THEN
         READ (LINC, ERR=100) FILVER
      ELSE
         READ (LINC, CFMT, ERR=100) FILVER
      ENDIF
      CALL CKCOMP (FILVER, LIST, NLIST, IND)
      IF (IND .LE. 0) THEN
         VOK = .FALSE.
      ELSE
         VOK = .TRUE.
      ENDIF
C
      IF (LBIN) THEN
         NREC = 2
         READ (LINC, ERR=100) PRVERS
         NREC = 3
         READ (LINC, ERR=100) PREC
         NREC = 4
         READ (LINC, ERR=100) KERR
      ELSE
         NREC = 2
         READ (LINC, CFMT, ERR=100) PRVERS
         NREC = 3
         READ (LINC, CFMT, ERR=100) PREC
         NREC = 4
         READ (LINC, LFMT, ERR=100) KERR
      ENDIF
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the Chemkin linkfile...',
     2      ' Check CHEMKIN INTERPRETER output for error conditions.'
         ENDIF
         IF (.NOT. VOK) WRITE (LOUT, '(/A)')
     1   ' Chemkin linkfile is incompatible with Chemkin-III Library'
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Precision of Chemkin linkfile does not agree with',
     2      ' precision of Chemkin library'
         ENDIF
         IFLAG = 20
         REWIND LINC
         RETURN
      ENDIF
C
      NREC = 5
      IF (LBIN) THEN
         READ (LINC, ERR=100) LENICK, LENRCK, LENCCK
      ELSE
         READ (LINC, IFMT, ERR=100) LENICK, LENRCK, LENCCK
      ENDIF
      REWIND LINC
C
      LENI = LENICK
      LENR = LENRCK
      LENC = LENCCK
      LI   = LENI
      LR   = LENR
      LC   = LENC
      RETURN
C
  100 CONTINUE
      IFLAG = NREC
      WRITE (LOUT, '(/A,/A,I5)')
     1   ' Error reading linkfile,',
     2   ' SUBROUTINE CKLEN record index #', IFLAG
      REWIND LINC
C
C     Generic Formats - limit lines to 132 characters
C
8001  FORMAT (10I12)
8002  FORMAT (1P,5E24.16)
8003  FORMAT (8A16)
C
C     end of SUBROUTINE CKLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKDLIM (STRING, DELIM, I1, I2)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDLIM
C  returns pointers into a character string of the first and
C  second occurrences of a particular character.
C
C  Arguments:
C  STRING - Character string.
C  DELIM  - Single character.
C  I1     - Integer scalar, location in STRING of first DELIM.
C  I2     - Integer scalar, location in STRING of second DELIM.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), DELIM*1
C
      I1 = INDEX(STRING, DELIM)
      I2 = I1 + INDEX(STRING(I1+1:), DELIM)
C
C     end of SUBROUTINE CKDLIM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKDTAB (STRING)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDTAB (STRING)
C  Replaces any tab character in a character string with one space.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), TAB*1
C
      TAB = CHAR(9)
   10 CONTINUE
      IND = INDEX(STRING,TAB)
      IF (IND .GT. 0) THEN
         STRING(IND:IND) = ' '
         GO TO 10
      ENDIF
C
C     end of SUBROUTINE CKDTAB
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      CHARACTER*(*) FUNCTION CKCHUP(ISTR, ILEN)
      CHARACTER*(*) ISTR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      CKCHUP = ' '
      CKCHUP = ISTR(1:ILEN)
      JJ = MIN (LEN(CKCHUP), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 05 N = 1,26
            IF (ISTR(J:J) .EQ. LCASE(N)) CKCHUP(J:J) = UCASE(N)
   05    CONTINUE
   10 CONTINUE
C
C     end of FUNCTION CKCHUP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      CHARACTER*(*) FUNCTION CKCHLO(ISTR, ILEN)
      CHARACTER*(*) ISTR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      CKCHLO = ' '
      CKCHLO = ISTR(1:ILEN)
      JJ = MIN (LEN(CKCHLO), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 05 N = 1,26
            IF (ISTR(J:J) .EQ. UCASE(N)) CKCHLO(J:J) = LCASE(N)
   05    CONTINUE
   10 CONTINUE
C
C     end of FUNCTION CKCHLO
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      INTEGER FUNCTION CKLKUP (ITEM, LIST, NLIST)
C
C  START PROLOGUE
C
C  INTEGER FUNCTION CKLKUP (ITEM, LIST, NLIST)
C
C      Looks up an item in an integer list. If an item is found,
C   it returns the first position of the item in the list. If an
C   item is not found, this routine returns the value 0.
C
C  INPUT
C
C  ITEM    - Integer scalar; Item to look up in the list
C  LIST(*) - Integer array;  List of entries
C  NLIST   - Integer scalar; Number of entries in the list
C
C  END PROLOGUE
C
      INTEGER   ITEM, LIST, NLIST, I
      DIMENSION  LIST(NLIST)
C
      CKLKUP = 0
      DO 10 I = NLIST, 1, -1
        IF (LIST(I) .EQ. ITEM) CKLKUP = I
   10 CONTINUE
C
C     end of FUNCTION CKLKUP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION CKFRCH (STR)
C
C  START PROLGUE
C
C  INTEGER FUNCTION CKFRCH (STR)
C
C  Returns the index of the first non-blank, non-tab character in
C  a string.
C
C  INPUT
C  STR   - Character string
C
C  END PROLOGUE
C
      CHARACTER STR*(*), TAB*1
      INTEGER ILEN, I
C
      ILEN = LEN(STR)
      TAB  = CHAR(9)
      CKFRCH = 0
      DO 10 I = 1, ILEN
         IF (STR(I:I).EQ.' ' .OR. STR(I:I).EQ.TAB) GO TO 10
         CKFRCH = I
         RETURN
   10 CONTINUE
C
C     end of FUNCTION CKFRCH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION CKLSCH (STR)
C
C  START PROLOGUE
C
C  INTEGER FUNCTION CKLSCH (STR)
C
C  Returns the index of the final non-blank, non-tab character in
C  a string.
C
C  INPUT
C  STR   - Character string
C
C  END PROLOGUE
C
      CHARACTER STR*(*), TAB*1, NUL*1
      INTEGER ILEN, I
C
      ILEN = LEN(STR)
      CKLSCH = 0
      TAB = CHAR(9)
      NUL = CHAR(0)
      DO 10 I = ILEN, 1, -1
         IF (STR(I:I).EQ.' ' .OR. STR(I:I).EQ.TAB .OR.
     1       STR(I:I).EQ.NUL) GO TO 10
         CKLSCH = I
         RETURN
   10 CONTINUE
C
C     end of FUNCTION CKLSCH
      RETURN
      END
C                                                                      C
      INTEGER FUNCTION CKSLEN (LINE)
C
C  BEGIN PROLOGUE
C
C  INTEGER FUNCTION CKSLEN (LINE)
C  Returns the effective length of a character string, i.e.,
C  the index of the last character before an exclamation mark (!)
C  indicating a comment.
C
C  INPUT
C  LINE     - Character string.
C
C  OUTPUT
C  CKSLEN   - Integer scalar, the effective length of LINE.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*)
      INTEGER CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
C
      IND = CKFRCH(LINE)
      IF (IND.EQ.0 .OR. LINE(IND:IND).EQ.'!') THEN
         CKSLEN = 0
      ELSE
         IND = INDEX(LINE,'!')
         IF (IND .GT. 0) THEN
            CKSLEN = CKLSCH(LINE(1:IND-1))
         ELSE
            CKSLEN = CKLSCH(LINE)
         ENDIF
      ENDIF
C
C     end of FUNCTION CKSLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXMIN (X, NN, XMIN, IMIN)
C
C  START PROLOGUE
C  Returns the minimum value in an array and its location in the array.
C
C  INPUT
C  X(*)      - Real array.
C  NN        - Integer scalar; size of X.
C  OUTPUT
C  XMIN      - Real scalar.
C  IMIN      - Integer scalar; location in X of XMIN.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(NN)
C
      IMIN = 1
      XMIN = X(IMIN)
      DO 10 N = 1, NN
         IF (X(N) .LT. XMIN) THEN
            XMIN = X(N)
            IMIN = N
         ENDIF
   10 CONTINUE
C
C     end of CKXMIN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXMAX (X, NN, XMAX, IMAX)
C
C  START PROLOGUE
C  Returns the maximum value in an array and its location in the array.
C
C  INPUT
C  X(*)      - Real array.
C  NN        - Integer scalar; size of X.
C  OUTPUT
C  XMAX      - Real scalar.
C  IMAX      - Integer scalar; location in X of XMAX.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(NN)
C
      IMAX = 1
      XMAX = X(IMAX)
      DO 10 N = 1, NN
         IF (X(N) .GT. XMAX) THEN
            XMAX = X(N)
            IMAX = N
         ENDIF
   10 CONTINUE
C
C     end of CKXMAX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWC (C, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWC (C, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given molar
C  concentrations;  see Eq. (5).
C
C  INPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION C(*), ICKWRK(*), RCKWRK(*)
C
      CTOT = 0.0
      WTM  = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
         WTM  = WTM  + C(K) * RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      WTM = WTM / CTOT
C
C     end of SUBROUTINKE CKMMWC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given mole
C  fractions;  see Eq. (4).
C
C  INPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
C
      WTM = 0.0
      DO 100 K = 1, NKK
         WTM = WTM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKMMWX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given mass
C  fractions;  see Eq. (3).
C
C  INPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
C
      SUMYOW=0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      WTM = 1.0 / SUMYOW
C
C     end of SUBROUTINE CKMMWY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMXTP (ICKWRK, MAXTP)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMXTP (ICKWRK, MAXTP)
C  Returns the maximum number of temperatures used in fitting the
C  thermodynamic properties of the species.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  MXTP      - Integer scalar, maximum number of temperatures used
C              to fit the thermodynamic properties of the species.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*)
C
      MAXTP = MXTP
C
C     end of SUBROUTINE CKMXTP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNCF  (MDIM, ICKWRK, RCKWRK, NCF)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNCF  (MDIM, ICKWRK, RCKWRK, NCF)
C  Returns the elemental composition of the species
C
C  INPUT
C  MDIM      - Integer scalar, first dimension of the matrix NCF;
C              MDIM must be equal to or greater than MM, the total
C              element count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NCF(*,*)  - Real matrix, the elemental composition of the species;
C              dimension at least MM for the first, the total element
C              count, and at least KK for the second, the total species
C              count.
C              NCF(M,K) is the quantity of the element M in species K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NCF(MDIM,*)
C
      I_KNCF = IcNC - 1
      DO 150 K = 1, NKK
         DO 100 M = 1, NMM
            I_KNCF = I_KNCF + 1
            NCF(M,K) = ICKWRK(I_KNCF)
  100    CONTINUE
150   CONTINUE
C
C     end of SUBROUTINE CKNCF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNORM (ARRAY, NN)
C
C START PROLOGUE
C
C SUBROUTINE CKNORM (ARRAY, NN)
C Utility to normalize the real members of an array.
C
C INPUT
C ARRAY(*)  - Real array.
C NN        - Integer scalar; the size of ARRAY.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      SUM = CKSUM (ARRAY, NN)
      IF (SUM .NE. 0.0) THEN
         DO 20 N = 1, NN
            ARRAY(N) = ARRAY(N) / SUM
   20    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKNORM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKSCAL (ARRAY, NN, SCAL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSCAL (ARRAY, NN)
C  Utility to scale the real members of an array.
C
C  INPUT
C  ARRAY(*)  - Real array.
C  NN        - Integer scalar; the size of ARRAY.
C  SCAL      - Real scalar; the multiplier for ARRAY members.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      DO 10 N = 1, NN
         ARRAY(N) = SCAL * ARRAY(N)
   10 CONTINUE
C
C     end of SUBROUTINE CKSCAL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C*****precision > double
      DOUBLE PRECISION FUNCTION CKSUM (ARRAY, NN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION CKSUM (ARRAY, NN)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      CKSUM = 0.0
      DO 10 N = 1, NN
         CKSUM = CKSUM + ARRAY(N)
   10 CONTINUE
C
C     end of FUNCTION CKSUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNPAR (LINE, NPAR, LOUT, IPAR, ISTART, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNPAR (LINE, NPAR, LOUT, IPAR, ISTART, KERR)
C  Searches a character string LINE from last to first character,
C  to create a substring IPAR containing NPAR blank-delimited numbers;
C  ISTART is the column of LINE containing IPAR. This allows format-
C  free input of combined alpha-numeric data.  For example,
C
C     input:  LINE*80   = "t1 t2 dt  300.0  3.0E3  50"
C             NPAR      = 3, the number of substrings requested
C             LOUT      = 6, a logical unit number on which to write
C                         diagnostic messages.
C     output: IPAR*80   = "300.0  3.0E3  50"
C             ISTART    = 13, the starting column in LINE of the
C                         NPAR substrings
C             KERR      = .FALSE.
C
C  INPUT
C  LINE      - Character string; length determined by calling routine.
C  NPAR      - Integer scalar, number of substrings expected.
C  LOUT      - Integer scalar, output unit for printed diagnostics.
C
C  OUTPUT
C  IPAR      - Character string, subset of LINE, containing only the
C              NPAR substrings.
C  ISTART    - Integer scalar, starting location in LINE of the NPAR
C              substrings.
C  KERR      - Logical, character length or syntax error flag.
C
C  END PROLOGUE
C
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), IPAR*(*)
      LOGICAL FOUND, KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
C----------Find Comment String (! signifies comment)
C
      ILEN = CKSLEN(LINE)
      KERR = .FALSE.
C
      IF (ILEN.GT.0) THEN
         FOUND = .FALSE.
         N = 0
         DO 40 I = ILEN, 1, -1
            IF (FOUND) THEN
               IF (LINE(I:I).EQ.' ') THEN
                  N = N+1
                  FOUND = .FALSE.
                  IF (N.EQ.NPAR) THEN
                     ISTART = I+1
                     L1 = ILEN - ISTART + 1
                     L2 = LEN(IPAR)
                     IF (L2 .GE. L1) THEN
                        IPAR = LINE(ISTART:ILEN)
                     ELSE
                        WRITE (LOUT,*)
     1               ' Error in CKNPAR...character length too small...'
                        KERR = .TRUE.
                     ENDIF
                     GO TO 100
                  ENDIF
               ENDIF
            ELSE
               IF (LINE(I:I).NE.' ') FOUND = .TRUE.
            ENDIF
   40    CONTINUE
      ENDIF
C
      WRITE (LOUT,*) ' Error in CKNPAR...',NPAR,' values not found...'
      KERR = .TRUE.
C
C     end of SUBROUTINE CKNPAR
  100 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNU   (KDIM, ICKWRK, RCKWRK, NUKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNU   (KDIM, ICKWRK, RCKWRK, NUKI)
C  Returns the stoichiometric coefficients of the reactions;
C  see Eq. (50).
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix NUKI;
C              KDIM must be greater than or equal to KK, the total
C              species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NUKI(*,*) - Integer matrix, stoichiometric coefficients of the
C              species in the reactions;  dimension at least KK for
C              the first, the total species count, and at least II
C              for the second, the total reaction count.
C              NUKI(K,I) is the stoichiometric coefficient of
C              species K in reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NUKI(KDIM,*)
C
      DO 100 I = 1, NII
         DO 50 K = 1, NKK
            NUKI(K,I) = 0
   50    CONTINUE
  100 CONTINUE
C
      I_NK = IcNK - MXSP
      I_NU = IcNU - MXSP
      DO 200 I = 1, NII
         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
         DO 150 N = 0, MXSP - 1
            K = ICKWRK(I_NK + N)
            IF (K .NE. 0) NUKI(K,I) = NUKI(K,I) + ICKWRK(I_NU + N)
  150    CONTINUE
200   CONTINUE
C
C     end of SUBROUTINE CKNU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNUF   (KDIM, ICKWRK, RCKWRK, NUFKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNUF   (KDIM, ICKWRK, RCKWRK, NUKI)
C  Returns the forward stoichiometric coefficients for reactions;
C  by definition, reactants' coefficients are negative;  see Eq. (50).
C  Contrast this subroutine with subroutine CKNU, which returns the
C  net stoichiometric coefficients for a reaction.
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix NUKI;
C              KDIM must be greater than or equal to KK, the total
C              species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NUFKI(*,*)- Integer matrix, stoichiometric coefficients of the
C              species in the forward direction of the reactions
C              (reactants only); dimension at least KK in the first,
C              the total species count, and at least II for the
C              second, the total reaction count.
C              NUKI(K,I) is the stoichiometric coefficient of
C              species K in forward direction of reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NUFKI(KDIM,*)
C
      DO 100 I = 1, NII
         DO 50 K = 1, NKK
            NUFKI(K,I) = 0
   50    CONTINUE
  100 CONTINUE
C
      I_NK = IcNK - MXSP
      I_NU = IcNU - MXSP
      NREAC = MXSP / 2
      DO 200 I = 1, NII
         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
         DO 150 N = 0, NREAC - 1
            K = ICKWRK(I_NK + N)
            IF (K .EQ. 0) GO TO 200
            NUFKI(K,I) = NUFKI(K,I) + ICKWRK(I_NU + N)
  150    CONTINUE
  200 CONTINUE
C
C     end of SUBROUTINE CKNUF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and molar concentrations;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*)
C
      CTTOT = 0.0
      DO 100 K = 1, NKK
         CTTOT = CTTOT + C(K) * T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
      P    = RCKWRK(NcRU) * CTTOT
C
C     end of SUBROUTINE CKPC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPHAZ (ICKWRK, RCKWRK, KPHASE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPHAZ (ICKWRK, RCKWRK, KPHASE)
C  Returns a set of flags indicating phases of the species
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  KPHASE(*) - Integer array, phases of the species;
C              dimension at least KK, the total species count.
C              KPHASE(K)=-1, species K is solid
C              KPHASE(K)= 0, species K is gaseous
C              KPHASE(K)=+1, species K is liquid
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), KPHASE(*)
C
      DO 100 K = 1, NKK
         KPHASE(K) = ICKWRK(IcPH + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKPHAZ
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, LC, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPNT (LSAVE, LOUT, NPOINT, VERS, PREC, LENI, LENR,
C                    LENC, KERR)
C  Reads from a file information about a Chemkin linkfile, and
C  pointers for work arrays.
C
C  INPUT
C  LSAVE     - Integer scalar, input unit for binary data file.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  NPOINT    - Integer scalar, total pointers count.
C  VERS      - Real scalar, version number of the Chemkin linkfile.
C  PREC      - Character string, machine precision of the linkfile.
C  LENI      - Integer scalar, length required for integer work array.
C  LENR      - Integer scalar, length required for real work array.
C  LENC      - Integer scalar, length required for character work array.
C  KERR      - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
C
C     Data about the machine dependent constants is carried in
C
      COMMON/MACH/SMALL,BIG,EXPARG
      LOGICAL KERR, IERR
      CHARACTER*16 PREC, FILVER, PRVERS, P, V
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
      EXPARG = LOG(BIG)
C
      KERR = .FALSE.
      READ (LSAVE, ERR=100)
     *                FILVER,   PREC,   LENI,   LENR,   LENC,
C
C     include file for CHEMKIN-III cklib.f, dated: March 1, 1966
C
C     Integer constants
C
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     4   MXORD, KEL, NKKI,
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
C
C     END include file for cklib.f
C
      NPOINT = 102
      V = FILVER
      P = PREC
      LI = LENI
      LR = LENR
      LC = LENC
      IERR = KERR
      RETURN
C
  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Chemkin linkfile data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      FILVER   = ' '
      PRVERS     = ' '
      V      = FILVER
      PREC   = ' '
      P      = PREC
C
C     end of SUBROUTINE CKPNT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPX   (RHO, T, X, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPX   (RHO, T, X, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and mole fractions;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      SUMXW = 0.0
      SUMXT = 0.0
      DO 100 K = 1, NKK
         SUMXW = SUMXW + X(K)*RCKWRK(NcWT + K - 1)
         SUMXT = SUMXT + X(K)*T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
      P = RHO * RCKWRK(NcRU) * SUMXT / SUMXW
C
C     end of SUBROUTINE CKPX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and mass fractions;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      SUMYOW = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUMYOW = SUMYOW + Y(K+1) * T(ICKWRK(IcKTF + K))
     1                            /RCKWRK(NcWT + K)
150   CONTINUE
      P = RHO * RCKWRK(NcRU) * SUMYOW
C
C     end of SUBROUTINE CKPY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQC   (T, C, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQC   (T, C, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given temperature(s)
C  and molar concentrations;  see Eqs. (51) and (58).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKRATX (T, C, ICKWRK(IcNU), ICKWRK(IcNK), RCKWRK(NcCO),
     2             ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcFT),
     3             ICKWRK(IcKF), RCKWRK(NcFL), ICKWRK(IcTB),
     3             ICKWRK(IcKN), RCKWRK(NcKT), ICKWRK(IcKT),
     4             RCKWRK(NcKF), RCKWRK(NcKR), RCKWRK(NcI1),
     5             RCKWRK(NcI2), RCKWRK(NcI3), ICKWRK(IcRNU),
     6             RCKWRK(NcRNU), ICKWRK(IcORD), ICKWRK(IcKOR),
     7             RCKWRK(NcKOR), ICKWRK(IcMO), ICKWRK(IcXS))
C
      NIM1 = NII - 1
      DO 100 I = 0, NIM1
         Q(I+1) = RCKWRK(NcI1 + I) - RCKWRK(NcI2 + I)
  100 CONTINUE
C
C     end of SUBROUTINE CKQC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQNUE (Q, T, P, C, ICKWRK, RCKWRK, XNUES, XNUEH)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQNUE (Q, T, P, C, ICKWRK, RCKWRK, XNUES, XNUEH)
C  Returns the electron momentum-transfer collision frequencies for
C  all species, using averages for unspecified species
C
C  INPUT
C  Q(*)      - Real array, rates of production of the reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm2
C  C(*)      - Real array, concentrations of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  OUTPUT
C  XNUES(*)  - Real array, momentum-transfer collision frequencies for
C              the electrons with the species;
C              dimension at least KK, the total species count.
C  XNUEH     - Real scalar, total momentum-transfer collision frequency
C              for the electrons
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, PERM0 = 8.854D-7, ECHRG=1.6022D-12,
     1           AVAG = 6.022D23, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, PERM0 = 8.854E-7, ECHRG=1.6022E-12,
C     1           AVAG = 6.022E23, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION Q(*), T(*), C(*), ICKWRK(*), RCKWRK(*), XNUES(*)
C
      DO 40 K = 1, NKK
         XNUES(K) = 0.0
 40   CONTINUE
C
      XNUEH = 0.0
      SUMQN = 0.0
      DO 100 I1 = 1, NMOM
         I = ICKWRK(IcMO + I1 - 1)
         K = ICKWRK(IcMK + I1 -1)
         XNUES(K) = Q(I) * C(K)*AVAG
         XNUEH = XNUEH + XNUES(K)
         SUMQN = SUMQN + Q(I)
 100  CONTINUE
      IF (NMOM.NE.0) THEN
         QAV = SUMQN / NMOM
      ELSE
         QAV = 0.0
      ENDIF
C
C    CALCULATE NUEI AND NUEE:
C
      EDEN = C(KEL)*AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      EMASS = RCKWRK(NcWT + KEL - 1)/AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      DO 200 KI = 1, NKKI
         K = ICKWRK(IcKI + KI - 1)
         IF (EDEN.GT.1.0) THEN
            TK = T(ICKWRK(IcKTF + K - 1))
            XDEN = C(K)*AVAG
            Z = ABS(FLOAT(ICKWRK(IcCH + K - 1)))
            XLAM = (PERM0 * BOLTZ * TE / ECHRG**2)**1.5
     1           * 12.0 * PI / (Z*SQRT(EDEN))
            XNUES(K) = XDEN * LOG(XLAM) * (4.*SQRT(2.*PI)/3.)
     1           * (EMASS/(BOLTZ*TK))**1.5
     2           * ((ECHRG**2)/(4.*PI*PERM0*EMASS))**2
         ELSE
            XNUES(K) = SMALL
         ENDIF
200   CONTINUE
      IF (EDEN.GT.1.0) THEN
         XNUES(KEL) = EDEN * (8./3.) * SQRT(PI) *
     1                (EMASS/(BOLTZ*TE))**1.5 *
     2            ((ECHRG**2)/(4*PI*PERM0*EMASS))**2 * LOG(XLAM)
      ELSE
         XNUES(KEL) = SMALL
      ENDIF
C
C    FILL IN ANY GAPS IN DATA PROVIDED BY AVERAGE NEUTRAL FREQ
C
      DO 110 K = 1, NKK
         IF (XNUES(K).LE.0.0) THEN
            XNUES(K) = QAV * C(K)*AVAG
            XNUEH = XNUEH + XNUES(K)
         ENDIF
110   CONTINUE
C
C     end of SUBROUTINE CKQNUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQNUI (Q, T, P, C, KDIM, ICKWRK, RCKWRK, XNUIK, XNUIM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQNUI (Q, T, P, C, KDIM, ICKWRK, RCKWRK, XNUIK, XNUIM)
C  Returns the ion momentum-transfer collision frequencies from
C  collision cross-sections for all species, using averages for
C  unspecified species.
C
C  INPUT
C  Q(*)      - Real array, rates of production for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm2.
C  C(*)      - Real array, concentrations of the mixtures;
C              dimension at least KK, the total species count.
C  KDIM      - Integer scalar, first dimension of the matrix XNUIK;
C              KDIM must be at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  XNUIK(*,*)- Real matrix, momentum-transfer collision frequencies
C              for the ions with all other species;
C              dimension at least KK for the first, the total
C              species count, and at least NKKI for the second, the
C              total ion count.
C  XNUIM(*)  - Real array, momentim-averaged momentum-transfer
C              collision frequencies for the ions;
C              dimension at least NKKI, the total ion count.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, PERM0 = 8.854D-7, ECHRG=1.6022D-12,
     1           AVAG = 6.022D23, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, PERM0 = 8.854E-7, ECHRG=1.6022E-12,
C     1           AVAG = 6.022E23, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION Q(*), T(*), C(*), ICKWRK(*), RCKWRK(*),
     1          XNUIK(KDIM,*), XNUIM(*)
C
      DO 10 KI = 1, NKKI
         DO 5 K = 1, NKK
            XNUIK(K,KI) = 0.0
 5       CONTINUE
         XNUIM(KI) = 0.0
 10   CONTINUE
C
      EDEN = C(KEL)*AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      BOLTZ = RCKWRK(NcRU)/AVAG
      DEBYE = SQRT (PERM0 * BOLTZ * TE / (EDEN * ECHRG**2))
      EMASS = RCKWRK(NcWT+KEL-1)/AVAG
C
      SUMQN = 0.0
      NEUT = 0
      DO 100 N = 1, NXSM
         I = ICKWRK(IcXS + N - 1)
         QXS = Q(I)
         K = ICKWRK(IcXK + N - 1)
         KIK = ICKWRK(IcXI + N - 1)
         DO 50 KI1 = 1, NKKI
            IF (KIK .EQ. ICKWRK(IcKI + KI1 - 1)) KI = KI1
 50      CONTINUE
         XMSS = RCKWRK(NcWT + KIK - 1)/AVAG
         TI = T(ICKWRK(IcKTF + KIK - 1))
         TK = T(ICKWRK(IcKTF + K - 1))
C         RMASS = (RCKWRK(NcWT + KIK - 1)*RCKWRK(NcWT + K - 1))
C     1              / (RCKWRK(NcWT+KIK-1) + RCKWRK(NcWT+K-1)) / AVAG
         XDEN = C(K)*AVAG
         XNUIK(K,KI) = XDEN * SQRT(3.*BOLTZ*TI/XMSS) * QXS
         SUMQN = SUMQN + QXS
         NEUT = NEUT + 1
         WTFAC = 2.*RCKWRK(NcWT + K - 1)
     1                / (RCKWRK(NcWT + KIK - 1)+RCKWRK(NcWT + K - 1))
         XNUIM(KI) = XNUIM(KI) + XNUIK(K,KI)*WTFAC
 100  CONTINUE
C
C  USE AVERAGE X-SECTION OF NEUTRAL SPECIES TO DETERMINE MISSING DATA
C
      QNAVG = SUMQN / NEUT
      DO 200 KI = 1, NKKI
         KIK = ICKWRK(IcKI + KI - 1)
         XMSS = RCKWRK(NcWT + KIK - 1)/AVAG
         TI = T(ICKWRK(IcKTF + KIK - 1))
         DO 150 K = 1, NKK
            IF (ICKWRK(IcCH+K-1).EQ.0 .AND. XNUIK(K,KI) .LE. 0.0) THEN
               XDEN = C(K)*AVAG
               XNUIK(K,KI) = XDEN * SQRT(3.*BOLTZ*TI/XMSS) * QNAVG
               WTFAC = 2.*RCKWRK(NcWT + K - 1)
     1                  / (RCKWRK(NcWT+KI-1)+RCKWRK(NcWT + K - 1))
               XNUIM(KI) = XNUIM(KI) + XNUIK(K,KI)*WTFAC
            ELSEIF (ICKWRK(IcCH + K - 1).NE.0 .AND. K.NE.KEL) THEN
               XDEN = C(K)*AVAG
               TK = T(ICKWRK(IcKTF + K - 1))
               B0FAC = ECHRG**2 / (12. * PI * PERM0 * BOLTZ * TK)
               B0BAR = ABS(ICKWRK(IcCH + K - 1))
     1                 * ABS(ICKWRK(IcCH + KIK - 1)) * B0FAC
               XLAM = DEBYE / B0BAR
               XNUIK(K,KI) = XDEN * (8./3.) * SQRT(PI) *
     1                      (EMASS/(BOLTZ*TK))**1.5 *
     2                      (ECHRG**2/(4*PI*PERM0*EMASS))**2 * LOG(XLAM)
            ELSEIF (K.EQ.KEL) THEN
               IF (EDEN.GT.1.) THEN
                  XIDEN = C(KIK)*AVAG
                  Z = ABS(ICKWRK(IcCH + KIK - 1))
                  XLAM = (PERM0 * BOLTZ * TE / (ECHRG**2))**1.5
     1                * 12.0 * PI / (Z * SQRT(EDEN))
                  XNUIK(KEL,KI) = XIDEN * LOG(XLAM)
     1                         * (4.0*SQRT(2.*PI)/3.0)
     2                         * (EMASS/(BOLTZ*TI))**1.5
     3                         * ((ECHRG**2)/(4.*PI*PERM0*EMASS))**2
               ELSE
                  XNUIK(KEL,KI) = SMALL
               ENDIF
            ENDIF
 150     CONTINUE
 200  CONTINUE
C
C     end of SUBROUTINE CKQNUI
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQXP  (P, T, X, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQXP  (P, T, X, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given pressure,
C  temperature(s) and mole fractions;  see Eqs. (51) and (58).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQYP  (P, T, Y, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQYP  (P, T, Y, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given pressure,
C  temperature(s) and mass fractions;  see Eqs. (51) and (58).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, Mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQYR  (RHO, T, Y, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQYR  (RHO, T, Y, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given mass density,
C  temperature(s) and mass fractions; see Eqs. (51) and (58).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKR2CH (RNUM, STR, I, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKR2CH (RNUM, STR, I, KERR)
C  Returns a character string representation of a real number
C  and the effective length of the string.
C
C  INPUT
C  RNUM      - Real scalar, to be converted to a string;
C              the maximum magnitude of RNUM is machine-dependent.
C
C  OUTPUT
C  STR      - Character string, left-justified representation of RNUM;
C             i.e., RNUM=  0.0      returns STR=" 0.00"
C                   RNUM= -10.5     returns STR="-1.05E+01"
C                   RNUM= 1.86E-100 returns in STR=" 1.86E-100"
C                   the minimum length of STR required is 5
C  I        - Integer scalar, total non-blank characters in RNUM.
C  KERR     - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STR*(*)
      LOGICAL KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      ILEN = LEN(STR)
      STR = ' '
      KERR = .FALSE.
      I = 0
C
      RDUM = ABS(RNUM)
      IF (ILEN .LT.10) THEN
         IF (RDUM.GE.9.995 .OR. RDUM.LT.0.995 .OR. ILEN.LT.5) THEN
            KERR = .TRUE.
            RETURN
         ENDIF
      ENDIF
C
      STR = ' 0.00'
      I = 5
      IF (RDUM .EQ. 0.0) RETURN
C
C     convert RDUM to a value between 1.0 and 10.0
C
      IF (RDUM.GT.BIG .OR. RDUM.LT. SMALL) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IE = 0
   10 CONTINUE
C
C     find a real number 10 < rnum <= 1
C
      IF (RDUM .GE. 10.0) THEN
         IE = IE - 1
         RDUM = RDUM / 10.0
         GO TO 10
C
      ELSEIF (RDUM .LT. 1.0) THEN
         IE = IE + 1
         RDUM = RDUM * 10.0
         GO TO 10
      ENDIF
C
      IF (NINT(100.0*RDUM) .GT. 100.0*RDUM) THEN
         RDUM = RDUM + .005
         IF (RDUM .GE. 10.0) GO TO 10
      ENDIF
      IVAL = INT(RDUM)
      IF (RNUM .LT. 0.0) THEN
         STR(1:) = '- .0'
      ELSE
         STR(1:) = '  .0'
      ENDIF
      CALL CKI2CH (IVAL, STR(2:2), LT, IERR)
      NREM = 100.0*RDUM - 100*IVAL
      IF (NREM .GE. 10) THEN
         CALL CKI2CH (NREM, STR(4:5), LT, IERR)
      ELSE
         CALL CKI2CH (NREM, STR(5:5), LT, IERR)
      ENDIF
C
      IF (ABS(IE) .NE. 0) THEN
         IF (IE .LT. 0) THEN
            STR(6:8) = 'E+0'
         ELSE
            STR(6:8) = 'E-0'
         ENDIF
         I = 9
         IF (ABS(IE) .GE. 10) THEN
            CALL CKI2CH (ABS(IE), STR(8:), LT, IERR)
            IF (ABS(IE) .GT. 99) I = 10
         ELSE
            CALL CKI2CH (ABS(IE), STR(9:9), LT, IERR)
         ENDIF
      ENDIF
      I = CKLSCH(STR)
C
C     end of SUBROUTINE CKR2CH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRAEX (I, RCKWRK, RA)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRAEX (I, RCKWRK, RA)*
C  Get/put the Pre-exponential coefficient of the Ith reaction
C
C  INPUT
C  I         - Integer scalar, reaction index;
C              I > 0 gets RA(I) from RCKWRK
C              I < 0 puts RA(I) into RCKWRK
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  If I < 1:
C  RA        - Real scalar, pre-exponential coefficient for reaction I.
C                 cgs units, mole-cm-sec-K
C
C  OUTPUT
C  If I > 1:
C  RA        - Real scalar, pre-exponential coefficient for reaction I.
C                 cgs units, mole-cm-sec-K
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*)
C
      NI = NcCO + (IABS(I)-1)*(NPAR+1)
      IF (I .GT. 0) THEN
         RA = RCKWRK(NI)
      ELSE
         RCKWRK(NI) = RA
      ENDIF
C
C     end of SUBROUTINE CKRAEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRAT  (RCKWRK, ICKWRK, II, KK, MAXSP, MAXTB, RU,
     1                   PATM, T, C, NSPEC, NU, NUNK, NPAR, PAR, NREV,
     2                   IREV, RPAR, NFAL, IFAL, IFOP, IFLO, KFAL, NFAR,
     3                   FPAR, NLAN, NLAR, ILAN, PLT, NRLT, IRLT, RPLT,
     4                   NTHB, ITHB, NTBS, AIK, NKTB, SMH, RKFT, RKRT,
     5                   RKF, RKR, EQK, CTB, NRNU, IRNU, RNU, NORD,
     6                   IORD, MXORD, KORD, RORD, NEIM, IEIM, IEIMT,
     7                   NJAN, NJAR, IJAN, PJAN, NFT1, NF1R, IFT1,
     8                   PF1, NMOM, IMOM, NXSM, IXSM, NTDE, ITDE, ITDK,
     9                   KTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRAT  (RCKWRK, ICKWRK, II, KK, MAXSP, MAXTB, RU,
C                     PATM, T, C, NSPEC, NU, NUNK, NPAR, PAR, NREV,
C                     IREV, RPAR, NFAL, IFAL, IFOP, IFLO, KFAL, NFAR,
C                     FPAR, NLAN, NLAR, ILAN, PLT, NRLT, IRLT, RPLT,
C                     NTHB, ITHB, NTBS, AIK, NKTB, SMH, RKFT, RKRT,
C                     RKF, RKR, EQK, CTB, NRNU, IRNU, RNU, NORD,
C                     IORD, MXORD, KORD, RORD, NEIM, IEIM, IEIMT,
C                     NJAN, NJAR, IJAN, PJAN, NFT1, NF1R, IFT1,
C                     PF1, NMOM, IMOM, NXSM, IXSM, NTDE, ITDE, ITDK,
C                     KTFL)**
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  II        - Integer scalar, total reaction count.
C  KK        - Integer scalar, total species count.
C  MAXSP     - Integer scalar, maximum number of species allowed in a
C              reaction; in the current formulation MAXSP=12.
C  MAXTB     - Integer scalar, maximum number of third bodies allowed
C              in a reaction.
C  RU        - Real scalar, universal gas constant.
C                 cgs units, 8.314510E7 ergs/(mole*K)
C  PATM      - Real scalar, pressure of one standard atmosphere.
C                 cgs units, 1.01325E6 dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                   cgs units, mole/cm**3
C  NSPEC(*)  - Integer array, count of participant species in reactions,
C              and the flag for reversibility of reactions;
C              dimension at least II, the total reaction count.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  NPAR      - Integer scalar, total number of parameters in the
C              Arrhenius rate expression for reactions;
C              in the current formulation NPAR=4.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  NREV      - Integer scalar, total number of reactions with
C              explicit reverse parameters.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C              IREV(N) is the reaction index for the Nth reaction
C              with explicit reverse parameters.
C  RPAR(*,*) - Real matrix,  reverse Arrhenius rate coefficients for
C              the NREV reactions; dimension at least NPAR for the
C              first, and at least NREV for the second.
C              RPAR(L,N) is the Lth coefficient for the Nth reaction
C              with explicit reverse rate coefficients.
C                 cgs units same as PAR(*,*)
C  NFAL      - Integer scalar, total number of pressure-dependent
C              reactions.
C  IFAL(*)   - Integer array,  reaction indices for the NFAL reactions;
C              dimension at least NFAL.
C              IFAL(N) = I, reaction I is a pressure-dependent reaction.
C  IFOP(*)   - Integer array, formulation type for the NFAL
C              reactions; dimension at least NFAL.
C  IFLO(*)   - Integer array, pressure-dependence option (unimolecular,
C              or chemically activated) for the NFAL reactions;
C              dimension at least NFAL.
C  KFAL(*)   - Integer array, Array of
C  NFAR      - Integer scalar, total number of additional parameters
C              allowed for reaction formulations.
C  FPAR(*,*) - Real matrix of additional parameters for the NFAL
C              reactions; dimension at least NFAR for the first, and
C              at least NFAL for the second.
C              FPAR(L,N) is the Lth parameter for the Nth pressure-
C              dependent reaction.
C  NLAN      - Integer scalar, total number of Landau-Teller reactions.
C  NLAR      - Integer scalar, number of additional parameters required
C              in a Landau-Teller rate expression;
C              in the current formulation NLAR=2.
C  ILAN(*)   - Integer array, reaction indices for the NLAN reactions;
C              dimension at least NLAN.
C              ILAN(N) is the reaction index for the Nth Landau-
C              Teller reaction.
C  PLT(*,*)  - Real matrix, the additional parameters for the NLAN
C              reactions; dimension at least NLAR for the first and at
C              least NLAN for the second.
C              PLAN(L,N) is the Lth parameter for the Nth
C              Landau-Teller reaction.
C  NRLT      - Integer scalar, total number of Landau-Teller reactions
C              with explicit reverse parameters.
C  IRLT(*)   - Integer array, reaction indices for the NRLT reactions;
C              dimension at least NRLT.
C              IRLT(N) is the reaction index for the Nth reaction
C              with Landau-Teller reverse parameters.
C  RPLT(*,*) - Real matrix, the additional rate parameters for the
C              NRLT reactions; dimension at least NLAR for the first
C              and at least NRLT for the second.
C              RPLT(L,N) is the Lth reverse parameter for the Nth
C              Landau-Teller reaction with reverse parameters.
C  NTHB      - Integer scalar, total number of third-body reactions.
C  ITHB(*)   - Integer array, reaction indices for the NTHB reactions;
C              dimension at least NTHB.
C              ITHB(N) is the reaction index for the Nth third-body
C              reaction.
C  NTBS(*)   - Integer array, total enhanced third-body count for a
C              third-body reaction; dimension at least NTHB.
C              NTBS(N) is the number of enhanced third bodies for
C              the Nth third-body reaction.
C  AIK(*,*)  - Real matrix, enhancement factors of third bodies
C              for the NTHB reactions; dimension MAXTB for the first,
C              the maximum number of enhancement factors, and NTHB
C              for the second.
C              AIK(L,N) is the enhancement factor for the Lth
C              enhanced third body in the Nth third-body reqction.
C  NKTB(*,*) - Integer matrix, species indices for the enhanced
C              third bodies in the NTHB reactions; dimension MAXTB
C              for the first and NTHB for the second.
C
C  NRNU      - Integer scalar, total number of reactions with real
C              stoichiometric coefficients.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  CTB(*)   -  Real array, concentration of third bodies for the
C              reactions; dimension at least II, the total reaction
C              count.
C              CTB(I) is the third-body concentration for reaction I.
C  NORD      - Integer scalar, total number of species changed-order
C              reaction.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  MXORD     - Integer scalar, maximum number of species change-orders
C              allowed in a reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  NEIM      - Integer scalar, total number of electron-impact
C              reactions.
C  IEIM(*)   - Integer array, reaction indices for the NEIM reactions;
C              dimension at least NEIM.
C              IEIM(N) is the reaction index for the Nth electron-
C              impact reaction.
C  IEIMT(*)  - Integer array, temperature-dependence indices for the
C              NEIM reactions; dimension at least NEIM.
C              IEIMT(N) is a pointer into a temperature array for
C              the Nth electron-impact reaction.
C  NJAN      - Integer scalar, total number of Janev, Langer, et al.
C              reactions.
C  NJAR      - Integer scalar, the number of additional rate
C              parameters required for Janev reactions.
C  IJAN(*)   - Integer array, reaction indices of the NJAN reactions;
C              dimension at least NJAN.
C              IJAN(N) is the reaction index for the Nth Janev et al.
C              reaction.
C  PJAN(*,*) - Real matrix, rate parameters for the NJAN reactions;
C              dimension at least NJAR for the first and NJAN for the
C              second.
C              PJAN(L,N) is the Lth parameter for the Nth Janev et al.
C              reaction.
C  NFT1      - Integer scalar, total number of fit-type reactions.
C  IFT1(*)   - Integer array, reaction indices for the NFT1 reactions;
C              dimension at least NFT1.
C              IFT1(N) is the reaction index for the Nth fit-type
C              reaction.
C  NF1R      - Integer scalar, number of additional rate parameters
C              required in a fit-type reactions.
C  PF1(*,*)  - Real matrix, the additional rate parameters for the
C              NFT1 reactions; dimension at least NF1R for the first
C              and at least NFT1 for the second.
C              PF1(L,N) is the Lth fit parameter for the Nth fit-type
C              reaction.
C  NMOM      - Integer scalar, total number of electron momentum-
C              transfer reactions.
C  IMOM(*)   - Integer array, reaction indices for the NMOM reactions;
C              dimension at least NMOM.
C              IMOM(N) is the reaction index for the Nth electron
C              momentum-transfer reaction.
C  NXSM      - Integer scalar, total number of ion momemtum-
C              transfer cross-section reactions.
C  IXSM(*)   - Integer array, reaction indices for the NXSM reactions;
C              dimension at least NXSM.
C  NTDE      - Integer scalar, total number of non-thermal-
C              equilibrium reactions.
C  ITDE(*)   - Integer array, reaction indices for the NTDE reactions;
C              dimension at least NTDE.
C              ITDE(N) is the reaction index for the Nth non-
C              thermal-equilibrium reaction.
C  KTFL(*)   - Integer array, indices into the temperature array,
C              for the KK species, as required by the NTDE reactions;
C              dimension at least KK.
C              KTFL(K) is the temperature array index for species K.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT; dimension at least KK.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units ,depend on the reaction
C  EQKC(*)   - Real array, equilibrium constants in concentration
C              units for reactions; dimension at least II.
C              EQKC(I) is the equilibrium constant for reaction I.
C                 cgs units, (mole/cm**3)**some power,
C                 depends on reaction
C  RKF(*)    - Real array
C              dimension at least II.
C  RKR(*)    - Real array
C              dimension at least II.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RCKWRK(*), ICKWRK(*), C(*), NSPEC(*), NU(MAXSP,*),
     1          NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),
     2          ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*),
     3          IFAL(*), IFOP(*), IFLO(*), KFAL(*), FPAR(NFAR,*),
     4          ITHB(*), NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*),
     5          SMH(*), RKFT(*), RKRT(*), RKF(*), RKR(*), EQK(*),
     6          CTB(*), IRNU(*), RNU(MAXSP,*), IORD(*), KORD(MXORD,*),
     7          RORD(MXORD,*), T(*), IEIM(*), IEIMT(*), IJAN(*),
     8          PJAN(NJAR,*), IFT1(*), PF1(NF1R,*), IMOM(*), IXSM(*),
     9          ITDE(*), ITDK(*), KTFL(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
     1             RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU, IEIM,
     2             IEIMT, IJAN, PJAN, IFT1, PF1, ITDE, ITDK, KTFL,
     3             RKFT, RKRT, EQK)
      CALL CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL, FPAR,
     1             ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF, RKR, CTB,
     2             IRNU, RNU, IORD, KORD, RORD, IMOM, IXSM)
C
C     end of SUBROUTINE CKRAT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
     1                   RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU,
     2                   IEIM, IEIMT, IJAN, PJAN, IFT1, PF1, ITDE,
     3                   ITDK, KTFL, RKFT, RKRT, EQK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
C                     RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU,
C                     IEIM, IEIMT, IJAN, PJAN, IFT1, PF1, ITDE,
C                     ITDK, KTFL, RKFT, RKRT, EQK)
C
C  INPUT
C
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  NSPEC(*)  - Integer array, total number of participant species
C              for reactions, and the flag for reversibility of
C              reactions; dimension at least II.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C              IREV(N) is the reaction index for the Nth reaction
C              with explicit reverse parameters.
C  RPAR(*,*) - Real matrix,  reverse Arrhenius rate coefficients for
C              the NREV reactions; dimension at least NPAR for the
C              first, and at least NREV for the second.
C              RPAR(L,N) is the Lth coefficient for the Nth reaction
C              reaction with explicit reverse rate coefficients.
C                 cgs units same as PAR(*,*)
C  ILAN(*)   - Integer array, reaction indices for the NLAN reactions;
C              dimension at least NLAN.
C              ILAN(N) is the reaction index for the Nth Landau-
C              Teller reaction.
C  PLT(*,*)  - Real matrix, the additional parameters for the NLAN
C              reactions; dimension at least NLAR for the first and at
C              least NLAN for the second.
C              PLAN(L,N) is the Lth parameter for the Nth
C              Landau-Teller reaction.
C  IRLT(*)   - Integer array, reaction indices for the NRLT reactions;
C              dimension at least NRLT.
C              IRLT(N) is the reaction index for the Nth reaction
C              with Landau-Teller reverse parameters.
C  RPLT(*,*) - Real matrix, the additional rate parameters for the
C              NRLT reactions; dimension at least NLAR for the first
C              and at least NRLT for the second.
C              RPLT(L,N) is the Lth reverse parameter for the Nth
C              Landau-Teller reaction with reverse parameters.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  IEIM(*)   - Integer array, reaction indices for the NEIM reactions;
C              dimension at least NEIM.
C              IEIM(N) is the reaction index for the Nth electron-
C              impact reaction.
C  IEIMT(*)  - Integer array, temperature-dependence indices for the
C              NEIM reactions; dimension at least NEIM.
C              IEIMT(N) is a pointer into a temperature array for
C              the Nth electron-impact reaction.
C  IJAN(*)   - Integer array, reaction indices of the NJAN reactions;
C              dimension at least NJAN.
C              IJAN(N) is the reaction index for the Nth Janev et al.
C              reaction.
C  PJAN(*,*) - Real matrix, rate parameters for the NJAN reactions;
C              dimension at least NJAR for the first and NJAN for the
C              second.
C              PJAN(L,N) is the Lth parameter for the Nth Janev et al.
C              reaction.
C  IFT1(*)   - Integer array, reaction indices for the NFT1 reactions;
C              dimension at least NFT1.
C              IFT1(N) is the reaction index for the Nth fit-type
C              reaction.
C  PF1(*,*)  - Real matrix, the additional rate parameters for the
C              NFT1 reactions; dimension at least NF1R for the first
C              and at least NFT1 for the second.
C              PF1(L,N) is the Lth fit parameter for the Nth fit-type
C              reaction.
C  ITDE(*)   - Integer array, reaction indices for the NTDE reactions;
C              dimension at least NTDE.
C              ITDE(N) is the reaction index for the Nth non-
C              thermal-equilibrium reaction.
C  ITDK(*)   - Integer array, special species for the NTDE reactions;
C              dimension at least NTDE.
C              ITDK(N) is the indentifying species, K, whose associated
C              temperature (KTFL(K)) should be used in calculating
C              quantities associated with the reaction.
C  KTFL(*)   - Integer array, indices into the temperature array,
C              for the KK species, as required by the NTDE reactions;
C              dimension at least KK.
C              KTFL(K) is the temperature array index for species K.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT; dimension at least KK.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units depend on the reaction
C  EQK(*)    - Real array, equilibrium constants in concentration
C              units for reactions; dimension at least II.
C              EQKC(I) is the equilibrium constant for reaction I.
C                 cgs units (mole/cm**3)**some power,
C                 depends on reaction
C               NOTE: EQK(I) as returned from this routine may be
C                     in error. Use CKEQ instead.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      include 'ckstrt.h'
C     Integer arrays
      DIMENSION ICKWRK(*), NSPEC(NII), NU(MXSP,NII), NUNK(MXSP,NII),
     1          IREV(NREV), ILAN(NLAN), IRLT(NRLT), IRNU(NRNU),
     2          IEIM(NEIM), IEIMT(NEIM), IJAN(NJAN), IFT1(NFT1),
     3          ITDE(NTDE), ITDK(NTDE), KTFL(*),
C     Real arrays
     4          RCKWRK(*), PAR(NPAR+1,NII), RPAR(NPAR+1,NREV),
     5          PLT(NLAR,NII), RPLT(NLAR,NRLT), SMH(NKK), RKFT(NII),
     6          RKRT(NII), RNU(MXSP,NRNU), T(*), EQK(NII),
     7          PJAN(NJAR,NJAN), PF1(NF1R,NFT1)
C
      COMMON /MACH/ SMALL, BIG, EXPARG
      INTEGER CKLKUP
      EXTERNAL CKLKUP, CKSMH
C
C     Find Gibbs/ RT for all species in the mechanism
C     Note: CKSMH takes a vector of temperatures, and the G/RT
C           value for each species will be evaluated at the
C           temperature associated with that species.
C
      CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
      TEMP = T(1)
      ALOGT = LOG(TEMP)
      TINV = 1.0/TEMP
      RU = RCKWRK(NcRU)
      PATM = RCKWRK(NcPA)
      PFAC = PATM / (RU * TEMP)
C
C.....Default calculation of the forward rate constant..................
C     Put it in a loop that vectorizes
C
      DO 30 I = 1, NII
        RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)*TINV)
   30 CONTINUE
C
C.....Default way to calculate the equilibrium constant.................
C     Loop over NII -> largest expense in the routine, so optimize)
C
      DO 50 I = 1, NII
C
C     First decide on whether the equilibrium constant needs to be
C     calculated. -> irreversible rxns don't need their equilibrium
C     constants calculated. REV rxns don't either, but a lookup
C     function on a generic loop over all gas reactions would be too
C     expensive.
C
       IF (NSPEC(I) .GT. 0) THEN
C
C         Initialize the net mole change number and the
C         the net DeltaG value with the contributions from
C         the first reactant and product
C
         NUSUMK = NU(1,I) + NU(7,I)
         SUMSMH = NU(1,I)*SMH(NUNK(1,I)) + NU(7,I)*SMH(NUNK(7,I))
         IF (NUNK(2,I) .NE. 0) THEN
           NUSUMK = NUSUMK + NU(2,I)
           SUMSMH = SUMSMH + NU(2,I)*SMH(NUNK(2,I))
           IF (NUNK(3,I) .NE. 0) THEN
             NUSUMK = NUSUMK + NU(3,I)
             SUMSMH = SUMSMH + NU(3,I)*SMH(NUNK(3,I))
             IF (NUNK(4,I) .NE. 0) THEN
               NUSUMK = NUSUMK + NU(4,I)
               SUMSMH = SUMSMH + NU(4,I)*SMH(NUNK(4,I))
               IF (NUNK(5,I) .NE. 0) THEN
                 NUSUMK = NUSUMK + NU(5,I)
                 SUMSMH = SUMSMH + NU(5,I)*SMH(NUNK(5,I))
                 IF (NUNK(6,I) .NE. 0) THEN
                   NUSUMK = NUSUMK + NU(6,I)
                   SUMSMH = SUMSMH + NU(6,I)*SMH(NUNK(6,I))
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
           NUSUMK = NUSUMK + NU(8,I)
           SUMSMH = SUMSMH + NU(8,I)*SMH(NUNK(8,I))
           IF (NUNK(9,I) .NE. 0) THEN
             NUSUMK = NUSUMK + NU(9,I)
             SUMSMH = SUMSMH + NU(9,I)*SMH(NUNK(9,I))
             IF (NUNK(10,I) .NE. 0) THEN
               NUSUMK = NUSUMK + NU(10,I)
               SUMSMH = SUMSMH + NU(10,I)*SMH(NUNK(10,I))
               IF (NUNK(11,I) .NE. 0) THEN
                 NUSUMK = NUSUMK + NU(11,I)
                 SUMSMH = SUMSMH + NU(11,I)*SMH(NUNK(11,I))
                 IF (NUNK(12,I) .NE. 0) THEN
                   NUSUMK = NUSUMK + NU(12,I)
                   SUMSMH = SUMSMH + NU(12,I)*SMH(NUNK(12,I))
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
C
C        Calculate the concentration equilibrium constant,
C        Protecting against overflow in the exponential
C
         IF (NUSUMK .NE. 0) THEN
           EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**NUSUMK)
         ELSE
           EQK(I) = EXP(MIN(SUMSMH,EXPARG))
         ENDIF
C
C        Calculate the reverse rate constant from the forward
C        rate constant and the equilibrium constant
C
          RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
        ELSE
C
C        Case or irreversible reactions
C
          RKRT(I) = 0.0
        ENDIF
   50 CONTINUE
C
C.....Fix-up look for rxn's with real-stoichiometry.....................
C
C     Must completely redo the calculation of the Gibbs free energy
C     of reaction because the stoichiometric coefficients have been
C     overridden.
C     There is also no need to do this loop if this real stoich
C     reaction is also a REV reaction.
C
      DO 70 N = 1, NRNU
         I = IRNU(N)
         ISREV = CKLKUP(I, IREV, NREV)
         IF (NSPEC(I).GT.0 .AND. ISREV.EQ.0) THEN
           RNUSUM = RNU(1,N) + RNU(7,N)
           SUMSMH = RNU(1,N)*SMH(NUNK(1,I)) + RNU(7,N)*SMH(NUNK(7,I))
           DO 60 L = 2, 6
              IF (NUNK(L,I) .EQ. 0) GO TO 61
                 SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
                 RNUSUM = RNUSUM + RNU(L,N)
   60      CONTINUE
   61      CONTINUE
           DO 62 L = 8, 12
             IF (NUNK(L,I) .EQ. 0) GO TO 63
                SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
                RNUSUM = RNUSUM + RNU(L,N)
   62      CONTINUE
   63      CONTINUE
           IF (RNUSUM .NE. 0.0) THEN
              EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**RNUSUM)
           ELSE
              EQK(I) = EXP(MIN(SUMSMH,EXPARG))
           ENDIF
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         END IF
   70 CONTINUE
C
C.....Rxns with explicitly defined reverse rate constants...............
C
C     Fix-up loop for reactions which have a defined reverse
C     rate constant. We will recalculate the reverse rate
C     constant in this case
C
      DO 90 N = 1, NREV
         I = IREV(N)
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)*TINV)
   90 CONTINUE
C
C.....Landau-Teller reactions...........................................
C
      DO 200 N = 1, NLAN
         I = ILAN(N)
C
C        modify K_f for reaction I
         TFAC = PLT(1,N)/TEMP**(1.0/3.0) + PLT(2,N)/TEMP**(2.0/3.0)
         RKFT(I) = RKFT(I) * EXP(TFAC)
C
C        Bail out for irreversible reactions
         IF (NSPEC(I) .LE. 0) GO TO 200
C
C        Lookup whether this is an Landau-Teller Rxn with an
C        explicitly given reverse reaction rate
C        - Note: If this is such a rxn, need to modify calculation
C                to add more parameters.
C                Even if it isn't such a reaction, we just changed
C                RKFT(I) above, so need to recalculate RKRT(I)
C
         ISRLT = CKLKUP(I, IRLT, NRLT)
         IF (ISRLT .EQ. 0) THEN
C           new K_r = K_eq / new K_f if no explicit reverse LT
C           parameters
            RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ELSE
C           new K_r = old K_r, modified by explicit parameters
            TFAC = RPLT(1,ISRLT)/TEMP**(1.0/3.0)
     1           + RPLT(2,ISRLT)/TEMP**(2.0/3.0)
            RKRT(I) = RKRT(I) * EXP(TFAC)
         ENDIF
  200 CONTINUE
C
C.....Electron-impact reactions.........................................
C
      DO 300 N = 1, NEIM
C
C        No change in rate constant unless temperature other
C        than first temperature is associated with this reaction.
C
         IF (IEIMT(N) .EQ. 1) GO TO 300
C
C        new K_f for reaction I is Arrhenius expression with new
C        temperature. Still check to see if TEMP is different than
C        T(1)
C
         TEMP = T(IEIMT(N))
         IF (TEMP .EQ. T(1)) GO TO 300
         I = IEIM(N)
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I) * LOG(TEMP) - PAR(3,I)/TEMP)
C
C          If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 300
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C             new K_r with explicit parameters and new temperature
            RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV)*LOG(TEMP)
     1              - RPAR(3,ISREV)/TEMP)
         ELSE
C
C            Fix concentration term in equilibrium constant, using the
C            previously calculated Gibbs Free energy part, which is
C            still good.
C
           PFAC2 = PATM / (RU*TEMP)
C          Does reaction have real stoichiometry?
           ISREAL = CKLKUP(I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
             L = ISREAL
             RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1              + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2              + RNU(11,L)+RNU(12,L)
             IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
             NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1              + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
             IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF
C
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
C
  300 CONTINUE
C
C.....Non-thermal-equilibrium, species-temperature-dependent reactions:.
C
      DO 400 N = 1, NTDE
C
C        Get the temperature index for the species which controls
C        which temperature to use in this reaction
C
         KTEMP = KTFL(ITDK(N))
C
C        No change in rate constant unless temperature other
C        than first temperature is associated with this reaction.
C
         IF (KTEMP .EQ. 1) GO TO 400
C
         TEMP = T(KTEMP)
         IF (TEMP .EQ. T(1)) GOTO 400
         I = ITDE(N)
C
C        Calculate a new K_f for reaction I, Arrhenius expression
C        with the new temperature
C
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 400
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C          Fix concentration term in equilibrium constant, using the
C          previously calculated Gibbs Free energy part, which is
C          still good.
C
           PFAC2 = PATM / (RU*TEMP)
C            real coefficients?
           ISREAL = CKLKUP(I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
             L = ISREAL
             RNUSUM= RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1             + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2             + RNU(11,L) +RNU(12,L)
             IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
             NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1              + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
             IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
  400 CONTINUE
C
C.....reactions using fit#1:  k = A * T^B * exp(v1/T+v2/T^2+v3/T^3...)..
C
      DO 500 N = 1, NFT1
         IF (N .EQ. 1) BLOG = LOG(BIG)
         I = IFT1(N)
C
C        Species temperature-dependent reaction?
C
         N2 = CKLKUP(I, ITDE, NTDE)
         IF (N2 .EQ. 0) THEN
           TEMP = T(1)
         ELSE
           TEMP = T(KTFL(ITDK(N2)))
         ENDIF
C
C        Electron-impact reaction?
C
         N2 = CKLKUP(I, IEIM, NEIM)
         IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C        Calculate a new K_f based on fit#1, and possibly a
C        different temperature
C
         SUMJ = 0.0
         DO 470 J = 1, NF1R
            ROOTJ = 1.0/J
            IF (TEMP.GE.BIG**ROOTJ) THEN
               TEMPPW = BIG
            ELSE
               TEMPPW = TEMP**J
            ENDIF
            IF (SUMJ .GE. BLOG) THEN
               SUMJ = BLOG
            ELSE
               SUMJ = SUMJ + PF1(J,N)/TEMPPW
            ENDIF
 470     CONTINUE
         SUMJ = MIN (SUMJ, BLOG)
         RKFT(I) = MIN(BIG, (PAR(1,I) * TEMP**PAR(2,I) * EXP(SUMJ)))
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 500
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C         Fix concentration term in equilibrium constant, using the
C         previously calculated Gibbs Free energy part, which is
C         still good. Only needs to be done, if we have a new
C         temperature
C
           IF (TEMP .NE. T(1)) THEN
C
              PFAC2 = PATM / (RU * TEMP)
C             real coefficients?
              ISREAL = CKLKUP(I, IRNU, NRNU)
              IF (ISREAL .GT. 0) THEN
                 L = ISREAL
                 RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1                  + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2                  + RNU(11,L)+RNU(12,L)
                 IF (RNUSUM .NE. 0.0) EQK(I)=EQK(I)*(PFAC2/PFAC)**RNUSUM
              ELSE
                 NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)
     1                  + NU(6,I)+NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)
     2                  + NU(11,I)+NU(12,I)
                 IF (NUSUMK .NE. 0) EQK(I)=EQK(I)* (PFAC2/PFAC)**NUSUMK
              ENDIF
           ENDIF
C
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
  500 CONTINUE
C
C.....jannev, langer, evans & post - type reactions.....................
C
      DO 600 N = 1, NJAN
         I = IJAN(N)
C
C        Specify the temperature to be used below by looking up
C        whether this rxn has a special temperature variable
C        associated with it. If it doesn't, use the default first
C        temperature variable.
C
         N2 = CKLKUP(I, ITDE, NTDE)
         IF (N2 .EQ. 0) THEN
           TEMP = T(1)
         ELSE
           TEMP = T(KTFL(ITDK(N2)))
         ENDIF
C
C        Lookup up whether this rxn is an electron impact rxn
C        If it is, override the temperature variable determined
C        in the previous loop.
C
         N2 = CKLKUP (I, IEIM, NEIM)
         IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C        Re-evaluate Arrhenius K_f, possibly different temperature,
C        then modify for jannev expression
C
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C        convert E- temperature to eV's
         TEV = TEMP / 11595.0
         SUMJ = 0.0
         DO 530 J = 1, NJAR
            SUMJ = SUMJ + PJAN(J,N) * (LOG(TEV))**(J-1)
  530    CONTINUE
         RKFT(I) =  MIN(BIG, RKFT(I) * EXP(SUMJ))
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 600
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C          Fix concentration term in equilibrium constant, using the
C          previously calculated Gibbs Free energy part, which is
C          still good. Only needs to be done, if we have a new
C          temperature
C
           IF (TEMP .NE. T(1)) THEN
C
             PFAC2 = PATM / (RU * TEMP)
C            real coefficients?
             ISREAL = CKLKUP(I, IRNU, NRNU)
             IF (ISREAL .GT. 0) THEN
               L = ISREAL
               RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1                + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2                + RNU(11,L)+RNU(12,L)
               IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
             ELSE
               NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1                + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)
     2                + NU(12,I)
               IF (NUSUMK.NE.0) EQK(I) = EQK(I) * (PFAC2/PFAC) ** NUSUMK
             ENDIF
           ENDIF
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
 600  CONTINUE
C
C     end of SUBROUTINE CKRATT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL,
     1                   FPAR, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,
     2                   RKR, CTB, IRNU, RNU, IORD, KORD, RORD, IMOM,
     3                   IXSM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL,
C                     FPAR, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,
C                     RKR, CTB, IRNU, RNU, IORD, KORD, RORD, IMOM,
C                     IXSM)
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  NSPEC(*)  - Integer array, total number of participant species
C              for reactions, and the flag for reversibility of
C              reactions; dimension at least II.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  IFAL(*)   - Integer array, reaction indices for the NFAL reactions;
C              dimension at least NFAL.
C              IFAL(N) = I, reaction I is a pressure-dependent reaction.
C  IFOP(*)   - Integer array, formulation type for the NFAL reactions;
C              dimension at least NFAL.
C  IFLO(*)   - Integer array, pressure-dependence type for the NFAL
C              reactions (unimolecular vs chemically activated);
C              dimension at least NFAL.
C  KFAL(*)   - Integer array, array of 3rd-body species indices for the
C              NFAL reactions; 0 to use total concentration of mixture,
C              K for concentration of species K.
C  FPAR(*,*) - Real matrix of parameters for the NFAL reactions;
C              dimension at least NFAR for the first, and at least NFAL
C              for the second.
C              FPAR(L,N) is the Lth parameter for the Nth pressure-
C              dependent reaction.
C  ITHB(*)   - Integer array, reaction indices for the NTHB reactions;
C              dimension at least NTHB.
C              ITHB(N) is the reaction index for the Nth third-body
C              reaction.
C  NTBS(*)   - Integer array, total number of enhanced third-bodies
C              in a third-body reaction; dimension at least NTHB.
C              NTBS(N) is the total enhanced third-body count for
C              the Nth third-body reaction.
C  AIK(*,*)  - Real matrix, enhancement factors of third bodies
C              for the NTHB reactions; dimension MAXTB for the first,
C              the maximum number of enhancement factors, and NTHB
C              for the second.
C              AIK(L,N) is the enhancement factor for the Lth
C              enhanced third body in the Nth third-body reqction.
C  NKTB(*,*) - Integer matrix, species indices for the enhanced
C              third bodies in the NTHB reactions; dimension MAXTB
C              for the first and NTHB for the second.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  CTB(*)   -  Real array, concentration of third bodies for the
C              reactions;
C              dimension at least II, the total reaction count.
C              CTB(I) is the third-body concentration for reaction I.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  IMOM(*)   - Integer array, reaction indices for the NMOM reactions;
C              dimension at least NMOM.
C              IMOM(N) is the reaction index for the Nth electron
C              momentum-transfer reaction.
C  IXSM(*)   - Integer array, reaction indices for the NXSM reactions;
C              dimension at least NXSM.
C  OUTPUT
C  RKF(*)-
C              dimension at least II, the total reaction count.
C  RKR(*)-
C              dimension at least II, the total reaction count.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      include 'ckstrt.h'
C     Integer arrays
      DIMENSION NU(MXSP,NII), NUNK(MXSP,NII),
     1          IFAL(NFAL), IFOP(NFAL), IFLO(NFAL), KFAL(NFAL),
     2          ITHB(NTHB), NTBS(NTHB), NKTB(MXTB,NTHB),
     3          IRNU(NRNU), IORD(NORD), KORD(MXORD, NORD),
     4          IMOM(NMOM), IXSM(NXSM),
C     Real arrays
     5          T(*), C(NKK), PAR(NPAR+1,NII), FPAR(NFAR,NII),
     6          AIK(MXTB,NTHB), RKFT(NII), RKRT(NII), RKF(NII),
     7          RKR(NII), CTB(NII), RNU(MXSP,NRNU), RORD(MXORD,NORD)
      COMMON /MACH/ SMALL,BIG,EXPARG
C
C     third-body reactions
C
      IF (NTHB .GT. 0) THEN
         CTOT = 0.0
         DO 10 K = 1, NKK
            CTOT = CTOT + C(K)
   10    CONTINUE
         DO 80 N = 1, NTHB
            I = ITHB(N)
            CTB(I) = CTOT
            DO 70 L = 1, NTBS(N)
               CTB(I) = CTB(I) + (AIK(L,N)-1.0)*C(NKTB(L,N))
   70       CONTINUE
   80    CONTINUE
      ENDIF
C
C     If pressure correction:
C
      IF (NFAL .GT. 0) THEN
         TK = T(1)
         ALOGT = LOG(TK)
C
         DO 90 N = 1, NFAL
C
            I = IFAL(N)
            IF (KFAL(N) .GT. 0) THEN
C              third-body species named
               CONC = C(KFAL(N))
            ELSE
C              no species named
               CONC = CTB(I)
               CTB(I) = 1.0
            ENDIF
C
            IF (IFLO(N) .EQ. 0) THEN
C              unimolecular reaction (RKFT is K_inf)
               PR = (FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/TK))
     1              * CONC / RKFT(I)
               PCOR = PR / (1.0 + PR)
            ELSE
C              chemically activated reaction (RKFT is K_zero)
               PR = RKFT(I) * CONC /
     1              (FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/TK))
               PCOR = 1.0 / (1.0 + PR)
            ENDIF
C
            IF (IFOP(N) .GT. 1) THEN
C
C              Fcentering (F^x)
C
               PRLOG = LOG10(MAX(PR,SMALL))
C
               IF (IFOP(N) .EQ. 2) THEN
C
C              8-PARAMETER SRI FORM
C
                  XP = 1.0/(1.0 + PRLOG**2)
                  FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/TK)
     1                   + EXP(-TK/FPAR(6,N))) **XP)
     2                  * FPAR(7,N) * TK**FPAR(8,N)
C
               ELSE
C
C              6-PARAMETER TROE FORM
C
                  FCENT = (1.0-FPAR(4,N)) * EXP(-TK/FPAR(5,N))
     1                  +       FPAR(4,N) * EXP(-TK/FPAR(6,N))
C
C              7-PARAMETER TROE FORM
C
                  IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/TK)
C
                  FCLOG = LOG10(MAX(FCENT,SMALL))
                  XN    = 0.75 - 1.27*FCLOG
                  CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
                  FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
                  FC = 10.0**FLOG
               ENDIF
               PCOR = FC * PCOR
            ENDIF
C
            RKFT(I) = RKFT(I) * PCOR
            RKRT(I) = RKRT(I) * PCOR
   90    CONTINUE
      ENDIF
C
C       Modify rate constants for third body collision reactions
C       now that we have taken care of the unimolecular reaction
C       special case
C
      DO 100 N = 1, NTHB
        I = ITHB(N)
        RKFT(I) = RKFT(I) * CTB(I)
        RKRT(I) = RKRT(I) * CTB(I)
  100 CONTINUE
C
C     Calculate the rate of each reaction from the rate
C     constant and product of the concentration raized
C     to their stoichiometric powers.
C      - NU(1,I) check is used, because real stoich. coefficient
C        rxns set NU(*,I)=0. They are fixed up in a later
C        loop so there is no need to initialize RKF and RKR.
C
      DO 150 I = 1, NII
         IF (NU(1,I) .EQ. 0) GO TO 150
C        4th parameter may be perturbation factor
         RKF(I) = RKFT(I)*C(NUNK(1,I))**IABS(NU(1,I)) * PAR(4,I)
         RKR(I) = RKRT(I)*C(NUNK(7,I))**NU(7,I) * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I)= RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
            IF (NUNK(3,I).NE.0) THEN
               RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
               IF (NUNK(4,I).NE.0) THEN
                  RKF(I) = RKF(I) * C(NUNK(4,I))**IABS(NU(4,I))
                  IF (NUNK(5,I).NE.0) THEN
                     RKF(I) = RKF(I) * C(NUNK(5,I))**IABS(NU(5,I))
                     IF (NUNK(6,I) .NE. 0)
     1                  RKF(I) = RKF(I) * C(NUNK(6,I))**IABS(NU(6,I))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I)*C(NUNK(8,I))**NU(8,I)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I)*C(NUNK(9,I))**NU(9,I)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I)*C(NUNK(10,I))**NU(10,I)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I)*C(NUNK(11,I))**NU(11,I)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I)*C(NUNK(12,I))**NU(12,I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  150 CONTINUE
C
      DO 160 N = 1, NRNU
C        For real coefficients, must ensure non-negative concentrations
         I = IRNU(N)
         RKF(I) = RKFT(I) * MAX(ZERO,C(NUNK(1,I))) ** ABS(RNU(1,N)) 
     1                    * PAR(4,I)
         RKR(I) = RKRT(I) * MAX(ZERO,C(NUNK(7,I))) ** RNU(7,N) 
     1            * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(2,I))) ** ABS(RNU(2,N))
            IF (NUNK(3,I) .NE. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(3,I))) ** ABS(RNU(3,N))
               IF (NUNK(4,I) .NE. 0) THEN
                  RKF(I) = RKF(I) * 
     1                     MAX(ZERO,C(NUNK(4,I))) ** ABS(RNU(4,N))
                  IF (NUNK(5,I) .NE. 0) THEN
                     RKF(I) = RKF(I) * 
     1                        MAX(ZERO,C(NUNK(5,I))) ** ABS(RNU(5,N))
                     IF (NUNK(6,I) .NE. 0) THEN
                       RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(6,I))) 
     1                          ** ABS(RNU(6,N))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(8,I))) ** RNU(8,N)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(9,I))) ** RNU(9,N)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(10,I))) ** RNU(10,N)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I) * 
     1                        MAX(ZERO,C(NUNK(11,I))) ** RNU(11,N)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I) * 
     1                           MAX(ZERO,C(NUNK(12,I))) ** RNU(12,N)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  160 CONTINUE
C
C       Handle reactions with different rxn orders than their
C       stoichiometry
C
      DO 200 N = 1, NORD
         I = IORD(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = RKRT(I) * PAR(4,I)
         DO 190 L = 1, MXORD
            NK = KORD(L,N)
            IF (NK .LT. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,C(-NK)) ** RORD(L,N)
            ELSEIF (NK .GT. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,C(NK)) ** RORD(L,N)
            ENDIF
  190    CONTINUE
  200 CONTINUE
C
C RETURN PURE COLLISION FREQUENCIES FOR E MOMENTUM-TRANSFER (NOT RATES)
C
      DO 300 N = 1, NMOM
         I = IMOM(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = 0.0
 300  CONTINUE
C
C RETURN PURE COLLISION X-SECTIONS FOR ION MOMENTUM-TRANSFER (NOT RATES)
C
      DO 350 N = 1, NXSM
         I = IXSM(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = 0.0
 350  CONTINUE
C
C     end of SUBROUTINE CKRATX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRCXP  (P, T, X, ICKWRK, RCKWRK, RCFT, RCRT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRCXP  (P, T, X, ICKWRK, RCKWRK, RCFT, RCRT)
C     Returns the forward and reverse rate constants for all reactions
C     given pressure, temperature and mole fractions;
C     see Eqs. (51) and (58).  Note this subroutine will calculate
C     a value for the reverse rate constant irrespective of
C     whether the reaction was deemed reversible in the interpretor
C     file.  Also note that the concentration of third bodies
C     for third body reactions is included in the returned rate
C     constant.  The units for the rate constant will depend
C     on the number of reactants.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T(*)   - Temperature.
C                   cgs units - K
C                   Data type - real array
C     X      - Mole fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension X(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     RCFT      - Rate constant for the forward reaction.
C                   cgs units - mole/(cm**3*sec) *
C                           (cm**3/mole)**(sum of forward stoich. coeff)
C                   Data type - real array
C                   Dimension RCFT(*) at least II, the total number
C                   of reactions.
C     RCRT      - Rate constant for the forward reaction.
C                   cgs units - mole/(cm**3*sec) *
C                           (cm**3/mole)**(sum of reverse stoich. coeff)
C                   Data type - real array
C                   Dimension RCRT(*) at least II, the total number
C                   of reactions.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), X(NKK), RCFT(NII), RCRT(NII), T(*)
C      CHARACTER*64 VERSN
C      DATA VERSN/'@(#) ckrcxp.f 1.5 12/22/95'/
C
      EXTERNAL CKRATT, CKRTCN
C
C  Replace the normal NSPEC vector, located at ICKWRK(IcNS),
C  with an all positive one. Set the NREAC vector to be negative
C  to flag this. The NREAC vector is never used!
C    This way CKRATT will calculate the reverse rate constant for
C  reactions that are defined as irreversible in the mechanism
C
      NIM1 = NII - 1
      DO 10 I = 0, NIM1
        IF (ICKWRK(IcNS + I) .LT. 0) THEN
           ICKWRK(IcNS + I) = - ICKWRK(IcNS + I)
           ICKWRK(IcNR + I) = - ICKWRK(IcNR + I)
        ENDIF
 10   CONTINUE
C
C  Call the normal temperature dependent part of the rate constant
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     1             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     2             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     3             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     4             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     5             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     6             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     7             ICKWRK(IcTK), ICKWRK(IcKTF), RCFT, RCRT,
     8             RCKWRK(NcI1))
C
C  Flip the normal NSPEC vector, located at ICKWRK(IcNS),
C  back to its previous state
C
      DO 20 I = 0, NIM1
        IF (ICKWRK(IcNR + I) .LT. 0) THEN
           ICKWRK(IcNS + I) = - ICKWRK(IcNS + I)
           ICKWRK(IcNR + I) = - ICKWRK(IcNR + I)
        ENDIF
 20   CONTINUE
C
C  Calculate the concentrations of all species
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
C  Fix up the rate constants due to fall-off and third-body
C  effects
C
      CALL CKRTCN(RCKWRK, ICKWRK, NII, NKK, MXTB, RCKWRK(NcRU),
     1            RCKWRK(NcPA), T, RCKWRK(NcK1), NPAR+1, RCKWRK(NcCO),
     2            NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF),
     3            NFAR, RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),
     4            RCKWRK(NcKT), ICKWRK(IcKT), RCFT, RCRT, RCKWRK(NcI4))
C
C     end of SUBROUTINE CKRCXP
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKRTCN (RCKWRK, ICKWRK, II, KK, MAXTB, RU, PATM, T, C,
     1                   NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR,
     2                   NTHB, ITHB, NTBS, AIK, NKTB, RCFT, RCRT, CTB)
C
C  START PROLOGUE
C
C   This subroutine modifies the forward and reverse rate constants
C obtained from CKRATT to account for those parts of CKRATX that
C don't involve multiplications of the concentrations of reactants
C or products. This specifically includes third body effects and
C and fall-off effects.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MACH/ SMALL, BIG, EXPARG
      DIMENSION ICKWRK(*), RCKWRK(*), C(KK), PAR(NPAR,II),
     1          FPAR(NFAR,NFAL), AIK(MAXTB,NTHB), RCFT(II), RCRT(II),
     2          CTB(II), IFAL(NFAL), IFOP(NFAL), KFAL(NFAL),
     3          NTBS(NTHB), ITHB(NTHB), NKTB(MAXTB,NTHB)
C
C         PRELIMINARIES
C
      ALOGT = LOG(T)
      AINVT = 1.0/T
C
C  Third-body reactions -
C       Note that fall-off reactions with no special fall-off species
C       are alwys declared as third body reactions. That way CTB(I)
C       will be defined correctly in the fall-off reaction section
C       below.
C
      IF (NTHB .GT. 0) THEN
C
C         Find the total concentration
C
         CTOT = 0.0
         DO 30 K = 1, KK
           CTOT = CTOT + C(K)
   30    CONTINUE
         DO 80 N = 1, NTHB
           CTB(ITHB(N)) = CTOT
           DO 40 L = 1, NTBS(N)
             CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.0)*C(NKTB(L,N))
   40      CONTINUE
   80   CONTINUE
      ENDIF
C
C  Corrections for fall-off reactions
C
      DO 90 N = 1, NFAL
C
C         Store the reaction number
C
        I = IFAL(N)
C
C         Store the special species, if there is one
C
        K = KFAL(N)
C
C         Calculate the low pressure reaction rate
C
        RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)*AINVT)
C
C         Find the concentration of the third body
C
        IF (K .EQ. 0) THEN
           PR = RKLOW * CTB(I) / RCFT(I)
        ELSE
           PR = RKLOW * C(K)   / RCFT(I)
        ENDIF
C
C         Unitize the effective concentration for fall-off
C         reactions. We don't want to multiply the rxn rate
C         constant by an extra concentration factor
C
        CTB(I) = 1.0
C
C         This is the Lindemann form , i.e., IFOP(N) = 1
C
        PCOR = PR / (1.0 + PR)
C
        IF (IFOP(N) .GT. 1) THEN
            PRLOG = LOG10(MAX(PR,SMALL))
            IF (IFOP(N) .EQ. 2) THEN
C
C              SRI FORM
C
               XP = 1.0/(1.0 + PRLOG**2)
               FC = ( ( FPAR(4,N)*EXP(-FPAR(5,N)*AINVT)
     $            + EXP(-T/FPAR(6,N)) )**XP ) * FPAR(7,N) * T**FPAR(8,N)
C
            ELSE
C
C              6-PARAMETER TROE FORM
C
               FCENT = (1.0-FPAR(4,N)) * EXP(-T/FPAR(5,N))
     $               + FPAR(4,N) * EXP(-T/FPAR(6,N))
C
C              7-PARAMETER TROE FORM
C
               IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)*AINVT)
C
               FCLOG = LOG10(MAX(FCENT,SMALL))
               XN    = 0.75 - 1.27*FCLOG
               CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
               FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
               FC = 10.0**FLOG
            ENDIF
            PCOR = PCOR * FC
         ENDIF
C
C           Correct both the forward and reverse rate constant
C           by the k/k_inf, PCOR
C
         RCFT(I) = RCFT(I) * PCOR
         RCRT(I) = RCRT(I) * PCOR
C
   90 CONTINUE
C
C Multiply the rate constant by the third body factor
C
      DO 140 N = 1, NTHB
         I = ITHB(N)
         RCFT(I) = RCFT(I)*CTB(I)
         RCRT(I) = RCRT(I)*CTB(I)
  140 CONTINUE
C
C Multiply the rate constant by the perturbation factor.
C
      DO 150 I = 1, II
         RCFT(I) = RCFT(I)*PAR(4,I)
         RCRT(I) = RCRT(I)*PAR(4,I)
  150 CONTINUE
C
C     end of SUBROUTINE CKRTCN
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKRDEX (I, RCKWRK, RD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRDEX (I, RCKWRK, RD)*
C  Get/put the perturbation factor of the Ith reaction
C
C  INPUT
C  I         - Integer scalar, reaction index;
C              I > 0 gets RD(I) from RCKWRK
C              I < 0 puts RD(I) into RCKWRK
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  If I < 1:
C  RD        - Real scalar, perturbation factor for reaction I;
C              cgs units, mole-cm-sec-K.
C
C  OUTPUT
C  If I > 1:
C  RD        - Real scalar, perturbation factor for reaction I;
C              cgs units, mole-cm-sec-K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ONE = 1.0)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*)
C
      NI = NcCO + (IABS(I)-1)*(NPAR+1) + NPAR
      IF (I .GT. 0) THEN
         RD = RCKWRK(NI)
      ELSE
C
C          Assign the perturbation factor
C
         RCKWRK(NI) = RD
      ENDIF
C
C     end of SUBROUTINE CKRDEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHEX (K, RCKWRK, A6)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHEX (K, RCKWRK, A6)
C
C  Returns an array of the sixth thermodynamic polynomial
C  coefficients for a species, or changes their value,
C  depending on the sign of K.
C
C  INPUT
C  K         - Integer scalar, species index;
C              K>0 gets A6(*) from RCKWRK,
C              K<0 puts A6(*) into RCKWRK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  A6(*)     - Real array, the 6th thermodynamic polynomial
C              coefficients for species K, over the number
C              of fit temperature ranges; dimension at least (MXTP-1),
C              where MXTP is the maximum number of temperatures used
C              to divide the thermodynamic fits.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*), A6(*)
C
      DO 100 L = 1, MXTP-1
         NA6 = NCAA + (L-1)*NCP2 + (IABS(K)-1)*NCP2T + NCP
         IF (K .GT. 0) THEN
            A6(L) = RCKWRK(NA6)
         ELSE
            RCKWRK(NA6) = A6(L)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE CKRHEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOC (P, T, C, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOC (P, T, C, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and molar concentrations;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*)
C
      RHO  = 0.0
      DO 100 K = 1, NKK
         RHO = RHO + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKRHOC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and mole fractions;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      SUMT = 0.0
      SUMW = 0.0
      DO 100 K = 1, NKK
         SUMW = SUMW + X(K)*RCKWRK(NcWT + K - 1)
         SUMT = SUMT + X(K)*T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
C
      RHO = P * SUMW / (SUMT * RCKWRK(NcRU))
C
C     end of SUBROUTINE CKRHOX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and mass fractions;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      SUM = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUM = SUM + Y(K+1) * T(ICKWRK(IcKTF + K)) / RCKWRK(NcWT + K)
150   CONTINUE
      RHO = P / (SUM * RCKWRK(NcRU))
C
C     end of SUBROUTINE CKRHOY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
C  Returns universal gas constants and the pressure of one standard
C  atmosphere
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RU        - Real scalar, universal gas constant.
C                 cgs units, 8.314510E7 ergs/(mole*K)
C  RUC       - Real scalar, universal gas constant used only in
C              conjuction with activation energy.
C                 preferred units, RU / 4.184 cal/(mole*K)
C  PA        - Real scalar, pressure of one standard atmosphere.
C                 cgs units, 1.01325E6 dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      RU  = RCKWRK(NcRU)
      RUC = RCKWRK(NcRC)
      PA  = RCKWRK(NcPA)
C
C     end of SUBROUTINE CKRP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C  Writes to a binary file information about a Chemkin linkfile,
C  pointers for the Chemkin Library, and Chemkin work arrays.
C
C  INPUT
C  LOUT      - Integer scalar, formatted output file unit number.
C  LSAVE     - Integer scalar, binary output file unit number.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
C
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER*(*) CCKWRK(*)
      CHARACTER*16 FILVER, PRVERS, PREC
      LOGICAL KERR
C
      NPOINT = 85
      WRITE (LSAVE, ERR=999)
     *                FILVER, PREC, LENI, LENR, LENC,
C
C     include file for CHEMKIN-III cklib.f, dated: March 1, 1966
C
C     Integer constants
C
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     4   MXORD, KEL, NKKI,
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
C
C     END include file for cklib.f
C
C
      WRITE (LSAVE, ERR=999) (ICKWRK(L), L = 1, LENI)
      WRITE (LSAVE, ERR=999) (RCKWRK(L), L = 1, LENR)
      WRITE (LSAVE, ERR=999) (CCKWRK(L), L = 1, LENC)
      RETURN
C
  999 CONTINUE
      WRITE (LOUT,*)' Error writing Chemkin binary file information...'
      KERR = .TRUE.
C
C     end of SUBROUTINE CKSAVE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSBML (P, T, X, ICKWRK, RCKWRK, SBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSBML (P, T, X, ICKWRK, RCKWRK, SBML)*
C  Returns the mean entropy of the mixture in molar units given
C  pressure, temperature(s) and mole fractions; see Eq. (42).
C
C  INPUT
C  P         - Real scalar, pressure.
C                cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SBML      - Real scalar, mean entropy.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))
      SBML = 0.0
      DO 100 K = 1, NKK
         SBML = SBML + X(K) * ( RCKWRK(NcK1 + K - 1) -
     1          RCKWRK(NcRU)*LOG(MAX(X(K),SMALL)) - RLNP )
  100 CONTINUE
C
C     end of SUBROUTINE CKSBML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSBMS (P, T, Y, ICKWRK, RCKWRK, SBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSBMS (P, T, Y, ICKWRK, RCKWRK, SBMS)*
C  Returns the mean entropy of the mixture in mass units given pressure,
C  temperature(s) and mass fractions; see Eq. (43).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SBMS      - Real scalar, mean entropy.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKMMWY(Y, ICKWRK, RCKWRK, WTM)
      RLNP = RCKWRK(NcRU) * LOG (P / RCKWRK(NcPA))
C
      SUM = 0.0
      RU = RCKWRK(NcRU)
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         SUM = SUM + RCKWRK(NcK2 + K) *
     1             ( RCKWRK(NcK1 + K) - RU *
     3               LOG(MAX(RCKWRK(NcK2 + K),SMALL)) - RLNP)
  100 CONTINUE
      SBMS = SUM / WTM
C
C     end of SUBROUTINE CKSBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)*
C  Returns the array of entropies minus enthalpies for species.
C  It is normally not called directly by the user.
C
C  INPUT
C  T(*)      - Real array, temepratures;
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for species,
C              SMH(K) = S(K)/R - H(K)/RT;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SMH(*)
      SAVE TN1, TNLOG, TNHALF, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TNHALF = TN1 / 2
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 6
         TN3 = TN3 / 12
         TN4 = TN4 / 20
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TKHALF = TK1 / 2
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TKHALF = TNHALF
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMH(K+1) = RCKWRK(NA1) * (TKLOG - 1.0)
     1            + RCKWRK(NA1+1) * TKHALF + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            - RCKWRK(NA1 + NCP1 - 1) / TK1
     4            + RCKWRK(NA1 + NCP2 - 1)
 250  CONTINUE
C
C     end of SUBROUTINE CKSMH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSML  (T, ICKWRK, RCKWRK, SML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSML  (T, ICKWRK, RCKWRK, SML)
C  Returns the standard state entropies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SML(*)    - Real array, standard state entropies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SML(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SML(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            + RCKWRK(NA1 + NCP2 - 1))
250   CONTINUE
C
C     end of SUBROUTINE CKSML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSMS  (T, ICKWRK, RCKWRK, SMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMS  (T, ICKWRK, RCKWRK, SMS)
C  Returns the standard state entropies in mass units; see Eq. (28).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SMS(*)    - Real array, standard state entropies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SMS(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMS(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     1            + RCKWRK(NA1 + NCP2 - 1) ) / RCKWRK(NcWT + K)
  250 CONTINUE
C
C     end of SUBROUTINE CKSMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,
     1                   RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,
C                     RVAL, KERR)
C  Search a character string, LINE, for (1) a character substring which
C  may also appear in an array of character substrings KRAY, and
C  (2) some number of character substrings representing numbers.
C  In the case of (1), if the character substring appears in KRAY,
C  KNUM is its index position.
C  In the case of (2), the character substrings are converted to
C  NVAL real numbers and stored in RVAL, until NEXP are converted.
C
C  This allows format-free input of combined alpha-numeric data.
C  For example, the subroutine might be called to find a Chemkin
C  species index and convert the other substrings to real values:
C
C     input:  LINE    = "N2  1.2"
C             NEXP    = 1, the number of values expected
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
C             NN      = 9, the number of entries in KRAY(*)
C     output: KNUM    = 3, the index number of the substring in
C                       KRAY(*) which corresponds to the first
C                       substring in LINE
C             NVAL    = 1, the number of values found in LINE
C                       following the first substring
C             RVAL(*) = 1.200E+00, the substring converted to a number
C             KERR    = .FALSE.
C  INPUT
C  LINE      - Character string; length depends on calling routine.
C  NEXP      - Integer scalar, number of values to be found in LINE.
C              If NEXP < 0, then IABS(NEXP) values are expected, but
C              it is not an error condition if less values are found.
C  LOUT      - Integer scalar, formatted output file unit.
C  KRAY(*)   - Character string array.
C  NN        - Integer scalar, total number of character strings
C              in KRAY.
C
C  OUTPUT
C  KNUM      - Integer scalar, index of character string in KRAY
C              which corresponds to the first substring in LINE.
C  NVAL      - Integer scalar, count of real values found in LINE.
C  RVAL(*)   - Real array, real values found in LINE; dimension at least
C              NEXP.
C  KERR      - Logical, syntax or dimensioning error flag;
C              corresponding string not found, or total of
C              values found is not the number of values expected,
C              will result in KERR = .TRUE.
C
C  END PROLOGUE
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), KRAY(*)*(*), ISTR*80
      DIMENSION RVAL(*)
      LOGICAL KERR, IERR
      INTEGER CKFRCH, CKLSCH, CKSLEN
      EXTERNAL CKFRCH, CKLSCH, CKSLEN
C
      NVAL = 0
      KERR = .FALSE.
      ILEN = CKSLEN(LINE)
      IF (ILEN .LE. 0) RETURN
C
      I1 = CKFRCH(LINE(1:ILEN))
      I3 = INDEX(LINE(I1:ILEN),' ')
      IF (I3 .EQ. 0) I3 = ILEN - I1 + 1
      I2 = I1 + I3
      ISTR = ' '
      ISTR = LINE(I1:I2-1)
C
      CALL CKCOMP (ISTR, KRAY, NN, KNUM)
      IF (KNUM.EQ.0) THEN
         LT = MAX (CKLSCH(ISTR), 1)
         WRITE (LOUT,'(A)')
     1   ' Error in CKSNUM...'//ISTR(1:LT)//' not found...'
         KERR = .TRUE.
      ENDIF
C
      ISTR = ' '
      ISTR = LINE(I2:ILEN)
      IF (NEXP .NE. 0)
     1      CALL CKXNUM (ISTR, NEXP, LOUT, NVAL, RVAL, IERR)
C
C     end of SUBROUTINE CKSNUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSOR  (T, ICKWRK, RCKWRK, SOR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSOR  (T, ICKWRK, RCKWRK, SOR)
C  Returns the nondimensional entropies;  see Eq. (21).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SOR(*)    - Real array, nondimensional entropies for species;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SOR(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SOR(K+1) = RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            + RCKWRK(NA1 + NCP2 - 1)
250   CONTINUE
C
C     end of SUBROUTINE CKSOR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
C  Returns an array of substrings in a character string with blanks
C  or tabs as delimiters
C
C  INPUT
C  LINE      - Character string; length determined by calling routine.
C  LOUT      - Integer scalar, formatted output file unit.
C  NDIM      - Integer scalar, dimension of a character string array.
C
C  OUTPUT
C  SUB(*)    - Character string array, the character substrings of
C              LINE; dimension SUB at least NDIM.
C  NFOUND    - Integer scalar, count of substrings found in LINE.
C  KERR      - Logical, error flag; dimensioning errors will result in
C              KERR = .TRUE.
C
C  END PROLOGUE
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER SUB(*)*(*), LINE*(*)
      LOGICAL KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
      NFOUND = 0
      ILEN = LEN(SUB(1))
C
      IEND = 0
      KERR = .FALSE.
   25 CONTINUE
C
      ISTART = IEND + 1
      DO 100 L = ISTART, CKSLEN(LINE)
C
         IF (LINE(L:L) .NE. ' ' .AND. LINE(L:L).NE.CHAR(9)) THEN
            IEND   = INDEX(LINE(L:), ' ')
            IF (IEND .EQ. 0) THEN
               IEND = CKSLEN(LINE)
            ELSE
               IEND = L + IEND - 1
            ENDIF
            IF (IEND-L+1 .GT. ILEN) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...substring too long'
               KERR = .TRUE.
            ELSEIF (NFOUND+1 .GT. NDIM) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...NDIM too small'
               KERR = .TRUE.
            ELSE
               NFOUND = NFOUND + 1
               SUB(NFOUND) = LINE(L:IEND)
            ENDIF
            GO TO 25
         ENDIF
C
  100 CONTINUE
C
C     end of SUBROUTINE CKSUBS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYME (CCKWRK, LOUT, ENAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYME (CCKWRK, LOUT, ENAME, KERR)*
C  Returns the character strings of element names.
C
C  INPUT
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C  LOUT      - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  ENAME(*)  - Character string array, element names; dimension at
C              least MM, the total element count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      CHARACTER CCKWRK(*)*(*), ENAME(*)*(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(ENAME(1))
      DO 150 M = 1, NMM
         LT = CKLSCH(CCKWRK(IcMM+M-1))
         ENAME(M) = ' '
         IF (LT .LE. ILEN) THEN
            ENAME(M) = CCKWRK(IcMM+M-1)
         ELSE
            WRITE (LOUT,'(A)')
     1      ' Error in CKSYME...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
C
C     end of SUBROUTINE CKSYME
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT, ISTR,
     1                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYMR (I, ICKWRK, RCKWRK, CCKWRK, LT, ISTR, KERR)*
C  Returns a character string which describes the Ith reaction,
C  and the effective length of the character string.
C
C  INPUT
C  I         - Integer scalar, reaction index.
C  LOUT      - Integer scalar, formatted output file unit.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C
C  OUTPUT
C  ISTR      - Character string, description of reaction I.
C  LT        - Integer scalar, number of non-blank characters in ISTR.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
C
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER CCKWRK(*)*(*), ISTR*(*), IDUM*80
      LOGICAL KERR, IERR
      INTEGER CKLKUP, CKLSCH
      EXTERNAL CKLKUP, CKLSCH
C
      ISTR = ' '
      ILEN = LEN(ISTR)
      KERR = .FALSE.
C
      ISRNU = CKLKUP (I, ICKWRK(IcRNU), NRNU)
      ISFAL = CKLKUP (I, ICKWRK(IcFL),  NFAL)
      ISTHB = CKLKUP (I, ICKWRK(IcTB),  NTHB)
      ISWL  = CKLKUP (I, ICKWRK(IcWL),  NWL)
C
      KFAL = 0
      IF (ISFAL.GT.0) KFAL = ICKWRK(IcKF + ISFAL - 1)
      IF (ISWL.GT.0)  WL   = RCKWRK(NcWL + ISWL  - 1)
C
      NSTART = 1
      NEND = MXSP/2
   50 CONTINUE
      N1 = NSTART
      N2 = NEND
      DO 100 N = N1, N2
         K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
C
         IF (K .GT. 0) THEN
C           append this species
            IF (N .GT. N1) THEN
C              need '+'
               NEXT = CKLSCH(ISTR)+1
               IF (NEXT .GT. ILEN) GO TO 100
               ISTR(NEXT:NEXT) = '+'
            ENDIF
            IDUM = ' '
            IF (ISRNU .EQ. 0) THEN
C              integer stoichiometric coefficient
               NU = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            ELSE
C              real stoichiometric coefficient
               RNU = ABS(RCKWRK(NcRNU + (ISRNU-1)*MXSP + N - 1))
C              integer part of coefficient
               NU = RNU
               IF (RNU-NU .NE. 0.0) THEN
C                 RNU has decimal part > 0, else can treat
C                 as whole number
                  NU = 0
C                 convert RNU to character string
                  CALL CKR2CH (RNU, IDUM, L, IERR)
                  IF (IERR) GO TO 200
                  IND = INDEX(IDUM,'.')
C                 restrict decimal part to 3 digits
                  L = MIN (L, IND+3)
               ENDIF
            ENDIF
            IF (NU .GT. 1) THEN
C              convert NU to character string
               CALL CKI2CH (NU, IDUM, L, IERR)
               IF (IERR) GO TO 200
            ENDIF
C
            IF (IDUM .NE. ' ') THEN
C              coefficient > 1 added to string
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+L .GT. ILEN) GO TO 400
               ISTR(NEXT:) = IDUM
            ENDIF
C
C           species name
            NEXT = CKLSCH(ISTR) + 1
            L = CKLSCH(CCKWRK(IcKK + K - 1))
            IF (NEXT+L .GT. ILEN) GO TO 400
            ISTR(NEXT:) = CCKWRK(IcKK + K - 1)
         ENDIF
C
         IF (K.EQ.0 .OR. N.EQ.N2) THEN
C           last product or reactant - supplemental info
            IF (ISWL .GT. 0) THEN
C              radiation enhancement factor
               IF (N .EQ. MXSP/2 .AND. WL.LT.0 .OR.
     1             N .EQ. MXSP   .AND. WL.GT.0) THEN
                   NEXT = CKLSCH(ISTR) + 1
                   IF (NEXT+3 .GT. ILEN) GO TO 400
                   ISTR(NEXT:) = '+HV'
               ENDIF
            ENDIF
C
            IF (ISFAL .GT. 0) THEN
C              pressure-dependence third-body species
               IDUM = '(+'
               IF (KFAL .GT. 0) THEN
                  IDUM(3:) = CCKWRK(IcKK + KFAL - 1)
               ELSE
                  IDUM(3:) = 'M'
               ENDIF
               NEXT = CKLSCH(IDUM) + 1
               IDUM(NEXT:) = ')'
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+CKLSCH(IDUM)-1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = IDUM
            ELSEIF (ISTHB .GT. 0) THEN
C              third-body reaction
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '+M'
            ENDIF
C
C           quit at end of products
            IF (N .GT. MXSP/2) GO TO 150
C
C           delimeter at end of reactants
            NEXT = CKLSCH(ISTR) + 1
            IF (ICKWRK(IcNS+I-1) .LT. 0) THEN
               IF (NEXT+1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '=>'
            ELSE
               IF (NEXT+2 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '<=>'
            ENDIF
            NSTART = NEND + 1
            NEND = MXSP
            GO TO 50
         ENDIF
  100 CONTINUE
C
  150 CONTINUE
      LT = CKLSCH(ISTR)
      RETURN
C
  200 CONTINUE
      WRITE (LOUT, 300)
  300 FORMAT (' Error in CKSYMR...character string length too small')
      ISTR = ' '
      LT = 0
      KERR = .TRUE.
      RETURN
C
  400 CONTINUE
      WRITE (LOUT, 500)
  500 FORMAT (' Syntax error in CKSYMR...')
      ISTR = ' '
      LT = 0
      KERR = .TRUE.
C
C     end of SUBROUTINE CKSYMR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)*
C  Returns the character strings of species names
C
C  INPUT
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENRCK.
C  LOUT      - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  KNAME(*)  - Character string array, species names;
C              dimension at least KK, the total species count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      CHARACTER CCKWRK(*)*(*), KNAME(*)*(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(KNAME(1))
      DO 150 K = 1, NKK
         LT = CKLSCH(CCKWRK(IcKK + K - 1))
         KNAME(K) = ' '
         IF (LT .LE. ILEN) THEN
            KNAME(K) = CCKWRK(IcKK + K - 1)
         ELSE
            WRITE (LOUT,*)
     1      ' Error in CKSYM...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
C
C     end of SUBROUTINE CKSYMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKTHB  (KDIM, ICKWRK, RCKWRK, AKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKTHB  (KDIM, ICKWRK, RCKWRK, AKI)
C  Returns matrix of enhanced third body coefficients; see Eq. (58).
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix AKI;
C              KDIM must be at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  AKI(*,*)  - Real matrix, enhanced third body efficiencies of the
C              species in the reactions;
C              dimension at least KK for first, the total species count,
C              and at least II for the second, the total reaction count.
C              AKI(K,I) is the enhanced efficiency of species K in
C              reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), AKI(KDIM,*)
C
      DO 150 I = 1, NII
         DO 100 K = 1, NKK
            AKI(K,I) = 1.0
  100    CONTINUE
  150 CONTINUE
C
      I_KTHB = IcKT
      I_KFAC = NcKT
      DO 250 N = 0, NTHB - 1
         I = ICKWRK(IcTB + N)
         DO 200 L = 0, ICKWRK(IcKN + N) - 1
            K        = ICKWRK(I_KTHB + L)
            AKI(K,I) = RCKWRK(I_KFAC + L)
  200    CONTINUE
         I_KTHB = I_KTHB + MXTB
         I_KFAC = I_KFAC + MXTB
  250 CONTINUE
C
C     end of SUBROUTINE CKTHB
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUBML (T, X, ICKWRK, RCKWRK, UBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUBML (T, X, ICKWRK, RCKWRK, UBML)
C  Returns the mean internal energy of the mixture in molar units;
C  see Eq. (39).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UBML      - Real scalar, mean internal energy.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      UBML = 0.0
      DO 100 K = 1, NKK
         UBML = UBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKUBM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUBMS (T, Y, ICKWRK, RCKWRK, UBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUBMS (T, Y, ICKWRK, RCKWRK, UBMS)
C  Returns the mean internal energy of the mixture in mass units;
C  see Eq. (40).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UBMS      - Real scalar, mean internal energy.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKUMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      UBMS = 0.0
      DO 100 K = 1, NKK
         UBMS = UBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKUBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUML  (T, ICKWRK, RCKWRK, UML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUML  (T, ICKWRK, RCKWRK, UML)
C  Returns the internal energies in molar units;  see Eq. (23).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UML(*)    - Real array, internal energies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), UML(*)
C
      CALL CKHML (T, ICKWRK, RCKWRK, UML)
      RU = RCKWRK(NcRU)
      DO 150 K = 1, NKK
         UML(K) = UML(K) - RU * T(ICKWRK(IcKTF + K - 1))
150   CONTINUE
C
C     end of SUBROUTINE CKUML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUMS  (T, ICKWRK, RCKWRK, UMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUMS  (T, ICKWRK, RCKWRK, UMS)
C  Returns the internal energies in mass units;  see Eq. (30).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UMS(*)    - Real array, internal energies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), UMS(*)
C
      CALL CKHMS (T, ICKWRK, RCKWRK, UMS)
      RU = RCKWRK(NcRU)
      DO 150 K = 1, NKK
         TEMP = T(ICKWRK(IcKTF + K - 1))
         UMS(K) = UMS(K) - TEMP*RU/RCKWRK(NcWT + K - 1)
150   CONTINUE
C
C     end of SUBROUTINE CKUMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWC   (T, C, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWC   (T, C, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given
C  temperature(s) and molar concentrations;  see Eq. (49).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKRATX (T, C, ICKWRK(IcNU), ICKWRK(IcNK), RCKWRK(NcCO),
     2             ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcFT),
     3             ICKWRK(IcKF), RCKWRK(NcFL), ICKWRK(IcTB),
     3             ICKWRK(IcKN), RCKWRK(NcKT), ICKWRK(IcKT),
     4             RCKWRK(NcKF), RCKWRK(NcKR), RCKWRK(NcI1),
     5             RCKWRK(NcI2), RCKWRK(NcI3), ICKWRK(IcRNU),
     6             RCKWRK(NcRNU), ICKWRK(IcORD), ICKWRK(IcKOR),
     7             RCKWRK(NcKOR), ICKWRK(IcMO), ICKWRK(IcXS))
      CALL CKDOT (RCKWRK(NcI1), RCKWRK(NcI2), ICKWRK, RCKWRK,
     1            RCKWRK(NcK1), RCKWRK(NcK2))
C
      NKM1 = NKK - 1
      DO 25 K = 0, NKM1
         WDOT(K+1) = RCKWRK(NcK1 + K) - RCKWRK(NcK2 + K)
   25 CONTINUE
C
C     end of SUBROUTINE CKWC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWL   (ICKWRK, RCKWRK, WL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWL   (ICKWRK, RCKWRK, WL)
C  Returns a set of flags providing information on the wave length
C  of photon radiation
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WL(*)     - Real array, radiation wavelengths for reactions;
C              dimension at least II, total reaction count.
C                 cgs units, angstrom.
C              WL(I)= 0.  reaction I does not have radiation as
C                         either a reactant or product
C              WL(I)=-A   reaction I has radiation of wavelength A
C                         as a reactant
C              WL(I)=+A   reaction I has radiation of wavelength A
C                         as a product
C              If A = 1.0 then no wavelength information was given;
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), WL(*)
C
      DO 100 I = 1, NII
         WL(I) = 0.0
  100 CONTINUE
      DO 150 N = 1, NWL
         WL(ICKWRK(IcWL+N-1)) = RCKWRK(NcWL+N-1)
  150 CONTINUE
C
C     end of SUBROUTINE CKWL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C  Returns the molecular weights of the species
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WT(*)     - Real array, molecular weights of the species;
C              dimension at least KK, the total species count.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), WT(*)
C
      DO 100 K = 1, NKK
         WT(K) = RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKWT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWXP  (P, T, X, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWXP  (P, T, X, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mole fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWXR  (RHO, T, X, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWXR  (RHO, T, X, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given mass
C  density, temperature(s) and mole fractions;  see Eq. (49).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWXR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mass fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYPK  (P, T, Y, RKFT, RKRT, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYPK  (P, T, Y, RKFT, RKRT, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mass fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), RKFT(*), RKRT(*), ICKWRK(*), RCKWRK(*),
     1          WDOT(*)
C
      DO 25 I = 1, NII
         RCKWRK(NcKF + I - 1) = RKFT(I)
         RCKWRK(NcKR + I - 1) = RKRT(I)
   25 CONTINUE
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKRATX (T, RCKWRK(NcK1), ICKWRK(IcNU), ICKWRK(IcNK),
     2             RCKWRK(NcCO), ICKWRK(IcFL), ICKWRK(IcFO),
     3             ICKWRK(IcFT), ICKWRK(IcKF), RCKWRK(NcFL),
     3             ICKWRK(IcTB), ICKWRK(IcKN), RCKWRK(NcKT),
     4             ICKWRK(IcKT), RCKWRK(NcKF), RCKWRK(NcKR),
     5             RCKWRK(NcI1), RCKWRK(NcI2), RCKWRK(NcI3),
     6             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcORD),
     7             ICKWRK(IcKOR), RCKWRK(NcKOR), ICKWRK(IcMO),
     8             ICKWRK(IcXS))
      CALL CKDOT (RCKWRK(NcI1), RCKWRK(NcI2), ICKWRK, RCKWRK,
     1            RCKWRK(NcK1), RCKWRK(NcK2))
C
      NKM1 = NKK - 1
      DO 50 K = 0, NKM1
         WDOT(K+1) = RCKWRK(NcK1 + K) - RCKWRK(NcK2 + K)
   50 CONTINUE
C
C     end of SUBROUTINE CKWYPK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given mass
C  density, temperature and mass fractions;  see Eq. (49).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature;
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C  Searches a character string, LINE, for NEXP space-delimited
C  substrings representing numbers, until NVAL real values are
C  converted and stored in the array, RVAL.
C  This allows format-free input of numerical data.  For example:
C
C     input:  LINE    = " 0.170E+14 0 47780.0"
C             NEXP    = 3, the number of values requested
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C     output: NVAL    = 3, the number of values found
C             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
C             KERR    = .FALSE.
C
C  INPUT
C  LINE   - Character string, length established by calling program.
C  NEXP   - Integer scalar, number of real values to be found in LINE;
C           If NEXP < 0 then IABS(NEXP) values are expected, but
C           it is not an error condition if fewer values are found.
C  LOUT   - Integer scalar, output unit for printed diagnostics.
C
C  OUTPUT
C  NVAL   - Integer scalar, count of real values found in LINE.
C  RVAL   - Real array, values converted from characters in LINE;
C           dimension at least NEXP.
C  KERR   - Logical, syntax or dimensioning error flag.
C
C  END PROLOGUE
C
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), ITEMP*80
      DIMENSION RVAL(*), RTEMP(80)
      LOGICAL KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
C----------Find Comment String (! signifies comment)
C
      ILEN = CKSLEN(LINE)
      NVAL = 0
      KERR = .FALSE.
C
      IF (ILEN .LE. 0) RETURN
      IF (ILEN .GT. 80) THEN
         WRITE (LOUT,*)     ' Error in CKXNUM...line length > 80 '
         WRITE (LOUT,'(A)') LINE
         KERR = .TRUE.
         RETURN
      ENDIF
C
      ITEMP = LINE(1:ILEN)
      CALL CKDTAB (ITEMP)
      IF (NEXP .LT. 0) THEN
         CALL CKPARR (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
      ELSE
         CALL CKPARR (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
         IF (IERR .EQ. 1) THEN
            WRITE (LOUT, *)    ' Syntax errors in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
         ELSEIF (NVAL .NE. NEXP) THEN
            WRITE (LOUT,*) ' Error in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
            WRITE (LOUT,*) NEXP,' values expected, ',
     1                     NVAL,' values found.'
         ENDIF
      ENDIF
      IF (NVAL .LE. IABS(NEXP)) THEN
         DO 20 N = 1, NVAL
            RVAL(N) = RTEMP(N)
   20    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKXNUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given pressure, temperature(s)
C  and mole fractions;  see Eq. (10).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUMXT = 0.0
      DO 100 K = 1, NKK
         SUMXT = SUMXT + X(K)*T(ICKWRK(IcKTF + K - 1))
 100  CONTINUE
      PRUT = P/(RCKWRK(NcRU)*SUMXT)
      DO 150 K = 1, NKK
         C(K) = X(K)*PRUT
150   CONTINUE
C
C     end of SUBROUTINE CKXTCP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTCR (RHO, T, X, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTCR (RHO, T, X, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given mass density, temperature(s),
C  and mole fractions;  see Eq. (11).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      RHOW = RHO / SUM
      DO 200 K = 1, NKK
         C(K) = X(K)*RHOW
200   CONTINUE
C
C     end of SUBROUTINE CKXTCR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
C  Returns the mass fractions given mole fractions; see Eq. (9).
C
C  INPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), Y(*)
C
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K) * RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      DO 200 K = 1, NKK
         Y(K) = X(K) * RCKWRK(NcWT + K - 1) / SUM
200   CONTINUE
C
C     end of SUBROUTINE CKXTY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given pressure, temperature(s)
C  and mass fractions;  see Eq. (7).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUMYOW = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUMYOW = SUMYOW + Y(K+1)*T(ICKWRK(IcKTF + K))
     1                          /RCKWRK(NcWT + K)
150   CONTINUE
      SUMYOW = SUMYOW*RCKWRK(NcRU)
      DO 200 K = 1, NKK
         C(K) = P*Y(K)/(SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
C
C     end of SUBROUTINE CKYTCP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTCR (RHO,T, Y, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCR (RHO,T, Y, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given mass density, temperature(s),
C  and mass fractions;  see Eq. (8).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), C(*)
C
      DO 150 K = 1, NKK
         C(K) = RHO * Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
C
C     end of SUBROUTINE CKYTCR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
C  Returns the mole fractions given mass fractions;  see Eq. (6).
C
C  INPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*), X(*)
C
      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      DO 200 K = 1, NKK
         X(K) = Y(K) / (SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
C
C     end of SUBROUTINE CKYTX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PKINDX (ICKWRK, KELECT, KKION)
C
C  START PROLOGUE
C
C  SUBROUTINE PKINDX (ICKWRK, KELECT, KKION)
C  Returns plasma indices for the particular reaction mechanism.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  KELECT    - Integer scalar, species array index for the electron.
C  KKION     - Integer scalar, total ion count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*)
C
      KKION  = NKKI
      KELECT = KEL
C
C     end of SUBROUTINE PKINDX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
C   BEGIN PROLOGUE  CKPARI
C   REFER TO  IPGETI
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851725   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses integer variables from a character variable.  Called
C            by IPGETI, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  CKPARI may be used for parsing an input record that contains integer
C  values, but was read into a character variable instead of directly
C  into integer variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by CKPARI to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of IVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1,   2,,40000   , ,60'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set IVAL(1) = 1
C     (2) set IVAL(2) = 2
C     (3) leave IVAL(3) unchanged
C     (4) set IVAL(4) = 40000
C     (5) leave IVAL(5) unchanged
C     (6) set IVAL(6) = 60
C
C   CKPARI will print diagnostics on the default output device, if
C   desired.
C
C   CKPARI is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume IVAL = (0, 0, 0) and NEXPEC = 3 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2 ,   3 45 '         (2, 3, 45)                0       3
C  '2.15,,3'               (2, 0, 3)                 1       0
C  '3X, 25, 2'             (0, 0, 0)                 1       0
C  '10000'                 (10000, 0, 0)             2       1
C
C      Assume IVAL = (0, 0, 0, 0) and NEXPEC = -4 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1, 2'                  (1, 2)                    0       2
C  ',,37  400'             (0, 0, 37, 400)           0       4
C  ' 1,,-3,,5'             (1, 0, -3, 0)             3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  IVAL (I,O) - the integer value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to IVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  CKFRCH,CKLSCH
C   END PROLOGUE  CKPARI
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION IVAL(*)
      CHARACTER *8 FMT(14)
      LOGICAL OKINCR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C   FIRST EXECUTABLE STATEMENT  CKPARI
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = CKLSCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set false when a space follows
C--- an integer value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR .OR. NC .EQ. IE) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into integer array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      IES = NC - 1
      NCH = IES - IBS + 1
      DATA FMT/' (I1)', ' (I2)', ' (I3)', ' (I4)', ' (I5)',
     1   ' (I6)', ' (I7)', ' (I8)', ' (I9)', '(I10)',
     2   '(I11)', '(I12)', '(I13)', '(I14)'/
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) IVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .NE. 0 .AND. ICARD .GE. 0) THEN
         IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
         IF (IERR .EQ. 1)
     1   WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
         IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
         IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      ENDIF
C
C     end of SUBROUTINE CKPARI
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
C   BEGIN PROLOGUE  CKPARR
C   REFER TO  IPGETR
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851625   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses real variables from a character variable.  Called
C            by IPGETR, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  CKPARR may be used for parsing an input record that contains real
C  values, but was read into a character variable instead of directly
C  into real variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by CKPARR to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of RVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1.,   2,,4.e-5   , ,6.e-6'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set RVAL(1) = 1.0
C     (2) set RVAL(2) = 2.0
C     (3) leave RVAL(3) unchanged
C     (4) set RVAL(4) = 4.0E-05
C     (5) leave RVAL(5) unchanged
C     (6) set RVAL(6) = 6.0E-06
C
C   CKPARR will print diagnostics on the default output device, if
C   desired.
C
C   CKPARR is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
C  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
C  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
C  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
C  '1.0'                   (1.0, 0.0, 0.0)           2       1
C
C      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1.,2.'                 (1.0, 2.0)                0       2
C  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
C  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  RVAL (I,O) - the real value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to RVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  CKLSCH
C   END PROLOGUE  CKPARR
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION RVAL(*)
      CHARACTER *8 FMT(22)
      LOGICAL OKINCR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C   FIRST EXECUTABLE STATEMENT  CKPARR
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = CKLSCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set negative when a space follows
C--- a real value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into real array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      DATA FMT/     ' (E1.0)', ' (E2.0)', ' (E3.0)', ' (E4.0)',
     1   ' (E5.0)', ' (E6.0)', ' (E7.0)', ' (E8.0)', ' (E9.0)',
     2   '(E10.0)', '(E11.0)', '(E12.0)', '(E13.0)', '(E14.0)',
     3   '(E15.0)', '(E16.0)', '(E17.0)', '(E18.0)', '(E19.0)',
     4   '(E20.0)', '(E21.0)', '(E22.0)'/
      IES = NC - 1
      NCH = IES - IBS + 1
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) RVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .NE. 0 .AND. ICARD .GE. 0) THEN
         IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
         IF (IERR .EQ. 1)
     1   WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
         IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
         IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      ENDIF
C
C     end of SUBROUTINE CKPARR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C     OBSOLETE ROUTINES...
C
      CHARACTER*(*) FUNCTION UPCASE (STR, ILEN)
      CHARACTER STR*(*), CKCHUP*128
      INTEGER ILEN
      EXTERNAL CKCHUP
      UPCASE = CKCHUP (STR, ILEN)
      RETURN
      END
C
      CHARACTER*(*) FUNCTION LOCASE (STR, ILEN)
      CHARACTER STR*(*), CKCHLO*128
      INTEGER ILEN
      EXTERNAL CKCHLO
      LOCASE = CKCHLO (STR, ILEN)
      RETURN
      END
C
      INTEGER FUNCTION IPPLEN (STR)
      CHARACTER*(*) STR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
      IPPLEN = CKSLEN(STR)
      RETURN
      END
C
      INTEGER FUNCTION ILASCH (STR)
      CHARACTER*(*) STR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
      ILASCH = CKLSCH(STR)
      RETURN
      END
C
      INTEGER FUNCTION IFIRCH (STR)
      CHARACTER*(*) STR
      INTEGER CKFRCH
      EXTERNAL CKFRCH
      IFIRCH = CKFRCH(STR)
      RETURN
      END
C
      SUBROUTINE IPPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
C*****precision > double
      DOUBLE PRECISION RVAL(*)
C*****END precision > double
C*****precision > single
C      REAL RVAL(*)
C*****END precision > single
      CHARACTER*(*) STRING
      INTEGER ICARD, NEXPEC, NFOUND, IERR, LOUT
      CALL CKPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
      RETURN
      END
C
      SUBROUTINE IPPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
      INTEGER ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT
      CHARACTER*(*) STRING
      DIMENSION IVAL(*)
      CALL CKPARI (STRING,ICARD,NEXPEC,IVAL,NFOUND,IERR,LOUT)
      RETURN
      END
