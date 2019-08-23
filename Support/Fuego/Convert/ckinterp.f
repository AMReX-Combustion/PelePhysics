C    CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:32 $
C
      SUBROUTINE CKINTP (MDIM, KDIM, IDIM, MAXTP, MAXSP, MAXTB,
     1                   LIN, LOUT, LTHRM, LINCK, MAXORD,
     2                   LIWORK, I, LRWORK, R, LCWORK, C, L)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NPAR=3, NCP=5, NFIT=NCP+2, NLAR=2, NFAR=8, NJAR=9,
     1           NF1R=4)
C     NPAR     - required number of reaction Arrhenius parameters
C     NCP      - number of polynomial coefficients in
C                species fits to specific heat
C     NFIT     - total number of species thermodynamic polynomial
C                fit coefficients
C     NTR      - number of temperature ranges for the fit
C     NLAR     - required number of Teller-Landau rate parameters
C     NFAR     - maximum number of pres-dependency rate parameters
C     NJAR     - required number of Jannev-Langer rate parameters
C     NF1R     - required number of FIT#1 rate parameters
C
      DIMENSION I(LIWORK), R(LRWORK)
      CHARACTER*16 C(LCWORK), PRVERS, PREC, FILVER, PRDATE
      LOGICAL L(KDIM), KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C      Version and Precision common blocks
      COMMON /CKVERS/ PRVERS, PREC, FILVER
C
C     Write to standard output information about the ckinterp program
      PRVERS = '6.15'
      PRDATE = '98/03/03'
      FILVER = '1.0'
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      WRITE  (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
     1   ' CHEMKIN-III GAS-PHASE MECHANISM INTERPRETER:',
     2    PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
     3    PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
     4   ' Copyright 1995, Sandia Corporation.',
     5' The U.S. Government retains a limited license in this software.'
C
      NTR = MAXTP - 1
      KERR = .FALSE.
C     integer arrays
      IKNCF  = 1
      IKPHSE = IKNCF + MDIM*KDIM
      IKCHRG = IKPHSE + KDIM
      INT    = IKCHRG + KDIM
      INSPEC = INT    + KDIM
      INREAC = INSPEC + IDIM
      INU    = INREAC + IDIM
      INUNK  = INU    + MAXSP*IDIM
      IIDUP  = INUNK  + MAXSP*IDIM
      IIREV  = IIDUP  + IDIM
      IILAN  = IIREV  + IDIM
      IIRLT  = IILAN  + IDIM
      IIWL   = IIRLT  + IDIM
      IIFAL  = IIWL   + IDIM
      IIFOP  = IIFAL  + IDIM
      IIFLO  = IIFOP  + IDIM
      IKFAL  = IIFLO  + IDIM
      IITHB  = IKFAL  + IDIM
      INTBS  = IITHB  + IDIM
      INKTB  = INTBS  + IDIM
      IIRNU  = INKTB  + MAXTB*IDIM
      IIORD  = IIRNU  + IDIM
      IKORD  = IIORD  + IDIM
      IKION  = IKORD  + MAXORD*IDIM
      IIEIM  = IKION  + KDIM
      IIEIMT = IIEIM  + IDIM
      IIJAN  = IIEIMT + IDIM
      IIFT1  = IIJAN  + IDIM
      IIEXC  = IIFT1  + IDIM
      IIMOM  = IIEXC  + IDIM
      IKMOM  = IIMOM  + IDIM
      IIXSM  = IKMOM  + IDIM
      IIXSK  = IIXSM  + IDIM
      IIXSI  = IIXSK  + IDIM
      IITDE  = IIXSI  + IDIM
      IITDK  = IITDE  + IDIM
      IITOT  = IITDK  + IDIM - 1
      IF (LIWORK .LT. IITOT) THEN
         WRITE (LOUT, *)
     1   'CKINTP ERROR: IWORK array needs to be at least ',IITOT
         KERR = .TRUE.
      ENDIF
C     Real arrays
      IAWT  = 1
      IWT   = IAWT + MDIM
      IA    = IWT  + KDIM
      IT    = IA   + NFIT*NTR*KDIM
      IPAR  = IT   + MAXTP*KDIM
      IRPAR = IPAR + NPAR*IDIM
      IPLAN = IRPAR+ NPAR*IDIM
      IRLAN = IPLAN+ NLAR*IDIM
      IWL   = IRLAN+ NLAR*IDIM
      IPFAL = IWL  + IDIM
      IAIK  = IPFAL+ NFAR*IDIM
      IRNU  = IAIK + MAXTB*IDIM
      IRORD = IRNU + MAXSP*IDIM
      IPJAN = IRORD+ MAXORD*IDIM
      IPFT1 = IPJAN+ NJAR*IDIM
      IPEXC = IPFT1+ NF1R*IDIM
      NTOT  = IPEXC+ IDIM - 1
      IF (LRWORK .LT. NTOT) THEN
         WRITE (LOUT, *)
     1   'CKINTP ERROR: RWORK array needs to be at least ',NTOT
         KERR = .TRUE.
      ENDIF
C     Character arrays
      ICKNAM = 1
      ICENAM = ICKNAM + KDIM
      ICTOT  = ICENAM + MDIM - 1
      IF (LCWORK .LT. ICTOT) THEN
         WRITE (LOUT, *)
     1   'CKINTP ERROR: CWORK array needs to be at least ',ICTOT
         KERR = .TRUE.
      ENDIF
C
      IF (KERR) RETURN
      CALL CKMECH (MDIM, KDIM, IDIM, NPAR, NCP, NFIT, MAXTP, NTR,
     1     MAXSP, MAXTB, NLAR, NFAR, LIN, LOUT, LTHRM, LINCK, MAXORD,
     2     NJAR, NF1R, I(IKNCF), I(IKPHSE), I(IKCHRG), I(INT),
     3     I(INSPEC), I(INREAC), I(INU), I(INUNK), I(IIDUP), I(IIREV),
     4     I(IILAN), I(IIRLT), I(IIWL), I(IIFAL), I(IIFOP), I(IIFLO),
     5     I(IKFAL), I(IITHB), I(INTBS), I(INKTB), I(IIRNU), I(IIORD),
     6     I(IKORD), I(IKION), I(IIEIM), I(IIEIMT), I(IIJAN), I(IIFT1),
     7     I(IIEXC), I(IIMOM), I(IKMOM), I(IIXSM), I(IIXSK), I(IIXSI),
     8     I(IITDE), I(IITDK), R(IAWT), R(IWT), R(IA), R(IT), R(IPAR),
     9     R(IRPAR), R(IPLAN), R(IRLAN), R(IWL), R(IPFAL), R(IAIK),
     *     R(IRNU), R(IRORD), R(IPJAN), R(IPFT1), R(IPEXC), C(ICKNAM),
     1     C(ICENAM), L)
      RETURN
      END
C
      SUBROUTINE CKMECH
     1(MDIM, KDIM, IDIM, NPAR, NCP, NFIT, MAXTP, NTR, MAXSP,
     2 MAXTB, NLAR, NFAR, LIN, LOUT, LTHRM, LINCK, MAXORD,
     3 NJAR, NF1R, KNCF, KPHSE, KCHRG, NT, NSPEC, NREAC, NU,
     4 NUNK, IDUP, IREV, ILAN, IRLT, IWL, IFAL, IFOP, IFLO, KFAL,
     5 ITHB, NTBS, NKTB, IRNU, IORD, KORD, KION, IEIM, IEIMT, IJAN,
     6 IFT1, IEXC, IMOM, KMOM, IXSM, IXSK, IXSI, ITDE, ITDK, AWT,
     7 WT, A, T, PAR, RPAR, PLAN, RLAN, WL, PFAL, AIK, RNU, RORD,
     8 PJAN, PFT1, PEXC, KNAME, ENAME, ITHRM)
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (CKMIN = 0.001)
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN,
     2               NFT1, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     3               MELECT, KELECT, NKION
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM),  KPHSE(KDIM),      KCHRG(KDIM),
     $          NT(KDIM),         NSPEC(IDIM),      NREAC(IDIM),
     $          NU(MAXSP,IDIM),   NUNK(MAXSP,IDIM), IDUP(IDIM),
     $          IREV(IDIM),       ILAN(IDIM),       IRLT(IDIM),
     $          IWL(IDIM),        IFAL(IDIM),       IFOP(IDIM),
     1          IFLO(IDIM),
     $          KFAL(IDIM),       ITHB(IDIM),       NTBS(IDIM),
     $          NKTB(MAXTB,IDIM), IRNU(IDIM),       IORD(IDIM),
     $          KORD(MAXORD,IDIM),KION(KDIM),       IEIM(IDIM),
     $          IEIMT(IDIM),      IJAN(IDIM),       IFT1(IDIM),
     $          IEXC(IDIM),       IMOM(IDIM),       KMOM(IDIM),
     $          IXSM(IDIM),       IXSK(IDIM),       IXSI(IDIM),
     $          ITDE(IDIM),       ITDK(IDIM)
C
C     Real arrays
      DIMENSION AWT(MDIM),        WT(KDIM),        A(NFIT,NTR,KDIM),
     $          T(MAXTP,KDIM),    PAR(NPAR,IDIM),   RPAR(NPAR,IDIM),
     $          PLAN(NLAR,IDIM),  RLAN(NLAR,IDIM),  WL(IDIM),
     $          PFAL(NFAR,IDIM),  AIK(MAXTB,IDIM),  RNU(MAXSP,IDIM),
     $          RORD(MAXORD,IDIM),PJAN(NJAR,IDIM),  PFT1(NF1R,IDIM),
     $          PEXC(IDIM)
C
C     Character arrays
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM), PRVERS, PREC, FILVER
      COMMON /CKVERS/ PRVERS, PREC, FILVER
C
C     Logical array and error flag
      LOGICAL ITHRM(KDIM), KERR
C
C
C     Description of COMMON /CKINT/ variables:
C     LENICK - minimum size for a CHEMKIN integer workspace array
C     LENRCK - minimum size for a CHEMKIN real workspace array
C     LENCCK - minimum size for a CHEMKIN character workspace array
C     MM     - count of mechanism elements
C     KK     - count of mechanism species
C     II     - count of mechanism reactions
C     NREV   - reactions with explicit reverse Arrhenius parameters
C     NFAL   - reactions with pres dependency
C     NTHB   - reactions with third-bodies
C     NLAN   - Landau-Teller reactions
C     NRLT   - Landau-Teller reactions with explicit reverse parameters
C     NWL    - wavelength-enhanced reactions
C     NCHRG  - count of ionic species
C     NEIM   - electron-impact reactions
C     NJAN   - Jannev-Langer formulation reactions
C     NFT1   - fit#1-type reactions
C     NEXC   - electron excitation reactions
C     NMOM   - electron-momentum reactions
C     NXSM   - ion-momentum reactions
C     NTDE   - non-thermal temp dependency reactions
C     NRNU   - real-stoichiometry reactions
C     NORD   - change-order reactions
C     MELECT - element index of the electron
C     KELECT - species index of the electron species
C     NKION  - non-electron ionic species
C
C     initialize arrays
      CALL CKSET
     1     (MDIM, KDIM, IDIM, NPAR, NCP, NFIT, MAXTP, NTR, MAXSP,
     1      MAXTB, NLAR, NFAR, MAXORD, KNAME, ENAME, AWT, KNCF, WT,
     2      KPHSE, KCHRG, A, T, NT, NSPEC, NREAC, NU, NUNK, PAR,
     3      IDUP, IREV, RPAR, ILAN, PLAN, IRLT, RLAN, IWL, WL, IFAL,
     4      IFOP, IFLO, KFAL, PFAL, ITHB, NTBS, AIK, NKTB, IRNU,
     5      RNU, IORD, KORD, RORD, KION, NJAR, NF1R, IEIM, IEIMT,
     6      IJAN, PJAN, IFT1, PFT1, IEXC, PEXC, IMOM, KMOM, IXSM,
     7      IXSK, IXSI, ITDE, ITDK, ITHRM)
C
C     mechanism interpretation
      KERR = .FALSE.
      CALL CKKEY
     1     (LIN, LOUT, LTHRM, MDIM, KDIM, IDIM, NPAR, NCP, NFIT,
     1      MAXTP, NTR, MAXSP, MAXTB, NLAR, NFAR, MAXORD, KNAME,
     2      ENAME, AWT, KNCF, WT, KPHSE, KCHRG, A, T, NT, NSPEC,
     3      NREAC, NU, NUNK, PAR, IDUP, IREV, RPAR, ILAN, PLAN,
     3      IRLT, RLAN, IWL, WL, IFAL, IFOP, IFLO, KFAL, PFAL,
     4      ITHB, NTBS, AIK, NKTB, IRNU, RNU, IORD, KORD, RORD,
     5      KION, NJAR, NF1R, IEIM, IEIMT, IJAN, PJAN, IFT1, PFT1,
     7      IEXC, PEXC, IMOM, KMOM, IXSM, IXSK, IXSI, ITDE, ITDK,
     8      ITHRM, CKMIN, KERR)
C
C     set LENICK, LENRCK, LENCCK
      MXTB = 0
      DO 111 N = 1, NTHB
         MXTB = MAX (MXTB, NTBS(N))
  111 CONTINUE
      MXTK = 0
      DO 112 K = 1, KK
         MXTK = MAX (MXTK, NT(K))
  112 CONTINUE
      CALL CKSIZE (NPAR, MAXSP, MXTB, MXTK, MXTK-1, NFIT,
     1             NFAR, NLAR, NJAR, NF1R, MAXORD)
C
C     Write linkfile out
C*****linkfile (gas) > binary
C      CALL CKBIN
C*****END linkfile (gas) > binary
C*****linkfile (gas) > ascii
      CALL CKFORM
C*****END linkfile (gas) > ascii
     1   (LINCK, LOUT, KERR, CKMIN, MDIM, KDIM, IDIM, NFIT, NTR, MAXSP,
     2    MAXTB, MAXTP, NCP, NPAR, NLAR, NFAR, NJAR, MAXORD, NF1R,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, T, A, KION,
     4    NSPEC, NREAC, NU, NUNK, PAR, IREV, RPAR, IFAL, IFOP, IFLO,
     5    KFAL, PFAL, ITHB, NTBS, NKTB, AIK, ILAN, PLAN, IRLT, RLAN,
     6    IWL, WL, IEIM, IEIMT, IJAN, PJAN, IFT1, PFT1, IEXC, PEXC,
     7    IMOM, KMOM, IXSM, IXSI, IXSK, ITDE, ITDK, IRNU, RNU, IORD,
     8    KORD, RORD)
C
      WRITE (LOUT, '(/A,3(/A,I6))')
     1      ' WORKING SPACE REQUIREMENTS ARE',
     2      '    INTEGER:   ',LENICK,
     3      '    REAL:      ',LENRCK,
     4      '    CHARACTER: ',LENCCK
C
C     end of SUBROUTINE CKMECH
 100  CONTINUE
      RETURN
      END
C                                                                      
      SUBROUTINE CKIABS
C
C     V.6.15 98/03/03 (E. Meeks)
C     1. Fix action#103: Remove unused variable CKLSCH in  CKMECH.
C     V.6.14 97/11/11 (E. Meeks)
C     1. Fix bug#128: Correct call to CKPARR at line 2791 to look for
C        3 parameters rather than 2 for chemically activated reactions
C     V.6.13 97/08/18
C     1. Fix bug#074b:  in CKTHRM check NVAL>0 to avoid error in 
C        IELEM=VALUE(1). (F. Rupley)
C     V.6.12 97/08/04 
C     1. Fix bug #033: missing line setting LRNU=.TRUE. just befor
C        loop 111 in CKREAC. (E. Meeks)
C     V.6.11 97/07/23 F. Rupley
C     1. Fix bug #012, extra spaces in logical operators; 
C        in CKCHAR, change (I1 .EQ . 1) to (I1 .EQ. 1). (F. Rupley)
C     V.6.10 97/07/22 F. Rupley
C     1. CKTHRM, comment out IF (KERR) RETURN in order to
C        continue through thermo data in spite of previous error
C     V.6.9 97/04/15
C     1. logical function SKFILE is used in CKKEY to check for
C        existence of opened thermodynamics file;
C        CHEMKIN Committee bugfix #001
C     V.6.8 97/03/11
C     1. correct element balance checking in CKBAL
C     V.6.7 97/03/05
C     1. IFMT should be RFMT in CKFORM to write RORD array
C     V.6.6 97/03/02
C     1. correct order-change bugs in CKAUXL (add NSPEC), and
C        CALL CKNORD (NSPEC and NREAC arguments were reversed,
C        indexing in loop 100 is N, not NKORD)
C     V.6.5 97/03/01
C     1. new main "driver" program to set up arrays,
C        to satisfy Reaction Design requirement of providing object
C        files instead of source code - move OPENs to driver
C     V.6.4 97/02/11
C     1. new SUBROUTINE CKIABS to store version comments
C     V.6.3, 97/01/25
C     1. allow more than 3 fit temperatures for THERMO data
C     V.6.2, 96/12/18
C     1. fix so that MOME or XSM will CALL CKAUXL
C     V.6.1, 96/12/11
C     1. several bugfixes
C     V.6.0, 96/11/11
C     1. moved some routines into cklib.f so that now ckinterp needs
C        to be linked to the library
C     V.5.9, 96/09/04
C     1. reduced number of continuation lines in CKSIZE, as exceeded
C        limits of the SUN/Solaris2.3/SunPro Fortran3.0
C     VERSION 5.8, 96/08/13
C     1. correct CKREAC indexing in searching for reaction-string
C         pressure-dependence species
C     2. change CKDUP to search for differences, rather than sameness
C     VERSION 5.7, 96/08/06
C     1. OPEN therm.dat as "OLD" and initialize ISTR=' '
C     2. CALL CKNEIM arguments were reversed
C     3. CALL CKUNIT superfluous arguments deleted
C     3. CALL CKUNIT superfluous arguments deleted
C     VERSION 5.6, 96/08/05
C     1. bugfixes for CKNMOM, CKFORM, CKBIN, CPMOME and CPXSMI
C        per Ellen Meeks.
C     VERSION 5.5, 96/05/23
C     1. initial sccs version
C     VERSION 5.4 (F. Rupley, May 2, 1996)
C     1. declare LRNU logical in CKAUXL, KERR in CKKEY, CKMOME,
C        CPXSMI
C     VERSION 5.3 (F. Rupley, May 1, 1996)
C     1. use MDIM, KDIM, and IDIM for dimensions in CKBIN and CKFORM,
C        as opposed to previously-used MM, KK and II (or NTHB, etc),
C        which might be zero
C     VERSION 5.2 (F. Rupley, April 26, 1996)
C     1. delete second occurrence of KRTB=-1 and KPTB=-1, which would
C        cause KFAL(I) to be set incorrectly, thereby affecting rate
C        calculations in fall-off (pressure-dependent) reactions
C     2. several occurrences of CALL CKINUM had (...KORD(L,NORD)...,
C        which should have been KORD(1,NORD)
C     3. CKNTHB call list needs II
C     VERSION 5.1 (F. Rupley, April 14, 1996)
C     CAUTION...THIS WORK IS STILL IN PROGRESS:
C     1. declare KERR logical in CKBIN
C     2. REWIND LTHRM in CKTHRM
C     3. ascii linkfile named "chem.asc"
C     VERSION 5.0 (initial CHEMKIN-III gas-phase interpreter)
C     1.  binary/ascii linkfile options
C     2.  LIN=5 for mechanism, now allows use of redirect "<" input.
C     3.  LOUT=6 for formatted output, now allows use of redirect '>'
C         output.
C     4.  "fall-off" reaction additional logic for the case where
C         Arrhenius parameters are for the low-pressure rate and
C         HIGH keyword parameters are for high-pressure rate,
C         in addition to the original formulation with LOW
C         keyword for the opposite case; this required additional
C         array IFLO(*) for LOW vs HIGH cases.
C     5.  reorganized logic in most of the program and subroutines,
C         added separate "CKN..." subroutines for auxiliary reaction
C         option processing, added CKSET for initialization.
C     6.  added PROLOG comment sections for subroutines.
C     CHANGES FOR VERSION 4.1 (8/14/95 E. Meeks)
C     1.  Added KEL, the electron species index; NKION, total
C         number of ions, and KION(NKION), species indices for
C         the NKION ions to linkfile.
C     2.  Added temperature dependent reaction auxiliary keyword, TDEP,
C         where TDEP/species name/ indicates the species upon whose
C         temperature the reaction depends.
C     3.  Added several additional auxiliary keywords to describe
C         plasma processes:  MOME// for electron momentum-transfer
C         collision frequency[/s]; XSMI// for ion momentum-transfer
C         collision cross section [cm2]; EXCI/energy loss in eV/
C         for excitation-only reactions that are used to track electron
C         energy losses, without tracking all excited states
C     CHANGES FOR VERSION 4.0 (2/27/95 F. Rupley)
C     1.  Change character index "(:" to "(1:"
C     CHANGES FOR VERSION 3.9 (8/2/94 H. Moffat)
C     1.  Changed gas constant and Avrog numbers to 1986 CODATA
C         recommendations. RUC is now compatible with RU up to
C         machine precision.
C     2.  Fixed error in CPREAC, in that conversion factors were
C         only being calculated at single precision values on
C         double precision workstations.
C     3.  Reduced default lengths of KMAX and IMAX to something
C         reasonable, executable decreased from 6 Meg to 1 Meg.
C     CHANGES FOR VERSION 3.8 (6/13/94 F. Rupley per H. Moffat)
C     1.  Changed gas constant and Avrog numbers to 1986 CODATA
C         recommendations.  RUC is now compatible with RU up to
C         machine precision.
C     2.  Fixed error in CPREAC, in that conversion facors were only
C         being calculated at single precision values.
C     CHANGES FOR VERSION 3.6c (6/3/94 F. Rupley per H. Moffat)
C     1.  add ERR= or END= logic and error messages for chem.inp
C         and therm.dat input
C     2.  Allow comment lines (!) in thermodynamic data
C     CHANGES FOR VERSION 3.6b (5/20/94 F. Rupley per E. Meeks)
C     1.  Incorporate plasma options
C     CHANGES FOR VERSION 3.6 (4/29/94 F. Rupley)
C     1.  Cannot change RORD if reaction is irreversible.
C     CHANGES FOR VERSION 3.5 (4/19/94 F. Rupley)
C     1.  Fix bug with index NUNK(N) for CKBAL and CKRBAL.
C     CHANGES FOR VERSION 3.4 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 3.3 (1/26/94 F. Rupley per R. Kee)
C     1.  Real stoichometric coefficients used in a supplemental way;
C         NRNU total number of reactions using real stoichometry,
C         IRNU array of reaction indices, RNU real coefficients.
C     CHANGES FOR VERSION 3.2 (11/11/93 F. Rupley per T.U.Delft)
C     1.  Ensure that SUBROUTINE CKUNIT does not check for units beyond
C         end of LINE.
C     CHANGES FOR VERSION 3.1 (2/24/93 F. Rupley per C. Westbrook,LLNL)
C     1.  Problem in CKREAC for species starting with "M", where
C         "+M" is signal for third-body.
C     CHANGES FOR VERSION 3.0 (4/13/92 F. Rupley per M. Coltrin)
C     1.  Correct logic in CKDUP, add argument to call list.
C     CHANGES FOR VERSION 2.9 (2/24/92 F. Rupley)
C     1.  Check that reverse (REV) parameters were given when
C         RTL reverse Teller-Landauer parameters are given.
C     2.  Add 2*II to length of real work space
C     CHANGES FOR VERSION 2.8
C     1.  Change output format to print all 16 characters for
C         a species name.
C     CHANGES FOR VERSION 2.7
C     1.  Two otherwise duplicate reactions are unique if one
C         is a third body reaction and the other not.
C     CHANGES FOR VERSION 2.6
C     1.  LENRCK lengthened by II+NREV to reflect additional
C         work space needed by CKRAT for a 4th parameter
C         (perturbation factor).
C     CHANGES FOR VERSION 2.5
C     1.  Need to get TLO,THI,TMID from database BEFORE reading
C         user's THERMO data (unless THERMO ALL option is used)
C     CHANGES FOR VERSION 2.4
C     1.  Additional checking of TLO,TMID,THI for species -
C         a) set initial values at -1.
C         b) if user has not provided a TLO,TMID, or THI, use
C            values provided by THERMO.DAT.
C         c) check that TLO < THI, TLO <= TMID <= THI
C      CHANGES FOR VERSION 2.3
C     1.  In CKPREAC, error correction of 10/18/90 (above, V2.1).
C     CHANGES FOR VERSION 2.2
C     1.  11/14/90 (F. Rupley per M. Coltrin):
C         Initialize variable NCHRG
C     CHANGES FOR VERSION 2.1
C     1.  10/18/90 (F. Rupley):
C         Error in scaling pre-exponential constants RPAR(3,*)
C         where REV is declared, and FPAL(3,*) for fall-off reactions,
C         as RPAR(3,II)*EFAC should read RPAR(3,NREV), and
C            FPAL(3,II)*EFAC should read FPAL(3,NFAL).
C         This error was introduced in CKINTERP.15 during refinement
C         Dof units conversion routines.
C     2.  Subroutine CKDUP modified to recognize that two reactions
C         may be duplicate except for a third-body species in a
C         fall-off reaction.
C     CHANGES FOR VERSION 2.0
C     1.  Error in UPCASE could cause interpreter to ignore some
C         keywords.
C     CHANGES FOR VERSION 1.9
C     1.  First record of binary file now consists of a character
C         string version, precision, and logical error flag
C     CHANGES FOR VERSION 1.8
C     1.  Change Subroutine CKUNIT to parse LINE instead of SUB(*)
C         in order to correct misinterpretation of unit strings
C         with slashes.
C     CHANGES FROM VERSION 1.7
C     1.  Further correction of molecules conversion for fall-off
C         and third-body reactions
C     CHANGES FROM VERSION 1.5
C     1.  Correct molecules to moles unit conversion
C     2.  Correct UPCASE to avoid dimensioning errors
C     CHANGES FROM VERSION 1.4
C     1.  Modify OPEN statements
C     CHANGES FROM VERSION 1.3
C     1.  Add "unix" change blocks
C     CHANGES FROM VERSION 1.2
C     1.  Reaction delimiters are now "=" or "<=>" if reversible,
C                                            " =>" if irreversible.
C     2.  Fixed an error with IFIRCH(LINE) in IPPLEN
C     CHANGES FROM VERSION 1.1
C     1.  Changed CHARACTER*100 to CHARACTER*80
C     2.  Added THERMO "ALL" option
C     3.  Write LENICK, LENRCK, LENCCK to binary file
C     4.  Allow reaction species to end in '=' or '-'
C     5.  Allow real values of elemental composition in THERMO cards
C     6.  Allow upper/lower case input
C     CHANGES FROM VERSION 1.0
C     1.  Changed from REAL*8 to DOUBLE PRECISION
C
      RETURN
      END
C
      SUBROUTINE CKAUXL (KDIM, IDIM, SUB, NSUB, KNAME, LOUT, MAXSP,
     1                   NPAR, NREAC, NSPEC, ITHB, NTBS, MAXTB, NKTB,
     2                   AIK,
     2                   IFAL, IDUP, NFAR, PFAL, IFOP, IFLO, ILAN, NLAR,
     3                   PLAN, IREV, RPAR, IRLT, RLAN, IWL, WL, KERR,
     4                   IORD, MAXORD, KORD, RORD, NUNK, NU, IRNU, RNU,
     5                   IEIM, IEIMT, IJAN, NJAR, PJAN, IFT1, NF1R,
     6                   PFT1, IEXC, PEXC, IMOM, IXSM, ITDE, ITDK,
     7                   KCHRG, AXUNIT, EXUNIT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKAUXL
C  processes character strings in lines which follow a reaction
C  to specify additional options;
C  strings are in the form of a "KEYword" (3 characters are usually
C  sufficient, but more may be used to better describe the option),
C  and may be followed by slash(/)-delimited arguments:
C
C  DUPlicate           - reaction is allowed to be duplicated.
C  MOM                 - electron momentum transfer, CALL CKNMOM.
C  XSM                 - ion momentum transfer,      CALL CKNXSM.
C  KNAME/val1/         - enhanced third-body input,  CALL CKNTHB.
C                        *required "+M" in reaction.
C  LOW/val(n), n=1,3/  - low-pressure fall-off,      CALL CKNFAL.
C  HIGH/val(n),n=1,3/  - high-pressure dependency,   CALL CKNFAL.
C  TROe/val(n), n=1,3/ - 3-parameter TROE option,    CALL CKNFAL.
C  TROe/val(n), n=1,4/ - 4-parameter TROE option,    CALL CKNFAL.
C  SRI/val(n), n=1,3/  - SRI formulation,            CALL CKNFAL.
C                        *LOW, TROE, SRI required
C                        "(+M)" or "(+KNAME)" in reaction.
C  LT/val1 val2/       - Landau-Teller formulation,  CALL CKNLAN.
C  REV/val(n), n=1,3/  - explicit reverse rate,      CALL CKNREV.
C  RLT/val(n), n=1,2/  - Landau-Teller reverse rate, CALL CKNRLT.
C  EIM/val1/           - electron-impact reaction,   CALL CKNEIM.
C  TDEp/val1/          - non-thermal temp dependency,CALL CKNTDE.
C  HV/val1/            - radiation-enhancement,      CALL CKNWL.
C                        *required "+HV" in reaction.
C  FORd/KNAME val1/    - species forward order,      CALL CKNORD.
C  RORd/KNAME val1/    - species reverse order,      CALL CKNORD.
C  JAN/val(n), n=1,9/  - Jannev, et al. rate,        CALL CKNJAN.
C  FIT1/val(n), n=1,4/ - fit #1-type rate,           CALL CKNFT1.
C  EXC/val1/           - excitation/energy loss/,    CALL CKNEXC.
C  UNIt/string/        - set rate parameter units,   CALL CKUNIT.
C
C  Arguments:
C  SUB(*)   - Character string array;
C             SUB(*) would contain an option's keyword and parameters.
C  NSUB     - Integer scalar, character string count.
C  KNAME(*) - Character string array, species names.
C  LOUT     - Integer scalar, formatted output file unit number.
C  MAXSP    - Integer scalar, maximum reaction species.
C  NPAR     - Integer scalar, required number of Arrhenius params.
C  NREAC(*) - Integer array, reactions count of reactants only.
C  NSPEC(*) - Integer array, reactions count of reactants+products.
C  ITHB(*)  - Integer array, third-body reaction indices.
C  NTBS(*)  - Integer array, NTHB reactions third-body counts.
C  MAXTB    - Integer scalar, maximum reaction third-bodies.
C  NKTB(*,*)- Integer matrix, reactions third-body species indices.
C  AIK(*,*) - Real matrix, third-body enhancement factors.
C  IFAL(*)  - Integer array, pres dependency reaction indices.
C  IDUP(*)  - Integer array, reaction duplication flag,
C             -1, reaction is allowed to be duplicated,
C              0, not.
C  NFAR     - Integer scalar, maximum pres dependency params.
C  PFAL(*,*)- Real matrix, pres dependency parameters.
C  IFOP(*)  - Integer array, pres dependency formulation,
C             0, none
C             1, Lindemann
C             2, SRI
C             3, 3-parameter TROE
C             4, 4-parameter TROE
C  IFLO(*)  - Integer array, description of pres dependency.
C             0, Arrhenius rate is high-pressure limit,
C                rate-modify the low-pressure region.
C             1, Arrhenius rate is low-pressure limit,
C                rate-modify the high-pressure region.
C  ILAN(*)  - Integer array, Landau-Teller reaction indices.
C  NLAR     - Integer scalar, required number Landau-Teller params.
C  PLAN(*,*)- Real matrix, Landau-Teller parameters.
C  IREV(*)  - Integer array, explicit reverse parameter rxn indices.
C  RPAR(*,*)- Real matrix, explicit reverse parameters.
C  IRLT(*)  - Integer array, Landau-Teller reaction indices.
C  RLAN(*,*)- Real matrix, reverse Landau-Teller parameters.
C  IWL(*)   - Integer array, radiation-enhanced reaction indices.
C  WL(*)    - Real array, radiation wavelengths.
C  KERR     - Logical error flag.
C  IORD(*)  - Integer array, changed-order species indices.
C  MAXORD   - Integer scalar, maximum reaction change-orders.
C  KORD(*,*)- Integer matrix, reaction change-order species indices,
C             < 0, change forward order for species KORD
C             > 0, change reverse order for species KORD
C  RORD(*,*)- Real matrix, reaction change-order species values.
C  NUNK(*)  - Integer array, reaction species indices.
C  NU(*)    - Integer array, reaction stoichiometric coefficients.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  IEIM(*)  - Integer array, the NEIM reaction indices.
C  IEIMT(*) - Integer array, NEIM temp dependency flags.
C  IJAN(*)  - Integer array, Jannev-Langer reaction indices.
C  NJAR     - Integer scalar, required number of Jannev params.
C  PJAN(*,*)- Real matrix, Jannev reaction parameters.
C  IFT1(*)  - Integer array, Fit#1 reaction indices.
C  NF1R     - Integer scalar, required number of FIT1 parameters.
C  PFT1(*,*)- Real matrix, FIT1 parameters.
C  IEXC(*)  - Integer array, the NEXC reaction indices.
C  PEXC(*)  - Real array, the NEXT reactions energy loss (eV).
C  IMOM(*)  - Integer array, momentum-transfer reaction indices.
C  IXSM(*)  - Integer array, cross-section reaction indices.
C  ITDE(*)  - Integer array, non-thermal reaction indices.
C  ITDK(*)  - Integer array, temp dependency species indices.
C  KCHRG(*) - Integer array, species ionic charges.
C  AXUNIT   - Character string, reaction A units description.
C  EXUNIT   - Character string, reaction E units description.
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
C----------------------------------------------------------------------C
C
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Integer arrays
      DIMENSION ITHB(IDIM), NTBS(IDIM), NKTB(MAXTB,IDIM), IDUP(IDIM),
     1          IFAL(IDIM), IFOP(IDIM), IFLO(IDIM), ILAN(IDIM),
     2          IREV(IDIM), IRLT(IDIM), IWL(IDIM), IORD(IDIM),
     3          KORD(MAXORD,IDIM), NUNK(MAXSP), NU(MAXSP), IRNU(IDIM),
     4          KCHRG(KDIM), IEIM(IDIM), IEIMT(IDIM), IJAN(IDIM),
     5          IFT1(IDIM), IEXC(IDIM), IMOM(IDIM), IXSM(IDIM),
     6          ITDE(IDIM), ITDK(IDIM)
C     Real arrays
      DIMENSION AIK(MAXTB,IDIM), PFAL(NFAR,IDIM), PLAN(NLAR,IDIM),
     1          RPAR(NPAR,IDIM), RLAN(NLAR,IDIM), WL(IDIM),
     2          RORD(MAXORD,IDIM), RNU(MAXSP,IDIM), PJAN(NJAR,IDIM),
     3          PFT1(NF1R,IDIM), PEXC(IDIM)
C
      CHARACTER SUB(NSUB)*80, KNAME(KDIM)*16, AXUNIT*(*), EXUNIT*(*),
     1          RSTR*80, IUNITS*80, KEY*4, CKCHUP*4
C
      LOGICAL KERR, LLAN, LRLT, LTHB, LRNU, LFAL, LTRO, LSRI, LWL,
     1        LREV, LORD, LEIM, LJAN, LFT1, LEXC, LMOM, LXSM, LTDE
C
      INTEGER CKLSCH
      EXTERNAL CKLSCH, CKCHUP
C
      DO 500 N = 1, NSUB
         ILEN = CKLSCH(SUB(N))
         IF (ILEN .LE. 0) GO TO 500
C        check status of this reaction so far
         LTHB = II .EQ. ITHB(MAX(NTHB,1))
         LRNU = II .EQ. IRNU(MAX(NRNU,1))
         LORD = II .EQ. IORD(MAX(NORD,1))
         LFAL = II .EQ. IFAL(MAX(NFAL,1))
         LTRO = LFAL .AND. IFOP(MAX(NFAL,1)).GT.2
         LSRI = LFAL .AND. IFOP(MAX(NFAL,1)).EQ.2
         LWL  = II .EQ. IWL (MAX(NWL, 1))
         LREV = II .EQ. IREV(MAX(NREV,1))
         LLAN = II .EQ. ILAN(MAX(NLAN,1))
         LRLT = II .EQ. IRLT(MAX(NRLT,1))
         LEIM = II .EQ. IEIM(MAX(NEIM,1))
         LTDE = II .EQ. ITDE(MAX(NTDE,1))
         LJAN = II .EQ. IJAN(MAX(NJAN,1))
         LFT1 = II .EQ. IFT1(MAX(NFT1,1))
         LEXC = II .EQ. IEXC(MAX(NEXC,1))
         LMOM = II .EQ. IMOM(MAX(NMOM,1))
         LXSM = II .EQ. IXSM(MAX(NXSM,1))
C
         KEY = ' '
         KEY = CKCHUP(SUB(N), 4)
         ILEN = CKLSCH(SUB(N))
C
         IF (KEY(1:3) .EQ. 'DUP') THEN
            IDUP(II) = -1
            WRITE (LOUT, 4000)
            GO TO 500
C
         ELSEIF (KEY(1:3) .EQ. 'MOM') THEN
C
C           electron momentum transfer collision frequency
            CALL CKNMOM (LOUT, KELECT, IDIM, II, LMOM, NMOM,
     1                   IMOM, KERR)
            GO TO 500
C
         ELSEIF (KEY(1:3) .EQ. 'XSM') THEN
C
C           momentum transfer cross section for ions
            CALL CKNXSM (LOUT, IDIM, II, NKION, LXSM, NXSM, IXSM, KERR)
            GO TO 500
        ENDIF
C
C        other options need parameters
         CALL CKDLIM (SUB(N), '/', I1, I2)
         IF (I2 .LE. I1) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A,A)')
     1         'Error...',KEY,' requires /-delimited arguments...'
            GO TO 500
         ENDIF
C
C        SUB is "KEY / RSTR /" , where RSTR contains any number
C        of parameters or arguments
         RSTR = ' '
         RSTR = SUB(N)(I1+1:I2-1)
C
         IF (KEY(1:3).EQ.'LOW' .OR. KEY(1:3).EQ.'HIG' .OR.
     1       KEY(1:3).EQ.'TRO' .OR. KEY(1:3).EQ.'SRI') THEN
C
C           pres dependency reaction
            CALL CKNFAL (LOUT, IDIM, II, KEY, SUB(N), RSTR, NFAL, IFAL,
     1                   IFOP, IFLO, NFAR, PFAL, LFAL, LLAN, LRLT,
     2                   LREV, LTRO, LSRI, KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'REV') THEN
C
C           reverse Arrhenius parameters
            CALL CKNREV (LOUT, IDIM, II, RSTR, LFAL, LREV, NSPEC, NPAR,
     1                   NREV, IREV, RPAR, KERR)
C
         ELSEIF (KEY(1:2) .EQ. 'LT') THEN
C
C           Landau-Teller parameters
            CALL CKNLAN (LOUT, IDIM, II, RSTR, LLAN, LFAL, NLAN, NLAR,
     1                   ILAN, PLAN, KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'RLT') THEN
C
C           reverse Landau-Teller parameters
            CALL CKNRLT (LOUT, IDIM, II, RSTR, LFAL, LRLT, NSPEC, NLAR,
     1                   NRLT, IRLT, RLAN, KERR)
C
         ELSEIF (KEY(1:2) .EQ. 'HV') THEN
C
C           radiation wavelength enhancement factor
            CALL CKNWL (LOUT, IDIM, II, RSTR, LWL, NWL, IWL, WL, KERR)
C
         ELSEIF (KEY.EQ.'FORD' .OR. KEY.EQ.'RORD') THEN
C
C           change-order species and paramter
            CALL CKNORD (LOUT, KDIM, IDIM, KEY, RSTR, II, KK, KNAME,
     1                   NORD, MAXORD, IORD, KORD, RORD, MAXSP, NUNK,
     2                   NU, NSPEC, NREAC, LORD, LRNU, NRNU, RNU,
     3                   KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'EIM') THEN
C
C           electron impact reaction
            CALL CKNEIM (LOUT, IDIM, II, RSTR, LEIM, LTHB, NEIM, IEIM,
     $                   IEIMT, KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'TDE') THEN
C
C           non-thermal-equilibrium temp dependency
            CALL CKNTDE (LOUT, KDIM, IDIM, II, RSTR, LTDE, LTHB, KNAME,
     1                   KK, NTDE, ITDE, ITDK, KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'JAN') THEN
C
C           Jannev, Langer, Evans & Post type reaction
            CALL CKNJAN (LOUT, IDIM, II, RSTR, LJAN, NJAN, IJAN, NJAR,
     1                   PJAN, KERR)
C
         ELSEIF (KEY(1:4) .EQ. 'FIT1') THEN
C
C           miscellaneous fit #1: k = A * T^B * exp[SUM(Vn/T^n)]
            CALL CKNFT1 (LOUT, IDIM, II, RSTR, LFT1, NFT1, IFT1, NF1R,
     1                   PFT1, KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'EXC') THEN
C
C           excitation-only reaction description (for energy loss)
            CALL CKNEXC (LOUT, IDIM, II, RSTR, LEXC, NEXC, IEXC, PEXC,
     1                   KERR)
C
         ELSEIF (KEY(1:3) .EQ. 'UNI') THEN
C
C           units conversion for reaction parameters is required:
C           get parameters from user's input string
            IUNITS = ' '
            CALL CKUNIT (IUNITS, AXUNIT, EXUNIT, IUNITS)
            WRITE (LOUT, 300) IUNITS
C
         ELSE
C
C           only thing left at this point is enhanced third bodies,
C           in which case instead of just KEY, need full species name
            CALL CKNTHB (LOUT, II, KDIM, IDIM, SUB(N)(1:I1-1), SUB(N),
     1                   RSTR, KNAME, KK, LTHB, LEIM, LTDE, NTHB, ITHB,
     2                   NTBS, NKTB, MAXTB, AIK, KERR)
         ENDIF
  500 CONTINUE
C
  300 FORMAT (6X,'Units for this reaction are: ', A)
 4000 FORMAT (6X,'Declared duplicate reaction...')
C
C     end of SUBROUTINE CKAUXL
      RETURN
      END
C                                                                      C
      SUBROUTINE CKBAL (LOUT, MDIM, KDIM, IDIM, MAXSP, NU, NUNK, KCHRG,
     1                  IRNU, RNU, KNCF, CKMIN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKBAL
C  checks a reaction to ensure that its reactants' mass and
C  electronic charge equals its product's mass and electronic charge.
C
C  Arguments:
C  LOUT      - Integer scalar, formatted output file unit number.
C  II        - Integer scalar, reaction count, and index of current
C              reaction.
C  MDIM      - Integer scalar, maximum number of elements.
C  MM        - Integer scalar, element count.
C  MAXSP     - Integer scalar, maximum number of species allowed
C               in a reaction.
C  NU(*)     - Integer array, species stoichiometric coefficients for
C              this reaction.
C  NUNK(*)   - Integer array, species indices for this reaction.
C  KCHRG(*)  - Integer array, species electronic charges.
C  NRNU      - Integer scalar, real stoichiometry reaction count.
C  IRNU(*)   - Integer array, real stoichiometry reaction indices.
C  RNU(*,*)  - Real matrix, real stoichiometric coefficients.
C  KNCF(*,*) - Integer matrix, species elemental composition.
C  CKMIN     - Real scalar, an error tolerance for balancing.
C  KERR      - Logical error flag.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN,
     2               NFT1, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     3               MELECT, KELECT, NKION
C
      DIMENSION NU(MAXSP), NUNK(MAXSP), IRNU(IDIM), RNU(MAXSP,IDIM),
     1          KNCF(MDIM,KDIM), KCHRG(KDIM)
      LOGICAL KERR, LRNU, IERR
C
      LRNU = IRNU(MAX(NRNU,1)) .EQ. II
C
C     element/mass balance of reaction
      IERR = .FALSE.
      DO 60 M = 1, MM
         SUMM = 0.0
         DO 50 N = 1, MAXSP
C           index number for this species
            NK = NUNK(N)
            IF (NK .NE. 0) THEN
C              stoichiometric coefficient for this species
               IF (LRNU) THEN
                  COEF = RNU(N, NRNU)
               ELSE
                  COEF = NU(N)
               ENDIF
               SUMM = SUMM + COEF * KNCF(M, NK)
            ENDIF
   50    CONTINUE
         IF (ABS(SUMM) .GT. CKMIN) IERR = .TRUE.
   60 CONTINUE
      IF (IERR) THEN  
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...reaction does not balance in elements/mass...'
      ENDIF
C
C     charge balance of reaction
      SUMCH = 0.0
      DO 70 N = 1, MAXSP
C        index number for this species
         NK = NUNK(N)
         IF (NK .NE. 0) THEN
C           stoichiometric coefficient for this species
            IF (LRNU) THEN
               COEF = RNU(N,NRNU)
            ELSE
               COEF = NU(N)
            ENDIF
C           sum of charge
            SUMCH = SUMCH + COEF * KCHRG(NK)
         ENDIF
   70 CONTINUE
      IF (ABS(SUMCH) .GT. CKMIN) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...reaction does not balance in electronic charge...'
      ENDIF
C
C     end of SUBROUTINE CKBAL
      RETURN
      END
C                                                                      C
      SUBROUTINE CKCHAR (SUB, NSUB, NDIM, STRAY, RAY, NN, KERR, LOUT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCHAR
C  copies character strings in an array into another
C  character string array, or converts the string to a real value and
C  copies into a real array; used to store element and atomic weight
C  data, species names, etc.  CKCHAR does not allow duplicate character-
C  strings to be entered into the array, and CKCHAR does not allow
C  the array to be over-filled.
C
C  Arguments:
C  SUB(*)   - Character string array.
C  NSUB     - Integer scalar, number of character strings in SUB.
C  NDIM     - Integer scalar, dimension of STRAY and RAY.
C  STRAY(*) - Character string array.
C  RAY(*)   - Real array.
C  NN       - Integer scalar, count of non-blank STRAY.
C  KERR     - Logical error flag.
C  LOUT     - Integer scalar, formatted output file unit number.
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
      DIMENSION RAY(NDIM)
      CHARACTER SUB(NSUB)*(*), STRAY(NDIM)*(*), ISTR*80, CKCHUP*4
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH, CKCHUP
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF ( CKCHUP(SUB(N), 3) .EQ. 'END') RETURN
         ISTR = ' '
         I1 = INDEX(SUB(N),'/')
         IF (I1 .EQ. 1) THEN
            KERR = .TRUE.
            WRITE (LOUT, 130) SUB(N)(1:CKLSCH(SUB(N)))
         ELSE
            IF (I1 .LE. 0) THEN
               ISTR = SUB(N)
            ELSE
               ISTR = SUB(N)(1:I1-1)
            ENDIF
            CALL CKCOMP (ISTR, STRAY, NN, INUM)
C
            IF (INUM .GT. 0) THEN
               WRITE (LOUT, 100) SUB(N)(1:CKLSCH(SUB(N)))
            ELSE
               IF (NN .LT. NDIM) THEN
                  IF (ISTR(ILEN+1:) .NE. ' ') THEN
                     WRITE (LOUT, 120) SUB(N)(1:CKLSCH(SUB(N)))
                     KERR = .TRUE.
                  ELSE
                     NN = NN + 1
                     STRAY(NN) = ' '
                     STRAY(NN) = ISTR(1:ILEN)
                     IF (I1 .GT. 0) THEN
                        I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                        ISTR = ' '
                        ISTR = SUB(N)(I1+1:I2-1)
                        CALL CKPARR (ISTR, 1, 1, RAY(NN),NVAL,IER,LOUT)
                        KERR = KERR .OR. (IER.NE.0)
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (LOUT, 110) SUB(N)(1:CKLSCH(SUB(N)))
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
      RETURN
C
  100 FORMAT (6X,'Warning...duplicate array element ignored...',A)
  110 FORMAT (6X,'Error...character array size too small for  ...',A)
  120 FORMAT (6X,'Error...character array element name too long...',A)
  130 FORMAT (6X,'Error...misplaced value...',A)
      END
C                                                                      C
      SUBROUTINE CKDUP (LOUT, IDIM, II, MAXSP, NSPEC, NREAC, NU, NUNK,
     1                  NFAL, IFAL, KFAL, IDUP, NTHB, ITHB, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDUP
C  checks reaction II against the (II-1) previous reactions for
C  duplication.
C
C  Arguments:
C  II      - Integer scalar, reaction count and index of
C            current reaction.
C  MAXSP   - Integer scalar, maximum number of species allowed in a
C            reaction.
C  NSPEC(*)- Integer array, reactions reactants+species count.
C  NREAC(*)- Integer array, reactions reactants only count.
C  NU(*,*) - Integer matrix, reactions stoichiometric coefficients.
C  NUNK(*,*)-Integer matrix, reactions species indices.
C  NFAL    - Integer scalar, count of pres dependency reactions.
C  IFAL(*) - Integer array, pressure-dependent reaction indices.
C  KFAL(*) - Integer array, third-body species indices.
C  IDUP(*) - Integer array, flag to allow duplication of reactions.
C  NTHB    - Integer scalar, count of 3rd-body type reactions.
C  ITHB(*) - Integer array, third-body reaction indices.
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
      DIMENSION NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     1          NUNK(MAXSP,IDIM), IFAL(IDIM), KFAL(IDIM),
     2          IDUP(IDIM), ITHB(IDIM)
      LOGICAL KERR
C
      ISAME = 0
      KHALF = MAXSP/2
C
C     II is the current reaction
C     NRI is II's reactant species count
C     NPI is II's product species count
C     LFAL is zero, or II's fall-off reaction index
C     LTHB is zero, or II's 3-rd body reaction index
      NRI = NREAC(II)
      NPI = ABS(NSPEC(II)) - NREAC(II)
      CALL CKINUM (II, IFAL, NFAL, LFAL)
      CALL CKINUM (II, ITHB, NTHB, LTHB)
C
C     check II against the previous reactions
      DO 500 J = 1, II-1
C
C        JFAL is zero, or J's fall-off reaction index
         CALL CKINUM (J, IFAL, NFAL, JFAL)
C        II is not like J if its pressure-dependency is different
         IF (JFAL.NE.0 .AND. LFAL.NE.0) THEN
            IF (KFAL(JFAL) .NE. KFAL(LFAL)) GO TO 500
         ENDIF
C
C        JTHB is zero, or J's 3rd-body reaction index
C        II is not like J if its 3rd-body presence is different
         CALL CKINUM (J, ITHB, NTHB, JTHB)
         IF ((JTHB.NE.0 .AND. LTHB.EQ.0) .OR.
     1       (JTHB.EQ.0 .AND. LTHB.NE.0)) GO TO 500
C
C        NRJ is J's reactant species count
C        NPJ is J's product species count
         NRJ = NREAC(J)
         NPJ = ABS(NSPEC(J)) - NREAC(J)
C
C        II is not like J if a species count is different
         IF (NRI.NE.NRJ .OR. NPI.NE.NPJ) GO TO 100
C
C        does J have same reactants as II
         DO 20 N = 1, NRI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL CKINUM (KI, NUNK(1,J), NRI, KJ)
C              II is not like J if KI reactant not in J,
C              or if it's coefficient is different
               IF (KJ .LE. 0) GO TO 100
               IF (NU(KJ,J) .NE. NU(N,II)) GO TO 100
            ENDIF
   20    CONTINUE
C
C        reactants are same, does J have same products as II
         DO 25 N = KHALF+1, KHALF + NPI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL CKINUM (KI, NUNK(KHALF+1,J), NPI, KJ)
               KJ = KHALF + KJ
C              II is not like J if KI product not in J,
C              or if it's coefficient is different
               IF (KJ.LE.KHALF .OR. NU(KJ,J).NE.NU(N,II)) GO TO 100
            ENDIF
   25    CONTINUE
C
C        same products, reactants, coefficients, pres dependency, and
C        third-body relationship
C
         ISAME = J
         GO TO 600
C
  100    CONTINUE
C
C        II is different from J in forward direction; check reverse,
C        if II has same number of reactants as J has products,
C        and same number of products as J has reactants
C
         IF (NPI.NE.NRJ .OR. NPJ.NE.NRI) GO TO 500
C
C        check I reactants against J products
         DO 30 N = 1, NRI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL CKINUM (KI, NUNK(KHALF+1,J), NPJ, KJ)
               KJ = KHALF + KJ
C              II is not like J if KI reactant not in J products,
C              or if it's coefficient is different
               IF (KJ.LE.KHALF .OR. NU(N,II).NE.-NU(KJ,J)) GO TO 500
            ENDIF
   30    CONTINUE
C
C        check I products against J reactants
         DO 35 N = KHALF+1, KHALF + NPI
            KI = NUNK(N,II)
            IF (KI .NE. 0) THEN
               CALL CKINUM (KI, NUNK(1,J), NRJ, KJ)
C              II is not like J if KI product not in J reactants,
C              or if it's coefficient is different
               IF (KJ .LE. 0) GO TO 500
               IF (-NU(N,II) .NE. NU(KJ,J)) GO TO 500
            ENDIF
   35    CONTINUE
C
C        J products same as II reactants, and
C        J reactants same as II products
         IF (NSPEC(J).LT.0 .AND. NSPEC(II).LT.0) THEN
C           only OK if J and II are both irreversible
         ELSE
            ISAME = J
            GO TO 600
         ENDIF
C
  500 CONTINUE
C
  600 CONTINUE
C
      IF (ISAME .EQ. 0) RETURN
C
C     Reaction #ISAME is a duplicate
      IF (IDUP(ISAME).NE.0 .AND. IDUP(II).NE.0) THEN
C        Reaction #ISAME is a legal duplicate
         IDUP(ISAME) = ABS(IDUP(ISAME))
         IDUP(II) = ABS(IDUP(II))
         RETURN
      ENDIF
C
      KERR = .TRUE.
      WRITE (LOUT, 1050) ISAME
 1050 FORMAT (6X,'Error...undeclared duplicate to reaction number ',I3)
C
C     end of SUBROUTINE CKCHAR
      RETURN
      END
C
      SUBROUTINE CKFORM
     1   (LINCK, LOUT, KERR, CKMIN, MDIM, KDIM, IDIM, NFIT, NTR, MAXSP,
     2    MAXTB, MAXTP, NCP, NPAR, NLAR, NFAR, NJAR, MAXORD, NF1R,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, T, A, KION,
     4    NSPEC, NREAC, NU, NUNK, PAR, IREV, RPAR, IFAL, IFOP, IFLO,
     5    KFAL, PFAL, ITHB, NTBS, NKTB, AIK, ILAN, PLAN, IRLT, RLAN,
     6    IWL, WL, IEIM, IEIMT, IJAN, PJAN, IFT1, PFT1, IEXC, PEXC,
     7    IMOM, KMOM, IXSM, IXSI, IXSK, ITDE, ITDK, IRNU, RNU, IORD,
     8    KORD, RORD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKFORM
C  writes mechanism information into a formatted linkfile.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
      LOGICAL       KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C     Character arrays
      CHARACTER*16  ENAME(MDIM), KNAME(KDIM), CFMT, IFMT, LFMT, RFMT,
     1              PRVERS, PREC, FILVER
C     Version and Precision common blocks
      COMMON /CKVERS/ PRVERS, PREC, FILVER
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KCHRG(KDIM), NT(KDIM), KPHSE(KDIM),
     1          KION(KDIM), NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     2          NUNK(MAXSP,IDIM), IREV(IDIM), IFAL(IDIM), IFOP(IDIM),
     3          IFLO(IDIM), KFAL(IDIM), ITHB(IDIM), NTBS(IDIM),
     4          NKTB(MAXTB,IDIM), ILAN(IDIM), IRLT(IDIM), IWL(IDIM),
     5          IEIM(IDIM), IEIMT(IDIM), IJAN(IDIM), IFT1(IDIM),
     6          IEXC(IDIM), IMOM(IDIM), KMOM(IDIM), IXSM(IDIM),
     7          IXSI(IDIM), IXSK(IDIM), ITDE(IDIM), ITDK(IDIM),
     8          IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM)
C    Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), T(MAXTP,KDIM),  A(NFIT,NTR,KDIM),
     1          PAR(NPAR,IDIM), RPAR(NPAR,IDIM), PFAL(NFAR,IDIM),
     2          AIK(MAXTB,IDIM), PLAN(NLAR,IDIM), RLAN(NLAR,IDIM),
     3          WL(IDIM), PJAN(NJAR,IDIM), PFT1(NF1R,IDIM),
     4          PEXC(IDIM), RNU(MAXSP,IDIM), RORD(MAXORD,IDIM)
C
C======================================================================
C            Write Out Linkfile in ascii format
C
C      Character data - A16
C      Integer data   - I12
C      Logical data   - L8
C      real*8 data    - 1PE24.16
C
C======================================================================
C
C     Open linkfile as an ascii file,
C     rewriting over an existing file
C     Linkile version (should be 1.0)
      WRITE (LINCK, CFMT, ERR=5003) FILVER
C
C     Ckinterp program version
      WRITE (LINCK, CFMT, ERR=5003) PRVERS
C
C     Precision
      WRITE (LINCK, CFMT, ERR=5003) PREC
C
C     Whether ckinterp successfully completed or not
      WRITE (LINCK, LFMT, ERR=5001)  KERR
C
      IF (KERR) THEN
         WRITE (LOUT, '(//A,A)')
     1   ' WARNING...THERE IS AN ERROR IN THE LINKING FILE,',
     2   '           DUE TO MECHANISM INPUT ERRORS.'
         RETURN
      ENDIF
C
C     Work space sizes
      WRITE (LINCK, IFMT)  LENICK, LENRCK, LENCCK
C
C     Fixed parameterization sizes:
C     MAXSP = Maximum number of reactants and products in the
C             stiochiometric matrix for each reaction
C     MAXTB = Maximum number of third body effects defined
C             for each reaction
C     MAXTP = Maximum number of temperature regions + 1
C             used in fits to thermodynamic functions
C     NCP   = Number of coefficients in the polynomial fit
C             to constant pressure heat capacity
C     NPAR  = Maximum number of parameters in the Arrhenious
C             coefficient parameterization of a reaction
C     NLAR  = Max numbers of parameters in a Landau-Teller rxn
C     NFAR  = Maximum number of parameters in the
C             parameterization of pres dependency effects
C     NJAR  = Maximum number of parameters in the
C             parameterization of Jannev, Langer,
C             Evans & Post reactions
C     MAXORD = Maximum number of species whose order in a
C              reaction may be variable
C     NFIR  = Maximum number of parameters in the
C             parameterization of fit#1 reactions
C
      MXTB = 0
      DO 111 N = 1, NTHB
         MXTB = MAX (MXTB, NTBS(N))
  111 CONTINUE
      MXTK = 0
      DO 112 K = 1, KK
         MXTK = MAX (MXTK, NT(K))
  112 CONTINUE
      WRITE (LINCK, IFMT, ERR=5001) MAXSP, MXTB, MXTK, NCP, NPAR,
     1                             NLAR, NFAR, NJAR, MAXORD, NF1R
C
C     Problem specifications:
C     MM    = element count
C     KK    = species count
C     II    = gas-phase reaction count
C     NREV  = count of rxns with reverse rate constants supplied
C     NFAL  = count of rxns with pres dependency parameters
C     NTHB  = count of rxns with third body coefficients
C     NLAN  = count of rxns with Landau-Teller coefficients
C     NRLT  = count of rxns with reverse L-T coefficients
C     NWL   = count of rxns with radiation wavelenghts
C     NCHRG = count of species with KCHRG <> 0
C     NEIM  = count of rxns with electron impact
C     NJAN  = count of Jannev, Langer, Evans & Post Rxns
C     NFT1  = count of rxns using fit #1
C     NEXC  = count of excitation reactions
C     NMOM  = count of rxns with momentum transfer
C     NXSM  = count of rxns with MT XS for ions
C     NTDE  = count of rxns with non-thermal T dependency
C     NRNU  = count of rxns with real stoich.
C     NORD  = count of rxns with order dependency different
C             than the stoichiometry
C     KELECT= ?
C     NKION =  ?
C
      WRITE (LINCK, IFMT, ERR=5001) MM, KK, II, NREV, NFAL, NTHB,
     1                             NLAN, NRLT, NWL, NCHRG, NEIM,
     2                             NJAN, NFT1, NEXC, NMOM, NXSM,
     3                             NTDE, NRNU, NORD, KELECT, NKION
C
C     Elemental and charge balance tolerance parameter
      WRITE (LINCK, RFMT, ERR=5002) CKMIN
C
C     Element names
      WRITE (LINCK, CFMT, ERR=5003) (ENAME(M), M = 1, MM)
C
C     Elemental atomic weights
      WRITE (LINCK, RFMT, ERR=5002) (AWT(M), M = 1, MM)
C
C     Species names
      WRITE (LINCK, CFMT, ERR=5003) (KNAME(K), K=1,KK)
C
C     Species molecular weights
      WRITE (LINCK, RFMT, ERR=5002) (WT(K), K = 1, KK)
C
C     Species elemental compositions
      WRITE (LINCK, IFMT, ERR=5001) ((KNCF(M,K), M = 1, MM), K = 1, KK)
C
C     Species electronic charges
      WRITE (LINCK, IFMT, ERR=5001)  (KCHRG(K), K = 1, KK)
C
C     Number of temperatures regions used to fit thermo for each species
      WRITE (LINCK, IFMT, ERR=5001)  (NT(K), K = 1, KK)
C
C     An integer representing physical state of each species
      WRITE (LINCK, IFMT, ERR=5001)  (KPHSE(K), K = 1, KK)
C
C     Array of temperatures used to specify boundaries of the
C     temperature regions used in fits to thermodynamic functions
      WRITE (LINCK, RFMT, ERR=5002)  ((T(L,K), L = 1, MXTK), K = 1, KK)
C
C     Fits to thermodynamic functions for each species
C     NTR=MAXTP-1, NFIT=NCP+2
      WRITE (LINCK, RFMT, ERR=5002)
     1            (((A(N,L,K), N = 1, NFIT), L = 1, MXTK-1), K = 1,KK)
C
      IF (NKION .GT. 0) THEN
         WRITE (LINCK, IFMT, ERR=5001) NKION
         WRITE (LINCK, IFMT, ERR=5001) (KION(N), N = 1, NKION)
      ENDIF
C
      IF (II .GT. 0) THEN
C        GAS PHASE REACTION OUTPUT
C        count of species defined as reactants and products
         WRITE (LINCK, IFMT, ERR=5001) (NSPEC(I), I = 1, II)
C        count of reactants defined for each reaction
C        A negative number flags an irreversible reaction also
         WRITE (LINCK, IFMT, ERR=5001) (NREAC(I), I = 1, II)
C        stoichiometric information for each reaction
         WRITE (LINCK, IFMT, ERR=5001)
     1          ((NU(M,I), NUNK(M,I), M = 1, MAXSP), I = 1, II)
C        forward Arrhenius coefficients for each reaction
         WRITE (LINCK, RFMT, ERR=5002) ((PAR(N,I), N = 1, NPAR), I=1,II)
      ENDIF
C
C     SPECIAL REACTIONS:
C
      IF (NREV .GT. 0) THEN
C        Reactions with overriden reverse rate constants
         WRITE (LINCK, IFMT, ERR=5001) NREV
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IREV(N), N = 1, NREV)
C        Arrhenious reverse rate constant coefficients
         WRITE (LINCK, RFMT, ERR=5002) ((RPAR(L,N), L=1,NPAR), N=1,NREV)
      ENDIF
C
      IF (NFAL .GT. 0) THEN
C        Reactions with pres dependency parameterizations
         WRITE (LINCK, IFMT, ERR=5001) NFAL, NFAR
C        reaction indices, pres dependency parameterization type,
C        pres dependency flow type, third-body species if any
         WRITE (LINCK, IFMT, ERR=5001)
     1          (IFAL(N), IFOP(N), IFLO(N), KFAL(N), N = 1, NFAL)
C        pres dependency parameters of those reactions
         WRITE (LINCK, RFMT, ERR=5002) ((PFAL(L,N),L=1,NFAR), N=1,NFAL)
      ENDIF
C
      IF (NTHB .GT. 0) THEN
C        Reactions with third body effects
         WRITE (LINCK, IFMT, ERR=5001) NTHB, MXTB
C        reaction ID, ITHB, for each third-body effect
C        number of third body effects for this reaction, NTBS
         WRITE (LINCK, IFMT, ERR=5001) (ITHB(N), NTBS(N), N = 1, NTHB)
C        thirdy body species IDs
         WRITE (LINCK, IFMT, ERR=5001) ((NKTB(M,N),M=1,MXTB), N=1,NTHB)
C        third body enhancement coefficient
         WRITE (LINCK, RFMT, ERR=5002) ((AIK(M,N),M=1,MXTB), N=1,NTHB)
      ENDIF
C
      IF (NLAN .GT. 0) THEN
C        Reactions with Landau reaction rate constants
         WRITE (LINCK, IFMT, ERR=5001) NLAN, NLAR
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (ILAN(N), N = 1, NLAN)
C        parameterization of rate constants for those reactions
        WRITE (LINCK, RFMT, ERR=5002) ((PLAN(L,N),L=1,NLAR), N=1,NLAN)
      ENDIF
C
      IF (NRLT .GT. 0) THEN
C        Reactions with specified reverse rate parameters given by
C        Landau-Teller parameterization
         WRITE (LINCK, IFMT, ERR=5001) NRLT, NLAR
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IRLT(N),  N = 1, NRLT)
C        Landau-Teller reverse paramerization coefficients
         WRITE (LINCK, RFMT, ERR=5002) ((RLAN(L,N),L=1,NLAR), N=1,NRLT)
      ENDIF
C
      IF (NWL .GT. 0) THEN
C        Reactions with radiation wavelength specifications
         WRITE (LINCK, IFMT, ERR=5001) NWL
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IWL(N), N = 1, NWL)
C        parameterization of the physical effect of the wavelength
         WRITE (LINCK, RFMT, ERR=5002) (WL(N), N = 1, NWL)
      ENDIF
C
      IF (NEIM .GT. 0) THEN
C        Reactions with electron impact effects
         WRITE (LINCK, IFMT, ERR=5001) NEIM
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IEIM(N), N = 1, NEIM)
C        temp dependency flags
         WRITE (LINCK, IFMT, ERR=5001) (IEIMT(N), N = 1, NEIM)
      ENDIF
C
      IF (NJAN .GT. 0) THEN
C        Reactions with Jannev, Langer, Evans & Post rates
         WRITE (LINCK, IFMT, ERR=5001) NJAN, NJAR
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IJAN(N), N = 1, NJAN)
C        parameterizations of these reaction types
         WRITE (LINCK, RFMT, ERR=5002) ((PJAN(L,N),L=1,NJAR), N=1,NJAN)
      ENDIF
C
      IF (NFT1 .GT. 0) THEN
C        Reactions that use fit #1
         WRITE (LINCK, IFMT, ERR=5001) NFT1, NF1R
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IFT1(N), N = 1, NFT1)
C        parameterizations of these reaction types
         WRITE (LINCK, RFMT, ERR=5002) ((PFT1(L,N), L=1,NF1R), N=1,NFT1)
      ENDIF
C
      IF (NEXC .GT. 0) THEN
C        Excitation reactions
         WRITE (LINCK, IFMT, ERR=5001) NEXC
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IEXC(N), N = 1, NEXC)
C        energy loss per collision event in units of eV
         WRITE (LINCK, RFMT, ERR=5002) (PEXC(N), N = 1, NEXC)
      ENDIF
C
      IF (NMOM .GT. 0) THEN
C        Electron momentum reactions
         WRITE (LINCK, IFMT, ERR=5001) NMOM
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IMOM(N), N = 1, NMOM)
C        heavy species for each reaction
         WRITE (LINCK, IFMT, ERR=5001) (KMOM(N), N = 1, NMOM)
      ENDIF
C
      IF (NXSM .GT. 0) THEN
C        Ion momentum-transfer cross-section reactions
         WRITE (LINCK, IFMT, ERR=5001) NXSM
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IXSM(N), N = 1, NXSM)
C        the ion indices
         WRITE (LINCK, IFMT, ERR=5001) (IXSI(N), N = 1, NXSM)
C        the ion's collision partner's indices
         WRITE (LINCK, IFMT, ERR=5001) (IXSK(N), N = 1, NXSM)
      ENDIF
C
      IF (NTDE .GT. 0) THEN
C        non-thermal temperature dependent reactions
         WRITE (LINCK, IFMT, ERR=5001) NTDE
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (ITDE(N), N = 1, NTDE)
C        species number for which to use its temperature
         WRITE (LINCK, IFMT, ERR=5001) (ITDK(N),  N = 1, NTDE)
      ENDIF
C
      IF (NRNU .GT. 0) THEN
C        Real stoichiometry reactions
         WRITE (LINCK, IFMT, ERR=5001) NRNU
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IRNU(N), N = 1, NRNU)
C        stiochiometric coefficient array
         WRITE (LINCK, RFMT, ERR=5002) ((RNU(L,N),L=1,MAXSP),N=1,NRNU)
      ENDIF
C
      IF (NORD .GT. 0) THEN
C        Modified species order reactions
         WRITE (LINCK, IFMT, ERR=5001) NORD, MAXORD
C        reaction indices
         WRITE (LINCK, IFMT, ERR=5001) (IORD(N), N = 1, NORD)
C        order dependency species indices
         WRITE (LINCK, IFMT, ERR=5001) ((KORD(L,N),L=1,MAXORD),N=1,NORD)
C        species order values
         WRITE (LINCK, RFMT, ERR=5002) ((RORD(L,N),L=1,MAXORD),N=1,NORD)
      ENDIF
C
C CLEANUP
C
      WRITE (LOUT, '(///A,/A,A,A)')
     1' NO ERRORS FOUND ON INPUT: ',
     2' ASCII Vers. ',FILVER(1:CKLSCH(FILVER)),
     3' CHEMKIN linkfile chem.asc written.'
      RETURN
C
C ERROR HANDLING
C
5001  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write ascii integer data to chem.asc'
      RETURN
5002  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write ascii real data to chem.asc'
      RETURN
5003  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write ascii character data to chem.asc'
      RETURN
C
      END
C
      SUBROUTINE CKBIN
     1   (LINCK, LOUT, KERR, CKMIN, MDIM, KDIM, IDIM, NFIT, NTR, MAXSP,
     2    MAXTB, MAXTP, NCP, NPAR, NLAR, NFAR, NJAR, MAXORD, NF1R,
     3    ENAME, AWT, KNAME, WT, KNCF, KCHRG, NT, KPHSE, T, A, KION,
     4    NSPEC, NREAC, NU, NUNK, PAR, IREV, RPAR, IFAL, IFOP, IFLO,
     5    KFAL, PFAL, ITHB, NTBS, NKTB, AIK, ILAN, PLAN, IRLT, RLAN,
     6    IWL, WL, IEIM, IEIMT, IJAN, PJAN, IFT1, PFT1, IEXC, PEXC,
     7    IMOM, KMOM, IXSM, IXSI, IXSK, ITDE, ITDK, IRNU, RNU, IORD,
     8    KORD, RORD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKBIN
C  writes mechanism information into a binary linkfile.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Character arrays
      CHARACTER*16  ENAME(MDIM), KNAME(KDIM), PRVERS, PREC, FILVER
C      Version and Precision common blocks
      COMMON /CKVERS/ PRVERS, PREC, FILVER
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KCHRG(KDIM), NT(KDIM), KPHSE(KDIM),
     1          KION(KDIM), NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     2          NUNK(MAXSP,IDIM), IREV(IDIM), IFAL(IDIM), IFOP(IDIM),
     3          IFLO(IDIM), KFAL(IDIM), ITHB(IDIM), NTBS(IDIM),
     4          NKTB(MAXTB,IDIM), ILAN(IDIM), IRLT(IDIM), IWL(IDIM),
     5          IEIM(IDIM), IEIMT(IDIM), IJAN(IDIM), IFT1(IDIM),
     6          IEXC(IDIM), IMOM(IDIM), KMOM(IDIM), IXSM(IDIM),
     7          IXSI(IDIM), IXSK(IDIM), ITDE(IDIM), ITDK(IDIM),
     8          IRNU(IDIM), IORD(IDIM), KORD(MAXORD,IDIM)
C    Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), T(MAXTP,KDIM),  A(NFIT,NTR,KDIM),
     1          PAR(NPAR,IDIM), RPAR(NPAR,IDIM), PFAL(NFAR,IDIM),
     2          AIK(MAXTB,IDIM), PLAN(NLAR,IDIM), RLAN(NLAR,IDIM),
     3          WL(IDIM), PJAN(NJAR,IDIM), PFT1(NF1R,IDIM),
     4          PEXC(IDIM), RNU(MAXSP,IDIM), RORD(MAXORD,IDIM)
C
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C======================================================================
C            Write Out Linkfile in binary format
C
C
C======================================================================
C
C     Linkfile version (should be 1.0)
      WRITE (LINCK, ERR=5003) FILVER
C
C     Ckinterp program version
      WRITE (LINCK, ERR=5003) PRVERS
C
C     Precision
      WRITE (LINCK, ERR=5003) PREC
C
C     Whether ckinterp successfully completed or not
      WRITE (LINCK, ERR=5001)  KERR
C
      IF (KERR) THEN
         WRITE (LOUT, '(//A,/A)')
     1   ' WARNING...THERE IS AN ERROR IN THE LINKING FILE,',
     2   '           DUE TO MECHANISM INPUT ERRORS.'
         RETURN
      ENDIF
C
C     Work space sizes
      WRITE (LINCK, ERR=5001)  LENICK, LENRCK, LENCCK
C
C     Fixed parameterization sizes:
C     MAXSP = Maximum number of reactants and products in the
C             stiochiometric matrix for each reaction
C     MAXTB = Maximum number of third body effects defined
C             for each reaction
C     MAXTP = Maximum number of temperature regions + 1
C             used in fits to thermodynamic functions
C     NCP   = Number of coefficients in the polynomial fit
C             to the constant pressure heat capacity
C     NPAR  = Maximum number of parameters in the Arrhenious
C             coefficient parameterization of a reaction
C     NLAR  = Max numbers of parameters in a Landau-Teller rxn
C     NFAR  = Maximum number of parameters in the
C             parameterization of pres dependency effects
C     NJAR  = Maximum number of parameters in the
C             parameterization of Jannev, Langer, Evans
C             & Post reactions
C     MAXORD = Maximum number of species whose order in a
C              reaction may be variable
C     NFIR  = Maximum number of parameters in the
C             parameterization of fit#1 reactions
C
      MXTB = 0
      DO 111 N = 1, NTHB
         MXTB = MAX (MXTB, NTBS(N))
  111 CONTINUE
      MXTK = 0
      DO 112 K = 1, KK
         MXTK = MAX (MXTK, NT(K))
  112 CONTINUE
      WRITE (LINCK, ERR=5001) MAXSP, MXTB, MXTK, NCP, NPAR, NLAR,
     1                       NFAR, NJAR, MAXORD, NF1R
C
C     Problem specifications:
C     MM    = element count
C     KK    = species count
C     II    = gas-phase reaction count
C     NREV  = count of rxns with reverse rate constants supplied
C     NFAL  = count of rxns with pres dependency parameters
C     NTHB  = count of rxns with third body coefficients
C     NLAN  = count of rxns with Landau-Teller coefficients
C     NRLT  = count of rxns with reverse L-T coefficients
C     NWL   = count of rxns with radiation wavelenghts
C     NCHRG = count of species with KCHRG <> 0
C     NEIM  = count of rxns with electron impact
C     NJAN  = count of Jannev, Langer, Evans & Post Rxns
C     NFT1  = count of rxns using fit #1
C     NEXC  = count of excitation reactions
C     NMOM  = count of rxns with momentum transfer
C     NXSM  = count of rxns with MT XS for ions
C     NTDE  = count of rxns with non-thermal T dependency
C     NRNU  = count of rxns with real stoich.
C     NORD  = count of rxns with order dependency different
C             than the stoichiometry
C     KELECT= ?
C     NKION =  ?
C
      WRITE (LINCK, ERR=5001) MM, KK, II, NREV, NFAL, NTHB, NLAN,
     1                       NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2                       NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     3                       KELECT, NKION
C
C     elemental and charge balance tolerance parameter
      WRITE (LINCK, ERR=5002) CKMIN
C
C     element names
      WRITE (LINCK, ERR=5003) (ENAME(M), M = 1, MM)
C
C     elemental atomic weights
      WRITE (LINCK, ERR=5002) (AWT(M), M = 1, MM)
C
C     species names
      WRITE (LINCK, ERR=5003) (KNAME(K), K=1,KK)
C
C     species molecular weights
      WRITE (LINCK, ERR=5002) (WT(K), K = 1, KK)
C
C     species elemental compositions
      WRITE (LINCK, ERR=5001) ((KNCF(M,K), M = 1, MM), K = 1, KK)
C
C     species electronic charges
      WRITE (LINCK, ERR=5001)  (KCHRG(K), K = 1, KK)
C
C     number of temperatures regions used to fit thermodynamics
C     for each species
      WRITE (LINCK, ERR=5001)  (NT(K), K = 1, KK)
C
C     an integer representing physical state for each species
      WRITE (LINCK, ERR=5001)  (KPHSE(K), K = 1, KK)
C
C     array of temperatures used to specify boundaries of
C     the temperature regions used in fits to thermodynic functions
      WRITE (LINCK, ERR=5002)  ((T(L,K), L = 1, MXTK), K = 1, KK)
C
C     fits to thermodynamic functions for each species;
C     NTR=MAXTP-1, NFIT=NCP+2
      WRITE (LINCK, ERR=5002)
     1            (((A(N,L,K), N = 1, NFIT), L = 1, MXTK-1), K = 1,KK)
C
C     integer species indices for ionic species
      IF (NKION .GT. 0) THEN
        WRITE (LINCK, ERR=5001) NKION
        WRITE (LINCK, ERR=5001) (KION(N), N = 1, NKION)
      ENDIF
C
      IF (II .LE. 0) THEN
         WRITE (LOUT, '(/A,/A)')
     1   ' WARNING...NO REACTIONS FOUND, ',
     2   ' LINKFILE HAS NO REACTION INFORMATION ON IT.'
         GO TO 500
      ENDIF
C
C     GAS PHASE REACTION OUTPUT
C     count of species defined as reactants and products
      WRITE (LINCK, ERR=5001) (NSPEC(I), I = 1, II)
C     count of reactants defined for each reaction, also,
C     a negative number flags an irreversible reaction
      WRITE (LINCK, ERR=5001) (NREAC(I), I = 1, II)
C     stoichiometric information for each reaction
      WRITE (LINCK, ERR=5001)
     1       ((NU(M,I), NUNK(M,I), M = 1, MAXSP), I = 1, II)
C     forward Arrhenius coefficients for each reaction
      WRITE (LINCK, ERR=5002) ((PAR(N,I), N = 1, NPAR), I = 1, II)
C
C     SPECIAL REACTIONS:
C
      IF (NREV .GT. 0) THEN
C        reactions with explicit reverse rate constants
         WRITE (LINCK, ERR=5001) NREV
C        reaction indices
         WRITE (LINCK, ERR=5001) (IREV(N), N = 1, NREV)
C        Arrhenious reverse rate constant coefficients
         WRITE (LINCK, ERR=5002) ((RPAR(L,N), L = 1,NPAR), N = 1, NREV)
      ENDIF
C
      IF (NFAL .GT. 0) THEN
C        reactions with pres dependency parameterizations
         WRITE (LINCK, ERR=5001) NFAL, NFAR
C        reaction indices, pres dependency parameterization type,
C        pres dependency flow type, third-body species if any
         WRITE (LINCK, ERR=5001)
     $          (IFAL(N), IFOP(N), IFLO(N), KFAL(N), N = 1, NFAL)
C        pres dependency parameters of those reactions
         WRITE (LINCK, ERR=5002) ((PFAL(L,N), L = 1, NFAR), N = 1, NFAL)
      ENDIF
C
      IF (NTHB .GT. 0) THEN
C        Reactions with third body effects
         WRITE (LINCK, ERR=5001) NTHB, MXTB
C        reaction indices, count of third body effects
         WRITE (LINCK, ERR=5001) (ITHB(N), NTBS(N), N = 1, NTHB)
C        third-body species indices
         WRITE (LINCK, ERR=5001) ((NKTB(M,N), M = 1, MXTB), N = 1, NTHB)
C        third body enhancement coefficient for
C        species identified in NKTB(M,N)
         WRITE (LINCK, ERR=5002) ((AIK(M,N), M = 1, MXTB), N = 1, NTHB)
      ENDIF
C
C     Reactions whose rate constants are specified with Landau reaction
C     rate constants:
      IF (NLAN .GT. 0) THEN
         WRITE (LINCK, ERR=5001) NLAN, NLAR
C        reaction indices
         WRITE (LINCK, ERR=5001) (ILAN(N), N = 1, NLAN)
C        parameterization of the rate constants for those reactions
         WRITE (LINCK, ERR=5002) ((PLAN(L,N), L = 1, NLAR), N = 1, NLAN)
      ENDIF
C
      IF (NRLT .GT. 0) THEN
C        Reactions with specified reverse rate parameters given by
C        Landau-Teller parameterization
         WRITE (LINCK, ERR=5001) NRLT, NLAR
C        reaction indices
         WRITE (LINCK, ERR=5001) (IRLT(N),  N = 1, NRLT)
C        Landau-Teller reverse rate paramerization coefficients
         WRITE (LINCK, ERR=5002) ((RLAN(L,N), L = 1, NLAR), N = 1, NRLT)
      ENDIF
C
      IF (NWL .GT. 0) THEN
C        Reactions with radiation wavelength specifications
         WRITE (LINCK, ERR=5001) NWL
C        reaction indices
         WRITE (LINCK, ERR=5001) (IWL(N), N = 1, NWL)
C        parameterization of the physical effect of the wavelength
         WRITE (LINCK, ERR=5002) (WL(N), N = 1, NWL)
      ENDIF
C
      IF (NEIM .GT. 0) THEN
C        Reactions with electron impact effects
         WRITE (LINCK, ERR=5001) NEIM
C        reaction indices
         WRITE (LINCK, ERR=5001) (IEIM(N), N = 1, NEIM)
C        temp dependency flags
         WRITE (LINCK, ERR=5001) (IEIMT(N), N = 1, NEIM)
      ENDIF
C
      IF (NJAN .GT. 0) THEN
C        Reactions with Jannev, Langer, Evans & Post rates
         WRITE (LINCK, ERR=5001) NJAN, NJAR
C        reaction indices
         WRITE (LINCK, ERR=5001) (IJAN(N), N = 1, NJAN)
C        parameterizations of these reaction types
         WRITE (LINCK, ERR=5002) ((PJAN(L,N), L = 1, NJAR), N = 1, NJAN)
      ENDIF
C
      IF (NFT1 .GT. 0) THEN
C        Reactions that use fit #1
         WRITE (LINCK, ERR=5001) NFT1, NF1R
C        reaction indices
         WRITE (LINCK, ERR=5001) (IFT1(N), N = 1, NFT1)
C        parameterizations of these reaction types
         WRITE (LINCK, ERR=5002) ((PFT1(L,N), L = 1, NF1R), N = 1, NFT1)
      ENDIF
C
      IF (NEXC .GT. 0) THEN
C        Excitation reactions
         WRITE (LINCK, ERR=5001) NEXC
C        reaction indices
         WRITE (LINCK, ERR=5001) (IEXC(N), N = 1, NEXC)
C        energy loss per collision event in units of eV
         WRITE (LINCK, ERR=5002) (PEXC(N), N = 1, NEXC)
      ENDIF
C
      IF (NMOM .GT. 0) THEN
C        Electron momentum reactions
         WRITE (LINCK, ERR=5001) NMOM
C        reaction indices
         WRITE (LINCK, ERR=5001) (IMOM(N), N = 1, NMOM)
C        heavy species for each reaction
         WRITE (LINCK, ERR=5001) (KMOM(N), N = 1, NMOM)
      ENDIF
C
      IF (NXSM .GT. 0) THEN
C        Ion momentum-transfer cross-section reactions
         WRITE (LINCK, ERR=5001) NXSM
C        reaction indices
         WRITE (LINCK, ERR=5001) (IXSM(N), N = 1, NXSM)
C        the ion indices
         WRITE (LINCK, ERR=5001) (IXSI(N), N = 1, NXSM)
C        the ion's collision partner's indices
         WRITE (LINCK, ERR=5001) (IXSK(N), N = 1, NXSM)
      ENDIF
C
      IF (NTDE .GT. 0) THEN
C        Non-thermal temperature dependent reactions
         WRITE (LINCK, ERR=5001) NTDE
C        reaction indices
         WRITE (LINCK, ERR=5001) (ITDE(N), N = 1, NTDE)
C        species number for which to use its temperature
         WRITE (LINCK, ERR=5001) (ITDK(N),  N = 1, NTDE)
      ENDIF
C
      IF (NRNU .GT. 0) THEN
C        Real stoichiometry reactions
         WRITE (LINCK, ERR=5001) NRNU
C        reaction indices
         WRITE (LINCK, ERR=5001) (IRNU(N), N = 1, NRNU)
C        stiochiometric coefficient array for that reaction
         WRITE (LINCK, ERR=5002) ((RNU(L,N), L = 1, MAXSP), N = 1, NRNU)
      ENDIF
C
      IF (NORD .GT. 0) THEN
C        Modified species order reactions
         WRITE (LINCK, ERR=5001) NORD, MAXORD
C        reaction indices
         WRITE (LINCK, ERR=5001) (IORD(N), N = 1, NORD)
C        order dependency species indices
         WRITE (LINCK, ERR=5001) ((KORD(L,N), L = 1, MAXORD),N=1,NORD)
C        species order values
C        assigned by KORD(L,N)
         WRITE (LINCK, ERR=5002) ((RORD(L,N), L = 1, MAXORD),N=1,NORD)
      ENDIF
C
  500 CONTINUE
C
C CLEANUP
C
      WRITE (LOUT, '(///A,/A,A,A)')
     1' NO ERRORS FOUND ON INPUT, ',
     2' BINARY Vers. ',FILVER(1:CKLSCH(FILVER)),
     3' CHEMKIN linkfile chem.bin written.'
      RETURN
C
C ERROR HANDLING
C
5001  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write binary integer data to chem.bin'
      RETURN
5002  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write binary real data to chem.bin'
      RETURN
5003  CONTINUE
      WRITE (LOUT, '(///A)')
     $   'ERROR: Failure to write binary character data to chem.bin'
      RETURN
      END
C                                                                      C
      SUBROUTINE CKINUM (I, IARRAY, NI, IND)
C
C  START PROLOGUE
C
C  SUBROUTINE CKINUM
C  searches integer array IARRAY of length NI for the
C  first occurrence of integer I, and returns its index IND.
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
      DIMENSION IARRAY(NI)
C
      IND = 0
      DO 50 N = 1, NI
         IF (I .EQ. IARRAY(N)) THEN
            IND = N
            RETURN
         ENDIF
   50 CONTINUE
C
C     end of SUBROUTINE CKINUM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKISUB (LINE, SUB, NSUB)
C
C  START PROLOGUE
C
C  SUBROUTINE CKISUB
C  creates an array of character-strings given one character
C  string using blanks or tabs for delimiters;
C  slash(/)-delimited characters are not added as separate members
C  of the array, but rather, appended to the current array member.
C
C  Arguments:
C  LINE    - Character string.
C  SUB(*)  - Character string array.
C  NSUB    - Integer scalar, count of character strings added to
C            SUB.
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
      CHARACTER SUB(*)*(*), LINE*(*)
      INTEGER CKFRCH, CKLSCH, CKSLEN
      EXTERNAL CKFRCH, CKLSCH, CKSLEN
C
      NSUB = 0
      CALL CKDTAB (LINE)
      IF (CKSLEN(LINE) .LE. 0) RETURN
C
      NSTART = CKFRCH(LINE)
      ILEN = CKLSCH(LINE)
C
   10 CONTINUE
      ISTART = NSTART
      NSUB = NSUB + 1
      SUB(NSUB) = ' '
C
      ILAST = INDEX(LINE(ISTART:),' ') - 1
      IF (ILAST .GT. 0) THEN
         ILAST = ISTART + ILAST - 1
      ELSE
         ILAST = ILEN
      ENDIF
      SUB(NSUB) = LINE(ISTART:ILAST)
      IF (ILAST .EQ. ILEN) GO TO 50
C
      NSTART = ILAST + CKFRCH(LINE(ILAST+1:))
C
C     Does SUB have any slashes?
C
      I1 = INDEX(SUB(NSUB),'/')
      IF (I1 .LE. 0) THEN
         IF (LINE(NSTART:NSTART) .NE. '/') GO TO 10
         NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
         IND = INDEX(SUB(NSUB),' ')
         SUB(NSUB)(IND:) = LINE(NSTART:NEND)
         IF (NEND .EQ. ILEN) GO TO 50
         NSTART = NEND + CKFRCH(LINE(NEND+1:))
         GO TO 10
      ENDIF
C
C     Does SUB have 2 slashes?
C
      I2 = INDEX(SUB(NSUB)(I1+1:),'/')
      IF (I2 .GT. 0) GO TO 10
C
      NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
      IND = INDEX(SUB(NSUB),' ') + 1
      SUB(NSUB)(IND:) = LINE(NSTART:NEND)
      IF (NEND .LT. ILEN) THEN
         NSTART = NEND + CKFRCH(LINE(NEND+1:))
         GO TO 10
      ENDIF
C
   50 CONTINUE
C
C     end of SUBROUTINE CKISUB
      RETURN
      END
C                                                                      C
      SUBROUTINE CKKEY (LIN, LOUT, LTHRM, MDIM, KDIM, IDIM, NPAR, NCP,
     1                  NFIT, MAXTP, NTR, MAXSP, MAXTB, NLAR, NFAR,
     2                  MAXORD, KNAME, ENAME, AWT, KNCF, WT, KPHSE,
     3                  KCHRG, A, T, NT, NSPEC, NREAC, NU, NUNK, PAR,
     4                  IDUP, IREV, RPAR, ILAN, PLAN, IRLT, RLAN, IWL,
     5                  WL, IFAL, IFOP, IFLO, KFAL, PFAL, ITHB, NTBS,
     6                  AIK, NKTB, IRNU, RNU, IORD, KORD, RORD, KION,
     7                  NJAR, NF1R, IEIM, IEIMT, IJAN, PJAN, IFT1, PFT1,
     8                  IEXC, PEXC, IMOM, KMOM, IXSM, IXSK, IXSI, ITDE,
     9                  ITDK, ITHRM, CKMIN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKEY
C  is the main parsing routine for character-string input
C  of a chemical mechanism.  Depending on certain character-string
C  'flags', character strings are further directed to other
C  subroutines for parsing and storage of mechanism information.
C
C  Required element input -
C     The word 'ELEMENTS' followed by a list of atomic element names,
C     terminated by the word 'END'; an element name may be followed
C     by a slash-enclosed atomic weight value.         (CALL CKCHAR)
C  Required species input -
C     The word 'SPECIES' followed by a list of molecular species names,
C     terminated by the word 'END'.                    (CALL CKCHAR)
C  Optional thermodynamic properties input -
C     The word 'THERMO' starts input of species thermodynamic data.
C     If 'THERMO' is followed on the same line by 'ALL', the input
C     line immediately following must have three values, TLO, TMID,
C     and THI, default settings for dividing temperatures of the
C     fits;  if 'THERMO' is not followed by 'ALL' (or if 'THERMO' is
C     not used), there must exist file 'therm.dat' for access to
C     thermodynamic data.  Formatted thermodynamic polynomial
C     coefficients and elemental composition data follows until 'END'.
C                                                      (CALL CKTHRM)
C  Optional reaction input -
C     The word 'REACTION' starts processing of loosely-formatted
C     reaction descriptions.  'REACTION' may be followed on the same
C     line by descriptions of reaction input parameter units:
C       'MOLES' - (default), pre-exponential units are moles-sec-K;
C       'MOLECULES' - pre-exponential units are molecules and
C                     will be converted to moles.
C       'KELVINS' - activation energies are Kelvins, else
C                   activation energies are converted to Kelvins;
C       'EVOLTS' - activation energies are electron volts
C       'CAL/MOLE' - (default), activation energies are cal/mole;
C       'KCAL/MOLE' - activation energies are Kcal/mole;
C       'JOULES/MOLE' - activation energies are joules/mole;
C       'KJOULES/MOLE' - activation energies are Kjoules/mole.
C                                                      (CALL CKUNIT)
C     Then, process a reaction line according to a set of rules,
C                                                      (CALL CKREAC)
C     process auxiliary reaction information on following lines,
C                                                      (CALL CKAUXL)
C     upon starting a new reaction, proof previous reaction
C                                                      (CALL CPREAC)
C
C  Arguments:
C  LIN     - Integer scalar, formatted input file unit number.
C  LOUT    - Integer scalar, formatted output file unit number.
C  LTHRM   - Integer scalar, thermodynamic data input file unit number.
C  MDIM    - Integer scalar, maximum number of elements allowed.
C  KDIM    - Integer scalar, maximum number of species allowed.
C  IDIM    - Integer scalar, maximum number of reactions allowed.
C  NPAR    - Integer scalar, required number of Arrhenius parameters.
C  NCP     - Integer scalar, number of CP fit coefficients.
C  NFIT   - Integer scalar, total number of thermodynamic fit coeff.
C  MAXTP   - Integer scalar, maximum number of temperatures used to
C            divide ranges of fit temperatures.
C  NTR     - Integer scalar, actual number of temperature ranges.
C  MAXSP   - Integer scalar, maximum number of species allowed to
C            participate in a reaction.
C  MAXTB   - Integer scalar, maximum number of enhanced third bodies
C            allowed in a reaction.
C  NLAR    - Integer scalar, required number of additional parameters
C            for Landau-Teller reaction formulation.
C  NFAR    - Integer scalar, maximum number of additional parameters
C            allowed for a pres dependency reaction formulation.
C  MAXORD  - Integer scalar, maximum number of species change-orders
C            allowed in a reaction.
C  KNAME(*)- Character string array, species names.
C  ENAME(*)- Character string array, element names.
C  AWT(*)  - Real array, element atomic weights.
C  KNCF(*,*)-Integer array, elemental composition of species.
C  WT(*)   - Real array, species molecular weights.
C  KPHSE(*)- Integer array, species phases.
C  KCHRG(*)- Integer array, species electronic charge.
C  A(*,*,*)- Real array, polynomial coefficients for species
C            thermodynamic properties over temperature ranges.
C  T(*,*)  - Real array, species thermodynamic fit temperatures.
C  NT(*)   - Integer array, number of species fit temperatures.
C  NSPEC(*)- Integer array, reactions reactants+products count.
C  NREAC(*)- Integer array, reactions reactant only count.
C  NU(*,*) - Integer matrix, reactions stoichiometric coefficients.
C  NUNK(*,*)-Integer matrix, participating species indices for
C            reactions.
C  PAR(*,*)- Real matrix, Arrhenius rate parameters for reactions.
C  IDUP(*) - Integer array, flag for duplication of reactions.
C  IREV(*) - Integer array, reaction indices for those with
C            explicit reverse parameters.
C  RPAR(*,*)-Real matrix, explicit reverse parameters if given.
C  ILAN(*) - Integer array, reaction indices for Landau-Teller
C            reactions.
C  PLAN(*,*)-Real matrix, additional Landau-Teller reaction
C            parameters, if given.
C  IRLT(*) - Integer array, reaction indices for Landau-Teller
C            reactions with explicit reverse Landau-Teller
C            parameters.
C  RLAN(*,*)-Real matrix, additional reverse Landau-Teller reaction
C            parameters, if given.
C  IWL(*)  - Integer array, reaction indices for radiation-
C            enhanced reactions.
C  WL(*)   - Real array, radiation enhancement factors if given.
C  IFAL(*) - Integer array, reaction indices for pres dependency
C            reactions.
C  IFOP(*) - Integer array, integer flags for type of pres
C            dependency rate formulation.
C  IFLO(*) - Integer array, integer flags for HIGH/LOW pres
C            dependency.
C  KFAL(*) - Integer array, species indices for pres dependency
C            third-body in pres dependency rates.
C  PFAL(*,*)-Real array, additional rate parameters for pres
C            dependency reactions.
C  ITHB(*)  -Integer array, reaction indices for third-body reactions.
C  NTBS(*)  -Integer array, reactions enhanced third body counts.
C  AIK(*,*) -Real matrix, species enhancement factors for third-body
C            reactions.
C  NKTB(*,*)-Integer matrix, species indices for enhanced third-
C            bodies.
C  IRNU(*)  -Integer array, reaction indices for those with real
C            stoichiometry.
C  RNU(*,*) -Real matrix, reactions real stoichiometri coefficients.
C  IORD(*)  -Integer array, reaction indices for those with changed
C            species orders.
C  KORD(*,*)-Integer matrix, species indices for changed-order
C            reactions.
C  RORD(*,*)-Real matrix, reactions change-order values.
C  KION(*)  -Integer array, ionic species indices.
C  NJAR     -Integer scalar, required number of Jannev, et al.
C            reaction parameters.
C  NF1R     -Integer scalar, required number of fit#1-type reaction
C            parameters.
C  IEIM(*)  -Integer array, electron impact reaction indices.
C  IEIMT(*) -Integer array, electron impact reaction temperature-
C            dependency flags.
C  IJAN(*)  -Integer array, Jannev et al. reaction indices.
C  PJAN(*,*)-Real matrix, Jannet et al. reaction parameters.
C  IFT1(*)  -Integer array, fit#1-type reaction indices.
C  PFT1(*,*)-Real matrix, fit#1-type reaction parameters.
C  IEXC(*)  -Integer array, excitation reaction indices.
C  PEXC(*)  -Real array, energy loss for excitation reactions.
C  IMOM(*)  -Integer array, electron momentum transfer reaction
C            indices.
C  KMOM(*)  -Integer array, heavy species indices for electron
C            momentum transfer reactions.
C  IXSM(*)  -Integer array, ion momentum transfer reaction indices.
C  IXSK(*)  -Integer array, non-ion partner species indices for
C            momentum transfer reactions.
C  IXSI(*)  -Integer array, ion species indices for momentum transfer
C            reactions.
C  ITDE(*)  -Integer array, non-thermal temperature dependent reaction
C            indices.
C  ITDK(*)  -Integer array, non-thermal temperature dependent species
C            indices.
C  ITHRM(*) -Logical array, flag of thermodynamic properties completion
C            for species.
C  KERR     -Logical error flag
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER(I-N)
C*****END precision > single
C
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KPHSE(KDIM), KCHRG(KDIM), NT(KDIM),
     1          NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     2          NUNK(MAXSP,IDIM), IDUP(IDIM), IREV(IDIM), ILAN(IDIM),
     3          IRLT(IDIM), IWL(IDIM), IFAL(IDIM), IFOP(IDIM),
     4          IFLO(IDIM), KFAL(IDIM), ITHB(IDIM), NTBS(IDIM),
     5          NKTB(MAXTB,IDIM), IRNU(IDIM), IORD(IDIM),
     6          KORD(MAXORD,IDIM), KION(KDIM), IEIM(IDIM),
     7          IEIMT(IDIM), IJAN(IDIM), IFT1(IDIM), IEXC(IDIM),
     8          IMOM(IDIM), KMOM(IDIM), IXSM(IDIM), IXSK(IDIM),
     9          IXSI(IDIM), ITDE(IDIM), ITDK(IDIM)
C     Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), A(NFIT,NTR,KDIM), T(MAXTP,KDIM),
     1          PAR(NPAR,IDIM), RPAR(NPAR,IDIM), PLAN(NLAR,IDIM),
     2          RLAN(NLAR,IDIM), WL(IDIM), PFAL(NFAR,IDIM),
     3          AIK(MAXTB,IDIM), RNU(MAXSP,IDIM), RORD(MAXORD,IDIM),
     4          PJAN(NJAR,IDIM), PFT1(NF1R,IDIM), PEXC(IDIM)
C
      LOGICAL ITHRM(KDIM), THERMO, KERR, CKFILE
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM)
      CHARACTER*4 KEY, AUNITS, EUNITS, AXUNIT, EXUNIT, CKCHUP, TEMP
      CHARACTER*80 SUB(80), LINE, ISTR, IUNITS
      INTEGER CKLSCH, CKSLEN
      EXTERNAL CKLSCH, CKSLEN, CKCHUP, CKFILE
C
C     initialize task setting
      ITASK = 0
C     default allows use of therm.dat thermodynamic coefficient file
      THERMO = .TRUE.
C
C     top of mechanism input file handling
  100 CONTINUE
      LINE = ' '
      READ (LIN,'(A)',END=5000) LINE
C     "de-tab" the input
      CALL CKDTAB (LINE)
C
  105 CONTINUE
      ILEN = CKSLEN(LINE)
      IF (ILEN .EQ. 0) THEN
C        blank, or commented input line
         GO TO 100
      ENDIF
C
C     space-delimited character-strings
      CALL CKISUB (LINE(1:ILEN), SUB, NSUB)
C
C     is there a "keyword?"
      KEY = ' '
      KEY = CKCHUP(SUB(1), 4)
C
      IF (KEY.EQ.'SPEC' .OR. KEY.EQ.'ELEM') THEN
C
         IF (KEY .EQ. 'ELEM') THEN
C           set element processing task flag
            ITASK = 1
         ELSE
C           set species processing task flag
            ITASK = 2
         ENDIF
C
C        is there more information on the line?
         IF (NSUB .EQ. 1) GO TO 100
C        some elements or species follow the keyword
         DO 25 N = 2, NSUB
            SUB(N-1) = ' '
            SUB(N-1) = SUB(N)
   25    CONTINUE
         NSUB = NSUB-1
C        process this data
         GO TO 555
C
      ELSEIF (KEY .EQ. 'THER') THEN
C
C        turn off any input option
         ITASK = 0
C
C        process user's thermodynamic data
         IF (NSUB.GT.1 .AND. CKCHUP(SUB(2),3).EQ.'ALL') THERMO=.FALSE.
C	
C        if not THERMO, don't use LTHRM at all, not even for
C        getting temperature ranges;
C        if THERMO, use LTHRM to get temperature ranges, then use
C        LIN to get thermodynamic data
         IF (THERMO) THEN
            IF (.NOT. CKFILE(LTHRM, 'FORM')) THEN
               WRITE (LOUT, 334)
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
C
         CALL CKTHRM (LIN, LTHRM, THERMO, MDIM, ENAME, AWT, KNAME,
     1                KNCF, KPHSE, KCHRG, WT, MAXTP, NT, NTR, T,
     2                NFIT, A, ITHRM, KERR, LOUT, LINE)
C
C        use last LINE read for finding next ITASK
         GO TO 105
C
      ELSEIF (KEY .EQ. 'REAC') THEN
C
         ITASK = 4
C        start of reactions; are units specified?
         CALL CKUNIT (LINE(1:ILEN), AUNITS, EUNITS, IUNITS)
         GO TO 100
      ELSEIF (KEY .EQ. 'END') THEN
         ITASK = 0
         GO TO 100
      ENDIF
C
  555 CONTINUE
C
C     to get here, NSUB>0 and have a process to complete
      IF (ITASK .EQ. 1) THEN
C
C        element data
         IF (MM .EQ. 0) THEN
C           print header for element data
            WRITE (LOUT, 200)
            WRITE (LOUT, 300)
            WRITE (LOUT, 200)
         ENDIF
C
         M1 = MM +1
         CALL CKCHAR (SUB, NSUB, MDIM, ENAME, AWT, MM, KERR, LOUT)
         DO 110 M = M1, MM
            IF (AWT(M) .LE. 0.0) AWT(M) = CKATOM(ENAME(M))
C           print element data
            WRITE (LOUT, 400) M,ENAME(M)(1:4),AWT(M)
            IF (ENAME(M).EQ.'e' .OR. ENAME(M).EQ.'E') MELECT=M
            IF (AWT(M) .LE. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1000) ENAME(M)
            ENDIF
  110    CONTINUE
C
      ELSEIF (ITASK .EQ. 2) THEN
C
C        process species data
         IF (KK .EQ. 0) THEN
C           finished printing element data
            WRITE (LOUT, 200)
         ENDIF
C
         CALL CKCHAR (SUB, NSUB, KDIM, KNAME, WT, KK, KERR, LOUT)
C        cannot print species data until thermodynamic data set
C
      ELSEIF (ITASK .EQ. 4) THEN
C
C        process reaction data
         IF (II .LE. 0) THEN
C
C           this is the first reaction
C           assume end of species; need thermo data if not THERMO ALL
            IF (THERMO) THEN
C              not THERMO ALL, so use LTHRM both for
C              temperature ranges and for thermo data
               IF (.NOT. CKFILE(LTHRM, 'FORM')) THEN
                  WRITE (LOUT, 334)
                  KERR = .TRUE.
                  RETURN
               ENDIF
               CALL CKTHRM (LTHRM, LTHRM, THERMO, MDIM, ENAME, AWT,
     1                      KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP,
     2                      NT, NTR, T, NFIT, A, ITHRM, KERR, LOUT,
     3                      ISTR)
C              done with thermodynamic data
               THERMO = .FALSE.
            ENDIF
C
C           in any THERMO case, check and print species
            CALL CKPRNT (MDIM, KDIM, MAXTP, ENAME, KNAME, WT, KPHSE,
     1                   KCHRG, NT, T, KNCF, ITHRM, LOUT, KERR, KION)
            WRITE (LOUT, 1800)
         ENDIF
C
         IND = 0
C        is this auxiliary reaction data?
         DO 120 N = 1, NSUB
            TEMP = ' '
            TEMP = CKCHUP (SUB(N), 3)
            IND = MAX(IND, INDEX(SUB(N),'/'), INDEX(TEMP,'DUP'),
     1                     INDEX(TEMP,'MOM'), INDEX(TEMP,'XSM'))
  120    CONTINUE
         IF (IND .GT. 0) THEN
C
C           process auxiliary reaction data
            CALL CKAUXL (KDIM, IDIM, SUB, NSUB, KNAME, LOUT, MAXSP,
     1                   NPAR, NREAC(II), NSPEC(II), ITHB, NTBS,
     2                   MAXTB, NKTB,
     2                   AIK, IFAL, IDUP, NFAR, PFAL, IFOP, IFLO,
     3                   ILAN, NLAR, PLAN, IREV, RPAR, IRLT, RLAN,
     4                   IWL, WL, KERR, IORD, MAXORD, KORD, RORD,
     5                   NUNK(1,II), NU(1,II), IRNU, RNU, IEIM,
     6                   IEIMT, IJAN, NJAR, PJAN, IFT1, NF1R, PFT1,
     7                   IEXC, PEXC, IMOM, IXSM, ITDE, ITDK, KCHRG,
     8                   AXUNIT, EXUNIT)
            GO TO 100
         ENDIF
C
C        process a reaction string
         IF (II .GE. IDIM) THEN
C           cannot add more reactions
            WRITE (LOUT, 1070)
            KERR = .TRUE.
            GO TO 100
         ENDIF
C
         IF (II .GT. 0) THEN
C
C           check previous reaction for completeness
            CALL CPREAC (MDIM, KDIM, IDIM, MAXSP, MAXORD,
     1                   NSPEC, NPAR, PAR, RPAR, AUNITS,
     2                   EUNITS, NREAC, NUNK, NU, KCHRG, KNCF,
     3                   IDUP, IFAL, KFAL, NFAR, PFAL, IFOP,
     4                   IFLO, IREV, ITHB, ILAN, IRLT, KERR,
     5                   LOUT, IRNU, RNU, IXSM, IXSI, IXSK,
     6                   IMOM, KMOM, IORD, KORD, RORD, CKMIN,
     7                   AXUNIT, EXUNIT)
         ENDIF
C
C        new reaction
         AXUNIT = ' '
         EXUNIT = ' '
         II = II + 1
         CALL CKREAC (KDIM, IDIM, LINE(1:ILEN), KNAME, LOUT,
     1                MAXSP, NSPEC(II), NREAC(II), NUNK(1,II),
     2                NU(1,II), NPAR, PAR(1,II), ITHB, IFAL,
     3                KFAL, IWL, WL, IRNU, RNU, KERR)
C
         GO TO 100
      ENDIF
C
      GO TO 100
C
 5000 CONTINUE
C
C     end of input
      IF (II .GT. 0) THEN
C
C         check final reaction for completeness and print
          CALL CPREAC (MDIM, KDIM, IDIM, MAXSP, MAXORD,
     1                 NSPEC, NPAR, PAR, RPAR, AUNITS,
     1                 EUNITS, NREAC, NUNK, NU, KCHRG,
     2                 KNCF, IDUP, IFAL, KFAL, NFAR, PFAL,
     3                 IFOP, IFLO, IREV, ITHB, ILAN, IRLT,
     4                 KERR, LOUT, IRNU, RNU, IXSM, IXSI,
     6                 IXSK, IMOM, KMOM, IORD, KORD,
     6                 RORD, CKMIN, AXUNIT, EXUNIT)
C
C        check reactions declared as duplicates
         DO 500 I = 1, II
            IF (IDUP(I) .LT. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1095) I
            ENDIF
  500    CONTINUE
C
         WRITE (LOUT, '(/1X,A)') ' NOTE: '//IUNITS(1:CKLSCH(IUNITS))
         RETURN
      ENDIF
C
C     there was no reaction data, make sure species data is complete
      IF (THERMO) THEN
C        user did not use THERMO ALL
         IF (.NOT. CKFILE(LTHRM, 'FORM')) THEN
            WRITE (LOUT, 334)
            KERR = .TRUE.
            RETURN
         ENDIF
         CALL CKTHRM (LTHRM, LTHRM, THERMO, MDIM, ENAME, AWT,
     2                 KNAME, KNCF, KPHSE, KCHRG, WT, MAXTP, NT,
     3                 NTR, T, NFIT, A, ITHRM, KERR,
     4                 LOUT, LINE)
      ENDIF
C
C     check and print summary of species data
      CALL CKPRNT (MDIM, KDIM, MAXTP, ENAME, KNAME, WT, KPHSE,
     1             KCHRG, NT, T, KNCF, ITHRM, LOUT, KERR, KION)
C
  200 FORMAT (26X,20('-'))
  300 FORMAT (26X,'ELEMENTS',5X,'ATOMIC',/26X,'CONSIDERED',3X,'WEIGHT')
  333 FORMAT (/6X,'Error...no TLO,TMID,THI given for THERMO ALL...'/)
  334 FORMAT (/6X,'Error accessing thermodynamic data file...'/)
  400 FORMAT (25X,I3,'. ',A4,G15.6)
C
 1000 FORMAT (6X,'Error...no atomic weight for element ',A)
 1070 FORMAT (6X,'Error...more than IDIM reactions...')
 1095 FORMAT (6X,'Error...no duplicate declared for reaction no.',I3)
 1800 FORMAT (///54X, '(k = A T**b exp(-E/RT))',/,
     1        6X,'REACTIONS CONSIDERED',30X,'A',8X,'b',8X,'E',/)
C
C     end of SUBROUTINE CKKEY
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNEIM (LOUT, IDIM, II, RSTR, LEIM, LTHB, NEIM, IEIM,
     1                   IEIMT, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNEIM
C  processes electron-impact reaction character-string input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LEIM     - Logical, .TRUE. if electron-impact data already used in
C             reaction II.
C  LTHB     - Logical, .TRUE. if II is third-body reaction.
C  NEIM     - Integer scalar, electron-impact type reaction count.
C  IEIM(*)  - Integer array, the EIM reaction indices.
C  IEIMT(*) - Integer array, NEIM temperature-dependency flags.
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
      CHARACTER*(*) RSTR
      DIMENSION IEIM(IDIM), IEIMT(IDIM)
      LOGICAL LEIM, LTHB, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LEIM) THEN
         KERR = .TRUE.
         WRITE (LOUT,'(6X,A)') 'Error...multiple EIM declarations...'
      ELSE
         NEIM = NEIM + 1
         IEIM(NEIM) = II
C
         IF (LTHB) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...EIM and 3rd-body reaction mutually exclusive...'
         ENDIF
      ENDIF
C
      CALL CKPARI (RSTR, 1, 1, IEIMT(NEIM), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.1) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...problem reading EIM parameter...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, '(6X,A,I5)')
     1   'Electron 3rd-body reaction; Temp. Dependence =', IEIMT(NEIM)
      ENDIF
C
C     end of SUBROUTINE CKNEIM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNEXC (LOUT, IDIM, II, RSTR, LEXC, NEXC, IEXC, PEXC,
     1                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNEXC
C  processes excitation-only reaction character-string input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LEXC     - Logical, .TRUE. if excitation input already used in
C             this reaction.
C  NEXC     - Integer scalar, excitation reaction count.
C  IEXC(*)  - Integer array, excitation reaction indices.
C  PEXC(*)  - Real array, NEXC energy losses (eV).
C  KERR     - Logical error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION IEXC(IDIM), PEXC(IDIM)
      LOGICAL LEXC, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LEXC) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple EXC declarations...'
      ELSE
         NEXC = NEXC + 1
         IEXC(NEXC) = II
      ENDIF
C
      CALL CKPARR (RSTR, 1, 1, PEXC(NEXC), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.1) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X, A, A)')
     1      'Error...problem finding EXC parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 130) PEXC(NEXC)
      ENDIF
C
  130 FORMAT (6X,'Excitation reaction, energy loss =',e10.3,' eV')
C
C     end of SUBROUTINE CKNEXC
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNFAL (LOUT, IDIM, II, KEY, SUB, RSTR, NFAL, IFAL,
     1                   IFOP, IFLO, NFAR, PFAL, LFAL, LLAN, LRLT, LREV,
     2                   LTRO, LSRI, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNFAL
C  processes pressure-dependency reaction auxiliary character-string
C  input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KEY      - Character string, pres dependency option type.
C  SUB      - Character string.
C  RSTR     - Character string, representation of real values.
C  NFAL     - Integer scalar, pres dependency reactions count.
C  IFOP(*)  - Integer array, pres dependency formulation type,
C             0, none
C             1, Lindemann
C             2, SRI
C             3, 3-parameter TROE
C             4, 4-parameter TROE
C  IFLO(*)  - Integer array, description of pressure dependency.
C             0, Arrhenius rate is high-pressure limit,
C                low-pressure region rate-modify.
C             1, Arrhenius rate is low-pressure limit,
C                high-pressure region rate-modify.
C  NFAR     - Integer scalar, maximum number of additional
C             pres dependency parameters allowed.
C  PFAL(*,*)- Real matrix, pres dependency parameters.
C  LFAL     - Logical, .TRUE. if main reaction string was properly
C             initialized for pres dependency options (requires "(+").
C  LLAN     - Logical, .TRUE. if Landau-Teller formulation used.
C             This option is not allowed in a pres dependency reaction.
C  LREV     - Logical, .TRUE. if explicit reverse parameters given.
C             This option is not allowed in a pres dependency reaction.
C  LTRO     - Logical, .TRUE. if TROE parameters have been given.
C  LSRI     - Logical, .TRUE. if SRI parameters have been given.
C  KERR     - Logical, error flag.
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
      CHARACTER*(*) KEY, SUB, RSTR
      DIMENSION IFAL(IDIM), IFOP(IDIM), IFLO(IDIM), PFAL(NFAR,IDIM)
      LOGICAL LFAL, LLAN, LRLT, LREV, LTRO, LSRI, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (.NOT. LFAL) THEN
C        pres dep option needed "(+M)" or "(+species)" in reaction
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...reaction does not meet requirements for ',KEY
         NFAL = NFAL + 1
         IFAL(NFAL) = II
      ENDIF
C
      IF (LLAN) THEN
C        Landau-Teller reaction
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A)')
     1      'Error...',KEY,' and LAN mutually exclusive...'
      ENDIF
C
      IF (LRLT) THEN
C        reverse Landau-Teller reaction
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A)')
     1      'Error...',KEY,' and RLT mutually exclusive...'
      ENDIF
C
      IF (LREV) THEN
C        reverse parameters given
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A)')
     1      'Error...',KEY,' and REV mutually exclusive...'
      ENDIF
C
      IF (KEY(1:3) .EQ. 'LOW') THEN
C
C        Lindemann parameters usually, unless modification
C        parameters given later
         IF (IFLO(NFAL) .EQ. 0) THEN
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A)')
     1      'Error...multiple LOW declarations...',SUB(1:ILEN)
            KERR = .TRUE.
         ELSEIF (IFLO(NFAL) .EQ. 1) THEN
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A)')
     1      'Error...LOW and HIGH mutually exclusive...',SUB(1:ILEN)
            KERR = .TRUE.
         ELSE
            IFLO(NFAL) = 0
C           no previous pres dependency option; set Lindemann
            IF (IFOP(NFAL) .LE. 0) IFOP(NFAL) = 1
         ENDIF
         CALL CKPARR (RSTR, 1, 3, PFAL(1,NFAL),NVAL,IER,LOUT)
         IF (IER.NE.0 .OR. NVAL.NE.3) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(RSTR)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...problem finding LOW parameters...',RSTR(1:ILEN)
         ELSE
            WRITE (LOUT, '(6X,A,3E13.5)')
     1         'Low pressure limit:',(PFAL(L,NFAL),L=1,3)
         ENDIF
C
      ELSEIF (KEY(1:3) .EQ. 'HIG') THEN
C
C        Lindemann parameters usually, unless  modification
C        parameters given later
         IF (IFLO(NFAL) .EQ. 0) THEN
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A)')
     1      'Error...HIGH and LOW mutually exclusive...',SUB(1:ILEN)
            KERR = .TRUE.
         ELSEIF (IFLO(NFAL) .EQ. 1) THEN
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A)')
     1      'Error...multiple HIGH declarations...',SUB(1:ILEN)
            KERR = .TRUE.
         ELSE
            IFLO(NFAL) = 1
C           no previous pres dependency option; set Lindemann
            IF (IFOP(NFAL) .LE. 0) IFOP(NFAL) = 1
         ENDIF
         CALL CKPARR (RSTR, 1, 3, PFAL(1,NFAL), NVAL, IER, LOUT)
         IF (IER.NE.0 .OR. NVAL.NE.3) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(RSTR)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...problem finding HIGH parameters...',RSTR(1:ILEN)
         ELSE
            WRITE (LOUT, '(6X,A,3E13.5)')
     1         'High pressure limit:',(PFAL(L,NFAL),L=1,3)
         ENDIF
C
      ELSEIF (KEY(1:4) .EQ. 'TROE') THEN
C        use TROE formulation
         IF (LTRO) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...multiple TROE declarations...',SUB(1:ILEN)
         ENDIF
         IF (LSRI) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...TROE and SRI mutually exclusive...',SUB(1:ILEN)
         ENDIF
C
C        have 3 LOW parameters, need 3 or 4 more
         CALL CKPARR (RSTR,1,-4,PFAL(4,NFAL),NVAL,IER,LOUT)
         IF (NVAL.LT.3 .OR. NVAL.GT.4 .OR. IER.NE.0) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(RSTR)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...problem finding TROE parameters...',RSTR(1:ILEN)
         ELSE
            IF (NVAL .EQ. 3) THEN
               IFOP(NFAL) = 3
               WRITE (LOUT, '(6X,A,3E13.5)')
     1         'TROE centering:    ', (PFAL(L,NFAL),L=4,6)
            ELSE
               IFOP(NFAL) = 4
               WRITE (LOUT, '(6X,A,4E13.5)')
     1         'TROE centering:    ', (PFAL(L,NFAL),L=4,7)
            ENDIF
         ENDIF
C
      ELSEIF (KEY(1:3) .EQ. 'SRI') THEN
C        use SRI formulation
         IF (LTRO) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...TROE and SRI mutually exclusive...',SUB(1:ILEN)
         ENDIF
         IF (LSRI) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...multiple SRI declarations...',SUB(1:ILEN)
         ENDIF
C
C        have 3 LOW parameters, need 5 additional
         CALL CKPARR (RSTR,1,-5,PFAL(4,NFAL),NVAL,IER,LOUT)
         IF (IER.NE.0 .OR. NVAL.LT.3 .OR. NVAL.EQ.4 .OR. NVAL.GT.5) THEN
            KERR = .TRUE.
            ILEN = CKLSCH(RSTR)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...problem finding SRI parameters...',RSTR(1:ILEN)
         ELSE
            IFOP(NFAL) = 2
            IF (NVAL .EQ. 3) THEN
C              user declared 3, use default values for others
               PFAL(7,NFAL) = 1.0
               PFAL(8,NFAL) = 0.0
               WRITE (LOUT, '(6X,A,3E13.5)')
     1         'SRI centering:     ', (PFAL(L,NFAL),L=4,6)
            ELSE
C              used declared all 5
               WRITE (LOUT, '(6X,A,3E13.5)')
     1         'SRI centering:     ', (PFAL(L,NFAL),L=4,8)
            ENDIF
         ENDIF
      ENDIF
C
C     end of SUBROUTINE CKNFAL
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNFT1 (LOUT, IDIM, II, RSTR, LFT1, NFT1, IFT1, NF1R,
     1                   PFT1, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNFT1
C  processes Fit #1-type reaction character-string input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum  number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real numbers.
C  LFT1     - Logical, .TRUE. if this reaction already declared FIT1.
C  NFT1     - Integer scalar, Fit#1 reactions count.
C  IFT1(*)  - Integer array, Fit#1 reaction indices.
C  NF1R     - Integer scalar, required number of FIT1 parameters.
C  PFT1(*,*)- Real matrix, FIT1 parameters.
C  KERR     - Logical, error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION IFT1(IDIM), PFT1(NF1R,IDIM)
      LOGICAL LFT1, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LFT1) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple FIT1 declarations...'
         RETURN
      ENDIF
C
      NFT1 = NFT1 + 1
      IFT1(NFT1) = II
C
      CALL CKPARR (RSTR, 1, NF1R, PFT1(1,NFT1), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NF1R) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X, A, A)')
     1      'Error...problem finding FIT1 parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 130) (PFT1(N,NFT1), N=1,NF1R)
      ENDIF
C
  130 FORMAT (6X,'Modified fit#1:  k= A * T^B * exp [SUM(Vn/T^n)]...',
     1       /6X,'Added parameters: ',E10.3,3(/27X,E10.3))
C
C     end of SUBROUTINE CKNFT1
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNJAN (LOUT, IDIM, II, RSTR, LJAN, NJAN, IJAN, NJAR,
     1                   PJAN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNJAN
C  processes Jannev et al. reaction character-string inputinput
C
C  Arguments:
C  LOUT     - Integer scalar, formatted input file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LJAN     - Logical, .TRUE. if JAN already used for this reaction.
C  NJAN     - Integer scalar, Jannev et al. reactions count.
C  IJAN(*)  - Integer array, Jannev et al. reaction indices.
C  NJAR     - Integer scalar, required number of Jannev parameters.
C  PJAN(*,*)- Real matrix, Jannev-Langer rate parameters.
C  KERR     - Logical error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION IJAN(IDIM), PJAN(NJAR,IDIM)
      LOGICAL LJAN, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LJAN) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple JAN declarations...'
      ELSE
         NJAN = NJAN + 1
         IJAN(NJAN) = II
      ENDIF
C
      CALL CKPARR (RSTR, 1, NJAR, PJAN(1,NJAN), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NJAR) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X, A, A)')
     1      'Error...problem finding JAN parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 130) (PJAN(N,NJAN), N=1,NJAR)
      ENDIF
C
  130 FORMAT (6X,'Jannev, Langer, Evans & Post type reaction:'
     1       /6X,'Coefficients: ',5(E10.3,1X)
     2      /23X,5(E10.3,1X))
C
C     end of SUBROUTINE CKNJAN
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNLAN (LOUT, IDIM, II, RSTR, LLAN, LFAL, NLAN, NLAR,
     1                   ILAN, PLAN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNLAN
C  processes Landau-Teller reaction character-string input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LLAN     - Logical, .TRUE. if LT already used in this reaction.
C  NLAN     - Integer scalar, Landau-Teller reactions count.
C  NLAR     - Integer scalar, number of additional Landau-Teller
C             parameters required.
C  ILAN(*)  - Integer array, Landau-Teller reaction indices.
C  PLAN(*,*)- Real matrix, Landau-Teller parameters.
C  KERR     - Logical error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION ILAN(IDIM), PLAN(NLAR,IDIM)
      LOGICAL LFAL, LLAN, KERR
C
      IF (LLAN) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple LT declarations...'
      ELSE
         NLAN = NLAN + 1
         ILAN(NLAN) = II
      ENDIF
C
      IF (LFAL) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1      'Error...pres-dependency and LT mutually exclusive...'
      ENDIF
C
      CALL CKPARR (RSTR,1,NLAR,PLAN(1,NLAN),NVAL,IER,LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NLAR) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1      'Error...problem reading LT parameters...',RSTR
      ELSE
         WRITE (LOUT, '(6X,2(A,E12.5))')
     1   'Landau-Teller Parameters: B=', PLAN(1,NLAN),
     2   ', C=',PLAN(2,NLAN)
      ENDIF
C
C     end of SUBROUTINE CKNLAN
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNMOM (LOUT, KELECT, IDIM, II, LMOM, NMOM, IMOM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNMOM
C  processes electron momentum-transfer reaction character-string input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  KELECT   - Integer scalar, electron species index.
C  LMOM     - Logical, .TRUE. if MOM option already used in reaction II.
C  NMOM     - Integer scalar, electron momentum-transfer collision
C             reactions count.
C  IMOM(*)  - Integer array, momentum-transfer reaction indices.
C  KERR     - Logical error flag.
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
      DIMENSION IMOM(IDIM)
      LOGICAL KERR, LMOM
C
      IF (LMOM) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple MOME declarations...'
      ELSE
         NMOM = NMOM + 1
         IMOM(NMOM) = II
         WRITE (LOUT, '(6X,A,A)')
     1   'Electron momentum-transfer collision frequency',
     2   ' per molecule [cm3/s]'
C
         IF (KELECT .LE. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1   'Error...species list does not contain an electron species...'
         ENDIF
C
      ENDIF
C
C     need exactly one electron species with coeff. approx. 1,
C     and  exactly one non-electron species with coeff. approx. 1,
C     but cannot check until CPREAC, when all auxiliary information
C     is processed, including FORD changes
C
C     end of SUBROUTINE CKNMOM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNORD (LOUT, KDIM, IDIM, KEY, RSTR, II, KK, KNAME,
     1                   NORD, MAXORD, IORD, KORD, RORD, MAXSP, NUNK,
     2                   NU, NSPEC, NREAC, LORD, LRNU, NRNU, RNU, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNORD
C  processes reaction change-of-order character-string information.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  KEY      - Character string, type of change-order.
C  RSTR     - Character string, representation of real value(s).
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  KK       - Integer scalar, species count.
C  KNAME(*) - Character string array, species names.
C  NORD     - Integer scalar, change-order reactions count.
C  MAXORD   - Integer scalar, maximum number of species change-
C             orders allowed in a change-order reaction.
C  IORD(*)  - Integer array, change-order species indices.
C  KORD(*,*)- Integer matrix, change-order species indices,
C             < 0, change forward order for species KORD
C             > 0, change reverse order for species KORD
C  RORD(*,*)- Real matrix, change-order values.
C  MAXSP    - Integer scalar, maximum number of species allowed in a
C             reaction.
C  NUNK(*)  - Integer array, reactions species indices.
C  NU(*)    - Integer array, reactions stoichiometric coefficients.
C  NSPEC    - Integer scalar, reactions reactants+products count.
C  LORD     - Logical, .TRUE. if II had previous change-order input.
C  LRNU     - Logical, .TRUE. if II has real stoichiometry.
C  NRNU     - Integer scalar, real stoichiometry reactions count.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  KERR     - Logical error flag.
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
      CHARACTER KNAME(KDIM)*16, KEY*(*), RSTR*(*)
      DIMENSION IORD(IDIM), NUNK(MAXSP), NU(MAXSP), KORD(MAXORD,IDIM),
     1          RORD(MAXORD,IDIM), RNU(MAXSP,IDIM)
      LOGICAL LORD, LFORD, LRORD, LRNU, KERR, IERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
      IERR = .FALSE.
C
C     LRORD is reverse change-of-order,
C     LFORD is forward reaction change-of-order
C
      LFORD = KEY(1:1).EQ.'F'
      LRORD = KEY(1:1).EQ.'R'

      IF (.NOT. LORD) THEN
C
C        need to initialize this change-of-order reaction;
C        increment the number of reactions with non-standard orders
         NORD = NORD + 1
C        add this reaction index to the list of NORD reactions
         IORD(NORD) = II
C
C        copy original stoichiometry
         NKORD = 0
         NKORD = ABS(NREAC)
         DO 100 N = 1, NKORD
C           reactants
            KORD(N, NORD) = -NUNK(N)
            IF (LRNU) THEN
               RORD(N, NORD) = ABS(RNU(N,NRNU))
            ELSE
               RORD(N, NORD) = IABS(NU(N))
            ENDIF
  100    CONTINUE
         DO 110 N = MAXSP/2+1, MAXSP
C           products
            IF (NUNK(N) .NE. 0) THEN
               NKORD = NKORD + 1
               KORD(NKORD, NORD) = NUNK(N)
               IF (LRNU) THEN
                  RORD(NKORD, NORD) = RNU(N,NRNU)
               ELSE
                  RORD(NKORD, NORD) = NU(N)
               ENDIF
            ENDIF
  110    CONTINUE
      ENDIF
C
C     need to split species name from change-order value
      KEND = 0
      ILEN = CKSLEN(RSTR)
      DO 55 I = 2, ILEN
         IF (KEND .EQ. 0) THEN
            IF (RSTR(I:I).EQ.' ' .AND. RSTR(I-1:I-1).NE.' ') KEND=I-1
            IF (RSTR(I:I).NE.' ' .AND. I.EQ.ILEN) KEND=I
         ENDIF
   55 CONTINUE
C
      KNUM = 0
      IF (KEND .GT. 0) CALL CKCOMP (RSTR(1:KEND), KNAME, KK, KNUM)
C
      IF (KNUM .EQ. 0) THEN
C        could not get species; issue error message and set flag
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A,A)')
     1      'Error...unrecognized ',KEY,' species...', RSTR(1:ILEN)
      ENDIF
C
      NVAL = 0
      IF (KEND .GT. 0)
     1   CALL CKPARR (RSTR(KEND+1:), 1, 1, VAL, NVAL, IER, LOUT)
C
      IF (IER.NE.0 .OR. NVAL.NE.1) THEN
C        could not get good order value for species
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A,A)')
     1   'Error...problem finding ',KEY,' parameter...', RSTR(1:ILEN)
      ENDIF
C
      IF (LRORD .AND. NSPEC.LT.0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...RORD incompatible with irreversible reaction...',
     2   RSTR(1:ILEN)
      ENDIF
C
      IF (KNUM .NE. 0) THEN
C        store species index; negative for reactants
         IF (LFORD) KNUM = -KNUM
C
         NFIX =0
         DO 200 N = 1, MAXORD
            NK = KORD(N, NORD)
            IF (NFIX.EQ.0 .AND. (NK.EQ.KNUM .OR. NK.EQ.0))
C              this is the place to put species
     1         NFIX = N
  200    CONTINUE
C
         IF (NFIX .EQ. 0) THEN
C           cannot store more species
            IERR = .TRUE.
            WRITE (LOUT, '(6X,A,I5,A)')
     1      'Error...more than ',MAXORD,' change-orders...',
     2      RSTR(1:ILEN)
         ENDIF
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
      ELSE
         KORD(NFIX, NORD) = KNUM
         RORD(NFIX, NORD) = VAL
         IF (LFORD) THEN
            WRITE (LOUT, '(6X, A, A, 1PE12.3)')
     1      'Forward order ', KNAME(ABS(KNUM)), VAL
         ELSE
            WRITE (LOUT, '(6X, A, A, 1PE12.3)')
     1      'Reverse order ', KNAME(KNUM), VAL
         ENDIF
      ENDIF
C
C     end of SUBROUTINE CKNORD
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNREV (LOUT, IDIM, II, RSTR, LFAL, LREV, NSPEC, NPAR,
     1                   NREV, IREV, RPAR, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNREV
C  processes reaction explicit reverse parameter character-string info.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LFAL     - Logical, .TRUE. if reaction has pres dependency.
C  LREV     - Logical, .TRUE. if reaction already has reverse
C             parameters.
C  NSPEC    - Integer scalar, reaction's count of reactants+products.
C  NPAR     - Integer scalar, required number of Arrhenius params.
C  NREV     - Integer scalar, count of reverse parameter reactions.
C  IREV(*)  - Integer array, explicit rev parameter reaction indices.
C  RPAR(*,*)- Real matrix, reverse Arrhenius parameters.
C  KERR     - Logical error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION IREV(IDIM), RPAR(NPAR,IDIM)
      LOGICAL LFAL, LREV, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LREV) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple REV declarations...'
      ELSE
         NREV = NREV + 1
         IREV(NREV) = II
         IF (NSPEC .LT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1   'Error...REV and irreversible reaction mutually exclusive...'
         ENDIF
      ENDIF
C
      IF (LFAL)  THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1      'Error...pres-dependency and REV mutually exclusive...'
      ENDIF
C
      CALL CKPARR (RSTR, 1, NPAR, RPAR(1,NREV), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NPAR) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X, A, A)')
     1      'Error...problem finding REV parameters...',RSTR(1:ILEN)
      ELSE
         WRITE (LOUT, 130) (RPAR(N,NREV), N=1,NPAR)
      ENDIF
C
  130 FORMAT (6X, 'Reverse Arrhenius coefficients:',
     1            T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
C
C     end of SUBROUTINE CKNREV
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNRLT (LOUT, IDIM, II, RSTR, LFAL, LRLT, NSPEC, NLAR,
     1                   NRLT, IRLT, RLAN, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNRLT
C  processes explicitly-declared reverse Landau-Teller character-string
C  input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LFAL     - Logical, .TRUE. if reaction II has pres dependency.
C  LRLT     - Logical, .TRUE. if reaction II has already specified
C             reverse Landau-Teller input.
C  NSPEC    - Integer scalar, reactions reactants+products count.
C  NLAR     - Integer scalar, number of additional Landau-Teller
C             parameters required.
C  NRLT     - Integer scalar, count of reactions with explicit
C             reverse Landau-Teller formulation parameters.
C  IRLT(*)  - Integer array, reverse Landau-Teller reaction indices.
C  RLAN(*,*)- Real matrix, reverse Landau-Teller parameters.
C  KERR     - Logical error flag.
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
      CHARACTER*(*) RSTR
      DIMENSION IRLT(IDIM), RLAN(NLAR,IDIM)
      LOGICAL LFAL, LRLT, KERR
C
      IF (LRLT) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple RLT declarations...'
      ELSE
         NRLT = NRLT + 1
         IRLT(NRLT) = II
         IF (NSPEC .LT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1   'Error...RLT and irreversible reaction mutually exclusive...'
         ENDIF
      ENDIF
C
      IF (LFAL)  THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1      'Error...RLT and pres-dependency mutually exclusive...'
      ENDIF
C
      CALL CKPARR (RSTR, 1, NLAR, RLAN(1,NRLT), NVAL, IER, LOUT)
      IF (IER.NE.0 .OR. NVAL.NE.NLAR) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...RLT problem finding parameters...',RSTR
      ELSE
         WRITE (LOUT, '(6X,2(A,E12.5))')
     1   'Reverse Landau-Teller parameters: B=',RLAN(1,NRLT),
     2   ', C=',RLAN(2,NRLT)
      ENDIF
C
C     end of SUBROUTINE CKNRLT
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNTDE (LOUT, KDIM, IDIM, II, RSTR, LTDE, LTHB, KNAME,
     1                   KK, NTDE, ITDE, ITDK, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNTDE
C  processes non-thermal temperature-dependent reaction character-string
C  input.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of current
C             reaction.
C  RSTR     - Character string, representation of real values.
C  LTDE     - Logical, .TRUE. if reaction II already declared
C             non-thermal temperature-dependency.
C  LTHB     - Logical, .TRUE. if reaction II is a third-body reaction.
C  KNAME(*) - Character string array, species names.
C  KK       - Integer scalar, species count.
C  NTDE     - Integer scalar, non-thermal temp dependency rxn count.
C  ITDE(*)  - Integer array, non-thermal reaction indices.
C  ITDK(*)  - Integer array, temperature-dependency species indices.
C  KERR     - Logical error flag.
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
      CHARACTER KNAME(KDIM)*16, RSTR*(*)
      DIMENSION ITDE(IDIM), ITDK(IDIM)
      LOGICAL LTDE, LTHB, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (LTDE) THEN
         KERR = .TRUE.
         WRITE (LOUT,'(6X,A)') 'Error...multiple TDEP declarations...'
      ELSE
         NTDE = NTDE + 1
         ITDE(NTDE) = II
      ENDIF
C
      IF (LTHB) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1      'Error...TDEP and 3rd-body reaction mutually exclusive...'
      ENDIF
C
      CALL CKCOMP (RSTR, KNAME, KK, KNUM)
      IF (KNUM .LE. 0) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(RSTR)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...unrecognized TDEP species...',RSTR(1:ILEN)
      ELSE
         ITDK(NTDE) = KNUM
         WRITE (LOUT, '(6X, A, A)')
     1'Non-thermal reaction, depends on Species Temp. of ',KNAME(KNUM)
      ENDIF
C
C     end of SUBROUTINE CKNTDE
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNTHB (LOUT, II,KDIM, IDIM, KEY, SUB, RSTR, KNAME, KK,
     1                   LTHB, LEIM, LTDE, NTHB, ITHB, NTBS, NKTB,
     2                   MAXTB, AIK, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNTHB
C  processes third-body reactant character-string input.
C
C  Arguments:
C  LOUT      - Integer scalar, formatted output file unit number.
C  KEY       - Character string, 3rd body species name only.
C  SUB       - Character string, 3rd body species name and
C              enhancement value.
C  RSTR      - Character string, enhancement value only.
C  KNAME(*)  - Character string array, species names.
C  KK        - Integer scalar, species count.
C  LTHB      - Logical, .TRUE. if reaction string contained required
C              "+M";
C              if LTHB=.FALSE., reaction is in error.
C  LEIM      - Logical, .TRUE. if electron-impact reaction has been
C              specified by EIM keyword;
C              if LEIM=.TRUE., reaction is in error.
C  LTDE      - Logical, .TRUE. if non-thermal-temperture dependency has
C              been specified by TDEP keyword;
C              if LTDE=.TRUE., reaction is in error.
C  NTHB      - Integer scalar, third-body reaction count.
C  NTBS(*)   - Integer scalar, enhanced third-body count for
C              third-body reactions.
C  NKTB(*,*) - Integer matrix, species indices of enhanced third-
C              bodies in third-body reactions.
C  MAXTB     - Integer scalar, maximum number of enhanced third-bodies
C              allowed in a reaction.
C  AIK(*,*)  - Real matrix, enhancement factors for third-bodies in
C              third-body reactions.
C  KERR      - Logical error flag.
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
      DIMENSION ITHB(IDIM), NTBS(IDIM), AIK(MAXTB,IDIM),
     1          NKTB(MAXTB,IDIM)
      CHARACTER KNAME(KDIM)*16, KEY*(*), SUB*(*), RSTR*(*)
      LOGICAL LTHB, LEIM, LTDE, KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      IF (.NOT. LTHB) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(SUB)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...3rd-body requires (+M) in reaction...',SUB(1:ILEN)
         NTHB = NTHB + 1
         ITHB(NTHB) = II
      ENDIF
C
      IF (NTBS(NTHB) .EQ. MAXTB) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(SUB)
         WRITE (LOUT, '(6X,A,I3,A)')
     1   'Error...more than ',MAXTB,', 3rd-bodies...',SUB(1:ILEN)
      ENDIF
C
      IF (LEIM) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(SUB)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...3rd-body and EIM mutually exclusive...',SUB(1:ILEN)
      ENDIF
C
      IF (LTDE) THEN
         KERR = .TRUE.
         ILEN = CKLSCH(SUB)
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...3rd-body and TDEP mutually exclusive...',SUB(1:ILEN)
      ENDIF
C
      CALL CKCOMP (KEY, KNAME, KK, KNUM)
      CALL CKPARR (RSTR, 1, 1, VAL, NVAL, IER, LOUT)
C
      IF (IER.NE.0 .OR. NVAL.NE.1 .OR. KNUM.LE.0) THEN
         KERR = .TRUE.
         IF (KNUM .LE. 0) THEN
            ILEN = CKLSCH(KEY)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...unrecognized 3rd-body species...',KEY(1:ILEN)
         ENDIF
         IF (IER.NE.0 .OR. NVAL.NE.1) THEN
            ILEN = CKLSCH(SUB)
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...problem finding 3rd-body parameter...',
     2      SUB(1:ILEN)
         ENDIF
      ELSE
         NTBS(NTHB) = NTBS(NTHB) + 1
         NKTB(NTBS(NTHB), NTHB) = KNUM
         AIK(NTBS(NTHB),NTHB) = VAL
         WRITE (LOUT, '(6X,A16,A,1PE12.3)')
     1      KNAME(KNUM),'Enhanced by ', AIK(NTBS(NTHB),NTHB)
      ENDIF
C
C     end of SUBROUTINE CKNTHB
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNWL (LOUT, IDIM, II, RSTR, LWL, NWL, IWL, WL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNWL
C  processes reaction radiation-enhancement character-string info.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index
C             of this reaction.
C  RSTR     - Character string.
C  LWL      - Logical flag, .TRUE. if IWL(NWL)=II (reaction string
C             had a "+HV").
C  NWL      - Integer scalar, count of radiation-enhanced reactions.
C  IWL(*)   - Integer array, radiation-enhanced reaction indices.
C  WL(*)    - Real array, radiation wavelengths.
C  KERR     - Logical error flag.
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
       DIMENSION IWL(IDIM), WL(IDIM)
       CHARACTER RSTR*(*)
       LOGICAL LWL, KERR
       INTEGER CKLSCH
       EXTERNAL CKLSCH
C
       IF (.NOT.LWL) THEN
          KERR = .TRUE.
          WRITE (LOUT, '(6X,A)') 'Error...HV not in reaction species...'
          NWL = NWL + 1
          IWL(NWL) = II
       ELSEIF (WL(NWL) .GE. 0.0) THEN
          KERR = .TRUE.
          WRITE (LOUT,'(6X,A)') 'Error...multiple HV declarations...'
       ENDIF
       CALL CKPARR (RSTR, 1, 1, VAL, NVAL, IER, LOUT)
       IF (IER.NE.0 .OR. NVAL.NE.1) THEN
          KERR = .TRUE.
          ILEN = CKLSCH(RSTR)
          WRITE (LOUT, '(6X,A,A)')
     1       'Error finding HV parameter...', RSTR(1:ILEN)
      ELSE
          WL(NWL) = WL(NWL)*VAL
          WRITE (LOUT, '(6X,A,F10.2)')
     1   'Radiation wavelength (A): ',ABS(WL(NWL))
      ENDIF
C
C     end of SUBROUTINE CKNWL
      RETURN
      END
C                                                                      C
      SUBROUTINE CKNXSM (LOUT, IDIM, II, NKION, LXSI, NXSM, IXSM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNXSM
C  processes ion momentum-transfer reaction character-string input.
C
C  Summary of XSMI ion momentum-transfer reaction description:
C     Multiple XSMI declaration generates a WARNING
C     XSMI reaction requires exactly one ion reactant and exactly
C     one non-ion reactant; this requires a check of FORD
C     species to enforce this rule.
C
C  Arguments:
C  LOUT     - Integer scalar, formatted output file unit number.
C  IDIM     - Integer scalar, maximum number of reactions.
C  II       - Integer scalar, reaction count, and index of
C             current reaction.
C  NKION    - Integer scalar, ionic species count.
C  LXSI     - Logical flag; .TRUE. if XSMI already defined for this
C             reaction.
C  NXSM     - Integer scalar, ion momentum-transfer reaction count.
C  IXSM(*)  - Integer array, ion momentum-transfer reaction indices.
C  KERR     - Logical error flag.
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
      DIMENSION IXSM(IDIM)
      LOGICAL KERR, LXSI
C
      IF (LXSI) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)') 'Error...multiple XSMI declarations...'
      ELSE
         IF (NKION .LE. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...species list does not have ionic species...'
         ENDIF
         NXSM = NXSM + 1
         IXSM(NXSM) = II
         WRITE (LOUT, 4020)
      ENDIF
C
C     need exactly one ion reactant with coef approx 1,
C     and  exactly one non-ion reactant with coef approx 1,
C     but cannot confirm until CPREAC, after all auxiliary information
 4020 FORMAT (6X,'Ion momentum-transfer cross-section [cm2]')
C
C     end of SUBROUTINE CKNXSM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKPRNT (MDIM, KDIM, MAXTP, ENAME, KNAME, WT, KPHSE,
     1                   KCHRG, NT, T, KNCF, ITHRM,
     2                   LOUT, KERR, KION)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPRNT
C  is called after element and species input has been processed,
C  to do further checking for completeness and to print formatted
C  output.
C
C  Arguments:
C  MDIM     - Integer scalar, maximum number of elemnts.
C  MAXTP    - Integer scalar, maximum number of fit temperatures.
C  KDIM     - Integer scalar, maximum number of species.
C  MM       - Integer scalar, element count.
C  ENAME(*) - Character string array, element names.
C  KK       - Integer scalar, species count.
C  KNAME(*) - Character string array, species names.
C  WT(*)    - Real array, species molecular weights.
C  KPHSE(*) - Integer array, species phases.
C  KCHRG(*) - Integer array, species ionic charges.
C  NT(*)    - Integer array, count of dividing temperatures for
C             thermodynamic properties for species.
C  T(*,*)   - Real matrix, species thermodynics fit temperatures.
C  KNCF(*,*)- Integer matrix, elemental composition of species.
C  ITHRM(*) - Logical array, species thermodynamic flag.
C  LOUT     - Integer scalar, formatted output file unit number.
C  KERR     - Logical error flag.
C  NCHRG    - Integer scalar, ionic species count.
C  KELECT   - Integer scalar, electron species index.
C  NKION    - Integer scalar, non-electron ion species count.
C  KION(*)  - Integer array, non-electron ion species indices.
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
C     Prints species interpreter output and checks for completeness.
C----------------------------------------------------------------------C
C
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
      DIMENSION WT(KDIM), KPHSE(KDIM), KCHRG(KDIM), T(MAXTP,KDIM),
     1          NT(KDIM), KNCF(MDIM,KDIM), KION(KDIM), IPLUS(10)
      LOGICAL KERR, ITHRM(KDIM)
      CHARACTER KNAME(KDIM)*16, ENAME(MDIM)*16, IPHSE(3)*1, INUM(10)*1
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C
      DATA IPHSE/'S','G','L'/,
     1     INUM/'0','1','2','3','4','5','6','7','8','9'/
C
      WRITE (LOUT, 400) (ENAME(M), M = 1, MM)
      WRITE (LOUT, 300)
C
      KELECT = 0
      NELECT = 0
      IF (MELECT .NE. 0) THEN
         DO 20 K = 1, KK
            IF (KNCF(MELECT,K) .NE. 1) GO TO 20
C           this species contains one electron element
            NELEM = 0
            DO 15 M = 1, MM
               IF (M.NE.MELECT) NELEM = NELEM + KNCF(M,K)
   15       CONTINUE
C           species has no elements other than one electron
            IF (NELEM .EQ. 0) THEN
               NELECT = NELECT + 1
               KELECT = K
               IF (NELECT .GT. 1) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, '(6X,A,A)')
     1            'Error...duplicate electron-only species...',
     2            KNAME(K)
               ENDIF
            ENDIF
   20    CONTINUE
      ENDIF
C
      DO 100 K = 1, KK
C
         WRITE (LOUT, 500) K, KNAME(K), IPHSE(KPHSE(K)+2), KCHRG(K),
     1                    WT(K), INT(T(1,K)), INT(T(NT(K),K)),
     2                   (KNCF(M,K),M=1,MM)
         DO 85 N = 2, NT(K)
            IF (T(N,K) .LT. T(N-1,K)) THEN
               WRITE (LOUT, 240)
               KERR = .TRUE.
            ENDIF
   85    CONTINUE
C
C        each species must have thermodynamic data
C
         IF (.NOT. ITHRM(K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, 200)
         ENDIF
C
C        a species cannot start with a number
C
         CALL CKCOMP (KNAME(K)(1:1), INUM, 10, I)
         IF (I .GT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, 210)
         ENDIF
C
C        if '+' sign is used in a species name,
C           examples of legal species symbols with + are:
C           OH(+)2, OH(+2), OH+, OH++, OH+++, OH(+), OH(++),
C           OH[+OH], OH2+, OH+2
C
C           examples of illegal species symbols with + are:
C           +OH        (symbol starts with a +, this will cause
C                       confusion in a reaction)
C           OH(+OH)    (symbol in parentheses is another species-
C                       this arrangement is reserved for a
C                       pres dependency reaction)
C           OH+OH      (plus delimits other species names, this
C                       will cause confusion in a reaction)
C
         NPLUS = 0
         ILAST = CKLSCH(KNAME(K))
         DO 50 N = 1, ILAST
            IF (KNAME(K)(N:N) .EQ. '+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = N
            ENDIF
   50    CONTINUE
         DO 60 N = 1, NPLUS
            I1 = IPLUS(N)
            IF (I1 .EQ. 1) THEN
               WRITE (LOUT, 220)
               KERR = .TRUE.
            ELSE
C
C              is there another species name in parentheses
C
               IF (KNAME(K)(I1-1:I1-1) .EQ. '(') THEN
                  I1 = I1 + 1
                  I2 = I1 + INDEX(KNAME(K)(I1:),')')-1
                  IF (I2 .GT. I1) THEN
                     CALL CKCOMP (KNAME(K)(I1:I2-1), KNAME, KK, KNUM)
                     IF (KNUM .GT. 0) THEN
                        WRITE (LOUT, 230)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
C
C              is there another species name after a +
C
               I1 = I1 + 1
               IF (N .LT. NPLUS) THEN
                  DO 55 L = N+1, NPLUS
                     I2 = IPLUS(L)
                     IF (I2 .GT. I1) THEN
                        CALL CKCOMP (KNAME(K)(I1:I2-1),KNAME,KK,KNUM)
                        IF (KNUM .GT. 0) THEN
                           WRITE (LOUT, 230)
                           KERR = .TRUE.
                        ENDIF
                     ENDIF
   55             CONTINUE
               ENDIF
C
               I2 = CKLSCH(KNAME(K))
               IF (I2 .GE. I1) THEN
                  CALL CKCOMP (KNAME(K)(I1:I2), KNAME, KK, KNUM)
                  IF (KNUM .GT. 0) THEN
                     WRITE (LOUT, 230)
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
   60    CONTINUE
C
         IF (KCHRG(K) .NE. 0) THEN
            NCHRG = NCHRG + 1
            IF (K .NE. KELECT) THEN
               NKION = NKION+1
               KION(NKION) = K
            ENDIF
         ENDIF
C
  100 CONTINUE
C
      WRITE (LOUT, 300)
C
  200 FORMAT (6X,'Error...no thermodynamic properties for species')
  210 FORMAT (6X,'Error...species starts with a number')
  220 FORMAT (6X,'Error...species starts with a plus')
  230 FORMAT (6X,'Error...illegal + in species name')
  240 FORMAT (6X,'Error...fit temperatures not in ascending order')
  300 FORMAT (1X,79('-'))
C
  400 FORMAT (1X,79('-'),/T27,'C',/T24,'P  H',/T24,'H  A',/T24,'A  R',
     1       /1X,'SPECIES',T24,'S  G',T30,'MOLECULAR',T41,'TEMPERATURE',
     2       T54,'ELEMENT COUNT',
     3       /1X,'CONSIDERED',T24,'E  E',T30,'WEIGHT',T41,'LOW',
     4       T48,'HIGH',T54,15(A3))
  500 FORMAT (1X,I3,'. ',A16,T24,A1,T26,I2,T29,F10.5,T39,I6,T46,I6,
     1       T53,15(I3))
C
C     end of SUBROUTINE CKPRNT
      RETURN
      END
C                                                                      C
      SUBROUTINE CKREAC (KDIM, IDIM, LINE, KNAME, LOUT, MAXSP, NSPEC,
     1                   NREAC, NUNK, NU, NPAR, PAR, ITHB, IFAL, KFAL,
     2                   IWL, WL, IRNU, RNU, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKREAC
C  processes the character string representing a gas-phase reaction.
C  It is expected that the character string will satisfy the following
C  conditions:
C
C  a)  the final 3 space-delimted substrings will represent the
C      Arrhenius rate parameters,
C  b)  the rest of the character string will contain a delimeter
C      '=' or '<=>' separating reactants from products in a
C      reversible reaction, or '=>' in an irreversible reaction,
C  c)  each product or reactant will be delimited by a '+' plus-sign,
C  d)  each product or reactant must have previously been listed in
C      the SPECIES section of the mechanism,
C  e)  each product or reactant may be prefixed by the character-string
C      representation of an integer or real value, its stoichiometric
C      coefficient, else its default value is 1,
C  f)  number of species in a reaction is limited to MAXSP,
C      number of reactants or products is limited to MAXSP/2.
C  g)  if a reactant or product is "M", reaction is considered a
C      third-body reaction; "M" must be both a reactant and a product.
C  h)  if "+reactant" or "+product" is parenthesis-enclosed, i.e.
C      "(+M)" or "(+KNAME)", reaction is considered both a
C      pres dependency reaction and a third-body reaction, and
C      the species must be both a reactant and a product.
C  i)  if a reactant or product is "HV", reaction is considered a
C      radiation-enhanced reaction; unless assigned a value with
C      auxiliary reaction information, the wavelength is -1.0 if
C      a reactant or +1.0 if a product.
C
C  Arguments:
C  LINE     - Character string.
C  II       - Integer scalar, reaction count and index of
C             current reaction.
C  KK       - Integer scalar, mechanism species count.
C  KNAME(*) - Character string array, species names.
C  LOUT     - Integer scalar, formatted output file unit number.
C  MAXSP    - Integer scalar, maximum number of reaction species.
C  NSPEC    - Integer scalar, reaction reactants+products count.
C  NR       - Integer scalar, reaction reactants only count.
C  NUNK(*)  - Integer array, reactions species indices.
C  NU(*)    - Integer array, reaction's stoichiometric coefficients.
C  NPAR     - Integer scalar, number of parameters required for
C             a reaction's Arrhenius rate expression.
C  PAR(*)   - Real array, reactions Arrhenius rate parameters.
C  NTHB     - Integer scalar, third-body reactions count.
C  ITHB(*)  - Integer array, third-body reaction indices.
C  NFAL     - Integer scalar, pres dependency reactions count.
C  IFAL(*)  - Integer array, pressure-dependency reaction indices.
C  KFAL(*) - Integer array, pres dependency 3rd body species indices.
C  NWL      - Integer scalar, radiation-enhanced reactions count.
C  IWL(*)   - Integer array, radiation-enhanced reaction indices.
C  WL(*)    - Real array, radiation wavelengths.
C  NRNU     - Integer scalar, real stoichiometry reactions count.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  KERR     - Logical error flag.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Integer arrays
      DIMENSION NUNK(MAXSP), NU(MAXSP), IFAL(IDIM), KFAL(IDIM),
     1          ITHB(IDIM), IWL(IDIM), IRNU(IDIM), IPLUS(20)
C     Real arrays
      DIMENSION PAR(NPAR), WL(IDIM), RNU(MAXSP,IDIM)
     1
      CHARACTER KNAME(KDIM)*16, LINE*(*), CKCHUP*4
      CHARACTER*80 ISTR, RSTR, PSTR, ISPEC, INAME
      LOGICAL KERR, IERR, LTHB, LFAL, LWL, LREV, LRNU
      INTEGER CKLSCH
      EXTERNAL CKLSCH, CKCHUP
C
      IERR = .FALSE.
      LTHB = .FALSE.
      LFAL = .FALSE.
      LWL  = .FALSE.
      LREV = .TRUE.
      LRNU = .FALSE.
      NSPEC = 0
C
      NWANT = NPAR
      ILAST = CKLSCH(LINE)
      DO 5 L = ILAST, 1, -1
         RSTR = ' '
         IF (LINE(L:).NE.' ' .AND. LINE(L-1:L-1).EQ.' ') THEN
            RSTR = LINE(L:)
            CALL CKPARR (RSTR, -1, 1, PAR(NWANT), NVAL, IER, LOUT)
            IF (IER .NE. 0) IERR = .TRUE.
C           done with these characters
            LINE(L:) = ' '
            NWANT = NWANT - 1
            IF (NWANT .EQ. 0) GO TO 6
         ENDIF
    5 CONTINUE
    6 CONTINUE
C
      INAME = ' '
      ILEN = 0
      ILAST = CKLSCH(LINE)
      DO 10 I = 1, ILAST
C        compressing reaction string to remove blanks
         IF (LINE(I:I) .NE. ' ') THEN
            ILEN = ILEN+1
            INAME(ILEN:ILEN) = LINE(I:I)
         ENDIF
   10 CONTINUE
C
C     finding delimeter to separate reaction string, product string
      I = INDEX(INAME,'<=>')
      PSTR = ' '
      RSTR = ' '
      IF (I .GT. 0) THEN
         RSTR = INAME(1:I-1)
         PSTR = INAME(I+3:CKLSCH(INAME))
      ELSE
         I = INDEX(INAME,'=>')
         IF (I .GT. 0) THEN
            RSTR = INAME(1:I-1)
            PSTR = INAME(I+2:ILEN)
            LREV = .FALSE.
         ELSE
            I = INDEX(INAME,'=')
            IF (I .GT. 0) THEN
               RSTR = INAME(1:I-1)
               PSTR = INAME(I+1:ILEN)
            ENDIF
         ENDIF
      ENDIF
C
      IF (ILEN.GE.45 .AND. I.GT.0) THEN
         WRITE (LOUT, 1900) II,RSTR, (PAR(N),N=1,NPAR)
         WRITE (LOUT, 1920) INAME(I:)
      ELSE
          WRITE (LOUT, 1900) II,INAME(1:45),(PAR(N),N=1,NPAR)
      ENDIF
C
      IF (NWANT .NE. 0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...problem finding Arrhenius parameters...'
      ENDIF
C
      IF (I .LE. 0) THEN
         WRITE (LOUT, '(6X,A)')
     1   'Error...did not find delimiter <=>, =>, or =...'
         IERR = .TRUE.
         RETURN
      ENDIF
C
      IF (INDEX(RSTR,'=').GT.0 .OR. INDEX(PSTR,'=').GT.0) THEN
         IERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...more than one reaction delimeter...'
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
C-----Is this a pres dependency reaction?
C     Find exactly one '(+M)' or '(+KNAME)' as a reactant
C
      KRTB = -1
C     KRTB negative to signify NO 3rd body at all
C     KRTB zero to signify "(+M)" reactant
C     KRTB positive to signify "(+species)" reactant
C
   45 CONTINUE
C
      I1 = INDEX(RSTR,'(+')
      I2 = I1 + INDEX(RSTR(I1+1:),')')
      IF (I1.GT.0 .AND. I2.GT.I1) THEN
         ISTR = ' '
         ISTR = RSTR(I1+2:I2-1)
         IF (LFAL) THEN
C           already had either a "(+reactant)" or "(+M)"
            IERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...multiple pres-dependency declaration...',
     2      RSTR(I1:I2)
         ELSE
C           flag this reaction for a pressure-dependency
            LFAL = .TRUE.
         ENDIF
         IF (ISTR .EQ. 'm' .OR. ISTR.EQ.'M') THEN
            IF (LTHB) THEN
C              already had a "+M" (either with or without parentheses)
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1        'Error...multiple "M" 3rd body reactant...',RSTR(I1:I2)
            ELSE
C              flag this reaction for a 3rd body reaction
               LTHB = .TRUE.
            ENDIF
            IF (KRTB .GT. 0) THEN
C              already had a "(+reactant)" (other than M)
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...conflicting pres-dependency reactant...',
     2         RSTR(I1:I2)
            ELSE
C              flag this pressure-dependency reactant as "M"
               KRTB = 0
            ENDIF
         ELSE
C           find which "(+reactant)"
            CALL CKCOMP (ISTR, KNAME, KK, KNUM)
            IF (KNUM .LE. 0) THEN
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error identifying pres-dependency 3rd body reactant...',
     2         RSTR(I1:I2)
            ELSE
C              flag this as the pressure-dependency reactant
               KRTB = KNUM
            ENDIF
         ENDIF
C
C        cancel this string out of reactant string
         ISTR = ' '
         ISTR = RSTR(I2+1:)
         RSTR(I1:) = ' '
         RSTR(I1:) = ISTR
C
C        go back to check rest of reactants
         GO TO 45
      ENDIF
C
C-----Is this a pres dependency reaction?
C     Find exactly one '(+M)' or '(+KNAME)' as a product
C
      KPTB = -1
   55 CONTINUE
C
      I1 = INDEX(PSTR,'(+')
      I2 = I1 + INDEX(PSTR(I1+1:),')')
      IF (I1.GT.0 .AND. I2.GT.I1) THEN
         ISTR = ' '
         ISTR = PSTR(I1+2:I2-1)
         IF (.NOT. LFAL) THEN
C           flag was not set by presence of "(+reactant)" or "(+M)"
            IERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)')
     1      'Error...pres-dependency not initialized in reactants...',
     2      PSTR(I1:I2)
            LFAL = .TRUE.
         ENDIF
         IF (ISTR .EQ. 'm' .OR. ISTR.EQ.'M') THEN
            IF (.NOT. LTHB) THEN
C              flag was not set for presence of "(+M)" reactant
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...no reactant M 3rd body...',PSTR(I1:I2)
               LTHB = .TRUE.
            ENDIF
            IF (KRTB .GT. 0) THEN
C              flag was set for "(+species)" reactant
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1 'Error...conflicts with reactant pres-dependency species...',
     2         PSTR(I1:I2)
            ENDIF
            KPTB = 0
         ELSE
            CALL CKCOMP (ISTR, KNAME, KK, KNUM)
            IF (KNUM .LE. 0) THEN
C              this is not a good "(+product)"
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error identifying pres-dependency 3rd body product...',
     2         PSTR(I1:I2)
            ELSE
               KPTB = KNUM
            ENDIF
         ENDIF
         ISTR = ' '
         ISTR = PSTR(I2+1:)
         PSTR(I1:) = ' '
         PSTR(I1:) = ISTR
         GO TO 55
      ENDIF
C
      IF (LFAL) THEN
C         pressure-dependency was declared
          IF (KRTB.GE.0 .AND. KPTB.GE.0 .AND. KRTB .NE. KPTB) THEN
C            KPTB species index does not match KRTB species index
             IERR = .TRUE.
             WRITE (LOUT, '(6X,A,A)') 'Error...',
     1       'Pres-dependency reactant does not match product...'
          ELSE
             IF (KRTB.LT.0) THEN
C               "(+M)" reactant changes default -1 to 0
                IERR = .TRUE.
                WRITE (LOUT, '(6X,A)')
     1          'Error...did not find pres-dependency reactant...'
             ENDIF
             IF (KPTB.LT.0) THEN
C               "(+M)" product changes default -1 to 0
                IERR = .TRUE.
                WRITE (LOUT, '(6X,A)')
     1          'Error...did not find pres-dependency product...'
             ENDIF
         ENDIF
      ENDIF
C
C----------Find reactants, products-------------------------
C
      NREAC = 0
      NPROD = 0
      DO 600 J = 1, 2
         ISTR = ' '
         IF (J .EQ. 1) THEN
            ISTR = RSTR
         ELSE
            ISTR = PSTR
         ENDIF
C
C        '+'-sign is delimeter between species
         NPLUS = 0
         ILAST = CKLSCH(ISTR)
         DO 100 N = 1, ILAST
            IF (N .EQ. ILAST) THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = ILAST+1
            ELSEIF
     1      (ISTR(N:N).EQ.'+' .AND. ISTR(N+1:N+1).NE.'+') THEN
C              index of the character after the end of a species
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = N
            ENDIF
  100    CONTINUE
C
C        loop over '+'delimited substrings
C
         DO 200 N = 1, NPLUS
C
            IERR = .FALSE.
            IF (N .EQ. 1) THEN
               ISTART = 1
            ELSE
               ISTART = IPLUS(N-1)+1
            ENDIF
C
            ISPEC = ' '
            ISPEC = ISTR(ISTART:IPLUS(N)-1)
            ILEN = CKLSCH(ISPEC)
C
            IF (ISPEC.EQ.'M' .OR. ISPEC.EQ.'m') THEN
               IF (LFAL) THEN
                  IERR = .TRUE.
                  WRITE (LOUT, '(6X,A,A)')
     1            'Error...pres-dependency 3rd body conflict...',
     2            ISPEC(1:ILEN)
               ELSE
                  IF (J.EQ.1 .AND. KRTB.EQ.0) THEN
                     IERR = .TRUE.
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...multiple M 3rd body reactant...',
     1               ISPEC(1:ILEN)
                  ENDIF
                  IF (J.EQ.2 .AND. KPTB.EQ.0) THEN
                     KERR = .TRUE.
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...multiple M 3rd body product...',
     2               ISPEC(1:ILEN)
                  ENDIF
                  IF (J.EQ.2 .AND. KRTB.LT.0) THEN
                     IERR = .TRUE.
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...no M 3rd body reactant...',ISPEC(1:ILEN)
                  ENDIF
                  IF (.NOT. LTHB) LTHB = .TRUE.
                  IF (J .EQ. 1) THEN
                     KRTB = 0
                  ELSE
                     KPTB = 0
                  ENDIF
               ENDIF
               KERR = KERR.OR.IERR
               GO TO 200
            ENDIF
C
            IF (CKCHUP(ISPEC,2) .EQ. 'HV') THEN
               LWL = .TRUE.
               IF (J .EQ. 1) THEN
                  WCOEF = -1.0
               ELSE
                  WCOEF = 1.0
               ENDIF
               GO TO 200
            ENDIF
C
            ISTART = 0
            ILAST = CKLSCH(ISPEC)
            DO 120 I = 1, ILAST
C              look for species match
               CALL CKCOMP (ISPEC(I:), KNAME, KK, KNUM)
C              if species found, go to find a coefficient
               IF (KNUM .GT. 0) THEN
                  ISTART = I
                  GO TO 125
               ENDIF
  120       CONTINUE
C
  125       CONTINUE
C
C           default integer coefficient
            ICOEF = 1
C           default real coefficient
            RCOEF = 1.0
C
            IER = 0
            NVAL = 1
C
C           species name starts at ISTART; if ISTART>1, coefficient
            IF (ISTART .GT. 1) THEN
               IND = INDEX(ISPEC(1:ISTART-1), '.')
               IF (IND .GT. 0) THEN
C                 real coefficient
                  CALL CKPARR (ISPEC(1:ISTART-1), -1, 1, RCOEF, NVAL,
     1                         IER, LOUT)
                  IF (.NOT. LRNU) THEN
C                    initialize this real-coefficient reaction
                     NRNU = NRNU + 1
                     IRNU(NRNU) = II
                     LRNU = .TRUE.
C                    convert any previous coefficients
                     DO 111 L = 1, MAXSP
                        RNU(L, NRNU) = NU(L)
                        NU(L) = 0
  111                CONTINUE
                  ENDIF
C
               ELSE
C                 integer coefficient
                  CALL CKPARI (ISPEC(1:ISTART-1), -1, 1, ICOEF, NVAL,
     1                          IER, LOUT)
                  RCOEF = ICOEF
               ENDIF
            ENDIF
C
            IF (IER.NE.0 .OR. NVAL.NE.1) THEN
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error finding stoichiometric coefficient...',
     2         ISPEC(1:CKLSCH(ISPEC))
            ENDIF
C
            IF (KNUM .EQ. 0) THEN
C              this was not a good species name
               IERR = .TRUE.
               WRITE (LOUT, '(6X,A,A)')
     1         'Error...undeclared species...', ISPEC(1:CKLSCH(ISPEC))
            ENDIF
C
            KERR = KERR.OR.IERR
            IF (IERR) GO TO 200
C
C           to get here, species=KNUM, coef=ICOEF and/or RCOEF
C           find place to put this species
            NPOS = 0
            IF (J .EQ. 1) THEN
C              check previous reactant species
               DO 221 L = 1, NREAC
                  IF (KNUM .EQ. NUNK(L)) NPOS = L
  221          CONTINUE
C
               IF (NPOS .EQ. 0) THEN
C                 this is a new reactant
                  IF (NREAC .EQ. MAXSP/2) THEN
                     IERR = .TRUE.
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...reactant array full, cannot add...',
     2               ISPEC(1:CKLSCH(ISPEC))
                  ELSE
                     NREAC = NREAC + 1
                     NPOS = NREAC
                  ENDIF
               ENDIF
            ELSE
C              check previous product species
               DO 222 L = MAXSP/2 + 1, MAXSP/2 + NPROD
                  IF (KNUM .EQ. NUNK(L)) NPOS = L
  222          CONTINUE
C
               IF (NPOS .EQ. 0) THEN
C                 this is a new product
                  IF (NPROD .EQ. MAXSP/2) THEN
                     IERR = .TRUE.
                     WRITE (LOUT, '(6X,A,A)')
     1               'Error...product array full, cannot add...',
     2               ISPEC(1:CKLSCH(ISPEC))
                  ELSE
                     NPROD = NPROD + 1
                     NPOS = MAXSP/2 + NPROD
                  ENDIF
               ENDIF
            ENDIF
C
            IF (NPOS .GT. 0) THEN
               NUNK(NPOS) = KNUM
               IF (LRNU) THEN
                  RNU(NPOS, NRNU) = RNU(NPOS, NRNU) + RCOEF
               ELSE
                  NU(NPOS) = NU(NPOS) +ICOEF
               ENDIF
            ENDIF
            KERR = KERR.OR.IERR
  200    CONTINUE
  600 CONTINUE
C
      IF (LRNU) THEN
         DO 650 N = 1, NREAC
            RNU(N,NRNU) = -RNU(N,NRNU)
  650    CONTINUE
      ELSE
         DO 700 N = 1, NREAC
            NU(N) = -NU(N)
  700    CONTINUE
      ENDIF
C
      NSPEC = NREAC + NPROD
      IF (.NOT. LREV) NSPEC = -NSPEC
      IF (LFAL) THEN
         NFAL = NFAL + 1
         IFAL(NFAL) = II
         KFAL(NFAL) = KRTB
      ENDIF
      IF (LTHB) THEN
         NTHB = NTHB + 1
         ITHB(NTHB) = II
      ENDIF
      IF (LWL) THEN
         NWL = NWL + 1
         IWL(NWL) = II
         WL(NWL) = WCOEF
      ENDIF
C
 1900 FORMAT (I4,'. ', A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 1920 FORMAT (6X,A)
C
C     end of SUBROUTINE CKREAC
      RETURN
      END
C                                                                      C
      SUBROUTINE CKSET (MDIM, KDIM, IDIM, NPAR, NCP, NFIT, MAXTP, NTR,
     1                  MAXSP, MAXTB, NLAR, NFAR, MAXORD, KNAME, ENAME,
     2                  AWT, KNCF, WT, KPHSE, KCHRG, A, T, NT, NSPEC,
     3                  NREAC, NU, NUNK, PAR, IDUP, IREV, RPAR, ILAN,
     4                  PLAN, IRLT, RLAN, IWL, WL, IFAL, IFOP, IFLO,
     5                  KFAL, PFAL, ITHB, NTBS, AIK, NKTB, IRNU, RNU,
     6                  IORD, KORD, RORD, KION, NJAR, NF1R, IEIM,
     7                  IEIMT, IJAN, PJAN, IFT1, PFT1, IEXC, PEXC,
     8                  IMOM, KMOM, IXSM, IXSK, IXSI, ITDE, ITDK, ITHRM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSET
C  sets initial values of arrays and matrices.
C
C  Arguments:
C  MDIM     - Integer scalar, maximum number of elements.
C  KDIM     - Integer scalar, maximum number of species.
C  IDIM     - Integer scalar, maximum number of reactions.
C  NPAR     - Integer scalar, required number of Arrhenius rate
C             parameters.
C  NCP      - Integer scalar, number of CP-fit polynomial coefficients.
C  NFIT    - Integer scalar, number of thermo polynomial coefficients.
C  MAXTP    - Integer scalar, maximum number of fit temperatures.
C  NTR      - Integer scalar, number of fit temperature ranges.
C  MAXSP    - Integer scalar, maximum reaction species allowed.
C  MAXTB    - Integer scalar, maximum reaction third-bodies allowed.
C  NLAR     - Integer scalar, required number of Landau-Teller params.
C  NFAR     - Integer scalar, maximum pres dependency parameters.
C  MAXORD   - Integer scalar, maximum reaction change-order species.
C  KNAME(*) - Character string array, species names.
C  ENAME(*) - Character string array, element names.
C  AWT(*)   - Real array, element atomic weights.
C  KNCF(*,*)- Integer matrix, species elemental composition.
C  WT(*)    - Real array, species molecular weights.
C  KPHSE(*) - Integer array, species physical state.
C  KCHRG(*) - Integer array, species ionic charges.
C  A(*,*,*) - Real 3-dimensional array, species thermodynics
C             polynomial coefficients.
C  T(*,*)   - Real matrix, species fit temperatures.
C  NT(*)    - Integer array, species number of fit temperatures.
C  NSPEC(*) - Integer array, reactions' species counts.
C  NREAC(*) - Integer array, reactions' reactants only counts.
C  NU(*,*)  - Integer matrix, reaction species' stoichiometric
C             coefficients.
C  NUNK(*,*)- Integer matrix, reaction species' indices.
C  PAR(*,*) - Real matrix, reaction rate parameters.
C  IDUP(*)  - Integer array, flag for reaction duplication.
C  IREV(*)  - Integer array, indices of reactions with explicit
C             reverse rate parameters.
C  RPAR(*,*)- Real matrix, explicit reverse rate parameters, if given.
C  ILAN(*)  - Integer array, Landau-Teller reaction indices.
C  PLAN(*,*)- Real matrix, Landau-Teller reaction parameters.
C  IRLT(*)  - Integer array, reaction indices, Landau-Teller
C             reactions with explicit reverse Landau-Teller
C             parameters.
C  RLAN(*,*)- Real matrix, Landau-Teller explicit reverse
C             parameters.
C  IWL(*)   - Integer array, radiation-enhanced reaction indices.
C  WL(*)    - Real array, radiation wavelengths.
C  IFAL(*)  - Integer array, pres dependency reaction indices.
C  IFOP(*)  - Integer array, pres dependency formulation flag.
C  IFLO(*)  - Integer array, pres dependency region flag.
C  KFAL(*)  - Integer array, pres dependency third-body species.
C  PFAL(*,*)- Real matrix, pres dependency additional rate
C             parameters.
C  ITHB(*)  - Integer array, third-body reaction indices.
C  NTBS(*)  - Integer array, enhanced third-body count for third-body
C             reactions.
C  AIK(*,*) - Real matrix, species enhancement factors for third-body
C             reactions.
C  NKTB(*,*)- Integer matrix, indices for third-body reaction
C             enhanced species.
C  IRNU(*)  - Integer array, real stoichiometry reaction indices.
C  RNU(*,*) - Real matrix, real stoichiometric coefficients.
C  IORD(*)  - Integer array, changed-order reaction indices.
C  KORD(*,*)- Integer matrix, changed-order species indices.
C---------------------------------
C  RORD(*,*)- Real matrix, changedpo
C  KION(*)
C  NJAR
C  NF1R
C  IEIM(*)
C  IEIMT(*)
C  IJAN(*)
C  PJAN(*,*)
C  IFT1(*)
C  PFT1(*,*)
C  IEXC(*)
C  PEXC(*)
C  IMOM(*)
C  KMOM(*)
C  IXSM(*)
C  IXSK(*)
C  IXSI(*)
C  ITDE(*)
C  ITDK(*)
C  ITHRM(*)
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Integer arrays
      DIMENSION KNCF(MDIM,KDIM), KPHSE(KDIM), KCHRG(KDIM), NT(KDIM),
     1          NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     2          NUNK(MAXSP,IDIM), IDUP(IDIM), IREV(IDIM), ILAN(IDIM),
     3          IRLT(IDIM), IWL(IDIM), IFAL(IDIM), IFOP(IDIM),
     4          IFLO(IDIM), KFAL(IDIM), ITHB(IDIM), NTBS(IDIM),
     5          NKTB(MAXTB,IDIM), IRNU(IDIM), IORD(IDIM),
     6          KORD(MAXORD,IDIM), KION(KDIM), IEIM(IDIM),
     7          IEIMT(IDIM), IJAN(IDIM), IFT1(IDIM), IEXC(IDIM),
     8          IMOM(IDIM), KMOM(IDIM), IXSM(IDIM), IXSK(IDIM),
     9          IXSI(IDIM), ITDE(IDIM), ITDK(IDIM)
C     Real arrays
      DIMENSION AWT(MDIM), WT(KDIM), A(NFIT,NTR,KDIM), T(MAXTP,KDIM),
     1          PAR(NPAR,IDIM), RPAR(NPAR,IDIM), PLAN(NLAR,IDIM),
     2          RLAN(NLAR,IDIM), WL(IDIM), PFAL(NFAR,IDIM),
     3          AIK(MAXTB,IDIM), RNU(MAXSP,IDIM), RORD(MAXORD,IDIM),
     4          PJAN(NJAR,IDIM), PFT1(NF1R,IDIM), PEXC(IDIM)
C
      LOGICAL ITHRM(KDIM)
      CHARACTER*16 KNAME(KDIM), ENAME(MDIM)
C
C
C     initialize counters
      LENICK = 0
      LENRCK = 0
      LENCCK = 0
      MM = 0
      KK = 0
      II = 0
      NREV = 0
      NFAL = 0
      NTHB = 0
      NLAN = 0
      NRLT = 0
      NWL = 0
      NCHRG = 0
      NEIM = 0
      NJAN = 0
      NFT1 = 0
      NEXC = 0
      NMOM = 0
      NXSM = 0
      NTDE = 0
      NRNU = 0
      NORD = 0
      KELECT = 0
      NKION = 0
C
C     initialize atomic weights and element names
C
      DO 5 M = 1, MDIM
         AWT(M) = 0.0
         ENAME(M) = ' '
    5 CONTINUE
C
      DO 100 K = 1, KDIM
         KNAME(K) = ' '
C        number of temperatures used to define polynomial fits,
C        i.e., TLOW, TMID, THIGH in the case of MAXTP=3
         NT(K)  = MAXTP
         KION(K)  = 0
C        initialize temperature bounds to a negative
C        number as a flag to indicate that they haven't been set
C
         DO 10 M = 1, MAXTP
            T(M,K) = -1.0
   10    CONTINUE
C
C        "A(m,n,k)" contains thermo polynomial coefficients;
C        first dimension runs over NFIT number of coefficients,
C        second dimension runs over temperature ranges,
C        third dimension is for species indices;
C        initialize values to zero
C
         DO 15 M = 1, NFIT
            DO 14 N = 1, NTR
               A(M, N, K) = 0.0
   14       CONTINUE
   15    CONTINUE
C
C        species electronic charge
         KCHRG(K) = 0
C
C        species physical states
         KPHSE(K) = 0
C
C        species molecular weight
         WT(K) = 0.0
C
C        elemental composition; first dimension runs over
C        elements, second is for species indices.
         DO 20 M = 1, MDIM
            KNCF(M,K) = 0
   20    CONTINUE
C
C        thermo fits have not yet been read for this species
         ITHRM(K) = .FALSE.
C
  100 CONTINUE
C
C     initialize information concerning reactions
      DO 300 I = 1, IDIM
C
C        count of reactants plus products in a reaction
C        (is set to a negative number if reaction is irreversible)
         NSPEC(I) = 0
C
C        count of reactants in reaction I
         NREAC(I) = 0
C
C        MAXSP is the maximum number of species per reaction
         DO 205 M = 1, MAXSP
C
C           stoichiometric coefficient of Mth reactant in
C           Ith reaction
            NU(M,I) = 0
C
C           species number corresponding to NU(M,I)
            NUNK(M,I) = 0
  205    CONTINUE
C
C        rate coefficients for reaction I
         DO 215 N = 1, NPAR
            PAR(N,I) = 0.0
  215    CONTINUE
C
C        reaction I is duplicate reaction
         IDUP(I) = 0
C
C        reaction I has declared reverse rate coefficients
         IREV(I) = 0
C
C        reverse rate coefficients for reaction I
         DO 222 N = 1, NPAR
            RPAR(N,I) = 0.0
  222    CONTINUE
C
C        reaction I has Landau-Teller coefficients
         ILAN(I) = 0
C        reaction I has reverse Landau-Teller coefficients
         IRLT(I) = 0
C
C        Landau-Teller parameters
         DO 225 N = 1, NLAR
            PLAN(N,I) = 0.0
            RLAN(N,I) = 0.0
  225    CONTINUE
C
C        reaction I has radiation wavelength
         IWL(I) = 0
C        radiation wavelength
         WL(I) = 0.0
C
C        reaction I is a pres dependency reaction
         IFAL(I) = 0
C        type of pres dependency for reaction I
         IFOP(I) = 0
C        high/low pres dependency flag
         IFLO(I) = -1
C        3rd body species index for pres dependency reaction I
         KFAL(I) = 0
C        pres dependency rate parameters for reaction I
         DO 240 N = 1, NFAR
            PFAL(N,I) = 0.0
  240    CONTINUE
C
C        reaction I is a 3rd-body reaction
         ITHB(I) = 0
C        count of enhanced 3rd bodies in reaction I
         NTBS(I) = 0
C        species index numbers and enhancement factors
         DO 250 N = 1, MAXTB
            NKTB(N,I) = 0
            AIK(N,I) = 0.0
  250    CONTINUE
C
C        reaction I is an xxxxxxxxxxx
         IEIM(I) = 0
         IEIMT(I) = 0
C
C        Jannev-xxxxxx
         IJAN(I) = 0
         DO 260 N = 1, NJAR
            PJAN(N,I) = 0.0
  260    CONTINUE
C
         IFT1(I) = 0
         DO 270 N = 1, NF1R
            PFT1(N,I) = 0.0
  270    CONTINUE
C
         IEXC(I) = 0
         PEXC(I) = 0.0
C
         IMOM(I) = 0
         KMOM(I) = 0
         IXSM(I) = 0
         IXSK(I) = 0
         IXSI(I) = 0
C
         ITDE(I) = 0
         ITDK(I) = 0
C
C        real stoichiometric coefficients are used for reaction I
         IRNU(I) = 0
C
C        real coefficients of species in reaction I
         DO 230 M = 1, MAXSP
            RNU(M,I) = 0.0
  230    CONTINUE
C
C        a non-standard reaction order dependency was specified for
C        reaction I
         IORD(I) = 0
         DO 235 M = 1, MAXORD
C
C           species number for which a non-standard reaction order
C           applies
            KORD(M,I) = 0
C
C           reaction order for that species
            RORD(M,I) = 0.0
  235    CONTINUE
C
  300 CONTINUE
C
C     end of SUBROUTINE CKSET
      RETURN
      END
C                                                                      C
      SUBROUTINE CKSIZE (NPAR, MAXSP, MAXTB, MAXTP, NTR, NFIT, NFAR,
     1                   NLAR, NJAR, NF1R, MAXORD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSIZE
C  computes the size required for CHEMKIN work arrays.
C
C  Arguments:
C  NPAR     - Integer scalar, required number of Arrhenius parameters.
C  MAXSP    - Integer scalar, maximum reaction species.
C  MAXTB    - Integer scalar, maximum reaction 3rd bodies.
C  MAXTP    - Integer scalar, maximum species fit temperatures.
C  NTR      - Integer scalar, number of species fit temperature ranges.
C  NFIT     - Integer scalar, number of species thermodynamic fit
C             polynomial coefficients.
C  NLAR     - Integer scalar, required number of Teller-Landau
C             reaction parameters.
C  NJAR     - Integer scalar, required number of Janev-Langer
C             reaction parameters.
C  NF1R     - Integer scalar, required number of fit#1 reaction
C             parameters.
C  MAXORD   - Integer scalar, maximum reaction species change-orders.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Size of integer work space in chemkin:
      LENICK = 1
C              (CKINIT flag)
     1       + KK*(2 + MM)
C              (species phases and charges, elemental composition)
     2       + NKION
C              (ionic species indices)
     3       + KK
C              (number of fit temperatures)
     4       + II*(2 + 2*MAXSP)
C              (nspecies, nreac, nu, nunk)
     5       + NLAN
C              (Landau-Teller reaction indices)
     6       + NRLT
C              (Reverse Landau-Teller reaction indices)
     7       + NFAL*4
C              (pres dependency rxn, species, high/low flag, rate type)
     8       + NTHB*(2 + MAXTB)
C              (3rd body reaction, species count, species indices)
     9       + NREV
C              (reverse rate reaction indices)
     *       + NWL
C              (radiation-enhanced reaction indices)
     1       + NEIM*2
C              (electron-impact reaction, temp dependency)
     2       + NTDE*2
C              (non-thermal-equilibrium reaction, species indices)
     2       + NJAN
C              (Jannev-et al. reaction indices)
     3       + NFT1
C              (Fit1 reaction indices)
     4       + NEXC
C              (excitation reaction indices)
     5       + NMOM*2
C              (ion-momentum transfer reaction, species partner indices)
     6       + NXSM*3
C              (Xsection reaction, ion, non-ion indices)
     7       + NRNU
C              (real stoichiometry reaction indices)
     8       + NORD*(1 + MAXORD)                  + KK*3
C              (change-order reactions, species,  species scratch space)
C
C     Size of character*16 work space in chemkin:
      LENCCK = MM + KK
C
C     Size of real*8 work space in chemkin:
      LENRCK = MM
C              (atomic weights)
     1       + KK*(1 + MAXTP + NFIT*NTR)
C              (species molecular weights, fit temperatures, thermo)
     2       + II*(NPAR+1)
C              (rate parameters)
     3       + NREV*(NPAR+1)
C              (reverse parameters)
     4       + NLAN*NLAR
C              (Landau-Teller parameters)
     5       + NRLT*NLAR
C              (reverse Landau-Teller parameters)
     6       + NFAL*NFAR
C              (pres dependency parameters)
     7       + NTHB*MAXTB
C              (3rd-body enhancement factors)
     8       + NWL
C              (radiation enhancement factors)
     9       + NJAN*NJAR
C              (Jannev et al. parameters)
     *       + NFT1*NF1R
C              (Fit1 parameters)
     1       + NEXC
C              (energy loss, excitation reactions)
     2       + NRNU*MAXSP
C              (real stoichiometric coefficients)
     3       + NORD*MAXORD
C              (changed-order species values)
     4       + 3
C              (gas constants)
     5       + II*2
C              (forward/reverse temperature-dependent rates)
     6       + KK*4
C              (species scratch space)
     7       + II*4
C              (reaction scratch space)
C
C     end of SUBROUTINE CKSIZE
      RETURN
      END
C
      SUBROUTINE CKTHRM (LUNIT, LTHRM, THERMO, MDIM, ENAME, AWT, KNAME,
     1                   KNCF, KPHSE, KCHRG, WT, MAXTP, NT, NTR,
     2                   T, NFIT, A, ITHRM, KERR, LOUT, ISTR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKTHRM
C  processes formatted thermodynamic properties inputs of species;
C  format -
C  Line 1: species name, optional comments, elemental composition,
C          phase, T(low), T(high), T(mid), additional elemental
C          composition, line number (col. 80);
C          format(A10,A14,4(A2,I3),A1,E10.0,E10.0,E8.0,(A2,I3),I1)
C  Line 2: coefficients a(1--5) for upper temperature range,
C          line number (col. 80); format (5(e15.0), I1)
C  Line 3: coefficients a(6--7) for upper temperature range,
C          coefficients a(1--3) for lower temperature range,
C          line number (col. 80); format (5(e15.0), I1)
C  Line 4: coefficients a(4--7) for lower temperature range,
C          line number (col. 80); format (4(e15.0), I1)
C
C  Arguments:
C  LUNIT     - Integer scalar, formatted input file unit number.
C  LTHRM     - Integer scalar, formatted input file unit number.
C  THERMO    - Logical, thermodynamic processing status flag.
C  MDIM      - Integer scalar, maximum number of elements.
C  ENAME(*)  - Character-string array, element names.
C  AWT(*)    - Real array, element atomic weights.
C  KNAME(*)  - Character-string array, species names.
C  KNCF(*,*) - Integer matrix, elemental composition of species.
C  KPHSE(*)  - Integer array, species physical states.
C  WT(*)     - Real array, species molecular weights.
C  MAXTP     - Integer scalar, maximum number of species fit temps.
C  NT(*)     - Integer array, count of species fit temperatures.
C  MAXTR     - Integer scalar, number of fit temperature ranges.
C  T(*,*)    - Real matrix, species fit temperatures.
C  NFIT      - Integer scalar, number of species fit coefficients.
C  A(*,*,*)  - Real three-dimensional array, species thermodynamic
C              polynomial coefficients.
C  ITHRM(*)  - Logical array, species thermodynamics status flag.
C  KERR      - Logical, error flag.
C  LOUT      - Integer scalar, formatted output file unit number.
C  ISTR      - Character string, the final string input to CKTHRM.
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
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
      DIMENSION WT(KK), NT(KK), T(MAXTP,KK), KPHSE(KK), KNCF(MDIM,*),
     1          KCHRG(KK), A(NFIT,NTR,KK), AWT(MM), VALUE(5)
C
C     local temperature array
      PARAMETER (MAXTMP=10)
      DIMENSION TFIT(MAXTMP)
      CHARACTER ENAME(MM)*16, KNAME(KK)*16, ISTR*(*), SUB(80)*80,
     1          LINE(15)*80, ELEM*2, CKCHUP*4, KEY*4
      LOGICAL KERR, ITHRM(KK), THERMO, SETTMP
      INTEGER CKSLEN
      EXTERNAL CKSLEN, CKCHUP
C
      IF (MM .LE. 0) THEN
         WRITE (LOUT, '(6X,A)')
     1   'Error...cannot use THERM until ELEMents have been input...'
         KERR = .TRUE.
      ENDIF
      IF (KK .LE. 0) THEN
         WRITE (LOUT, '(6X,A)')
     1   'Error...cannot use THERM until SPECies have been input...'
         KERR = .TRUE.
      ENDIF
C      IF (KERR) RETURN
C
      IF (THERMO) THEN
         REWIND LTHRM
         LTEMP = LTHRM
      ELSE
C        THERMO ALL option eliminates need for CHEMKIN thermo data
         LTEMP = LUNIT
      ENDIF
C
  300 CONTINUE
C     stepping through THERMO data
C     (even if no species, need to get through thermo lines)
      ISTR = ' '
      READ (LTEMP,'(A)',ERR=22222, END=400) ISTR
      IF (CKSLEN(ISTR) .LE. 0) GO TO 300
C
C     CHEMKIN data has the THERM line to contend with
      IND = MAX (INDEX(ISTR,'THERM'), INDEX(ISTR,'therm'))
      IF (IND .GT. 0) GO TO 300
C
C     first non-blank line after THERM has NTFIT fit temperatures
C     (LTEMP here could be either user's input or CHEMKIN data file);
C     NTFIT is the number of temperatures (nominally 3) on the line,
C     must be <= MAXTP set by calling program, and in ascending order.
C     The species data which follows must then have 7 x NTFIT
C     coefficients, up to 5 per line, from HIGHEST to LOWEST
C     temperature ranges.
C     (To change the requirements temporarily for a specific species,
C     add a line TEMPS value(1) ... value(n), after Line 1.)
C
C     initial array of fit temperatures
      DO 5 N = 1, MAXTMP
         TFIT(N) = -1.0
    5 CONTINUE
      CALL CKDTAB(ISTR)
C
C     fit temperatures to be used as defaults
      CALL CKPARR (ISTR, -1, -MAXTMP, TFIT, NTFIT, IER, LOUT)
      IF (NTFIT .GT. MAXTP) THEN
         WRITE (LOUT, '(6X,A,I3)')
     1   'Error...must increase MAXTP parameter to at least ',NTFIT
         KERR = .TRUE.
         NTFIT = MAXTP
      ENDIF
      IF (IER .NE. 0) THEN
         WRITE (LOUT, '(6X,A,A)')
     1   'Error reading default temperatures...',ISTR
         KERR = .TRUE.
      ENDIF
C
C     temperatures given must be in ascending order
      DO 10 N = 2, NTFIT
         IF (TFIT(N) .LE. TFIT(N-1)) THEN
            WRITE (LOUT, '(6X,A)')
     1      'Error...fit temperatures are not in ascending order'
            KERR = .TRUE.
         ENDIF
   10 CONTINUE
C
C     Species thermodynamic data from here on,
C     LUNIT could be either user's input or CHEMKIN data file.
C
   25 CONTINUE
C     initialize species data storage
      NLINES = 0
      KNUM = 0
      SETTMP = .FALSE.
C
   50 CONTINUE
      ISTR = ' '
      READ (LUNIT,'(A)', ERR=22222, END=400) ISTR
C     skip blank lines
      CALL CKDTAB (ISTR)
      ILEN = CKSLEN(ISTR)
      IF (ILEN .LE. 0) GO TO 50
C
      CALL CKISUB (ISTR(1:ILEN), SUB, NSUB)
      KEY = ' '
      KEY = CKCHUP(SUB(1), 4)
C
      IND = MAX (INDEX(KEY,'END'), INDEX(KEY,'REAC'),
     1           INDEX(KEY,'ELEM'),INDEX(KEY,'SPEC'))
      IF (IND .GT. 0) THEN
C        end of thermo data
         IF (NLINES .EQ. 0) RETURN
C        still need to process previous species, K
         BACKSPACE (LUNIT)
C        500 is the main species processing section
         GO TO 500
      ENDIF
C
      IF (ILEN.GE.80 .AND. ISTR(80:80) .EQ. '1') THEN
C        '1' in col. 80 starts a set of species data
C
         IF (KNUM .GT. 0) THEN
C           have previously-stored species data to process
C           before proceeding to check this species
            BACKSPACE (LUNIT)
            GO TO 500
         ENDIF
C        new species data; compare to those in mechanism
         CALL CKNCMP (SUB(1), KNAME, KK, KNUM, KTIMES)
C        this species is not used
         IF (KNUM .LE. 0) GO TO 50
C        this species already has thermo data
         IF (ITHRM(KNUM)) GO TO 25
C        initialize data for the new species
         NLINES = 1
         LINE(NLINES) = ISTR(1:ILEN)
C        look for more input for this species
         GO TO 50
      ENDIF
C
      IF (KNUM .GT. 0) THEN
C        continue storing data for this species
         NLINES = NLINES + 1
         LINE(NLINES) = ISTR(1:ILEN)
         IF (KEY .EQ. 'TEMP') SETTMP = .TRUE.
C        look for more input for this species
      ENDIF
C     next input
      GO TO 50
C
  500 CONTINUE
C
C     Line 1 Input:
C     a) get/store elemental composition;
C     b) use that to calculate/store molecule mass and charge,
C     c) get/store integer representation of molecule phase
C     d) get/store T(1), T(3), T(2) for this particular species,
C        (used only if present and if NTFIT=3)
C
      ICOL = 25
      DO 110 N = 1, 5
         ELEM = ' '
         ELEM = LINE(1)(ICOL:ICOL+1)
         CALL CKPARR (LINE(1)(ICOL+2:ICOL+4), -1, 1, VALUE, NVAL,
     1                IER, LOUT)
         IF (ELEM.NE.' ' .AND. NVAL.GT.0) THEN
            IELEM = INT (VALUE(1))
            IF (IELEM .NE. 0) THEN
               CALL CKCOMP (ELEM, ENAME, MM, M)
               IF (M .GT. 0) THEN
C                 composition
                  KNCF(M,KNUM) = KNCF(M,KNUM) + IELEM
C                 molecular weight
                  WT(KNUM) = WT(KNUM) + AWT(M) * VALUE(1)
C                 electronic charge
                  IF (M .EQ. MELECT) KCHRG(KNUM)=KCHRG(KNUM)-VALUE(1)
               ELSE
                  WRITE (LOUT, '(6X,A,A,A)')
     1            'Error...element ', ELEM, 'not declared for species ',
     2            KNAME(KNUM)(1:10)
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
         IF (N .LT. 4) THEN
            ICOL = ICOL + 5
         ELSE
            ICOL = 74
         ENDIF
  110 CONTINUE
      IF (CKCHUP(LINE(1)(45:),1) .EQ. 'L') KPHSE(KNUM)=1
      IF (CKCHUP(LINE(1)(45:),1) .EQ. 'S') KPHSE(KNUM)=-1
C
      IF (SETTMP) THEN
C        line 2 has substitute fit temperatures
         CALL CKISUB (LINE(2), SUB, NSUB)
         IF (NSUB-1 .GT. MAXTP) THEN
            WRITE (LOUT, '(6X,A,I3,A,A)')
     1      'Error...MAXTP must be increased at least to ',NSUB-1,
     2      'for species ', KNAME(KNUM)
            KERR = .TRUE.
C           can only use MAXTP data for this species
            SETTMP = .FALSE.
         ELSE
            NT(KNUM) = NSUB - 1
            DO 116 N = 2, NSUB
               CALL CKPARR (SUB(N), -1, 1, T(N-1,KNUM), NVAL, IER, LOUT)
  116       CONTINUE
         ENDIF
         LINE(2) = ' '
      ENDIF
C
      IF (.NOT. SETTMP) THEN
         NT(KNUM) = NTFIT
         DO 115 N = 1, NTFIT
            T(N,KNUM) = TFIT(N)
  115    CONTINUE
      ENDIF
C
      IF (NT(KNUM).LE.3 .AND. LINE(1)(46:73).NE.' ') THEN
C        Line 1 fields OK for standard 2-range case, or 1 range
         IF (LINE(1)(46:55) .NE. ' ') CALL CKPARR
     1      (LINE(1)(46:55), 0, 1, T(1,KNUM), NVAL, IER, LOUT)
         IF (LINE(1)(66:73) .NE. ' ') CALL CKPARR
     1      (LINE(1)(66:73), 0, 1, T(2,KNUM), NVAL, IER, LOUT)
         IF (LINE(1)(56:65) .NE. ' ') CALL CKPARR
     1      (LINE(1)(56:65), 0, 1, T(3,KNUM), NVAL, IER, LOUT)
      ENDIF
C
C     Lines 2..NLINES:
C     15-column polynomial coefficients, HIGHEST to LOWEST temperatures
C
      NA = 0
      NRANGE = NT(KNUM) - 1
      DO 80 N = 2, NLINES
         ICOL = 1
         DO 75 L = 1, 5
            IF (LINE(N)(ICOL:ICOL+14) .NE. ' ') THEN
               NA = NA + 1
               IF (NA .EQ. 8) THEN
C                 back up to previous temperature range 7-polynomial set
                  NA = 1
                  NRANGE = NRANGE - 1
                  IF (NRANGE .EQ. 0) GO TO 90
               ENDIF
               READ (LINE(N)(ICOL:ICOL+14), '(E15.8)')
     1               A(NA, NRANGE, KNUM)
            ENDIF
            ICOL = ICOL + 15
   75    CONTINUE
   80 CONTINUE
      IF (NRANGE .GT. 1) THEN
            WRITE (LOUT, '(6X,A,I3,A)')
     1      'Error...data represents less than ',NT(KNUM)-1,
     2      ' temperature ranges...'
         KERR = .TRUE.
      ENDIF
C
   90 CONTINUE
C     set thermo flag for this species
      ITHRM(KNUM) = .TRUE.
C     reset defaults
      GO TO 25
C
C     end of SUBROUTINE CKTHRM thermo input
C
  400 CONTINUE
C     still have last species to process?
      IF (NLINES .GT. 0) GO TO 500
      RETURN
C
22222 CONTINUE
      WRITE (LOUT,*) ' Error reading thermodynamic data...'
      KERR = .TRUE.
C
C     end of SUBROUTINE CKTHRM
      RETURN
      END
C
      SUBROUTINE CKUNIT (LINE, AUNITS, EUNITS, IUNITS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUNIT
C  processes a character string containing information about
C  units used in a reaction description.
C
C  Arguments:
C  LINE
C  AUNITS
C  EUNITS
C  IUNITS
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
      CHARACTER*(*) LINE, IUNITS, AUNITS, EUNITS
      CHARACTER*4 CKCHUP
      INTEGER CKLSCH
      EXTERNAL CKLSCH, CKCHUP
C
      AUNITS = ' '
      EUNITS = ' '
      IUNITS = ' '
      ILAST = CKLSCH(LINE)-3
      DO 85 N = 1, ILAST
         IND = CKLSCH(IUNITS)
         IF (EUNITS .EQ. ' ') THEN
            IF (CKCHUP(LINE(N:), 4)     .EQ. 'CAL/') THEN
               EUNITS = 'CAL/'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units cal/mole'
               ELSE
                  IUNITS(IND:) = ', E units cal/mole'
               ENDIF
            ELSEIF (CKCHUP(LINE(N:), 4) .EQ. 'KCAL') THEN
               EUNITS = 'KCAL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kcal/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kcal/mole'
               ENDIF
            ELSEIF (CKCHUP(LINE(N:), 4) .EQ. 'JOUL') THEN
               EUNITS = 'JOUL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Joules/mole'
               ELSE
                  IUNITS(IND:) = ', E units Joules/mole'
               ENDIF
            ELSEIF (CKCHUP(LINE(N:), 4) .EQ. 'KJOU') THEN
               EUNITS = 'KJOU'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kjoule/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kjoule/mole'
               ENDIF
            ELSEIF (CKCHUP(LINE(N:), 4) .EQ. 'KELV') THEN
               EUNITS = 'KELV'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kelvins'
               ELSE
                  IUNITS(IND:) = ', E units Kelvins'
               ENDIF
            ELSEIF (CKCHUP(LINE(N:), 4) .EQ. 'EVOL') THEN
               EUNITS = 'EVOL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units electron Volts'
               ELSE
                  IUNITS(IND:) = ', E units electron Volts'
               ENDIF
            ENDIF
         ENDIF
         IF (AUNITS .EQ. ' ') THEN
            IF (CKCHUP(LINE(N:), 4) .EQ. 'MOLE') THEN
               IF (N+4.LE.CKLSCH(LINE) .AND.
     1                    CKCHUP(LINE(N+4:),1).EQ.'C') THEN
C
                  AUNITS = 'MOLC'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units molecules'
                  ELSE
                      IUNITS(IND:) = ', A units molecules'
                  ENDIF
               ELSE
                  AUNITS = 'MOLE'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units mole-cm-sec-K'
                  ELSE
                     IUNITS(IND:) = ', A units mole-cm-sec-K'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   85 CONTINUE
C
      IF (AUNITS .EQ. ' ') THEN
         AUNITS = 'MOLE'
         IND = CKLSCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', A units mole-cm-sec-K'
         ELSE
            IUNITS(IND:) = ' A units mole-cm-sec-K'
         ENDIF
      ENDIF
C
      IF (EUNITS .EQ. ' ') THEN
         EUNITS = 'CAL/'
         IND = CKLSCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', E units cal/mole'
         ELSE
            IUNITS(IND:) = ' E units cal/mole'
         ENDIF
      ENDIF
C
C     end of SUBROUTINE CKUNIT
      RETURN
      END
C                                                                      C
      SUBROUTINE CPREAC (MDIM, KDIM, IDIM, MAXSP, MAXORD,
     1                   NSPEC, NPAR, PAR, RPAR, AUNITS, EUNITS,
     1                   NREAC, NUNK, NU, KCHRG, KNCF, IDUP, IFAL,
     2                   KFAL, NFAR, PFAL, IFOP, IFLO, IREV, ITHB, ILAN,
     3                   IRLT, KERR, LOUT, IRNU, RNU, IXSM, IXSI, IXSK,
     5                   IMOM, KMOM, IORD, KORD, RORD, CKMIN,
     5                   AXUNIT, EXUNIT)
C
C  START PROLOGUE
C
C  SUBROUTINE CPREAC
C  does final checking for completeness and correctness of
C  a reaction once all of its input descriptions are processed;
C  formatted information and diagnostics is then printed.
C
C  Arguments:
C  MDIM      - Integer scalar, maximum number of elements.
C  KDIM      - Integer scalar, maximum number of species.
C  IDIM      - Integer scalar, maximum number of reactions.
C  MAXSP     - Integer scalar, maximum number of reaction species.
C  MAXORD    - Integer scalar, reaction maximum change-orders.
C  NSPEC(*)  - Integer array, reaction count of reactants+products.
C  NPAR      - Integer scalar, required number of Arrhenius params.
C  PAR(*,*)  - Real matrix, reactions Arrhenius parameters.
C  RPAR(*,*) - Real matrix, explicit reverse Arrhenius parameters.
C  AUNITS    - Character string, A units expression for mechanism.
C  EUNITS    - Character string, E units expression for mechanism.
C  NREAC(*)  - Integer array, reaction count of reactants only.
C  NUNK(*,*) - Integer matrix, reactions species indices.
C  NU(*,*)   - Integer matrix, reactions stoichiometric coefficients.
C  KCHRG(*)  - Integer array, species electronic charges.
C  KNCF(*,*) - Integer matrix, elemental composition of species.
C  IDUP(*)   - Integer array, reaction duplication flags.
C  IFAL(*)   - Integer array, pres dependency reaction indices.
C  KFAL(*)   - Integer array, 3rd-bodies for pres dependency reactions.
C  NFAR      - Integer scalar, maximum number of pres dependency
C              parameters.
C  PFAL(*,*) - Real matrix, pres dependency parameters.
C  IFOP(*)   - Integer array, pres dependency types.
C  IFLO(*)   - Integer array, pres dependency options.
C  IREV(*)   - Integer array, explicit reverse parameter reaction
C              indices.
C  ITHB(*)   - Integer array, third-body reaction indices.
C  ILAN(*)   - Integer array, Landau-Teller reaction indices.
C  IRLT(*)   - Integer array,  reactions with explicit reverse
C              Landau-Teller parameters.
C  KERR      - Logical, error flag.
C  LOUT      - Integer scalar, formatted output file unit number.
C  IRNU(*)   - Integer array, real stoichiometry reaction indices.
C  RNU(*,*)  - Real matrix, real stoichiometric coefficients.
C  IXSM(*)   - Integer array, electron momentum transfer reactions.
C  IMOM(*)   - Integer array, ion momentum transfer reactions.
C  CKMIN     - Real scalar, minimum used for content/balance checks.
C  AXUNIT    - Character string, A units expression for reaction.
C  EXUNIT    - Character string, E units expression for reaction.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C     (Value of Avrogadro's Constant from 1986 CODATA
C      recommended values (1993 CRC)
C      J. Research National Bureal of Standards, 92, 95, 1987
C      6.0221367(39) mol-1 )
      PARAMETER (RU_JOUL = 8.314510D0, AVAG = 6.0221367D23, ONE=1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C      PARAMETER (RU_JOUL = 8.314510E0, AVAG = 6.0221367E23, ONE=1.0E0)
C*****END precision > single
C
      COMMON /CKINT/ LENICK, LENRCK, LENCCK, MM, KK, II, NREV, NFAL,
     1               NTHB, NLAN, NRLT, NWL, NCHRG, NEIM, NJAN, NFT1,
     2               NEXC, NMOM, NXSM, NTDE, NRNU, NORD, MELECT,
     3               KELECT, NKION
C
C     Integer arrays
      DIMENSION NSPEC(IDIM), NREAC(IDIM), NUNK(MAXSP,IDIM),
     1          NU(MAXSP,IDIM), KCHRG(KDIM), KNCF(MDIM,KDIM),
     2          IDUP(IDIM), IFAL(IDIM), KFAL(IDIM), IFOP(IDIM),
     3          IFLO(IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     4          IREV(IDIM), ITHB(IDIM), ILAN(IDIM), IRLT(IDIM),
     5          IRNU(IDIM), IXSM(IDIM), IXSI(IDIM), IXSK(IDIM),
     6          IMOM(IDIM), KMOM(IDIM)
C     Real arrays
      DIMENSION PAR(NPAR,IDIM), RPAR(NPAR,IDIM), PFAL(NFAR,IDIM),
     1          RORD(MAXORD,IDIM), RNU(MAXSP,IDIM)
C
      CHARACTER*(*) AUNITS, EUNITS, AXUNIT, EXUNIT
      CHARACTER*16  UNITSTR
      LOGICAL KERR, LREV, LLAN, LRLT, LTHB, LFAL, LMOM, LXSM
C
      LREV = II .EQ. IREV(MAX(NREV,1))
      LLAN = II .EQ. ILAN(MAX(NLAN,1))
      LRLT = II .EQ. IRLT(MAX(NRLT,1))
      LTHB = II .EQ. ITHB(MAX(NTHB,1))
      LFAL = II .EQ. IFAL(MAX(NFAL,1))
      LMOM = II .EQ. IMOM(MAX(NMOM,1))
      LXSM = II .EQ. IXSM(MAX(NXSM,1))
C
C     balance reaction for element/mass and electronic charge
      CALL CKBAL (LOUT, MDIM, KDIM, IDIM, MAXSP, NU(1,II), NUNK(1,II),
     1            KCHRG, IRNU, RNU, KNCF, CKMIN, KERR)
C
C     check for duplicated reactions
      CALL CKDUP (LOUT, IDIM, II, MAXSP, NSPEC, NREAC, NU, NUNK, NFAL,
     1            IFAL, KFAL, IDUP, NTHB, ITHB, KERR)
C
      IF (LFAL) THEN
         IF (IFOP(NFAL).LE.0) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...no pres-dependency option defined...'
         ELSE
            IF (IFLO(NFAL).LT.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, '(6X,A)')
     1      'Error...no LOW/HIGH pres-dependency parameters given...'
            ENDIF
         ENDIF
      ENDIF
C
      IF (LRLT) THEN
         IF (.NOT. LLAN) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)') 'Error...',
     1      'no forward Landau-Teller parameters given for RLT...'
         ENDIF
         IF (.NOT. LREV) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A)')
     1      'Error...no reverse Arrhenius parameters given for RLT...'
         ENDIF
      ELSE
         IF (LLAN .AND. LREV) THEN
            KERR = .TRUE.
            WRITE (LOUT, '(6X,A,A)') 'Error...',
     1      'no reverse Landau-Teller parameters given for REV...'
         ENDIF
      ENDIF
C
      IF (LMOM .AND. KELECT.GT.0) THEN
C        need exactly one electron reactant (KELECT) with coeff approx 1
C        and exactly one non-electron (KMOM) with coeff approx 1
C        after all change-orders
         CALL CPMOME (LOUT, IDIM, II, KK, KELECT, CKMIN, NREAC(II),
     1                NUNK(1,II), NU(1,II), MAXSP, NRNU, IRNU, RNU,
     2                NORD, IORD, MAXORD, KORD, RORD, KMOM(NMOM), KERR)
      ENDIF
C
      IF (LXSM .AND. NKION.GT.0) THEN
C        need exactly one ion reactant (IXSI) with coeff approx 1
C        and exactly one non-ion reactant (IXSK) with coeff approx 1
         CALL CPXSMI (LOUT, KDIM, IDIM, II, KK, KCHRG, CKMIN,
     1                NREAC(II), NUNK(1,II), NU(1,II), MAXSP, NRNU,
     2                IRNU, RNU, NORD, IORD, MAXORD, KORD, RORD,
     3                IXSI(NXSM), IXSK(NXSM), KERR)
      ENDIF
C
C     find appropriate conversion factor between the energy unit
C     specified by the user and Kelvin units
C
      UNITSTR = ' '
      IF (EXUNIT .NE. ' ') THEN
C        user has changedunits for this reaction by auxiliary keyword
         UNITSTR = EXUNIT
      ELSE
C        default units, or user has changed units for all reactions
         UNITSTR = EUNITS
      ENDIF
C
      EFAC = 1.0
C
      IF (UNITSTR .NE. 'KELV') THEN
         IF (EUNITS .EQ. 'EVOL') THEN
C           convert B and E from eV to Kelvin
            EFAC = 11595.
         ELSEIF (EUNITS .EQ. 'CAL/') THEN
C           convert E from cal/mole to Kelvin
            EFAC = 4.184  / RU_JOUL
         ELSEIF (EUNITS .EQ. 'KCAL') THEN
C           convert E from kcal/mole to Kelvin
            EFAC = 4184.0 / RU_JOUL
         ELSEIF (EUNITS .EQ. 'JOUL') THEN
C           convert E from Joules/mole to Kelvin
            EFAC = 1.00  / RU_JOUL
         ELSEIF (EUNITS .EQ. 'KJOU') THEN
C           convert E from Kjoules/mole to Kelvin
            EFAC = 1000.0 / RU_JOUL
         ENDIF
         PAR(3,II) = PAR(3,II) * EFAC
C
         IF (LREV) RPAR(3,NREV) = RPAR(3,NREV) * EFAC
         IF (LFAL) PFAL(3,NFAL) = PFAL(3,NFAL) * EFAC
      ENDIF
C
      UNITSTR = ' '
      IF (AXUNIT .NE. ' ') THEN
C        pre-exponential units for this reaction
         UNITSTR= AXUNIT
      ELSE
C        default units
         UNITSTR = AUNITS
      ENDIF
C
      IF (UNITSTR .NE. 'MOLC') RETURN
C
      NSTOR = 0
      NSTOP = 0
      DO 50 N = 1, MAXSP
         IF (NU(N,II) .LT. 0) THEN
C           sum of stoichiometric coefficients of reactants
            NSTOR = NSTOR + ABS(NU(N,II))
         ELSEIF (NU(N,II) .GT. 0) THEN
C           sum of stoichiometric coefficients of products
            NSTOP = NSTOP + NU(N,II)
         ENDIF
   50 CONTINUE
C
      IF (LFAL) THEN
C
C        pres dependency reaction, "(+M)" or "(+species name)" does not
C        count except in "LOW" or "HIGH" A-factor;
C        reverse-rate declarations are not allowed
C
         IF (NSTOR.GT.1) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
         NSTOR = NSTOR + 1
         IF (NSTOR.GT.1) PFAL(1,NFAL) = PFAL(1,NFAL)*AVAG**(NSTOR-1)
C
      ELSEIF (LTHB) THEN
C
C        third body reaction, "+M" counts as species in
C        forward and reverse A-factor conversion
C
         NSTOR = NSTOR + 1
         NSTOP = NSTOP + 1
         IF (NSTOR.GT.1) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
         IF (LREV .AND. NSTOP.GT.1)
     1          RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
C
      ELSEIF (LXSM) THEN
C
C        Don't change momentum-transfer x-sections -units are [cm2]
C
      ELSEIF (LMOM) THEN
C
C        Don't change e- momentum-transfer collision freq [cm3/s]
C
      ELSE
C
C        not third-body or pres dependency reaction, but may have
C        reverse rates.
C
         IF (NSTOR .GT. 1) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
         IF (LREV .AND. NSTOP.GT.1)
     1          RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
      ENDIF
C
C     end of SUBROUTINE CPREAC
      RETURN
      END
C
      SUBROUTINE CPMOME (LOUT, IDIM, II, KK, KELECT, CKMIN, NREAC,
     1                   NUNK, NU, MAXSP, NRNU, IRNU, RNU, NORD, IORD,
     2                   MAXORD, KORD, RORD, KMOM, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CPMOME
C  does the final checking of the electron momentum-transfer
C  reaction after all auxiliary input is processed.
C
C  Arguments:
C  LOUT
C  II
C  KK
C  KELECT
C  CKMIN
C  NREAC
C  NUNK(*)
C  NU(*)
C  MAXSP
C  NRNU
C  IRNU(*)
C  RNU(*,*)
C  NORD
C  IORD(*)
C  MAXORD
C  KORD(*,*)
C  RORD(*,*)
C  KMOM
C  KERR
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
      DIMENSION NUNK(MAXSP), NU(MAXSP), IRNU(IDIM), RNU(MAXSP,IDIM),
     1          IORD(IDIM), KORD(MAXORD,IDIM), RORD(MAXORD,IDIM)
      LOGICAL LRNU, LORD, KERR
C
C     need exactly one electron reactant (KELECT) with coeff approx 1
      NELECT = 0
      ECOEF = 0.0
C     and exactly one non-electron reactant (KMOM) with coeff approx 1
      NGAS = 0
      GCOEF = 0.0
C
      LORD = IORD(MAX(NORD,1)) .EQ. II
      LRNU = IRNU(MAX(NRNU,1)) .EQ. II
C
      DO 300 N = 1, NREAC
         IF (NUNK(N) .EQ. KELECT) THEN
C           the electron species is present
            NELECT = NELECT + 1
            ECOEF = ABS(NU(N))
            IF (LRNU) ECOEF = ABS(RNU(N,NRNU))
            IF (LORD) THEN
C              look for forward change-order
               CALL CKINUM(-KELECT, KORD(1,NORD), MAXORD, IND)
               IF (IND .GT. 0) ECOEF = RORD(IND,NORD)
            ENDIF
C           not enough of the electron species to consider
            IF (ECOEF .LT. 1.0-CKMIN) NELECT = NELECT-1
C
         ELSE
C           a non-electron species is present
            NGAS = NGAS + 1
            KMOM = NUNK(N)
            GCOEF = ABS(NU(N))
            IF (LRNU) GCOEF = ABS(RNU(N,NRNU))
C           is there a later coefficient for the species?
            IF (LORD) THEN
               CALL CKINUM(-KMOM, KORD(1,NORD), MAXORD, IND)
               IF (IND .GT. 0) GCOEF = RORD(IND,NORD)
            ENDIF
C           not enough of the non-electron species to consider
            IF (GCOEF .LT. 1.0-CKMIN) NGAS = NGAS - 1
         ENDIF
  300 CONTINUE
C
C     any other non-electron species involved?
C     (don't need to check electron species, as the mechanism
C      is supposed to have only one electron species)
      IF (LORD) THEN
         DO 260 K = 1, KK
            IF (K .EQ. KELECT .OR. K .EQ. KMOM) GO TO 260
C           K is another species
            CALL CKINUM (-K, KORD(1,NORD), MAXORD, IND)
C           IND is its location in the change-orders
            IF (IND.EQ.0 .OR. RORD(IND,NORD).LT.CKMIN) GO TO 260
C           K is present as forward change-order > CKMIN
            IF (KMOM.LE.0 .OR. GCOEF.LT.CKMIN) THEN
C              did not previously have a good gas species
               KMOM = ABS(KORD(IND,NORD))
               GCOEF = RORD(IND,NORD)
            ELSEIF (KMOM.GT.0 .AND. GCOEF.GT.CKMIN) THEN
C              too many non-zero gas species with order > CKMIN
               NGAS = NGAS + 1
            ENDIF
  260    CONTINUE
      ENDIF
C
      IF (NGAS.NE.1 .OR. KMOM.LE.0 .OR.
     1    (GCOEF.GT.1.0+CKMIN .OR. GCOEF.LT.1.0-CKMIN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A,A)') 'Error...',
     1   'MOME requires exactly one non-electron ',
     2   'with coefficient/order approx. 1...'
      ENDIF
      IF (NELECT.NE.1 .OR.
     1    (ECOEF.GT.1.0+CKMIN .OR. ECOEF.LT.1.0-CKMIN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A)')
     1   'Error...MOME requires exactly one electron participant ',
     2   'with stoichiometric coefficient/order approx. 1...'
      ENDIF
C
C     end of SUBROUTINE CPMOME
      RETURN
      END
C
      SUBROUTINE CPXSMI (LOUT, KDIM, IDIM, II, KK, KCHRG, CKMIN, NREAC,
     1                   NUNK, NU, MAXSP, NRNU, IRNU, RNU, NORD, IORD,
     2                   MAXORD, KORD, RORD, IXSI, IXSK, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CPXSMI
C  checks the ion momentum-transfer, cross-section reaction
C  after all of the reaction auxiliary data is processed.
C
C  Arguments:
C  LOUT
C  II
C  KK
C  KCHRG(*)
C  CKMIN
C  NREAC
C  NUNK(*)
C  NU(*)
C  MAXSP
C  NRNU
C  IRNU(*)
C  RNU(*,*)
C  NORD
C  IORD(*)
C  MAXORD
C  KORD(*,*)
C  RORD(*,*)
C  IXSI(*)
C  IXSK(*)
C  KERR
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
      DIMENSION KCHRG(KDIM), NUNK(MAXSP), NU(MAXSP), IRNU(IDIM),
     1          RNU(MAXSP,IDIM), IORD(IDIM), KORD(MAXORD,IDIM),
     2          RORD(MAXORD,IDIM)
      LOGICAL LRNU, LORD, KERR
C
C     need exactly one ion reactant (IXSI) with coeff approx 1
      NION = 0
      CION = 0.0
C     and exactly one non-ion reactant (IXSK) with coeff approx 1
      NGAS = 0
      GCOEF = 0.0
C
      LORD = IORD(MAX(NORD,1)) .EQ. II
      LRNU = IRNU(MAX(NRNU,1)) .EQ. II
C
      DO 350 N = 1, NREAC
         NK = NUNK(N)
         IF (KCHRG(NK) .NE. 0) THEN
C           this is an ion
            NION = NION + 1
            IXSI = NK
            CION = ABS(NU(N))
            IF (LRNU) CION = ABS(RNU(N,NRNU))
            IF (LORD) THEN
C              look for forward change-order
               CALL CKINUM(-NK, KORD(1,NORD), MAXORD, IND)
               IF (IND .GT. 0) CION = RORD(IND,NORD)
            ENDIF
C           not enough of the electron species to consider
            IF (CION .LT. 1.0-CKMIN) NION = NION - 1
C
         ELSE
C           a non-ion species is present
            NGAS = NGAS + 1
            IXSK = NK
            GCOEF = ABS(NU(N))
            IF (LRNU) GCOEF = ABS(RNU(N,NRNU))
C           is there a later coefficient for the species?
            IF (LORD) THEN
               CALL CKINUM(-NK, KORD(1,NORD), MAXORD, IND)
               IF (IND .GT. 0) GCOEF = RORD(IND,NORD)
            ENDIF
C           not enough of the non-ion species to consider
            IF (GCOEF .LT. 1.0-CKMIN) NGAS = NGAS - 1
         ENDIF
  350 CONTINUE
C
C     any other ion species involved?
      IF (LORD) THEN
         DO 360 K = 1, KK
            IF (K .EQ. IXSI .OR. K .EQ. IXSK) GO TO 360
C           K is another species
            CALL CKINUM(-K, KORD(1,NORD), MAXORD, IND)
            IF (IND.EQ.0 .OR. RORD(IND,NORD).LT.CKMIN) GO TO 360
C           K is present as a forward change-order > CKMIN
            IF (KCHRG(K) .NE. 0) THEN
C              ion species
               IF (IXSI.LE.0 .OR. CION.LT.CKMIN) THEN
C                 did not previously have a good ion species
                  IXSI = ABS(KORD(IND,NORD))
                  CION = RORD(IND,NORD)
               ELSEIF (IXSI.GT.0 .AND. CION.GT.CKMIN) THEN
C                 too many ion species with order > CKMIN
                  NGAS = NGAS + 1
               ENDIF
            ELSE
C              non-ion species
               IF (IXSK.LE.0 .OR. GCOEF.LT.CKMIN) THEN
C                 did not previously have a good non-ion species
                  IXSK = ABS(KORD(IND,NORD))
                  GCOEF = RORD(IND,NORD)
               ELSEIF (IXSK.GT.0 .AND. GCOEF.GT.CKMIN) THEN
C                 too many non-ion species with order > CKMIN
                  NGAS = NGAS + 1
               ENDIF
            ENDIF
  360    CONTINUE
      ENDIF
C
      IF (NGAS.NE.1 .OR.IXSK.LE.0 .OR.
     1    (GCOEF.GT.1.0+CKMIN .OR. GCOEF.LT.1.0-CKMIN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...XSMI requires exactly one non-ion participant ',
     2   'with stoichiometric coefficient/order approx. 1...'
      ENDIF
      IF (NION.NE.1 .OR.
     1    (CION.GT.1.0+CKMIN .OR. CION.LT.1.0-CKMIN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, '(6X,A,A)')
     1   'Error...XSMI requires exactly one ion participant ',
     2   'with stoichiometric coefficient/order approx. 1...'
      ENDIF
C
C     end of SUBROUTINE CPXSMI
      RETURN
      END
C
      LOGICAL FUNCTION CKFILE (LUNIT, LFORM)
      INTEGER LUNIT
      CHARACTER LFORM*(*), LTYPE*16, STR1*1
      LOGICAL LOPEN
C
C     INQUIRE returns logical OPENED, character LTYPE
      INQUIRE (LUNIT, OPENED=LOPEN, FORM=LTYPE)
      CKFILE = LOPEN .AND. (LFORM(1:1).EQ.LTYPE(1:1)) 
      IF (CKFILE) THEN
         IF (LTYPE(1:1) .EQ. 'U') THEN
            READ (LUNIT, ERR=500)
         ELSE
            READ (LUNIT, '(A)', ERR=500) STR1
         ENDIF
         BACKSPACE (LUNIT)
         RETURN
  500    CONTINUE
C        error reading file
         CKFILE = .FALSE.
      ENDIF
      RETURN
      END
