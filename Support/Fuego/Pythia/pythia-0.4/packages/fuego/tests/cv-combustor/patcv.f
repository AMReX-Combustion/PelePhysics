         PROGRAM CV
C
C 'CV' Constant Volume explosion STRUCTURE CALCULATION 
C
C J. E. Shepherd
C Graduate Aeronautical Laboratories, MS 105-50
C California Institute of Technology
C Pasadena, CA 91125  USA
C jeshep@galcit.caltech.edu
C
C --- recent history ------ 
c  1990  Chemkin II version for PC                  JES 02-05-90
C  1995 fixed density input bug, write t_explsn to
C        screen                                      "  04-06-95
C       changed dimensions, fixed NOSHK as default   "  04-22-95
C       Fixed file handling, added UDET output       "  04-25-95
C       Fixed bug in multiple run cases, P output    "  04-25-95 V1.35D
c  UNIX mods to date and time, tested with cklib V4.9
c  and ckinterp 3.8                                  JES10-09-96 V1.4D
c  1996 Increased NK to 100 (was 50), LENIWK to 6000
c       (was 3000), and LENWK to 6000 (was 2500)
c       for Unix                                    MJK 10-23-96
c  1997 Increased 10000 (was 6000), and LENWK to 10000 (was 6000)
c                                                   MJK 6-3-97
c
c  1998 Increased LENIWK to 20000, LENWK to 20000, lencwk to 1000 (was 500)
c	NK to 150
c       Line 377, first test condition changed from .gt. to .lt. to fix
c	90% peak time bug (was finding 90% peak *after* peak)
c       Eric Schultz 7/6/98
c  2000 V1.5D Changed NK=500, changed ifix to int in readcv,
c       changed common block to eliminate gcc warnings
c       changed dimensions and checked gcc 2.95.2 compilation
c       LENIWK = 40000, LENWK = 40000, lencwk = 2000   JES 6-3-00
c
c  2000 Now program can also take internal energy as input (mass unit)
c       to faciliate fixing E,rho for CV calculations.
c       (use readcv2.f) 
c       The output file is mildly reformatted for later parsing PH 8-31-00
c
C
C  References:
C  Shepherd - Notes on PC-Chemkin and ZND, undated.
C  Kee, Rupley, Miller "CHEMKIN-II", SAND89-8009, 1989.
C
C  Must be linked to DDEBDF, ZEROIN, CKLIB, READCV
C  Requires CHEMKIN linking file 'Chem.bin' to run.
C
C        THE INDEPENDENT VARIABLE IN THIS PROBLEM IS TIME.
C        THE DEPENDENT VARIABLES ARE THE PRESSURE,P, DENSITY,RHO,
C        AND THE MASS FRACTIONS, Y(K).  ALL DIMENSIONS ARE IN 
C        ORIGINAL CHEMKIN PSUEDO-CGS SYSTEM.
C
C        X=TIME
C        Z(1)=T
C        Z(K+1)=Y(K)
C
C I/O UNITS USED:
C lin     5       INPUT UNIT 
C LOUT    7       OUTPUT UNIT
C LPLOT   8       PLOT FILE (ASCII COLUMNAR DATA)
C LINK    25      CHEMKIN LINKING FILE
C
C DIMENSIONS FOR CHEMKIN AND INTEGRATOR ARE CURRENTLY SET UP FOR ABOUT
C 50 SPECIES.  THE DATA IN THE PARAMETER STATEMENT
C BELOW DETERMINES THE DIMENSIONS.  THE PARAMETERS ARE:
C       NK -    NUMBER OF SPECIES
C       NEL -    NUMBER OF ELEMENTS
C NOTE: CHANGE PARAMETER STATEMENTS IN SUBROUTINES ALSO!
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NEL = 6, NK = 500, NQ = 1+NK,
     &           SMALL = 1.D-100, ITMAX=5000)
      PARAMETER (LRW = 250 + NQ*(10 + NQ), LIW = 55 + NQ)
      PARAMETER (lin = 5, lout = 7, lplot = 8, LINK = 25, NCHAR = 76)
      PARAMETER (LENIWK = 40000, LENWK = 40000, lencwk = 2000)
      character cckwrk(lencwk)*16, ksym(nk)*16, case*80
      DIMENSION INFO(15),ELWRK(LRW),IELWRK(LIW)
      DIMENSION Z(NQ),ZP(NQ),Y(NK), timmax(itmax), dtmax(itmax)
      LOGICAL PLOT, PLOTE,KERR, IERR,DINC,noshock
      character filold*40,po*40,old*3
      character filnew*40,pn*40,unknown*7
      CHARACTER DATE*10, CLOCK*8
      COMMON/WRK/ ICKWRK(LENIWK),CKWRK(LENWK)
      COMMON/PARAM/RU,WT(NK),P1,RHO1,PSC,UVEL,SON,XM,UDET,X(NK),T,
     &             P,RHO,RPAR(4),ALPHA,Y,WTM,T1,KK
      EXTERNAL FUN,JAC,SJUMP
      DATA PLOTE/.FALSE./
      data iplot /6/, RTOL/1.D-05/, ATOL/1.D-07/,DINC/.FALSE./
      DATA DDIST/1D-04/,TMAX/1.D0/,noshock/.false./,udet/0.d0/
      data filnew/'cv.out'/,pn/'Output file name'/,unknown/'unknown'/
      data filold/'chem.bin'/,po/'Linking file name'/,old/'old'/
      alpha = 0.d0
c
c open output file, allows using new name or overwriting old cv.out
c
      call ipofil(LOUT,0,pn,unknown,filnew)
C
C Open linking file, this is a special open for unformatted files
C Linking file can have any name now and be in any directory.
c 'cklink.bin' is still default
c
      call ipofilu(link,0,po,old,filold)
C
C        INITIALIZE CHEMKIN
C
      CALL FGINIT (LENIWK, LENWK, lencwk, LINK, LOUT, ICKWRK, CKWRK,
     &              cckwrk)
      CALL FGINDX (ICKWRK, CKWRK, MM, KK, II, NFIT)
      if (kk.gt.nk) then
          write(lout,
     & '('' Species dimension too small..must be at least '',I4)') KK
          stop
      endif
      CALL FGSYMS (cckwrk, lout,ksym,ierr)
      if (ierr) kerr = .true.
      CALL FGwt(ICKWRK,CKWRK,WT)
      CALL FGRP (ICKWRK, CKWRK, RU, RUC, PATM)
C
C PRINT OUT HEADER AND SPECIES NAMES FOR THIS REACTION SET
C
C      call tim(CLOCK)
C      call dat(date)
      WRITE(LOUT ,2000) DATE,CLOCK
2000  format(//' CV: EXPLOSION STRUCTURE CALCULATION '/,
     &' VERSION 1.5D - LAST MOD ON 3 June 00'/,
     &' CALCULATION RUN ON ',A10,' AT ',A8//,
     &' SPECIES USED IN THIS PROBLEM')
      WRITE(LOUT, 2100) (KSYM(K)(:10),K=1,KK)
2100  FORMAT(7(1X,1A10))
C
      OPEN(UNIT=lplot,FILE='CV.plt')
      WRITE(lplot,'('' **** BEGINNING OF PLOT FILE*****'')')
      WRITE(LPLOT, 2000) DATE,CLOCK
      WRITE(LPLOT, 2101) 'sp= ', (KSYM(K)(:7),K=1,KK)
2101  FORMAT(1X,1A04,15(1X,1A10))
c      WRITE(LPLOT, 2100) (KSYM(K)(:10),K=1,KK)
C
C RETURN POINT FOR BEGINNING ANOTHER CASE
C
      ICASE=0
 10   CONTINUE
C
C INITIALIZE ALL MASS AND MOLE FRACTIONS, Y(K)AND X(K), TO ZERO
C
      DO 20 K=1,KK
C         Y(K)= small
C         X(K)= small
         Y(K)= 0d0
         X(K)= 0d0
 20   CONTINUE
C
C****************************************************
C READ IN INPUT PARAMETERS FOR NEXT CASE
C        T1 = INITIAL  TEMPERATURE (K)
C        P1 = INITIAL  PRESSURE (ATM)
c        rho1 = initial density
c        e1   = internal energy (mass unit)
c need two of three
C*****************************************************
C
      call readcv2(case,rho1,t1,p1,e1,
     &plot,iplot,x,kk,ksym,rtol,atol,DINC,DDIST,TMAX,udet,noshock)
      IF (PLOT) PLOTE = .TRUE.
      if (case.eq.'STOP') goto 5000
      write (6,'(''  Running...'')')
C
C INCREMENT CASE NUMBER
C
      ICASE=ICASE+1
      CALL FGXTY(X,ICKWRK,CKWRK,Y)
      CALL FGMMWX(x,ICKWRK,CKWRK,WTM)
c
c check to see if density, pressure or temperature has to be computed
c 
      if (T1.ne.0d0 .and. rho1.ne.0d0) then
         P1 = RHO1*T1*RU/wtm
      else if (T1.ne.0d0 .and. P1.ne.0d0) then
         P1 = P1*PATM
         RHO1 = P1*WTM/(T1*RU)
      else if (rho1.ne.0d0 .and. p1.ne.0d0) then
         P1 = P1*PATM
         T1 = P1*WTM/(RU*RHO1)
      else if (rho1.ne.0d0 .and. e1.ne.0d0) then
         CALL GETT(e1,Y,ICKWRK,CKWRK,T1)
         P1 = RHO1*T1*RU/wtm
      else if (P1.ne.0d0 .and. e1.ne.0d0) then
         CALL GETT(e1,Y,ICKWRK,CKWRK,T1)
         RHO1 = P1*WTM/(T1*RU)
      else
         write(lout,*) 'Wrong inputs (from patcv.f)!'
      endif
c
c Is shock velocity non-zero and NOSHK not set?
c
      if (.not.NOSHOCK.and.UDET.eq.0.d0) NOSHOCK=.TRUE.
C
C ECHO OUT INPUT
C
      WRITE(LOUT, 2200) icase,case,T1,P1/PATM,RHO1
      do jj = 1, kk
      if (x(jj).gt.small) then
        WRITE(LOUT,'(1X,A10,2X,1PG12.4)')KSYM(JJ)(:10),X(JJ)
      endif
      enddo
2200  format(//' CASE #  ',i3,5x,1A80//,' INITIAL CONDITIONS '/,
     &' TEMPERATURE (K)  ',T30,1pg12.4/,
     &' PRESSURE (ATM)  ',T30,G12.4/,
     &' DENSITY (G/CC)  ',T30,G12.4/,
     &' SPECIES MOLE FRACTIONS: ')
      xmeancv=0d0
      xmeanu =0d0
      CALL FGCVBS(T1,Y,ICKWRK,CKWRK,xmeancv)
      CALL FGUBMS(T1,Y,ICKWRK,CKWRK,xmeanu)
      CALL GETT(xmeanu,Y,ICKWRK,CKWRK,TT)
      write(lout,*) 'Mean CV:    ', xmeancv
      write(lout,*) 'Mean U :    ', xmeanu
      write(lout,*) 'gett T :    ', TT
      IF (.not.NOSHOCK) THEN
C
C DETERMINE THE shock PRESSURE AND TEMPERATURE FROM THE JUMP
C CONDITIONS USING THE RAYLEIGH LINE-HUGONIOT INTERSECTION METHOD
C
C SET UP THE ERROR TOLERANCES FOR THE ROOT SOLVER
C
      AEE=1.D-07
      REE=1.D-12
      XLOW=1.D0/40.D0
      XHIGH=1.D0/1.005D0
C
C NOTE THAT 1/X IS THE DENSITY RATIO ACROSS THE SHOCK, 5.5 IS A TYPICAL VALUE
C for cj detonation, eg, a strong shock.
C
      CALL ZEROIN(SJUMP,XLOW,XHIGH,REE,AEE,IFLAG)
      IF (IFLAG.NE.1) THEN
         WRITE(LOUT,'(/'' SHCK CNVRGNC PRBLM, IFLG = '',I2)') IFLAG
         write(lout,'('' Is shock speed greater than sound speed?'')')
      ENDIF
C
C COMPUTE THE TEMP AND PRESSURE FROM THE DENSITY RATIO
C
      RHO=RHO1/XLOW
      P=P1+RHO1*UDET*UDET*(1.D0-XLOW)
      CALL FGMMWX(X,ICKWRK,CKWRK,WTM)
      T=P*WTM/(RU*RHO)
      ulab = udet*(1.0d0-rho1/rho)
      write(lout,2210)t,p/patm,rho,udet,ulab
2210  format(//,' POST SHOCK CONDITIONS '/,
     &' TEMPERATURE (K)  ',T30,1pg12.4/,
     &' PRESSURE (ATM)  ',T30,G12.4/,
     &' DENSITY (G/CC)  ',T30,G12.4/,
     &' SHOCK VELOCITY (CM/S)   ',T30,G12.4/,
     &' LAB PARTICLE VELOCITY (CM/S)   ',T30,G12.4/)
      else
c no shock in front
      p = p1
      t = t1
      rho = rho1
      endif
c header for output file
      WRITE(LOUT,2300)
 2300 FORMAT(/' REACTION ZONE STRUCTURE:'//,
     &//' ITER ',t10,' TIME',T25,' TEMP',T40,' DT/dt ',T50,' P'/t10,
     &' (S)',T25,' (K)',T40,' (K/S)',T50,' (atm) ')
C
C DESCRIPTORS FOR PLOT FILE
C
      IF (PLOT) THEN
        WRITE(LPLOT, 2200) icase,case,T1,P1/PATM,RHO1
        do jj = 1, kk
        if (x(jj).gt.small) then
          WRITE(LPLOT,'(1X,A10,2X,1PG12.4)')KSYM(JJ)(:10),X(JJ)
        endif
        enddo
        WRITE(LPLOT, '(/'' REACTION ZONE STRUCTURE''/)')
        WRITE(LPLOT,'('' THE OUTPUT DATA COLUMNS ARE:'')')
        WRITE(LPLOT ,2300)
        WRITE(LPLOT ,
     & '('' THEN THE SPECIES MOLE FRACTIONS IN THIS ORDER:'')')
        WRITE(LPLOT ,2100)(KSYM(K)(:10),K=1,KK)
      ENDIF
C
C*****************************************
C     SET THE INTEGRATION CONTROL PARAMETERS FOR LSODE
C     (NOTE: DEBDF IS JUST A DRIVER FOR A MODIFIED VERSION OF LSODE)
C     SEE THE SLATEC NOTES FOR DETAILS
C
      NEQ=KK+3
      DO 40 K1=1,15
         INFO(K1)=0
 40   CONTINUE
C
C CHOSE A LARGE OUTPUT TIME AND INTERMEDIATE OUTPUT MODE
C
      INFO(3) = 1
      TIME = 0.0D0
      TOUT = tmax
C
C INIT REACTION ZONE LENGTHS
C
      DTDTMAX = -1.0D+16
      tpeak = 0.d0
C
C INITIAL CONDITIONS FOR Z VECTOR
C
      Z(1)=T
      DO 50 K=1,KK
         Z(K+1)=Y(K)
 50   CONTINUE
      ITER = 0
C
C CALCULATE INITIAL POINT AND OUTPUT
C
      CALL FUN(TIME,Z,ZP,CKWRK,ICKWRK)
      WRITE(LOUT,2450) ITER,TIME,Z(1),ZP(1),P/PATM
 2450 FORMAT(1X,I5,t10,1PG10.3,T25,G10.3,T40,G10.3,T50,G10.3)
C
C
C     CALL THE DIFFERENTIAL EQUATION SOLVER.
C
 60   CALL DDEBDF(FUN,NEQ,TIME,Z,TOUT,INFO,RTOL,ATOL,IDID,
     *            ELWRK,LRW,IELWRK,LIW,CKWRK,ICKWRK,JAC)
      ITER = ITER + 1
      timmax(iter) = time
      dtmax(iter) = zp(1)
C      WRITE(6, 2450) iter,time,Z(1),Zp(1),P/PATM
C
C CHECK THE NUMBER OF INTERMEDIATE OUTPUT STEPS, IF GREATER THAN ITMAX
C GO TO THE NEXT PROBLEM
C
      IF (ITER.GE.ITMAX) THEN
         WRITE(LOUT ,'(10X,'' GREATER THAN '',
     &         I5,'' STEPS AND NO ANSWER'')') ITMAX
          WRITE(LOUT ,'(10X,'' IDID = '',I4)') IDID
          CALL FUN(TIME,Z,ZP,CKWRK,ICKWRK)
          GOTO 1000
      ENDIF
C
C CHECK THE RETURN FROM THE INTEGRATION
C
      IF (IDID.EQ.-33) THEN
        WRITE(LOUT,'('' FATAL ERROR FROM DEBDF, IDID ='',I4)') IDID
        GOTO 1000
      ENDIF
C
C SUCCESSFUL STEP TOWARD FINISH
C CALCULATE THE TEMPERATURE DERIVATIVE 
C
      IF (IDID.GE.-2) THEN
        CALL FUN(TIME,Z,ZP,CKWRK,ICKWRK)
        DTDT = ZP(1)
        IF (zp(1).GT.DTDTMAX) THEN
           DTDTMAX = zp(1)
           tPEAK = time
        ENDIF
C
C CODE IS HAVING TROUBLE, KEEP ON GOING UNTIL IT CRASHES FATALLY!
C
      ELSE
        INFO(1) = 0
        GOTO 60 
      ENDIF
C
C OUTPUT PLOT DATA AT THIS POINT, DO THIS EVERY IPLOTth ITERATION
C
      IF (PLOT.AND.0.EQ.MOD(ITER,IPLOT)) THEN
CPAT             WRITE(LPLOT,2500)ICASE,time,Z(1),Zp(1),(X(K),K=1,KK)
             WRITE(LPLOT,*)     ICASE, time, Z(1), Zp(1)
C Enable  111 for normal output  222 for phi output (mol/kg) (for ildm) 
C111             WRITE(LPLOT,*) (X(K),K=1,KK)
C222             WRITE(LPLOT,*) (Y(K)/WT(K)/1d-3,K=1,KK)
             WRITE(LPLOT,*) (Y(K),K=1,KK)

      ENDIF
CPAT 2500  FORMAT(1P,I2,5(13(E12.5,:,',')/))
CPAT 2600  FORMAT(6(E12.5,:,' ')/))
      
C
C OUTPUT DATA 
C     
      if (0.EQ.MOD(ITER,1)) THEN
          WRITE(lout, 2450) iter,time,Z(1),Zp(1),P/PATM
      endif
C
C KEEP ON COMPUTING
C SET INFO(1) IF IDID= -1 (500 STEPS) OR -2 (RELAXED ATOL, RTOL)
C
      IF (IDID.LE.1) THEN
         IF (IDID.LT.0) INFO(1) = 1
         GOTO 60 
      ENDIF
C
C IF IDID GT 1, TIME LIMIT REACHED!  (TIME=TMAX)
C
      WRITE (LOUT, '('' TIME = TMAX, LAST POINT:'')')
C
C FINISHED - PRINT OUT  LAST POINT AND GO TO THE NEXT CASE
C
1000  WRITE (LOUT, 2450) iter,time,Z(1),Zp(1),P/PATM
      WRITE(LOUT ,'('' NUMBER OF INTEGRATOR STEPS '',I4/)') ITER
      WRITE(LOUT, 2600) (X(K),K=1,KK)
 2600 FORMAT(' FINAL SPECIES MOLE FRACTIONS :'/,8(5(1X,G12.4)/))
C
C************ WRITE TO SUMMARY FILE ***********************
C
c     find 10 and 90 % points
      do i = 1,iter
      if (timmax(i).lt.tpeak.and.dtmax(i).lt.0.1*dtdtmax) then
          t_s = timmax(i)
      endif
      if (timmax(i).lt.tpeak.and.dtmax(i).gt.0.1*dtdtmax) then
          t_f = timmax(i)
      endif
      enddo
      WRITE(LOUT, 2700) tpeak, t_s, t_f
      write(*, *) ' Temperature at t_max: ' , Z(1)
      WRITE(6, 2700) tpeak, t_s, t_f
2700  FORMAT(/' Time to peak  (s) : ',T40,1PG12.4/,
     &       /' time to 10% of peak  (s) : ',T40,1PG12.4/
     &       /' time to 90% of peak (s) : ',T40,1PG12.4/)
C
C*************** NEXT CASE **********************************
C
      IF (PLOT) WRITE(LPLOT,'('' **** END OF PLOT FILE ****'')')
      GOTO 10
C
C     ERROR RETURN AND EOF
C
5000    WRITE(LOUT,'('' END OF DATA ON INPUT, ZND FINISHED'')')
        IF (PLOTE) THEN
           CLOSE(UNIT=LPLOT,STATUS='KEEP')
        ELSE
           CLOSE(UNIT=LPLOT,STATUS='DELETE')
        ENDIF
        close (unit=lout)
        close (unit=link)
        END
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      subroutine gett(E1,Y1,ICKWRK,CKWRK,T1)
      implicit real*8 (a-h,o-z)
      parameter (maxiter = 20)
c
      dx = 0.001
      tol= 1d-5
      T1 = 1000.
      do 50 i = 1,maxiter
        CALL FGUBMS(T1,Y1,ICKWRK,CKWRK,el)
        if (abs(e1-el) .lt. tol) then
          goto 60
        end if
        CALL FGUBMS(T1+dx,Y1,ICKWRK,CKWRK,er)
        dedt = (er-el)/dx
        T1   = T1+(e1-el)/dedt
 50   continue
 60   return 
      end 
C
      SUBROUTINE FUN(TIME,Z,ZP,CKWRK,ICKWRK)
C*******************************************************************
C THIS ROUTINE EVALUATES THE DERIVATIVES OF TEMPERATURE AND SPECIES
C ********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NEL = 6, NK = 500, NQ = 1+NK,
     &           SMALL = 1.D-100)
      PARAMETER (LRW = 250 + NQ*(10 + NQ), LIW = 55 + NQ)
      PARAMETER (LENIWK = 40000, LENWK = 40000, lencwk = 2000)
      COMMON/PARAM/RU,WT(NK),P1,RHO1,PSC,UVEL,SON,XM,UDET,X(NK),T,
     &             P,RHO,RPAR(4),ALPHA,Y,WTM,T1,KK
      DIMENSION Z(*),ZP(*),SIG(NK),Y(NK)
      DIMENSION UMS(NK),WDOT(NK)
C
C     UNRAVEL PHYSICAL VARIABLES FROM Z
C
      T=Z(1)
      DO 100 K=1,KK
         Y(K)=Z(K+1)
 100  CONTINUE
      CALL FGYTX(Y,ICKWRK,CKWRK,X)
      CALL FGMMWY(Y,ICKWRK,CKWRK,WTM)
      P = RHO*T*RU/WTM

      CALL FGCVBS(T,Y,ICKWRK,CKWRK,CVBMS)
C
C CALCULATE REACTION RATES, 
C SPECIES ENTHALPIES
C
      CALL FGWYP(P,T,Y,ICKWRK,CKWRK,WDOT)
      CALL FGUMS(T,ICKWRK,CKWRK,UMS)
C
C CALCULATE DERIVATIVES OF TEMPERATURE and MASS FRACTIONS 
C
      SUM =0.0D0
      DO 200 K=1,KK
         ZP(k+1)=WDOT(K)*WT(K)/RHO
         SUM=SUM+(UMS(K)/(CVBMS))*ZP(1+K)
 200  CONTINUE
      zp(1) = -1.d0*sum
      RETURN
      END
C
      SUBROUTINE JAC(X,U,PD,NROWPD,CKWRK,ICKWRK)
      DIMENSION CKWRK(*),ICKWRK(*)
C
C DUMMY SUBROUTINE FOR DEBDF (USER CAN PUT JACOBIAN HERE)
C
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SJUMP(XX)
C
C THIS FUNCTION IS USED WITH THE ROOT SOLVER ZEROIN TO COMPUTE THE
C JUMP CONDITIONS FOR A FROZEN SHOCKS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NEL = 6, NK = 500, NQ = 1+NK, SMALL = 1.D-100)
      PARAMETER (LRW = 250 + NQ*(10 + NQ), LIW = 55 + NQ)
      PARAMETER (LENIWK = 40000, LENWK = 40000, lencwk = 2000)
      COMMON/PARAM/RU,WT(NK),P1,RHO1,PSC,UVEL,SON,XM,UDET,X(NK),T,
     &             P,RHO,RPAR(4),ALPHA,Y,WTM,T1,KK
      COMMON/WRK/ICKWRK(LENIWK),CKWRK(LENWK)
      DIMENSION Z(NQ),ZP(NQ),Y(NK)
C
C COMPUTE THE PRESSURE USING THE DETONATION VELOCITY AND THE
C MOMENTUM JUMP CONDITION WITH THE ASSUMED DENSITY RATIO 'XX'
C
      P=P1+RHO1*UDET*UDET*(1.D0-XX)
C
C USE ASSUMED DENSITY RATIO TO GET DENSITY BEHIND THE SHOCK
      RHO=RHO1/XX
C
C COMPUTE TEMPERATURE AND  ENTHALPY BEHIND
C THE SHOCK FROM THE PRESSURE COMPUTED WITH THE ASSUMED DENSITY RATIO
C
      CALL FGMMWX(X,ICKWRK,CKWRK,WTM)
      T=P*WTM/(RU*RHO)
      CALL FGHBMS(T,Y,ICKWRK,CKWRK,GH)
C
C COMPUTE THE TEMPERATURE AND ENTHALPY IN FRONT OF THE SHOCK
C
      T1=P1*WTM/(RU*RHO1)
      CALL FGHBMS(T1,Y,ICKWRK,CKWRK,GH1)
C
C COMPUTE THE MISMATCH IN THE ENERGY JUMP CONDITION USING THE QUANTITIES
C FOUND ABOVE. FZERO WILL DRIVE THIS MISMATCH BELOW THE SPECIFIED TOLERANCE
C TO DETERMINE THE CORRECT DENSITY RATIO 'XX'
C
      SJUMP=GH-GH1-UDET*UDET*(1.D0-XX*XX)/2.D0
      RETURN
      END
c
      SUBROUTINE IPOFILU(LUNIT, MAXLEN, PROMPT, STAT, FILNAM)
C***BEGIN PROLOGUE  IPOFIL
C***DATE WRITTEN   850626
C***REVISION DATE  850626
C***CATEGORY NO.  J4.,K2.
C***KEYWORDS  INPUT,INTERACTIVE,FILE OPEN
C***AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C***PURPOSE  Opens a file from interactive FORTRAN code.
C***DESCRIPTION
C
C-----------------------------------------------------------------------
c  IPOFIL is intended to make opening of files from FORTRAN codes
c  simpler for the programmer.  It prints a prompt at the screen
c  to accept a file name, then attempts to open a file of the desired
c  type.  If the file open fails, due to an improper name, type, or
c  other system-dependent reasons, a diagnostic will be printed,
c  and the user prompted again, until the open is successful.
c  The file name has a default value of the entry value, which the
c  user can select by entering a carriage return only in response
c  to the prompt.
c
c  arguments: (I=input, O=output)
c  ------------------------------
c  LUNIT (I) - logical unit number desired
C  MAXLEN (I) - max. # characters in file name (ignored if =0)
C  PROMPT (I) -  character string for prompting user
C  STAT (I) - file status (character variable)
C  FILNAM (I/O) - name of file (character variable)
c
c  IPOFIL is part of IOPAK, and is written in ANSI FORTRAN 77
C-----------------------------------------------------------------------
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ILASCH
C***END PROLOGUE IPOFIL
C
      CHARACTER* (*)FILNAM, PROMPT, STAT
      CHARACTER*40 LOCNAM
C
C***FIRST EXECUTABLE STATEMENT IPOFIL
      IOS = 0
      LIN = ILASCH(FILNAM)
      LPR = ILASCH(PROMPT)
80    CONTINUE
      PRINT '(/1X,A,A4,A,A2)' , PROMPT(1:LPR), '?  <' ,
     1   FILNAM(1:LIN), '> '
      READ '(A)' , LOCNAM
      LLOC = ILASCH(LOCNAM)
C
      IF (MAXLEN .NE. 0 .AND. LLOC .GT. MAXLEN) THEN
         PRINT '(A,I2)' ,
     1   '!!! ERROR - MAXIMUM # CHARACTERS IN NAME IS ' ,
     2   MAXLEN
         GO TO 80
      ENDIF
C
      IF (LOCNAM .NE. ' ') THEN
         FILNAM = LOCNAM(1:LLOC)
         LIN    = LLOC
      ENDIF
C
      OPEN (UNIT = LUNIT, FILE = FILNAM(1:LIN), STATUS = STAT,
     1   FORM='UNFORMATTED',ERR = 90, IOSTAT = IOS)
      RETURN
C
90    CONTINUE
      PRINT '(A,A,A,I3,A)' , ' !! ERROR IN OPENING FILE ' ,
     1   FILNAM(1:LIN), ' IOS =' , IOS, ' - TRY AGAIN'
      GO TO 80
C
      END
      subroutine tim(temp)
c  **** NOTE - you may need to edit this routine ****
c Machine specific call to get time       
c  MS Fortran 
c      integer*2 ihr,imin,isec,ihun
c UNIX
      integer iarray(3)
      character*8 temp
c  MS Fortran 
c      call gettim(ihr,imin,isec,ihun)
C UNIX
c      call itime(iarray)
c MS Fortran
c      write(temp,'(I2,'':'',I2,'':'',I2)')ihr,imin,isec
C UNIX
      write(temp,'(I2,'':'',I2,'':'',I2)')iarray(1),iarray(2),iarray(3)
      return
      end
c
      subroutine dat(temp)
c  **** NOTE - you may need to edit this routine ****
c Machine specific call to get date       
c  MS Fortran 
c      integer*2 iyr,imn,idy
C  UNIX
      integer iarray(3)
      character*10 temp
c  MS Fortran 
c      call getdat(iyr,imn,idy)
c  UNIX
c      call idate( iarray)
c MS Fortran
c      write(temp,'(I2,''-'',I2,''-'',I2)')idy,imn,iyr
C UNIX
      write(temp,'(I2,''-'',I2,''-'',I4)')iarray(1),iarray(2),iarray(3)
      return
      end
c


