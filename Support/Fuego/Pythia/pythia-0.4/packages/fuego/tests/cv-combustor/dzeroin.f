      SUBROUTINE ZEROIN(F,B,C,RE,AE,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2646
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980
C                   *************************
C                   *       ISSUED BY       *
C                   *  SANDIA LABORATORIES, *
C                   *   A PRIME CONTRACTOR  *
C                   ********     TO THE     *
C                          *  UNITED STATES *
C                          *   DEPARTMENT   *
C                          *       OF       *
C                          *     ENERGY     *
C      *********************  ---NOTICE---  *********************
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *
C      *               NOR ANY OF THEIR EMPLOYEES,              *
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *
C      *          **********    ACCURACY,   **********          *
C      *          *        *  COMPLETENESS  *        *          *
C      *          *        *  OR USEFULNESS *        *          *
C      *          *        *     OF ANY     *        *          *
C      *          *        *  INFORMATION,  *        *          *
C      *          *        *   APPARATUS,   *        *          *
C      *       ****        *     PRODUCT    *        ****       *
C      *       *           *   OR PROCESS   *           *       *
C      *       *           *   DISCLOSED,   *           *       *
C      *       *           *  OR REPRESENTS *           *       *
C      *       *          **    THAT ITS    **          *       *
C      *       *          **  USE WOULD NOT **          *       *
C      *********          **    INFRINGE    **          *********
C                         **    PRIVATELY   **
C                         **      OWNED     **
C                         **     RIGHTS.    **
C                         **                **
C                         **                **
C                         **                **
C                         ********************
C
C     BASED ON A METHOD BY T J DEKKER
C     WRITTEN BY L F SHAMPINE AND H A WATTS
C     MODIFIED FOR THE MATH LIBRARY BY C B BAILEY
C
C     ABSTRACT
C        ZEROIN SEARCHES FOR A ZERO OF A FUNCTION F(X) BETWEEN
C        THE GIVEN VALUES B AND C UNTIL THE WIDTH OF THE INTERVAL
C        (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE SPECIFIED BY
C        THE STOPPING CRITERION, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C        THE METHOD USED IS AN EFFICIENT COMBINATION OF BISECTION AND
C        THE SECANT RULE.  IN ORDER TO INSURE THAT ZEROIN WILL CONVERGE
C        TO A ZERO, THE USER SHOULD PICK VALUES FOR B AND C AT WHICH
C        THE FUNCTION DIFFERS IN SIGN.
C
C     DESCRIPTION OF ARGUMENTS
C     F,B,C,RE AND AE ARE INPUT PARAMETERS
C     B,C AND IFLAG ARE OUTPUT PARAMETERS
C        F     - NAME OF THE REAL VALUED EXTERNAL FUNCTION.  THIS NAME
C                MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING
C                PROGRAM.  F MUST BE A FUNCTION OF ONE REAL ARGUMENT.
C        B     - ONE END OF THE INTERVAL (B,C).  THE VALUE RETURNED FOR
C                B USUALLY IS THE BETTER APPROXIMATION TO A ZERO OF F.
C        C     - THE OTHER END OF THE INTERVAL (B,C)
C        RE    - RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
C                IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
C                THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
C        AE    - ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
C                GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
C                NONZERO VALUE SHOULD BE CHOSEN FOR AE.
C        IFLAG - A STATUS CODE.  USER MUST CHECK IFLAG AFTER EACH CALL.
C                CONTROL RETURNS TO THE USER FROM ZEROIN IN ALL CASES.
C                XERROR DOES NOT PROCESS DIAGNOSTICS IN THESE CASES.
C                 1 B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C), AND
C                   F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C                 2 F(B) = 0.  HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
C                   COLLAPSED TO THE REQUESTED TOLERANCE.
C                 3 B MAY BE NEAR A SINGULAR POINT OF F(X).
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE AND THE FUNCTION CHANGES SIGN IN (B,C) BUT
C                   F(X) INCREASED IN MAGNITUDE AS (B,C) COLLAPSED,I.E.
C                     ABS(F(B OUT)) .GT. MAX(ABS(F(B IN)),ABS(F(C IN)))
C                 4 NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
C                   INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
C                   THE USER MUST EXAMINE THIS CASE AND DECIDE WHETHER
C                   B IS NEAR A LOCAL MINIMUM OF F(X), OR B IS NEAR A
C                   ZERO OF EVEN MULTIPLICITY, OR NEITHER OF THESE.
C                 5 TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.
C
C     REFERENCES
C       1.  L F SHAMPINE AND H A WATTS, ZEROIN, A ROOT-SOLVING CODE,
C           SC-TM-70-631, SEPT 1970.
C       2.  T J DEKKER, FINDING A ZERO BY MEANS OF SUCCESSIVE LINEAR
C           INTERPOLATION, *CONSTRUCTIVE ASPECTS OF THE FUNDAMENTAL
C           THEOREM OF ALGEBRA*, EDITED BY B DEJON AND P HENRICI, 1969.
C
C
C     ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS
C     DEFINED HERE TO BE THE VALUE FOR THE IBM PC DOUBLE PRECISION
C
      ER = 2.0D-13
C
C     INITIALIZE
      RW=DMAX1(RE,ER)
      AW=DMAX1(AE,0.0D0)
      IC=0
      ACBS=DABS(B-C)
      A=C
      T=A
      FA=F(T)
      T=B
      FB=F(T)
      FC=FA
      KOUNT=2
      FX=DMAX1(DABS(FB),DABS(FC))
C
    1 IF (DABS(FC) .GE. DABS(FB)) GO TO 2
C     PERFORM INTERCHANGE
      A=B
      FA=FB
      B=C
      FB=FC
      C=A
      FC=FA
C
    2 IF (FB .EQ. 0.0D0) GO TO 11
      CMB=0.5D0*(C-B)
      ACMB=DABS(CMB)
      TOL=RW*DABS(B)+AW
C
C     TEST STOPPING CRITERION
      IF (ACMB .LE. TOL) GO TO 10
C
C     CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
C     WHERE WE ARRANGE P .GE. 0.
C     THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
      P=(B-A)*FB
      Q=FA-FB
      IF (P .GE. 0.0D0) GO TO 3
      P=-P
      Q=-Q
C
C     UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
C     IN THE SIZE OF OUR BOUNDING INTERVAL.
    3 A=B
      FA=FB
      IC=IC+1
      IF (IC .LT. 4) GO TO 4
      IF (8.0D0*ACMB .GE. ACBS) GO TO 6
      IC=0
      ACBS=ACMB
C
C     TEST FOR TOO SMALL A CHANGE
    4 IF (P .GT. DABS(Q)*TOL) GO TO 5
C
C     INCREMENT BY TOLERANCE
      B=B+DSIGN(TOL,CMB)
      GO TO 7
C
C     ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.0D0
    5 IF (P .GE. CMB*Q) GO TO 6
C
C     INTERPOLATE
      B=B+P/Q
      GO TO 7
C
    6 B=0.5D0*(C+B)
C     BISECT
C
C     HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
    7 T=B
      FB=F(T)
      IF (FB .EQ. 0.0D0) GO TO 11
C
C     DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
      IF (DSIGN(1.0D0,FB) .NE. DSIGN(1.0D0,FC)) GO TO 8
      C=A
      FC=FA
    8 KOUNT=KOUNT+1
      IF (KOUNT .GT. 500) GO TO 15
      GO TO 1
C
C
C     FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
C
   10 IF (DSIGN(1.0D0,FB) .EQ. DSIGN(1.0D0,FC)) GO TO 13
      IF (DABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN
   11 IFLAG = 2
      RETURN
   12 IFLAG = 3
      RETURN
   13 IFLAG = 4
      RETURN
   15 IFLAG = 5
      RETURN
      END

