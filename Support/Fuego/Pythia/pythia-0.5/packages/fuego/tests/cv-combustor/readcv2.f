      subroutine readcv2(case,rho1,t1,p1,e1,plot,iplot,x,kk,ksym,
     &rtol,atol,dinc,ddist,tmax,ushock,noshock)
      implicit real*8 (A-H,O-Z)
      parameter (kkey = 14)
c parses input data file to get values associated with key words.
      character line*80,filold*40,po*40,old*3,kray(kkey)*5,po2*40
      character case*80,yorn*1,ksym(*)*16
      dimension x(*),val(10)
      logical plot, open, kerr,ierr,dinc,noshock
      data filold/'cv.inp'/,po/'Input file name'/,old/'old'/
      data po2/'Error in input file, ignore'/, yorn/'Y'/,open/.FALSE./
C define key words
      data kray(1)/'CASE'/,kray(2)/'RHO1'/,kray(3)/'P1'/,kray(4)/'T1'/,
     & kray(5)/'PLOT'/,kray(6)/'REAC'/,kray(7)/'END'/,kray(8)/'RTOL'/,
     & kray(9)/'ATOL'/,kray(10)/'TMAX'/,kray(11)/'USHK'/,
     & kray(12)/'NOSHK'/kray(13)/'STOP'/,kray(14)/'E1'/
c input file unit numbers
      data lufi/4/, lout/6/, small/1.d-100/
c initialize distance and plot default
      dinc = .false.
      plot = .FALSE.
      noshock = .false.
      ushock = 0.0d0
      rho1 = 0.0d0
      T1 = 0.0d0
      P1 = 0.0d0
      e1 = 0.0d0
      TMAX = 1.d0
      DDIST = 1.d0
c open input data file
      if (.not.open) then
         call ipofil(lufi,0,po,old,filold)
         open = .TRUE.
      endif
c read it until end
10    read(lufi,'(A)',end=100)line
cd    write(6,'('' line :'',1x,A)') line
c size it up, skip blank lines
      ilast = ilasch(line)
cd      write (6,'(A, i5)') ' line ends at ',ilast
      if (ilast.eq.0) then
cd         write (6,'(A)') ' Blank line'
         goto 10
      endif
c first nonblank
      istart = ifirch(line)  
cd      write (6,'(A, i5)') ' line starts  at',istart
c comment card ?
      if (line(istart:istart).eq.'!') then
cd         write (6,'(A)') ' Comment line'
        goto 10
      endif
c trailing comment ?
      do 15 i  = ilast,istart, -1
        if (line(i:i).eq.'!') then
           ilast = i-1
cd           write (6,'(A)') ' Found a trailing comment'
           goto 17
        endif
15    continue
c isolate the first substring
17    do 20 i=istart, ilast
         if (line(I:I).eq.' '.or.line(I:I).eq.',') then
            keyend = i-1
cd            write (6,'(A, i5)') ' key ends at',keyend
            goto 30
         endif
20    continue
c entire line is a keyword
       keyend = ilast
cd            write (6,'(A)') ' keyword only on this line'
c compare with keyword list
30     key = 0
       do 40 k=1,kkey
          if (line(istart:keyend).eq.kray(k)) then
             key = k
cd           write (6,'(A, i5)') ' key ',key
             goto 50
          endif
40     continue
50     continue
       if (key.eq.0) then
c error! - couldn't find the keyword in the list!
          write(6,'('' Error in input - keyword not found'')')
          write(6,'(1x,A)')line
          goto 10
       endif
c process keywords
       key1 = keyend + 2
       if (key.eq.1) then
          case = line(key1:ilast)
          write(lout, '(1x,A5,1X,80A)')kray(key),case
          goto 10
       elseif (key.eq.2) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             rho1 = xx
          write(lout, '(1x,A5,1X,1pg14.6)') kray(key),rho1
             goto 10
          endif
       elseif (key.eq.3) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             p1 = xx
          write(lout, '(1x,A5,1X,1Pg14.6)')kray(key),p1
             goto 10
          endif
       elseif (key.eq.4) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             t1 = xx
          write(lout, '(1x,A5,1X,1Pg14.6)')kray(key),t1
             goto 10
          endif
       elseif (key.eq.14) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             e1 = xx
          write(lout, '(1x,A5,1X,1Pg14.6)')kray(key),e1
             goto 10
          endif
       elseif (key.eq.5) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
             if (.not.ierr) then
             iplot = int(xx)
             write(lout, '(1x,A5,1X,1I3)')kray(key),iplot
             plot = .true.
             endif
             goto 10
       elseif (key.eq.6) then
          write(lout, '(1x,A5)') kray(key)
c read reactants
   60 CONTINUE
      READ  (lufi,'(1A80)',END=70)line
cd      write (lout,'(1x, 1a80)')line
      ILEN = INDEX (line, '!')
      IF (ILEN .EQ. 1) GO TO 60
C
      ILEN = ILEN - 1
      IF (ILEN .LE. 0) ILEN = LEN(line)
      IF (INDEX(line(:ILEN), 'END') .EQ. 0) THEN
         IF (line(:ILEN) .NE. ' ') THEN
            CALL CKSNUM (line(:ILEN), 1, LOUT, ksym, KK, knum,
     1                   NVAL, val, IERR)
            IF (IERR) THEN
               WRITE (LOUT,*) ' Error reading moles...'
               KERR = .TRUE.
            ELSE
               x(knum) = val(1)
         WRITE (lout,'(1X,1A10,'' X = '',G12.4)') ksym(knum)(:10),
     &               x(knum)
            ENDIF
         ENDIF
         GO TO 60
      ENDIF
C
   70 CONTINUE
C
C        NORMALIZE THE MOLE FRACTIONS
C
      XTOT=0.0D0
      DO 80 K=1,KK
      x(k) = max(x(k),small)
      XTOT = XTOT + x(K)
   80 CONTINUE
      DO 90 K=1,KK
      x(K) = x(K)/XTOT
   90 CONTINUE
      GOTO 10
       elseif (key.eq.7) then
          write(lout, '(1x,A5)') kray(key)
c end of this case
          return
       elseif (key.eq.8) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             rtol = xx
          write(lout, '(1x,A5,1X,1pg14.6)') kray(key),rtol
             goto 10
          endif
       elseif (key.eq.9) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             atol = xx
          write(lout, '(1x,A5,1X,1pg14.6)') kray(key),atol
             goto 10
          endif
       elseif (key.eq.10) then
          call ckxnum(line(key1:ilast),1,6,nf,val,ierr)
          if (.not.ierr) then
             tmax = val(1)
          write(lout, '(1x,A5,(2x,1pg14.6))') kray(key),tmax
          dinc=.true.
             goto 10
          endif
       elseif (key.eq.11) then
          call ckxnum(line(key1:ilast),1,6,nf,xx,ierr)
          if (.not.ierr) then
             ushock = xx
             write(lout, '(1x,A5,1x,1pg14.6)') kray(key),ushock
             goto 10
          endif
       elseif (key.eq.12) then
             write(lout, '(1x,A5)') kray(key)
             NOSHOCK = .true.
             goto 10
       elseif(key.eq.13) then
          write(lout, '(1x,A5)') kray(key)
c end of input
          goto 100
       endif
c drop through point for all errors
      call ipyorn(po2,yorn)
      if (yorn.eq.'Y') then
        goto 10
      else
        close(unit = lufi)
        stop
      endif
100   continue
c   all keywords have been found
      close(unit = lufi)
      case = 'STOP'
      return
      end
C
      SUBROUTINE IPYORN(PROMPT, YORN)
C***BEGIN PROLOGUE  IPYORN
C***DATE WRITTEN   850626
C***REVISION DATE  850626
C***CATEGORY NO.  J4.
C***KEYWORDS  INPUT,INTERACTIVE
C***AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C***PURPOSE  Prints yes or no question at terminal, accepts
c            user response.
C***DESCRIPTION
C-----------------------------------------------------------------
c  IPYORN is designed to handle "yes or no" questions often
c  used to control the flow of interactive codes.  It prints
c  a prompt followed by a question mark and default response,
c  reads the user response, and continually diagnoses any errors
c  until the response is valid.  Valid responses are any string
c  beginning with a 'y' or an 'n', either upper case or lower, or
c  a carriage return.  The carriage return is used to select the
c  default value, which is the entry value of YORN.
c
c  arguments: (I=input, O=output)
c  ------------------------------
c  PROMPT (I) - The prompting string. (Character variable)
c  YORN (I/O) - The user response, converted to upper case
c               for return
c
c  IPYORN is part of IOPAK, and is written in ANSI FORTRAN 77
C-----------------------------------------------------------------
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IFIRCH,ILASCH
C***END PROLOGUE IPYORN
C
      CHARACTER* (*)PROMPT, YORN
      CHARACTER*10 ANS, SUB
      CHARACTER*80 PRPLUS
C
C***FIRST EXECUTABLE STATEMENT IPYORN
100   CONTINUE
      IFIRS = IFIRCH(PROMPT)
      ILAS = ILASCH(PROMPT)
      PRPLUS = PROMPT(IFIRS:ILAS) // '? <' // YORN(1:1) // '>  '
      PRINT '(1X,A)', PRPLUS
      READ '(A)', ANS
C
      IFIRS = IFIRCH(ANS)
      ILAS = ILASCH(ANS)
c
c   check for default input of carriage return
c
      IF (IFIRS .EQ. 0) THEN
C
C            . . . convert input value to upper case, if necessary,
c                  for return
C
         IF (YORN .EQ. 'n') YORN = 'N'
         IF (YORN .EQ. 'y') YORN = 'Y'
         RETURN
      ENDIF
C
c   form substring of significant characters, and check first
c   character only
c
      SUB = ANS(IFIRS:IFIRS)
C
      IF (SUB .EQ. 'Y' .OR. SUB .EQ. 'y') THEN
         YORN = 'Y'
      ELSEIF (SUB .EQ. 'N' .OR. SUB .EQ. 'n') THEN
         YORN = 'N'
      ELSE
         PRINT *,
     1   ' !!! ERROR . . ACCEPTABLE RESPONSES ARE: YES, NO,'
         PRINT *,
     1   ' Y, N (UPPER OR lower CASE), OR A CARRIAGE RETURN.'
         GO TO 100
      ENDIF
C
      PRINT *, ' '
C
      END
      SUBROUTINE IPOFIL(LUNIT, MAXLEN, PROMPT, STAT, FILNAM)
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
     1   ERR = 90, IOSTAT = IOS)
      RETURN
C
90    CONTINUE
      PRINT '(A,A,A,I3,A)' , ' !! ERROR IN OPENING FILE ' ,
     1   FILNAM(1:LIN), ' IOS =' , IOS, ' - TRY AGAIN'
      GO TO 80
C
      END
