      PROGRAM CKDRIV
      implicit none
      CHARACTER*72 CHEMKIN_input, THERMO_input, TRANLIB_input, log_file,
     &     CHEMKIN_linking_file, TRANLIB_c_file, file_list
C
      integer LIN, LINCK, LOUT, LTHRM, LINMC, LTRAN, LNAME
      PARAMETER (LIN = 34, LINCK = 35, LOUT = 36, LTHRM = 37,
     &     LINMC = 38, LTRAN = 39, LNAME=40)
C
      integer IDIM, KDIM, MDIM, MAXORD, MAXSP, MAXTB, MAXTP
      PARAMETER
     +  (IDIM = 1000,
     +   KDIM = 200,
     +   MDIM = 20,
     +   MAXORD = 10,
     +   MAXSP = 12,
     +   MAXTB = 20,
     +   MAXTP = 3)
C
      integer LCWORK, LIWORK, LLWORK, LRWORK
      PARAMETER
     +  (LCWORK = KDIM + MDIM,
     +   LIWORK = IDIM * (27 + MAXORD + 2 * MAXSP + MAXTB)
     +      + KDIM * (4 + MDIM),
     +   LLWORK = KDIM,
     +   LRWORK = IDIM * (33 + MAXORD + MAXSP + MAXTB)
     +      + 2 * KDIM * (4 * MAXTP - 3) + MDIM)
C
      integer LCMCWK, LIMCWK, LLMCWK, LRMCWK
      PARAMETER
     +  (LCMCWK = LCWORK + KDIM,
     +   LIMCWK = LIWORK + 3 * KDIM,
     +   LLMCWK = LLWORK,
     +   LRMCWK = LRWORK + 300 + KDIM * (48 + 4 * KDIM + 8 * MAXTP))


      CHARACTER*16 C(LCMCWK)
      INTEGER I(LIMCWK)
      LOGICAL L(LLMCWK)
      DOUBLE PRECISION R(LRMCWK)
C
      namelist / files / CHEMKIN_input, THERMO_input, TRANLIB_input, log_file,
     &     CHEMKIN_linking_file, TRANLIB_c_file

      call getarg(1,file_list)
      OPEN (LNAME, STATUS = 'OLD', FORM = 'FORMATTED', FILE = trim(file_list))
      read(LNAME,files)
      close(LNAME)

      OPEN (LIN, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(CHEMKIN_input))
      OPEN (LTHRM, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(THERMO_input))
      OPEN (LTRAN, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(TRANLIB_input))
      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(CHEMKIN_linking_file))
      OPEN (LINMC, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(TRANLIB_c_file))
      OPEN (LOUT, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(log_file))
C
      CALL CKINTP
     +  (MDIM, KDIM, IDIM, MAXTP, MAXSP, MAXTB, LIN, LOUT, LTHRM, LINCK,
     +   MAXORD, LIWORK, I, LRWORK, R, LCWORK, C, L)

      CALL TRANFT
     +  (LINMC, LINCK, LTRAN, MAXTP, LIMCWK, LRMCWK, LCMCWK, I, R,
     +   C)


      CLOSE (LIN)
      CLOSE (LTRAN)
      CLOSE (LTHRM)
      CLOSE (LINCK)
      CLOSE (LINMC)
      CLOSE (LOUT)

      END
