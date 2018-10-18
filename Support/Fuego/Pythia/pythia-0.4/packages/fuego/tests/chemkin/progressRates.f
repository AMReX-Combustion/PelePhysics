      PROGRAM rates
c
c This program reads in a mechanism file, set the molar
c concentrations of some species to some values and 
c compute the equilibrium constants and rate of progress
c at T=1000K
C
      implicit  double precision (a-h,o-z)
      parameter (leniwk = 40000)
      parameter (lenwk  = 40000)
      parameter (lencwk = 2000)
      parameter (nmax = 500)
      parameter (link = 5)
      parameter (lout = 7)
      dimension ickwrk(leniwk), rckwrk(lenwk)
      character cckwrk(lencwk)*16
      dimension c(nmax), q(nmax), eq_c(nmax)
      
      write(*,*) 'Opening chem.bin... better be there :> '
      open(unit=link, file='chem.bin', form='unformatted')
      write(*,*) 'Opening chemkin log file: ck.log'
      open(unit=lout, file='ck.log')
      
      call ckinit(leniwk, lenwk, lencwk, link, lout, ickwrk, 
     &            rckwrk, cckwrk)
      call ckindx(ickwrk, ckwrk, mm, ns, nr, nfit)

c     set up x here
      do i = 1, ns
         c(i) = 1.0d0
      end do

      do i = 1, nr
         q(i) = 0.0d0
         eq_c(i) = 0.0d0
      end do

      tkelv  = 1000d0

      call ckqc(tkelv, c, ickwrk, rckwrk, q)
      call ckeqc(tkelv, c, ickwrk, rckwrk, eq_c)

      write(*,*) 'rates = ( '
      do i = 1, nr
         write(*,'(A4E15.7A2E15.7A2)') '   (', q(i), ',', eq_c(i), '),'
      end do
      write(*,*) '     )'
      END

c 
c END OF FILE
