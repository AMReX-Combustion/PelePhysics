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
      dimension c(nmax), wdot(nmax)
      
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

      tkelv  = 1000d0

      call ckwc(tkelv, c, ickwrk, rckwrk, wdot)

      write(*,*) 'productionRates = ( '
      do i = 1, ns
         write(*,'(1x,E15.7A2)') wdot(i), ','
      end do
      write(*,*) '     )'
      END

c 
c END OF FILE
