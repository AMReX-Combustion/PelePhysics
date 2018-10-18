      subroutine burncv(dt, phi, phihist, thist, rho, 
     &                   ene, nsp,itmax, nsteps)
c
c dt:       Will integrate ODE by dt
c phi:      array of size nsp. (initial condition)
c           will be replaced by final state on exit
c phihist:  Matrix (nsp,*). Will hold solution vector
c thist:    array of size >= itmax. Will hold output t
c rho,ene:  density and energy invariants of problem
c nsp:      number of species
c itmax:    maximum number of iterations 
c
      implicit real*8 (a-h,o-z)
      parameter (nk = 500)
      parameter (lrw = 250 + nk*(10 + nk), liw = 55 + nk)
      parameter (rtol = 1.0d-7, atol = 1.0d-8)
      dimension phi(*), phihist(nsp, *), thist(*)
      dimension info(15),elwrk(lrw),ielwrk(liw)
      dimension ckwrk(1), ickwrk(1)
      double precision xrho, xene
      common/burncvcommon/ xrho, xene
      external cvvf, cvvfjac
c
C*****************************************
c     set the integration control parameters for ddebdf
c
      do k1=1,15
         info(k1)=0
      end do
c     
c intermediate output mode and jac
      info(3) = 1
      info(5) = 1
c
c set variables in common block
      xrho = rho
      xene = ene
c
      time = 0.0d0
      thist(1) = time
c
      do k = 1, nsp
         phihist(k, 1) = phi(k)
      end do
c
      do i = 2, itmax
        call ddebdf(cvvf,nsp,time,phi,dt,info,rtol,atol,idid,
     *            elwrk,lrw,ielwrk,liw,ckwrk,ickwrk,cvvfjac)
        nsteps = i
        thist(i) = time
        do k = 1, nsp
           phihist(k, i) = phi(k)
        end do
c
c check the return from the integration
c
        if (idid.eq.-33) then
          write(lout,'(''Fatal error in ddebdf, idid ='',i4)') idid
          goto 999
        endif
c
c check to see whether integration reached tmax
        if (idid.gt.1) then
           goto 999
        endif

      end do
 999  continue
c
      return
      end
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      subroutine cvvf(time,phi,phidot,rckwrk,ickwrk)
c*********************************************************************
c RHS of CV Combuster in implicit form
c phi is specific mole numbers: [phi] = [molar conc] * 1000 / rho
c ********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter (nsp = 500)
      dimension phi(*),phidot(*), rckwrk(*), ickwrk(*)
      dimension y(nsp), wdot(nsp),tmparr(nsp)
      double precision xrho, xene
      common/burncvcommon/ xrho, xene
c
      call fgindx(ickwrk, rckwrk, idum1, ns, idum2, idum3)
      call fephity(phi, ickwrk, rckwrk, y)
      call feeytt( xene, y, ickwrk, rckwrk, tem)
      call fgpy( xrho, tem, y, ickwrk, rckwrk, pre)
      call fgwyp(pre, tem, y, ickwrk, rckwrk, wdot)
c     compute phidot
      do k = 1,ns
         phidot(k) = wdot(k) / (xrho/1.0d3)
      end do
      return
      end
c
      subroutine cvvfjac(time,phi,dphi,nrowpd,rckwrk,ickwrk)
c
c simple finite difference jacobian. Use this or set
c info(5) = 0 for ddebdf's fd jacobian.
c
      parameter (nsp = 500)
      implicit real*8 (a-h,o-z)
      dimension phi(*),rckwrk(*),ickwrk(*)
      dimension dphi(nrowpd,*)
      dimension fvec(nsp), fv(nsp)
      double precision xrho, xene
      common/burncvcommon/ xrho, xene
c
      call cvvf(time, phi, fvec, rckwrk, ickwrk)
      call fgindx(ickwrk,rckwrk, idum1, ns, idum2, idum3)
      ep = 1d-6
      do j = 1,ns
        temp = phi(j)
        h = ep*abs(temp)
        if (h .eq. 0) then
            h = ep
        end if
        phi(j) = temp + h
        h = phi(j) - temp
        call cvvf(time, phi, fv, rckwrk, ickwrk)
        phi(j) = temp
        do k = 1,ns
           dphi(k,j) = (fv(k)-fvec(k))/h 
        end do
      end do
c     
      return
      end
c
c
c END OF FILE  burncv.f
