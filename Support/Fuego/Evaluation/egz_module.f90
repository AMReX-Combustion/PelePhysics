module egz_module
  
  use fuego_chemistry

  implicit none

  private

  double precision, parameter :: Ru = 8.314d7
  double precision, parameter :: Patmos = 1.01325d6

  logical, save :: use_bulk_visc = .true.

  integer, parameter :: nfit=7
  integer, save :: iflag = -1
  integer, save :: np = -1
  integer, save :: ns, no
  double precision, allocatable, save :: wt(:), iwt(:), eps(:), sig(:), dip(:), pol(:), zrot(:)
  integer, allocatable, save :: nlin(:) 
  double precision, allocatable, save :: cfe(:,:), cfl(:,:), cfd(:,:,:)
  double precision, allocatable, save :: eps2(:,:)
  double precision, allocatable, save :: fita(:,:,:)
  double precision, parameter :: fita0(nfit) = (/ &
       .1106910525D+01, -.7065517161D-02, -.1671975393D-01, .1188708609D-01, &
       .7569367323D-03, -.1313998345D-02,  .1720853282D-03 /)

  !
  ! dimension(np,ns)
  !
  double precision, allocatable, save :: xtr(:,:), ytr(:,:), aux(:,:)
  double precision, allocatable, save :: cxi(:,:), cint(:,:)
  !
  ! dimension(np)
  !
  double precision, allocatable, save :: sumtr(:), wwtr(:)
  !
  double precision, allocatable, save :: dlt(:,:)

  !
  ! dimension(np,ns)
  !
  double precision, allocatable, save :: beta(:,:), eta(:,:), etalg(:,:), &
       rn(:,:), an(:,:), zn(:,:), dmi(:,:)
  ! 
  ! dimension(np,ns,ns)
  !
  double precision, allocatable, save :: G(:,:,:), bin(:,:,:), A(:,:,:)

  !$omp threadprivate(xtr,ytr,aux,cxi,cint,sumtr,wwtr,dlt,beta,eta,etalg)
  !$omp threadprivate(rn,an,zn,dmi,G,bin,A,np)

  public :: iflag
  public :: egz_init, egz_close, EGZINI, EGZPAR, EGZE1, EGZE3, EGZK1, EGZK3, EGZL1, EGZVR1
  ! egz_init and egz_close should be called outside OMP PARALLEL,
  ! whereas others are inside

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine egz_init(use_bulk_visc_in)
    logical, intent(in) :: use_bulk_visc_in

    use_bulk_visc = use_bulk_visc_in

    if (use_bulk_visc) then
       iflag = 5
    else
       iflag = 3
    end if
    
    call egtransetKK(ns)
    call egtransetNO(no)
    
    allocate(wt(ns))
    allocate(iwt(ns))
    allocate(eps(ns))
    allocate(sig(ns))
    allocate(dip(ns))
    allocate(pol(ns))
    allocate(zrot(ns))
    allocate(nlin(ns))
       
    allocate(cfe(no,ns))
    allocate(cfl(no,ns))
    allocate(cfd(no,ns,ns))
    
    call egtransetWT(wt)
    iwt = 1.d0/wt
    call egtransetEPS(eps)
    call egtransetSIG(sig)
    call egtransetDIP(dip)
    call egtransetPOL(pol)
    call egtransetZROT(zrot)
    call egtransetNLIN(nlin)
    call egtransetCOFETA(cfe)
    call egtransetCOFLAM(cfl)
    call egtransetCOFD(cfd)
    
    allocate(eps2(ns,ns))
    allocate(fita(nfit,ns,ns))
    
    call LEVEPS()

    call EGZABC(fita,fita0)

  end subroutine egz_init

  ! This subroutine should be called outside OMP PARALLEL
  subroutine egz_close()
    if (allocated(wt)) deallocate(wt)
    if (allocated(iwt)) deallocate(iwt)
    if (allocated(eps)) deallocate(eps)
    if (allocated(sig)) deallocate(sig)
    if (allocated(dip)) deallocate(dip)
    if (allocated(pol)) deallocate(pol)
    if (allocated(zrot)) deallocate(zrot)
    if (allocated(nlin)) deallocate(nlin)
    if (allocated(cfe)) deallocate(cfe)
    if (allocated(cfl)) deallocate(cfl)
    if (allocated(cfd)) deallocate(cfd)
    if (allocated(eps2)) deallocate(eps2)    
    if (allocated(fita)) deallocate(fita)
    !$omp parallel
    call egz_close_np()
    !$omp end parallel
  end subroutine egz_close


  ! This subroutine can be called inside OMP PARALLEL
  subroutine EGZINI(np_in)
    integer, intent(in) :: np_in
    logical, save :: first_call = .true.
    !$omp threadprivate(first_call)

    if (first_call) then
       np = np_in
    end if

    if (first_call .or. np.ne.np_in) then

       call egz_close_np()

       np = np_in

       allocate(xtr(np,ns))
       allocate(ytr(np,ns))
       allocate(aux(np,ns))
       allocate(sumtr(np))
       allocate(wwtr(np))
       allocate(dlt(np,6))

       allocate(beta(np,ns))
       allocate(eta(np,ns))
       allocate(etalg(np,ns))
       allocate(rn(np,ns))
       allocate(an(np,ns))
       allocate(zn(np,ns))
       allocate(dmi(np,ns))

       if (iflag .gt. 1) then
          allocate(bin(np,ns,ns))
       end if
       if (iflag.eq.3 .or. iflag.eq.5) then
          allocate(G(np,ns,ns))
          allocate(A(np,ns,ns))
       end if

       if (iflag > 3) then
          allocate(cxi(np,ns))
          allocate(cint(np,ns))
       end if
       
    end if

    first_call = .false.

  end subroutine EGZINI


  ! This subroutine can be called inside OMP PARALLEL
  subroutine egz_close_np()
    if (allocated(xtr)) deallocate(xtr)
    if (allocated(ytr)) deallocate(ytr)
    if (allocated(aux)) deallocate(aux)
    if (allocated(cxi)) deallocate(cxi)
    if (allocated(cint)) deallocate(cint)
    if (allocated(sumtr)) deallocate(sumtr)
    if (allocated(wwtr)) deallocate(wwtr)
    if (allocated(dlt)) deallocate(dlt)
    if (allocated(beta)) deallocate(beta)
    if (allocated(eta)) deallocate(eta)
    if (allocated(etalg)) deallocate(etalg)
    if (allocated(rn)) deallocate(rn)
    if (allocated(an)) deallocate(an)
    if (allocated(zn)) deallocate(zn)
    if (allocated(dmi)) deallocate(dmi)
    if (allocated(G)) deallocate(G)
    if (allocated(bin)) deallocate(bin)
    if (allocated(A)) deallocate(A)
  end subroutine egz_close_np


  subroutine LEVEPS()
    double precision, parameter :: pi = 3.1415926535D0, &
         fac = 1.0D-12, dipmin = 1.0D-20, boltz = 1.38056D-16
    integer :: j, k
    double precision :: rooteps(ns)
    do j=1,ns
       rooteps(j) = sqrt(EPS(j))
    end do
    do j=1,ns
       !DEC$ IVDEP
       do k=1,j
          IF((DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!                K IS POLAR, J IS NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.0D0 + 0.25D0*(POL(J)/SIG(J)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP(K)**2/(EPS(K)*SIG(K)**3)) * &
                  rooteps(k)/rooteps(j)
          ELSE IF((DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN)) THEN
!-----------------------------------------------------------------------
!             J IS POLAR, K IS NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.0D0 + 0.25D0*(POL(K)/SIG(K)**3) * &
                  (FAC/BOLTZ) * &
                  (DIP(J)**2/(EPS(J)*SIG(J)**3)) * &
                  rooteps(j)/rooteps(k)
          ELSE
!-----------------------------------------------------------------------
!              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
!-----------------------------------------------------------------------
             eps2(K,J) = 1.d0
          ENDIF
          eps2(K,J) = log(rooteps(j)*rooteps(k)* eps2(K,J)*eps2(K,J))
       end do
    end do
    do j=1,ns
       do k=j+1,ns
          eps2(k,j) = eps2(j,k)
       end do
    end do
  end subroutine LEVEPS


  subroutine egzabc(FA, FA0)
    double precision, intent(in) :: FA0(nfit)
    double precision, intent(out) :: FA(nfit,ns,ns)
    integer i,j,k,l,m,mm
    double precision :: SUMA, prod
    DO J = 1, NS
       DO I = J, NS
          do m = 1, nfit
             SUMA = 0.0D0
             mm   = m - 1
             do k = mm, nfit-1
                prod = 1.0d0
                do l = 1, k-mm
                   prod = prod * (-eps2(i,j)) * dble(mm+l) / dble(l)
                enddo
                SUMA = SUMA + FA0(k+1) * PROD
             enddo
             FA(m,I,J) = SUMA
          enddo
       ENDDO
    ENDDO
    DO J = 1, NS
       DO I = 1, J-1
          do m = 1, nfit
             FA(m,I,J) = FA(m,J,I)
          enddo
       ENDDO
    ENDDO
  end subroutine egzabc


  ! This subroutine can be called inside OMP PARALLEL
  subroutine EGZPAR(T, X, cpms)
    double precision, intent(in) :: T(np), X(np,ns)
    double precision, intent(in), optional :: cpms(np,ns)

    integer :: i, n
    double precision :: aaa(np)
    double precision, parameter :: sss = 1.d-16

    call LZPAR(T, cpms)

!-----------------------------------------------------------------------
!     Add a small constant to the mole and mass fractions
!-----------------------------------------------------------------------
      sumtr = 0.d0
      do n=1,ns
         do i=1,np
            sumtr(i) = sumtr(i) + X(i,n)
         end do
      enddo

      do i=1,np
         aaa(i)  = sumtr(i) / dble(ns)
      end do

      wwtr = 0.0d0
      do n=1,ns
         do i=1,np
            xtr(i,n) = X(i,n) + sss*(aaa(i) - X(i,n))
            wwtr(i) = wwtr(i) + xtr(i,n) * wt(n)
         end do
      end do

!-----------------------------------------------------------------------
!     AUX(i) = \sum_{j .ne. i} YTR(j)
!-----------------------------------------------------------------------
      aaa = 0.d0
      do n=1,ns
         do i=1,np
            ytr(i,n) = xtr(i,n) * wt(n) / wwtr(i)
            aaa(i) = aaa(i) + ytr(i,n)
         end do
      end do
      do n=1,ns
         do i=1,np
            aux(i,n) = aaa(i) - ytr(i,n)
         end do
      end do

  end subroutine EGZPAR


  subroutine LZPAR(T, cpms)
    double precision, intent(in) :: T(np)
    double precision, intent(in), optional :: cpms(np,ns)
    integer :: i, m, n
    double precision :: tmp(np), crot(np) 
    double precision :: wru, dr, sqdr, dr32, aaaa1, dd, sqdd, dd32, bbbb
    double precision, parameter :: PI1=1.d0/3.1415926535D0, PI32O2=2.7842D+00, &
         P2O4P2=4.4674D+00, PI32=5.5683D+00

    !DEC$ SIMD
    do i=1,np
       dlt(i,1) = log(T(i))
       dlt(i,2) = dlt(i,1) * dlt(i,1)
       dlt(i,3) = dlt(i,2) * dlt(i,1)
       dlt(i,4) = dlt(i,3) * dlt(i,1)
       dlt(i,5) = dlt(i,4) * dlt(i,1)
       dlt(i,6) = dlt(i,5) * dlt(i,1)
    end do

    do n=1,ns
       !DEC$ SIMD
       do i=1,np
          etalg(i,n) = cfe(1,n) + cfe(2,n)*dlt(i,1) + cfe(3,n)*dlt(i,2) + cfe(4,n)*dlt(i,3)
          eta(i,n) = exp(etalg(i,n))
       end do
    end do

    if (iflag .le. 1) return

    do n=1,ns
       do m=1,n-1
          !DEC$ SIMD
          do i=1,np
             tmp(i) = -(cfd(1,m,n)+cfd(2,m,n)*dlt(i,1)+cfd(3,m,n)*dlt(i,2) &
                  + cfd(4,m,n)*dlt(i,3))
          end do
          !DEC$ SIMD
          do i=1,np
             bin(i,m,n) = exp(tmp(i))
             bin(i,n,m) = bin(i,m,n)
          end do
       end do
       do i=1,np
          bin(i,n,n) = 0.d0
       end do
    end do

    if (iflag .le. 2) return

    if (iflag.eq.3 .or. iflag.eq.5) then
       do n=1,ns
          do m=1,n-1
             do i=1,np
                A(i,m,n) = fita(1,m,n) + fita(2,m,n)*dlt(i,1) + fita(3,m,n)*dlt(i,2) &
                     + fita(4,m,n)*dlt(i,3) + fita(5,m,n)*dlt(i,4) &
                     + fita(6,m,n)*dlt(i,5) + fita(7,m,n)*dlt(i,6)
                A(i,n,m) = A(i,m,n)
             end do
          end do
          do i=1,np
             A(i,n,n) = fita(1,n,n) + fita(2,n,n)*dlt(i,1) + fita(3,n,n)*dlt(i,2) &
                  + fita(4,n,n)*dlt(i,3) + fita(5,n,n)*dlt(i,4) &
                  + fita(6,n,n)*dlt(i,5) + fita(7,n,n)*dlt(i,6)
          end do
       end do
    end if

    if (iflag .eq. 3) return

!-----------------------------------------------------------------------
!         COMPUTE PARKER CORRECTION FOR ZROT
!         AND ALSO THE ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
!-----------------------------------------------------------------------
    do n=1,ns
       select case(nlin(n))
       case (0)
          do i=1,np
             crot(i) = 0.d0
             cint(i,n) = 0.d0
          end do
       case (1)
          wru = wt(n) / Ru
          do i=1,np
             crot(i) = 1.d0
             cint(i,n) = cpms(i,n) * wru - 2.50d0
          end do
       case (2)
          wru = wt(n) / Ru
          do i=1,np
             crot(i) = 1.5d0
             cint(i,n) = cpms(i,n) * wru - 2.50d0
          end do
       case default
          print *, "EFZ: wrong value in nlin"
          stop
       end select
       
       dr = eps(n) / 298.d0
       sqdr = sqrt(dr)
       dr32 = sqdr*dr
       aaaa1 = 1.d0/((1.0d0 + PI32O2*sqdr + P2O4P2*dr + PI32*dr32) * max(1.0d0, zrot(n)))

       do i=1,np
          dd = eps(n) / T(i)
          sqdd = sqrt(dd)
          dd32 = sqdd*dd
          bbbb = (1.0d0 + PI32O2*sqdd + P2O4P2*dd + PI32*dd32) 
          cxi(i,n) = crot(i) * PI1 * bbbb * aaaa1
       end do
    end do

    return
  end subroutine LZPAR


  ! shear viscosity
  subroutine EGZE1(alpha, X, mu)
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: mu(np)

    integer :: i, n
    double precision :: alpha1

    mu = 0.d0
    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*etalg(i,n)
          end do
       end do
       do i=1,np
          mu(i) = exp(mu(i))
       end do
    else if (alpha .eq. 1.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*eta(i,n)
          end do
       end do
    else if (alpha .eq. -1.d0) then
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)/eta(i,n)
          end do
       end do
       do i=1,np
          mu(i) = 1.d0/mu(i)
       end do
    else
       do n=1,ns
          do i=1,np
             mu(i) = mu(i) + X(i,n)*exp(alpha*etalg(i,n))
          end do
       end do
       alpha1 = 1.d0/alpha
       do i=1,np
          mu(i) = mu(i)**alpha1
       end do       
    end if

  end subroutine EGZE1


  ! shear viscosity
  subroutine EGZE3(T, mu)
    double precision, intent(in) :: T (np)
    double precision, intent(out) :: mu(np)

    integer :: i, n

    call EGZEMH(T)
    
    rn = beta

    call EGZCG1(1)

    mu = 0.d0
    do n=1,ns
       do i=1,np
          mu(i) = mu(i) + an(i,n) * beta(i,n)
       end do
    end do
  end subroutine EGZE3

  subroutine EGZEMH(T)
    double precision, intent(in) :: T (np)
    integer :: i, m, n
    double precision :: FAC(np), CCC, wtfac, wtmn, wtnm, aaa

    do i = 1, np
       ! EVALUATE THE MATRIX H
       ! Note:    FAC * BIN = 2 W_{ij} / \eta_{ij} / A_{ij}
       FAC(i) = (6.0D0 * RU / ( 5.0D0 * PATMOS )) * T(i)
    end do

    CCC = 5.0D0 / 3.0D0

    do n=1,ns
       do i=1,np
          G(i,n,n) = xtr(i,n) * xtr(i,n) / eta(i,n)
          ! EVALUATE THE RHS BETA
          beta(i,n) = xtr(i,n)
       end do
    end do

    do n=1,ns
       do m=1,n-1
          wtfac = 1.d0/(wt(m) + wt(n))
          wtmn = wt(m)*iwt(n)
          wtnm = wt(n)*iwt(m)
          !DEC$ SIMD PRIVATE(aaa)
          do i=1,np
             aaa = bin(i,m,n) * xtr(i,n) * xtr(i,m) * FAC(i) * wtfac
             G(i,m,m) = G(i,m,m) + aaa*(A(i,m,n)*wtnm + CCC)
             G(i,m,n) = aaa * (A(i,m,n) - CCC)
             G(i,n,m) = G(i,m,n)
             G(i,n,n) = G(i,n,n) + aaa*(A(i,m,n)*wtmn + CCC)
          end do
       end do
    end do

  end subroutine EGZEMH


  subroutine EGZCG1(itmax)
    integer, intent(in) :: itmax

    integer :: niter, i, n
    double precision :: betan(np), aaa(np), bbb(np), ccc(np), temp(np,ns)

    do i=1,np
       aaa(i) = 0.d0
       betan(i) = 0.d0
    end do

    do n = 1, ns
       do i=1,np
          an(i,n) = 0.0d0
          zn(i,n) = 0.0d0
          dmi(i,n) = 1.0D0 / G(i,n,n)
          aaa(i) = aaa(i) + dmi(i,n) * rn(i,n)*rn(i,n)
       end do
    enddo

    do niter=1, itmax
       do n=1, ns
          do i = 1, np
             zn(i,n) = dmi(i,n)*rn(i,n) + betan(i)*zn(i,n)
          enddo
       end do

       CALL EGZAXS(G, zn, temp)

       bbb = 0.d0
       do n=1,ns
          do i=1,np
             bbb(i) = bbb(i) + zn(i,n) * temp(i,n)
          end do
       end do

       do n=1,ns
          do i=1,np
             an(i,n) = an(i,n) + aaa(i)/bbb(i)*zn(i,n)
             rn(i,n) = rn(i,n) - aaa(i)/bbb(i)*temp(i,n)
          end do
       end do

       if (niter .eq. itmax) exit

       ccc = 0.d0
       do n=1,ns
          do i=1,np
             ccc(i) = ccc(i) + dmi(i,n) * rn(i,n)*rn(i,n)
          end do
       end do

       do i=1, np
          betan(i) = ccc(i) / aaa(i)
          aaa(i) = ccc(i)
       end do

    end do

  end subroutine EGZCG1

  subroutine EGZAXS(AA, X, B) ! B = AA.X, AA is symmetric.
    double precision, intent(in) :: AA(np,ns,ns)
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: B(np,ns)
    integer :: i, m, n
    B = 0.d0
    do n=1,ns
       B(:,n) = 0.d0
       do m=1,ns
          do i=1,np
             B(i,n) = B(i,n) + AA(i,m,n) * X(i,m)
          end do
       end do
    end do
  end subroutine EGZAXS


  ! volume viscosity
  subroutine EGZK1(alpha, X, VV)
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: VV(np)

    integer :: i, n
    double precision :: sxp(np), ccc, vvv, alpha1

    sxp = 0.d0
    do n=1,ns
       do i=1,np
          if (cxi(i,n) .ne. 0.d0) then
             sxp(i) = sxp(i) + X(i,n)
          end if
       end do
    end do
    do i=1,np
       sxp(i) = 1.d0/sxp(i)
       VV(i) = 0.d0
    end do

    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = etalg(i,n) + log ( 0.25d0*ccc*ccc/cxi(i,n) )
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
       do i=1,np
          VV(i) = exp(VV(i))
       end do
    else if (alpha .eq. 1.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = eta(i,n) * 0.25d0*ccc*ccc/cxi(i,n) 
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
    else if (alpha .eq. -1.d0) then
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = cxi(i,n) / (eta(i,n) * 0.25d0*ccc*ccc)
                VV(i) = VV(i) + sxp(i) * X(i,n) * vvv
             end if
          end do
       end do
       do i=1,np
          VV(i) = 1.d0/VV(i)
       end do       
    else
       do n=1,ns
          do i=1,np
             if (cxi(i,n) .ne. 0.d0) then
                ccc = cint(i,n) / ( 1.5d0 + cint(i,n) )
                vvv = etalg(i,n) + log ( 0.25d0*ccc*ccc/cxi(i,n) )
                VV(i) = VV(i) + sxp(i) * X(i,n) * exp(alpha*vvv)
             end if
          end do
       end do
       alpha1 = 1.d0/alpha
       do i=1,np
          VV(i) = VV(i)**alpha1
       end do       
    end if

  end subroutine EGZK1

  ! volume viscosity
  subroutine EGZK3(T, VV)
    double precision, intent(in) :: T(np)
    double precision, intent(out) :: VV(np)
    
    integer :: i, m, n
    double precision :: ccc(np), wtfac, bb
    double precision, parameter :: denfac = 2.4d0 * Ru / Patmos

    if (.not. use_bulk_visc) then
       VV = 0.d0
       return
    end if

    ccc = 0.0d0
    do n=1,ns
       do i=1,np
          ccc(i) = ccc(i) + xtr(i,n)*cint(i,n)
       end do
    end do

    do n=1,ns
       do i=1,np
          G(i,n,n) = 4.d0*cxi(i,n)/eta(i,n)*xtr(i,n)*xtr(i,n)
          beta(i,n) = -xtr(i,n) * cint(i,n) / (ccc(i)+1.5d0)          
       end do
    end do

    do n=1,ns
       do m=1,n-1
          wtfac = iwt(m) + iwt(n)
          !DEC$ SIMD PRIVATE(bb)
          do i=1,np
             bb = xtr(i,m)*xtr(i,n)*bin(i,m,n)*denfac*T(i)*A(i,m,n)*wtfac
             G(i,m,m) = G(i,m,m) + bb*cxi(i,m)
!             G(i,n,m) = 0.d0
!             G(i,m,n) = 0.d0
             G(i,n,n) = G(i,n,n) + bb*cxi(i,n)
          end do
       end do
    end do

    VV = 0.d0
    do n=1,ns
       if (cxi(1,n) .eq. 0.d0) then
          do i=1,np
             VV(i) = VV(i) + beta(i,n) * beta(i,n) 
          end do
       else
          do i=1,np
             VV(i) = VV(i) + beta(i,n) * beta(i,n) / G(i,n,n)
          end do
       end if
    end do
  end subroutine EGZK3


  ! therml conductivity
  subroutine EGZL1(alpha, X, con)
    double precision, intent(in) :: alpha
    double precision, intent(in) :: X(np,ns)
    double precision, intent(out) :: con(np)

    integer :: i, n
    double precision :: asum(np), alpha1

    asum = 0.d0
    if (alpha .eq. 0.d0) then
       do n=1,ns
          do i=1,np
             asum(i) = asum(i) + X(i,n)*(cfl(1,n) + cfl(2,n)*dlt(i,1) &
                  + cfl(3,n)*dlt(i,2) + cfl(4,n)*dlt(i,3))
          end do
       end do
       do i=1,np
          con(i) = exp(asum(i)) 
       end do
    else
       alpha1 = 1.d0 / alpha
       do n=1,ns
          !DEC$ SIMD
          do i=1,np
             asum(i) = asum(i) + X(i,n)*exp(alpha*(cfl(1,n) + cfl(2,n)*dlt(i,1) &
                  + cfl(3,n)*dlt(i,2) + cfl(4,n)*dlt(i,3)))
          end do
       end do
       if (alpha .eq. 1.d0) then
          do i=1,np
             con(i) = asum(i)
          end do
       else if (alpha .eq. -1.d0) then
          do i=1,np
             con(i) = 1.d0/asum(i)
          end do
       else
          !DEC$ SIMD
          do i=1,np
             con(i) = asum(i)**alpha1
          end do
       end if
    end if
    
  end subroutine EGZL1


  ! rho * flux diffusion coefficients
  subroutine EGZVR1(T, D)
    double precision, intent(in) :: T(np)
    double precision, intent(out) :: D(np,ns)

    integer :: i, m, n
    double precision :: fac(np)

    D = 0.d0
    do n=1,ns
       do m=1,ns
          !DEC$ SIMD
          do i=1,np
             D(i,m) = D(i,m) + xtr(i,n)*bin(i,m,n)
          end do
       end do
    end do

    do i=1,np
       fac(i) = (Patmos/Ru) / T(i)
    end do

    do n=1,ns
       do i=1,np
          D(i,n) = wt(n) * fac(i) * aux(i,n) / D(i,n)
       end do
    end do

  end subroutine EGZVR1

end module egz_module

