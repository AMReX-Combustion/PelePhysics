module actual_transport_module

  use amrex_fort_module, only : amrex_real
  use eos_type_module
  use transport_type_module
  use chemistry_module, only : Ru

  implicit none

  character (len=64) :: transport_name = "smp"
  logical, save, private :: smp_initialized = .false.
  integer, save, private :: npts_smp = 0
  real(amrex_real), save, allocatable :: Tloc(:), Yloc(:,:), Xloc(:,:), rholoc(:), wbar(:)
  real(amrex_real), save, allocatable :: muloc(:,:), lamloc(:,:), dbinloc(:,:,:)
  real(amrex_real), save, allocatable :: xiloc(:,:)
  real(amrex_real), save, allocatable :: logT(:,:)

  real(amrex_real), allocatable, save :: wt(:), iwt(:), eps(:), sig(:), dip(:), pol(:), zrot(:)
  integer, allocatable, save :: nlin(:)

  real(amrex_real), allocatable, save :: fitmu(:,:),fitlam(:,:),fitdbin(:,:,:)
  integer, save :: nfit 
  integer, parameter::  norder = 3
  integer, parameter::  iter = 1
  real(amrex_real), parameter ::  trace = 1.d-15

  logical, parameter :: use_bulk_viscosity = .true.

  real(amrex_real), save :: A_cst, B_cst, C_cst, D_cst, E_cst, F_cst, G_cst, H_cst, S_cst, W_cst
  real(amrex_real), save, allocatable :: Afac(:,:)
  real(amrex_real), save, allocatable :: Bfac(:,:)
  real(amrex_real), save, allocatable :: Sigmaij(:,:)
  real(amrex_real), parameter :: Avna = 6.022140857e23
  real(amrex_real), parameter :: PI = 4 * atan (1.0)
  real(amrex_real), save, allocatable :: Kappai(:)

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine actual_transport_init

    implicit none

    integer :: n,i,j
    integer :: iH2O
    
    call egtransetKK(nspec)
    call egtransetNO(nfit)

    if(.not.allocated(wt)) then

       allocate(wt(nspec))
       allocate(iwt(nspec))
       allocate(eps(nspec))
       allocate(sig(nspec))
       allocate(dip(nspec))
       allocate(pol(nspec))
       allocate(zrot(nspec))
       allocate(nlin(nspec))
       
       allocate(fitmu(nfit,nspec))
       allocate(fitlam(nfit,nspec))
       allocate(fitdbin(nfit,nspec,nspec))

       allocate(Sigmaij(nspec,nspec))

       allocate(Kappai(nspec))
       
    end if

    call egtransetWT(wt)
    call egtransetEPS(eps)
    call egtransetSIG(sig)
    call egtransetDIP(dip)
    call egtransetPOL(pol)
    call egtransetZROT(zrot)
    call egtransetNLIN(nlin)
    call egtransetCOFETA(fitmu)
    call egtransetCOFLAM(fitlam)
    call egtransetCOFD(fitdbin)

    iwt = 1.d0/wt

    !   do n=1,nspec

    !      write(6,*)"coeflam",fitlam(:,n)
    !      write(6,*)"coefeta",fitmu(:,n)
    !   do i=1,nspec
    !      write(6,*)"coefdebin",fitdbin(:,i,n)
    !   enddo
    !   enddo
    !   stop

    smp_initialized = .true.

    if (.not.allocated(Afac)) then
       allocate(Afac(10,4))
       allocate(Bfac(7,4))
    end if

    Afac(1,1) = 6.32402d0;     Afac(1,2) = 50.41190d0;  Afac(1,3) = -51.68010d0; Afac(1,4) = 1189.020d0

    Afac(2,1) = 0.12102e-2;    Afac(2,2) = -0.11536e-2; Afac(2,3) = -0.62571e-2; Afac(2,4) = 0.37283d-1

    Afac(3,1) = 5.28346d0;     Afac(3,2) = 254.209d0;   Afac(3,3) = -168.48100d0; Afac(3,4) = 3898.27000d0

    Afac(4,1) = 6.62263d0;     Afac(4,2) = 38.09570d0;  Afac(4,3) = -8.46414d0; Afac(4,4) = 31.41780d0 

    Afac(5,1) = 19.74540d0;    Afac(5,2) = 7.63034d0;   Afac(5,3) = -14.35440d0; Afac(5,4) = 31.52670d0

    Afac(6,1) = -1.89992d0;    Afac(6,2) = -12.5367d0;  Afac(6,3) = 4.98529d0; Afac(6,4) = -18.15070d0

    Afac(7,1) = 24.27450d0;    Afac(7,2) = 3.44945d0;   Afac(7,3) = -11.29130d0; Afac(7,4) = 69.34660d0

    Afac(8,1) = 0.79716d0;     Afac(8,2) = 1.11764d0;   Afac(8,3) = 0.12348e-1; Afac(8,4) = -4.11661d0 

    Afac(9,1) = -0.23816d0;    Afac(9,2) = 0.67695e-1;  Afac(9,3) = -0.81630d0; Afac(9,4) = 4.02528d0

    Afac(10,1) = 0.68629e-1;   Afac(10,2) = 0.34793d0;  Afac(10,3) = 0.59256d0; Afac(10,4) = -0.72663d0

    A_cst = 1.16145d0;  B_cst = 0.14874d0
    C_cst = 0.52487d0;  D_cst = 0.77320d0
    E_cst = 2.16178d0;  F_cst = 2.43787d0
    G_cst = -6.435e-4;  H_cst = 7.27371d0
    S_cst = 18.0323d0;  W_cst = -0.76830d0

    Bfac(1,1) = 2.41657d0;   Bfac(1,2) = 0.74824d0;   Bfac(1,3) = -0.91858d0;  Bfac(1,4) = 121.721d0
    Bfac(2,1) = -0.50924d0;  Bfac(2,2) = -1.50936d0;  Bfac(2,3) = -49.99120d0; Bfac(2,4) = 69.9834d0
    Bfac(3,1) = 6.61069d0;   Bfac(3,2) = 5.62073d0;   Bfac(3,3) = 64.75990d0;  Bfac(3,4) = 27.0389d0
    Bfac(4,1) = 14.54250d0;  Bfac(4,2) = -8.91387d0;  Bfac(4,3) = -5.63794d0;  Bfac(4,4) = 74.3435d0
    Bfac(5,1) = 0.79274d0;   Bfac(5,2) = 0.82019d0;   Bfac(5,3) = -0.69369d0;  Bfac(5,4) = 6.31734d0
    Bfac(6,1) = -5.86340d0;  Bfac(6,2) = 12.80050d0;  Bfac(6,3) = 9.58926d0;   Bfac(6,4) = -65.5292d0
    Bfac(7,1) = 81.17100d0;  Bfac(7,2) = 114.15800d0; Bfac(7,3) = -60.84100d0; Bfac(7,4) = 466.7750d0

 
    ! Compute Sigma_ij
    do j = 1,nspec
       do i = 1,nspec
!!$         Sigmaij(i,j) = sqrt(sig(i)*sig(j))*1e-8  ! converted into cm
          Sigmaij(i,j) = 0.5d0*(sig(i)+sig(j))*1e-8  ! converted into cm
       end do
    end do

    Kappai(:) = 0.0d0
    iH2O = 3
    Kappai(iH2O) = 0.076d0
    
    
  end subroutine actual_transport_init


  subroutine actual_transport_close

    implicit none

    if( allocated(wt)) deallocate(wt)
    if( allocated(iwt)) deallocate(iwt)
    if( allocated(eps)) deallocate(eps)
    if( allocated(sig)) deallocate(sig)
    if( allocated(dip)) deallocate(dip)
    if( allocated(pol)) deallocate(pol)
    if( allocated(zrot)) deallocate(zrot)
    if( allocated(nlin)) deallocate(nlin)

    deallocate(fitmu)
    deallocate(fitlam)
    deallocate(fitdbin)

    smp_initialized = .false.

  end subroutine actual_transport_close


  subroutine build_internal(npts)
    integer, intent(in) :: npts

    if (npts_smp .ne. npts .and. npts.gt.0) then
       if (npts_smp .ne. 0) then
          call destroy_internal()
       endif
       allocate(Tloc(npts))
       allocate(rholoc(npts))
       allocate(Yloc(npts,nspec))
       allocate(Xloc(npts,nspec))
       allocate(logT(npts,norder))
       allocate(wbar(npts))

       allocate(muloc(npts,nspec) )
       allocate(lamloc(npts,nspec) )
       allocate(xiloc(npts,nspec) )
       allocate(dbinloc(npts,nspec,nspec) )
       npts_smp = npts
    endif

  end subroutine build_internal


  subroutine destroy_internal

    deallocate(Tloc)
    deallocate(rholoc)
    deallocate(Yloc)
    deallocate(Xloc)
    deallocate(wbar)
    deallocate(logT)
    deallocate(muloc)
    deallocate(lamloc)
    deallocate(xiloc)
    deallocate(dbinloc)

    npts_smp = 0

  end subroutine destroy_internal



  subroutine actual_transport(which, coeff)

    use amrex_error_module

    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff
    integer :: i, n, nn

    if (.not. smp_initialized) then
       call amrex_error('simple ::actual_transport called before initialized')
    endif

    if (npts_smp .ne. coeff%npts) then
       call build_internal(coeff%npts)
    endif

    do i=1,coeff%npts

       Yloc(i,:)  = coeff % eos_state(i) % massfrac(:)    
       Tloc(i)    = coeff % eos_state(i) % T
       rholoc(i)  = coeff % eos_state(i) % rho
       logT(i,1) = log(Tloc(i))
       logT(i,2) = logT(i,1)**2
       logT(i,3) = logT(i,1)*logT(i,2)

    enddo

    ! define trace mole and mass fraction to get correct limiting
    ! behavior of species diffusion as species disappear

    do i = 1,coeff%npts

       !  just to not waste space
       Xloc(i,1) = sum(Yloc(i,:))
       wbar(i) = 0.d0

    enddo

    !  do i = 1,coeff%npts
    !       write(6,*)"yloc", Xloc(i,1),Yloc(i,:)
    !  enddo

    do n = 1, nspec
       do i = 1,coeff%npts

          Yloc(i,n) = Yloc(i,n) + trace*(Xloc(i,1)/dble(nspec)-Yloc(i,n)) 

       enddo
    enddo

    do n = 1, nspec
       do i = 1,coeff%npts

          wbar(i) = wbar(i) + Yloc(i,n)*iwt(n)

       enddo
    enddo

    do i = 1,coeff%npts

       wbar(i) = 1.d0/ wbar(i)

    enddo

    !  do i = 1,coeff%npts
    !       write(6,*)"new yloc", sum(Yloc(i,:)),Yloc(i,:)
    !  enddo

    !  write(6,*)"wt ",wt
    !  write(6,*)"iwt",wt
    !  write(6,*)"wbar",wbar


    do n = 1, nspec
       do i = 1,coeff%npts

          Xloc(i,n) = Yloc(i,n)*wbar(i)*iwt(n)

       enddo
    enddo

    !  do i = 1,coeff%npts
    !       write(6,*)"xloc",sum(Xloc(i,:)), Xloc(i,:)
    !  enddo
    !  do i = 1,coeff%npts
    !       write(6,*)"temp",Tloc(i),logT(i,1)
    !  enddo



    if (which % wtr_get_mu) then
       do n=1,nspec
          do i=1,coeff%npts

             muloc(i,n) = fitmu(1,n)+fitmu(2,n)*logT(i,1)+ fitmu(3,n)*logT(i,2)  &
                  + fitmu(4,n)*logT(i,3)
             muloc(i,n) = exp(muloc(i,n))

          enddo
       enddo

       coeff % mu(:) = 0.d0

       do n=1,nspec
          do i=1,coeff%npts

             coeff%mu(i) = coeff%mu(i)+ Xloc(i,n)*muloc(i,n)**6.d0

          enddo
       enddo


       do i=1,coeff%npts
          coeff%mu(i) = coeff%mu(i)**(1.d0/6.d0)
!         write(6,*)" mu",  coeff % eos_state(i) % T,coeff % mu(i)
       enddo


       !  assumption that we only get bulk viscosity is we are already getting shear viscosity
       if (which % wtr_get_xi) then

          call comp_pure_bulk(coeff, coeff%npts)

          coeff % xi(:) = 0.d0

          do n=1,nspec
             do i=1,coeff%npts

                coeff%xi(i) = coeff%xi(i)+ Xloc(i,n)*xiloc(i,n)**0.75d0

             enddo
          enddo

          do i=1,coeff%npts
             coeff%xi(i) = coeff%xi(i)**(4.d0/3.d0)
          enddo


       endif

    endif

    if (which % wtr_get_lam) then

       do n=1,nspec
          do i=1,coeff%npts

             lamloc(i,n) = fitlam(1,n)+fitlam(2,n)*logT(i,1)+ fitlam(3,n)*logT(i,2)  &
                  + fitlam(4,n)*logT(i,3)
             lamloc(i,n) = exp(lamloc(i,n))

          enddo
       enddo

       coeff % lam(:) = 0.d0

       do n=1,nspec
          do i=1,coeff%npts

             coeff%lam(i) = coeff%lam(i)+ Xloc(i,n)*lamloc(i,n)**0.25d0

          enddo
       enddo


       do i=1,coeff%npts
          coeff%lam(i) = coeff%lam(i)**4

   !      write(6,*)" lam",  coeff % eos_state(i) % T,coeff % lam(i)
       enddo

    endif

    call nonIdeal_Chung(coeff,which)
 !    do i=1,coeff%npts
 !        write(6,*)"  mu after",  coeff % eos_state(i) % T,coeff % mu(i)
 !        write(6,*)" lam after",  coeff % eos_state(i) % T,coeff % lam(i)
 !    enddo


    if (which % wtr_get_Ddiag .or. which % wtr_get_Dmat) then

       do n=1,nspec
          do nn=1,n-1
             do i=1,coeff%npts

                dbinloc(i,n,nn) = fitdbin(1,n,nn)+fitdbin(2,n,nn)*logT(i,1)   &
                     + fitdbin(3,n,nn)*logT(i,2)+ fitdbin(4,n,nn)*logT(i,3)
                dbinloc(i,n,nn) = exp(dbinloc(i,n,nn))

                dbinloc(i,nn,n) = dbinloc(i,n,nn)

             enddo
          enddo
          do i=1,coeff%npts
             dbinloc(i,n,n) = 0.d0
          enddo
       enddo

       call nonIdeal_Binary_Diff(coeff)

       if (which % wtr_get_Ddiag ) then

          call mixture(coeff,coeff%npts)

          do n=1,nspec
             do i=1,coeff%npts

                coeff%Ddiag(i,n) = rholoc(i)*coeff%Ddiag(i,n)

             enddo
          enddo

       elseif (.false.) then

          call matrix(coeff, coeff%npts)

          do n=1,nspec
             do nn=1,nspec
                do i=1,coeff%npts

                   coeff%Dmat(i,nn,n) = rholoc(i)*coeff%Dmat(i,nn,n)

                enddo
             enddo
          enddo
       

       endif
    endif

  end subroutine actual_transport


  subroutine comp_pure_bulk(coeff, npts)

    implicit none

    type (trv_t), intent(inout) :: coeff
    integer :: npts

    real(amrex_real) :: cvk(npts,nspec), cvkint(npts,nspec), cvkrot(npts,nspec)
    real(amrex_real) :: rwrk 
    real(amrex_real) :: FofT(npts,nspec), Fnorm(nspec), epskoverT

    real(amrex_real), parameter :: pi = 3.141592653589793238d0


    integer n,i,j,k,iwrk

    do n=1,npts
       call ckcvms( coeff % eos_state(n) % T, iwrk, rwrk, coeff % eos_state(n)%cvi )
    enddo

    do i=1,nspec
       do n=1,npts

          if(nlin(i) .eq.0)then

             cvkint(n,i) = 0.d0
             cvkrot(n,i) = 0.d0

          elseif(nlin(i).eq.1)then

             cvkint(n,i) =( coeff% eos_state(n)% cvi(i) ) * wt(i)/Ru - 1.5d0
             cvkrot(n,i) = 1.d0

          else

             cvkint(n,i) =( coeff% eos_state(n)% cvi(i) ) * wt(i)/Ru - 1.5d0
             cvkrot(n,i) = 1.5d0

          endif

       enddo
    enddo

    do i = 1,nspec

       epskoverT = eps(i)/298.d0

       Fnorm(i) = 1.d0 + 0.5d0*pi**1.5*sqrt(epskoverT) + (2.d0+.5d0*pi**2)*epskoverT &
            +(pi*epskoverT)**1.5

    enddo

    do i = 1,nspec

       do n=1,npts

          epskoverT = eps(i)/ Tloc(n)

          FofT(n,i) = 1.d0 + 0.5d0*pi**1.5*sqrt(epskoverT) + (2.d0+.5d0*pi**2)*epskoverT &
               +(pi*epskoverT)**1.5

       enddo

    enddo

    do i=1,nspec

       if(nlin(i) .ne. 0)then

          do n=1,npts

             !   zrot/crot approximately zint / cint by assuming vibrational internal energy is small
             !   cvkrot is scaled by wk / Ru = mk / kb relative to standard specific cv

             xiloc(n,i) = 0.25d0*pi*(cvkint(n,i)/(cvkint(n,i)+1.5d0))**2* zrot(i)/cvkrot(n,i)*  &
                  Fnorm(i)/FofT(n,i) * muloc(n,i)

          enddo

       else

          !     no bulk viscosity for monotmic species

          do n=1,npts

             xiloc(n,i) = 0.d0

          enddo

       endif

    enddo


  end subroutine comp_pure_bulk

  subroutine mixture(coeff,npts)

    implicit none

    type (trv_t), intent(inout) :: coeff

    integer :: npts

    real(amrex_real) :: term1(npts),term2(npts)
    integer i,j,k,n


    do j = 1, nspec
       term1 = 0.d0
       term2 = 0.d0
       do k = 1, nspec
          if(k.ne.j) then
             do i = 1,npts

                term1(i) = term1(i) + Yloc(i,k)
                term2(i) = term2(i) + Xloc(i,k)/dbinloc(i,k,j)

             enddo
          endif
       enddo

       do i=1,npts
          coeff%Ddiag(i,j) = wt(j)* term1(i)/term2(i) / wbar(i)
       enddo

    enddo

  end subroutine mixture

  subroutine matrix(coeff,npts)

    implicit none

    type (trv_t), intent(inout) :: coeff

    integer :: npts

    integer ::  i, j, k, jj, n
    real(kind=8) ::  term1(npts), term2(npts)
    real(kind=8) :: D_tilde(1:npts,1:nspec,1:nspec), Di(1:npts,1:nspec), Diff_ij(1:npts,1:nspec,1:nspec)
    real(kind=8) :: Deltamat(1:npts,1:nspec,1:nspec), Zmat(1:npts,1:nspec,1:nspec)
    real(kind=8), dimension(1:npts,1:nspec,1:nspec) :: Pmat, Jmat
    real(kind=8), dimension(1:npts,1:nspec) :: Minv, Mmat
    real(kind=8), dimension(1:npts,1:nspec,1:nspec) :: PJ, matrix1, matrix2
    real(kind=8) :: scr(npts)


    ! Find Di matrix 
    do i = 1, nspec
       term1 = 0.0d0  
       term2 = 0.0d0  
       do j = 1, nspec
          if(j.ne.i) then
             do n=1,npts
                term1(n) = term1(n) + Yloc(n,j)
                term2(n) = term2(n) + Xloc(n,j)/dbinloc(n,i,j)
             enddo
          endif
       enddo
       do n=1,npts
          Di(n,i) = term1(n)/term2(n) 
       enddo
    enddo


    ! Compute Mmat and Minv
    do i = 1, nspec
       do n=1,npts

          Mmat(n,i) = Xloc(n,i)/Di(n,i)
          Minv(n,i) = Di(n,i)/Xloc(n,i)

       enddo
    enddo


    ! Compute P matrix
    Pmat = 0.0d0
    do i = 1, nspec
       do j = 1, nspec
          do n=1,npts
             Pmat(n,i,j) = - Yloc(n,j)
          enddo
          if(i.eq.j) then
             do n=1,npts
                Pmat(n,i,j) =  Pmat(n,i,j) + 1.0d0
             enddo
          endif
       enddo
    enddo


    ! Compute Deltamat
    Deltamat = 0.0d0
    do i = 1, nspec
       do j = 1, nspec
          if(i.eq.j) then
             term1 = 0.0d0
             do k = 1, nspec
                if(k.ne.i) then

                   do n=1,npts

                      term1(n) = term1(n) + Xloc(n,i)*Xloc(n,k)/dbinloc(n,i,k)

                   enddo

                endif
             enddo

             do n=1,npts
                Deltamat(n,i,i) = term1(n)
             enddo

          else
             do n=1,npts
                Deltamat(n,i,j) = -Xloc(n,i)*Xloc(n,j)/dbinloc(n,i,j)
             enddo
          endif
          do n=1,npts
             Zmat(n,i,j) = -Deltamat(n,i,j)
          enddo
       enddo
    enddo


    ! Compute Zmat
    do i = 1, nspec
       do n=1,npts
          Zmat(n,i,i) = Zmat(n,i,i) + Mmat(n,i)
       enddo
    enddo

    ! Compute Jmat
    Jmat = 0.d0
    do i = 1, nspec
       do j = 1, nspec
          do n=1,npts
             Jmat(n,i,j) = Minv(n,i)*Zmat(n,i,j)
          enddo
       enddo
    enddo

    ! Compute PJ
    PJ = 0.0d0
    do i = 1, nspec
       do j = 1, nspec
          do k = 1, nspec
             do n=1,npts
                PJ(n,i,j) = PJ(n,i,j) + Pmat(n,i,k)*Jmat(n,k,j)
             enddo
          enddo
       enddo
    enddo



    ! Compute P M^-1 Pt; store it in matrix2
    do i = 1, nspec
       do j = 1, nspec
          scr = 0.d0
          do k = 1, nspec
             do n=1,npts
                scr(n) = scr(n) + Pmat(n,i,k)*Minv(n,k)*Pmat(n,j,k)
             enddo
             ! notice the change in indices for Pmat to represent Pmat^t
          enddo
          do n=1,npts
             matrix2(n,i,j) = scr(n)
             Diff_ij(n,i,j) = scr(n)
          enddo
       enddo
    enddo

    if(iter.gt.0)then

       do jj = 1,iter


          !         matrix1=0
          do i = 1, nspec
             do j = 1, nspec
                scr = 0.d0
                do k = 1, nspec
                   do n=1,npts
                      scr(n) = scr(n) + PJ(n,i,k)*Diff_ij(n,k,j)
                   enddo
                enddo
                do n=1,npts
                   matrix1(n,i,j) = scr(n)+matrix2(n,i,j)
                enddo
             enddo
          enddo

          Diff_ij=matrix1

       enddo

    endif



    ! Compute D_tilde
    do i = 1, nspec
       do j = 1, nspec
          do n=1,npts
             coeff%Dmat(n,i,j) = Diff_ij(n,i,j)*Yloc(n,i)
          enddo
       enddo
    enddo

  end subroutine  matrix
  !==================================================!
  ! Compute non-Ideal viscosity and conductivity     !
  !         using Chung's method !                   !
  !==================================================!
  subroutine nonIdeal_Chung(coeff, which)
    implicit none
    type (wtr_t), intent(in   ) :: which
    type (trv_t), intent(inout) :: coeff
    real(amrex_real) :: mu
    integer ::  i,j,k
    real(amrex_real) :: sigma_M,sigma_M_3, Tstar, Epsilon_M
    real(amrex_real) :: Omega_M, DP_red_4, DP_m_4
    real(amrex_real) :: MW_m, Vcm, Tcm
    real(amrex_real) :: sigma_ij,epsilon_ij, MW_ij,Omega_ij
    real(amrex_real) :: T1, T2, T3,InvSigma1,InvSigma2,InvSigma3
    real(amrex_real) :: KappaM
    real(amrex_real) :: sumMWij,lambda_p,beta
    real(amrex_real) :: y,G1,G2,eta_P, H2
    real(amrex_real) :: A(10), B(7)
    real(amrex_real), parameter :: dpFactor = 297.2069113e6


    real(amrex_real)  :: scr

    do k = 1,coeff%npts

       sigma_M_3   = 0.0d0
       Epsilon_M   = 0.0d0
       Tstar       = 0.0d0
       Omega_M     = 0.0d0
       MW_m        = 0.0d0
       DP_m_4      = 0.0d0
       KappaM      = 0.0d0

       do j = 1,nspec
          do i = 1,nspec

             T1 = sqrt(sig(i)*sig(j))
             T2 = sqrt(sig(i)*sig(j))*T1
             T3 = T2 * T1

             sigma_M_3 = sigma_M_3 + Xloc(k,i) * Xloc(k,j)*T3

             epsilon_ij = sqrt(eps(i)*eps(j))

             Epsilon_M = Epsilon_M + Xloc(k,i) * Xloc(k,j)*epsilon_ij*T3

             Omega_ij = 0.5d0*omega(i) + 0.5d0*omega(j)

             Omega_M = Omega_M + Xloc(k,i) * Xloc(k,j)*Omega_ij*T3

             MW_ij = 2.0d0/(iwt(i)+iwt(j))

        !    sumMWij = Xloc(k,i) * Xloc(k,j)* epsilon_ij * T2 *sqrt(MW_ij)

        !    MW_m = MW_m  + sumMWij*sumMWij
             sumMWij = Xloc(k,i) * Xloc(k,j)* epsilon_ij * T2 *sqrt(MW_ij)

             MW_m = MW_m  + sumMWij

             DP_m_4 = DP_m_4 + Xloc(k,i) * Xloc(k,j)*dip(i)*dip(i)*dip(j)*dip(j)/(T3*epsilon_ij)


             KappaM = KappaM +  Xloc(k,i) * Xloc(k,j) * sqrt(Kappai(i)*Kappai(j))

          end do
       end do


       Mw_m = MW_m **2

       sigma_M = sigma_M_3**(1.0d0/3.0d0)

       InvSigma1 = 1.0d0/sigma_M
       InvSigma2 = InvSigma1*1.0d0/sigma_M
       InvSigma3 = 1.0d0/sigma_M_3

       Epsilon_M = Epsilon_M*InvSigma3

       Tstar = coeff%eos_state(k)%T/Epsilon_M  ! non-dimensional

       Tcm = 1.2593d0*Epsilon_M              ! K

       Vcm = 1.8887*sigma_M_3                ! cm3/mol

       Omega_M = Omega_M*InvSigma3          

       MW_m = MW_m*InvSigma3*InvSigma1/(Epsilon_M*Epsilon_M) ! g/mol

       DP_m_4 = (DP_m_4* sigma_M_3* Epsilon_M)

       DP_red_4 = dpFactor*DP_m_4/(Vcm*Tcm*Vcm*Tcm)

       y = Vcm * coeff%eos_state(k)%rho/(6.0d0 * wbar(k))

       do i = 1, 10
          A(i) = Afac(i,1) + Afac(i,2)*Omega_M + Afac(i,3)*DP_red_4 + Afac(i,4)*KappaM

       end do

       G1 = (1.0d0-0.5d0*y)/((1.0-y)*(1.0-y)*(1.0-y))
       
       G2 = (A(1)*(1.0d0-exp(-A(4)*y))/y + A(2)*G1*exp(A(5)*y) + A(3)*G1)/(A(1)*A(4)+A(2)+A(3))

       eta_P = (36.344e-6*sqrt(MW_m*Tcm)/(Vcm**(2.0d0/3.0d0)))*A(7)*y*y*G2*exp(A(8) + A(9)/Tstar + A(10)/(Tstar*Tstar))

       ! Viscosity
!      if (which % wtr_get_mu) then

!          scr = coeff%mu(k)

!          coeff%mu(k) = coeff%mu(k)*(1.0d0/G2 + A(6)*y) + eta_P   ! in CGS units

!          if(coeff%mu(k).lt.0)then
!                 write(6,*)"kappa", kappaM, kappai(3),xloc(k,3),xloc(k,3)**2*kappai(3)

!                  write(6,*)" bad state Rho, T", coeff%eos_state(k)%rho,coeff%eos_state(k)%T
!                  write(6,*)" bad state Rho,moleface ", Xloc(k,:)

!                  write(6,*)" in chung", scr,coeff%mu(k),G2,A(6),Y,eta_P
!                  write(6,*)" MW_m, Tcm. Vcm", MW_m, Tcm, Vcm
!                  write(6,*)" A(1), A(2),A(3), A(4), G1",A(1), A(2),A(3), A(4), G1
!                  write(6,*)"A(7), A(8),A(9), A(10), Tstar",A(7), A(8),A(9), A(10), Tstar
!                  write(6,*)"MW_m,Tcm,vcm.Dm_m_4, epsm, omegm,sigm",MW_m,Tcm,Vcm, Dp_red_4,Epsilon_M,Omega_M,Sigma_M
!          endif

!      endif

       do i = 1, 7
          B(i) = Bfac(i,1) + Bfac(i,2)*Omega_M + Bfac(i,3)*DP_red_4 + Bfac(i,4)*KappaM
       end do

       H2 = (B(1)*(1.0d0-exp(-B(4)*y))/y + B(2)*G1*exp(B(5)*y) + B(3)*G1)/(B(1)*B(4)+B(2)+B(3))

       lambda_p = 3.039e-4*sqrt(Tcm/MW_m)/(Vcm**(2.0d0/3.0d0))*B(7)*y*y*H2*sqrt(Tstar)  ! (cal/cm s K)

       lambda_p = lambda_p* 4.184e+7   ! erg/(cm s K)

       beta = (1.0d0/H2)+B(6)*y  ! Non-dimensional 

       ! Thermal conductivity
       if (which % wtr_get_lam) then

           coeff%lam(k) = coeff%lam(k)*beta + lambda_p  ! erg/(cm s K)
 
       endif

    end do

  end subroutine nonIdeal_Chung
  !==================================================!
  ! Compute non-Ideal Binary diffusion coefficient   !
  !==================================================!
  subroutine nonIdeal_Binary_Diff(coeff)
    implicit none
    type (trv_t), intent(inout) :: coeff
    integer :: i,j,k, n
    real(amrex_real) :: Upsilon(nspec,nspec)
    real(amrex_real), parameter :: Pst = 1013250.0d0


    do n = 1,coeff%npts


       do j = 1, nspec
          do i = 1, nspec

             Upsilon(i,j) = 0.0d0 
             do k = 1, nspec
                Upsilon(i,j)  = Upsilon(i,j) + iwt(k)*coeff%eos_state(n)%massfrac(k)* ( &
                     8.0d0*(Sigmaij(i,k)**3.0d0 + Sigmaij(j,k)**3.0d0)    &
                     -6.0d0*(Sigmaij(i,k)*Sigmaij(i,k) + Sigmaij(j,k)*Sigmaij(j,k))*Sigmaij(i,j)    &
                     -3.0d0*((Sigmaij(i,k)*Sigmaij(i,k) - Sigmaij(j,k)*Sigmaij(j,k))**2.0d0)/Sigmaij(i,j) &
                     + Sigmaij(i,j)*Sigmaij(i,j)*Sigmaij(i,j) )
             end do
             Upsilon(i,j)  = Upsilon(i,j)*coeff%eos_state(n)%rho*Avna*PI/12.0d0
             Upsilon(i,j)  = Upsilon(i,j) + 1.0d0

             dbinloc(n,i,j) = dbinloc(n,i,j)*Pst*wbar(n)/(Ru*coeff%eos_state(n)%T*coeff%eos_state(n)%rho)/Upsilon(i,j)
             
          end do
       end do
       
    end do


  end subroutine nonIdeal_Binary_Diff

end module actual_transport_module

