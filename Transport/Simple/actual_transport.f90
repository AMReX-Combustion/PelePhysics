module actual_transport_module

  use amrex_fort_module, only : amrex_real
  use eos_type_module
  use transport_type_module
  use chemistry_module, only : Ru

  implicit none

  character (len=64) :: transport_name = "smp"
  logical, save, private :: smp_initialized = .false.
  integer, save, private :: npts_smp = 0
  real(amrex_real), save, allocatable :: Tloc(:), Yloc(:,:), Xloc(:,:), rholoc(:), wbar(:), pscale(:)
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
  real(amrex_real), parameter ::  Patm = 1.01325d6

  logical, parameter :: use_bulk_viscosity = .true.

  !$omp threadprivate(npts_smp,Tloc,Yloc,Xloc,rholoc,wbar,pscale,muloc,lamloc,dbinloc,xiloc,logT)

contains

  ! This subroutine should be called outside OMP PARALLEL
  subroutine actual_transport_init

    implicit none
 
    call egtransetKK(nspec)
    call egtransetNO(nfit)

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

    smp_initialized = .true.

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
       allocate(pscale(npts))

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
    deallocate(pscale)
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

   do n = 1, nspec
   do i = 1,coeff%npts
       
     Xloc(i,n) = Yloc(i,n)*wbar(i)*iwt(n)

   enddo
   enddo

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
         enddo


    endif


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

        pscale = 0.d0
        do i=1,coeff%npts

            pscale(i) = Patm * wbar(i)/ (Ru * coeff % eos_state(i) % T * coeff % eos_state(i) % rho)

        enddo

       if (which % wtr_get_Ddiag ) then

           call mixture(coeff,coeff%npts)

           do n=1,nspec
           do i=1,coeff%npts

                 coeff%Ddiag(i,n) = rholoc(i)*pscale(i)*coeff%Ddiag(i,n)

           enddo
           enddo

       elseif (.false.) then

           call matrix(coeff, coeff%npts)

           do n=1,nspec
           do nn=1,nspec
           do i=1,coeff%npts

                 coeff%Dmat(i,nn,n) = rholoc(i)*pscale(i)*coeff%Dmat(i,nn,n)

           enddo
           enddo
           enddo

       endif
    endif

  end subroutine actual_transport


  subroutine comp_pure_bulk(coeff, npts)

  implicit none

  type (trv_t), intent(inout) :: coeff
  integer npts

  real(amrex_real) :: cvk(npts,nspec), cvkint(npts,nspec), cvkrot(npts,nspec)
  real(amrex_real) :: FofT(npts,nspec), Fnorm(nspec), epskoverT

  real(amrex_real), parameter :: pi = 3.141592653589793238d0
 

  integer n,i

  do n=1,npts
      call ckcvms( coeff % eos_state(n) % T, coeff % eos_state(n)%cvi )
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
  integer i,j,k

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
  real(kind=8) :: Di(1:npts,1:nspec), Diff_ij(1:npts,1:nspec,1:nspec)
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
end module actual_transport_module

