module fuego_chemistry

    use amrex_fort_module,    only : amrex_real
    use extern_probin_module, only : mwt_scalar

#include "mechanism.h"

    implicit none

    integer, parameter :: nelements   = NUM_ELEMENTS  ! number of elements
    integer, parameter :: nspecies    = NUM_SPECIES   ! number of species
    integer, parameter :: nreactions  = NUM_REACTIONS ! number of reactions

    integer :: naux 
    character (len=16), save, allocatable :: aux_names(:)

    logical, save :: chemistry_initialized = .false.
    logical, save :: network_initialized = .false.

    integer, parameter :: L_elem_name = 3  ! Each element name has at most 3 characters
    character*(L_elem_name), save :: elem_names(NUM_ELEMENTS)

    integer, parameter :: L_spec_name = 16 ! Each species name has at most 8 characters
    character*(L_spec_name), save :: spec_names(NUM_SPECIES)

    real(amrex_real), save :: molecular_weight(NUM_SPECIES), inv_mwt(NUM_SPECIES)

    real(amrex_real), save :: Ru, Ruc, Patm

    interface
        subroutine egtransetLENIMC(LENIMC) bind(c,name='egtransetLENIMC')
            integer, intent(inout) :: LENIMC
        end subroutine

        subroutine egtransetLENRMC(LENRMC) bind(c,name='egtransetLENRMC')
            integer, intent(inout) :: LENRMC
        end subroutine

        subroutine egtransetNO(NO) bind(c,name='egtransetNO')
            integer, intent(inout) :: NO
        end subroutine

        subroutine egtransetKK(KK) bind(c,name='egtransetKK')
            integer, intent(inout) :: KK
        end subroutine

        subroutine egtransetNLITE(NLITE) bind(c,name='egtransetNLITE')
            integer, intent(inout) :: NLITE
        end subroutine

        subroutine egtransetPATM(PATM) bind(c,name='egtransetPATM')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: PATM
        end subroutine

        subroutine egtransetWT(WT) bind(c,name='egtransetWT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: WT(*)
        end subroutine

        subroutine egtransetEPS(EPS) bind(c,name='egtransetEPS')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: EPS(*)
        end subroutine

        subroutine egtransetSIG(SIG) bind(c,name='egtransetSIG')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: SIG(*)
        end subroutine

        subroutine egtransetDIP(DIP) bind(c,name='egtransetDIP')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: DIP(*)
        end subroutine

        subroutine egtransetPOL(POL) bind(c,name='egtransetPOL')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) ::POL(*)
        end subroutine

        subroutine egtransetZROT(ZROT) bind(c,name='egtransetZROT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: ZROT(*)
        end subroutine

        subroutine egtransetNLIN(NLIN) bind(c,name='egtransetNLIN')
            use amrex_fort_module, only : amrex_real
            integer, intent(inout) :: NLIN(*)
        end subroutine

        subroutine egtransetCOFETA(COFETA) bind(c,name='egtransetCOFETA')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFETA(*)
        end subroutine

        subroutine egtransetCOFLAM(COFLAM) bind(c,name='egtransetCOFLAM')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFLAM(*)
        end subroutine

        subroutine egtransetCOFD(COFD) bind(c,name='egtransetCOFD')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFD(*)
        end subroutine

        subroutine egtransetKTDIF(KTDIF) bind(c,name='egtransetKTDIF')
            use amrex_fort_module, only : amrex_real
            integer, intent(inout) :: KTDIF(*)
        end subroutine

        subroutine get_t_given_eY(e,y,t,ierr) bind(c,name='GET_T_GIVEN_EY') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: e
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: t
            integer, intent(out )          :: ierr
        end subroutine

        subroutine get_t_given_hY(h,y,t,ierr) bind(c,name='GET_T_GIVEN_HY') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: h
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(inout ) :: t
            integer, intent(out )          :: ierr
        end subroutine

        subroutine ckxty(x,y) bind(c,name='CKXTY') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: x(*)
            real(amrex_real), intent(out ) :: y(*)
        end subroutine

        subroutine ckpy(rho,T,y,P) bind(c,name='CKPY') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: rho,T
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: P
        end subroutine

        subroutine ckrp(ru,ruc,pa) bind(c,name='CKRP') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout ) :: ru,ruc,pa
        end subroutine

        subroutine ckytx(y,x) bind(c,name='CKYTX')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: x(*)
        end subroutine

        subroutine ckcpms(T,cpms) bind(c,name='CKCPMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: T
            real(amrex_real), intent(out ) :: cpms(*)
        end subroutine

        subroutine ckcvms(T,cvms) bind(c,name='CKCVMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: T
            real(amrex_real), intent(out ) :: cvms(*)
        end subroutine

        subroutine ckhms(T,hms) bind(c,name='CKHMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: T
            real(amrex_real), intent(inout) :: hms(*)
        end subroutine

        subroutine cksms(T,sms) bind(c,name='CKSMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: T
            real(amrex_real), intent(inout) :: sms(*)
        end subroutine

        subroutine get_critparams(Tci,ai,bi,acentric_i) bind(c,name='GET_CRITPARAMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout   ) :: Tci(*)
            real(amrex_real), intent(inout   ) :: ai(*)
            real(amrex_real), intent(inout   ) :: bi(*)
            real(amrex_real), intent(inout   ) :: acentric_i(*)
        end subroutine

        subroutine ckhbms(T,y,hbms) bind(c,name='CKHBMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: T
            real(amrex_real), intent(in   ) :: y(*)
            real(amrex_real), intent(out)   :: hbms
        end subroutine

        subroutine ckums(T,ums) bind(c,name='CKUMS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: T
            real(amrex_real), intent(inout) :: ums(*)
        end subroutine

        subroutine ckcvbs(T,y,cvbs) bind(c,name='CKCVBS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: T
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: cvbs
        end subroutine

        subroutine ckcpbs(T,y,cpbs) bind(c,name='CKCPBS') 
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: T
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: cpbs
        end subroutine

        subroutine ckytcr(rho,T,y,c) bind(c,name='CKYTCR')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: rho,T
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: c(*)
        end subroutine

        subroutine ckytcp(P,T,y,c) bind(c,name='CKYTCP')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: P,T 
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: c(*)
        end subroutine

        subroutine ckrhoy(P,T,y,rho) bind(c,name='CKRHOY')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in  ) :: P,T
            real(amrex_real), intent(in  ) :: y(*)
            real(amrex_real), intent(out ) :: rho
        end subroutine

        subroutine vckhms(np,T,hms) bind(c,name='VCKHMS') 
            use amrex_fort_module, only : amrex_real
            integer, intent(in  )           :: np
            real(amrex_real), intent(in   ) :: T(*)
            real(amrex_real), intent(inout) :: hms(*)
        end subroutine

        subroutine vckytx(np,y,x) bind(c,name='VCKYTX') 
            use amrex_fort_module, only : amrex_real
            integer, intent(in  )           :: np
            real(amrex_real), intent(in   ) :: y(*)
            real(amrex_real), intent(out  ) :: x(*)
        end subroutine

        subroutine cksyme(kname,lenkname) bind(c,name='CKSYME') 
            use amrex_fort_module, only : amrex_real
            integer, intent(inout )         :: kname(*)
            integer, intent(in    )         :: lenkname
        end subroutine

        subroutine cksyms(kname,lenkname) bind(c,name='CKSYMS') 
            use amrex_fort_module, only : amrex_real
            integer, intent(inout )         :: kname(*)
            integer, intent(in    )         :: lenkname
        end subroutine

        subroutine ckinit() bind(c,name='CKINIT') 
        end subroutine

        subroutine ckfinalize() bind(c,name='CKFINALIZE') 
        end subroutine

        subroutine ckindx(mm,kk,ii,nfit) bind(c,name='CKINDX')
            integer, intent(inout )         :: mm,kk,ii,nfit
        end subroutine

        subroutine ckwt(wt) bind(c,name='CKWT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: wt(*)
        end subroutine

        subroutine ckawt(awt) bind(c,name='CKAWT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: awt(*)
        end subroutine

        subroutine ckwc(T,c,wdot) bind(c,name='CKWC')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: T
            real(amrex_real), intent(in   ) :: c(*)
            real(amrex_real), intent(inout) :: wdot(*)
        end subroutine

        subroutine ckmmwy(y,wtm) bind(c,name='CKMMWY')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(in   ) :: y(*)
            real(amrex_real), intent(inout) :: wtm
        end subroutine

        subroutine vckwyr(np,rho,t,y,wdot) bind(c,name='VCKWYR')
            use amrex_fort_module, only : amrex_real
            integer, intent(in   ) :: np
            real(amrex_real), intent(inout) :: rho(*)
            real(amrex_real), intent(in) :: t(*)
            real(amrex_real), intent(inout) :: y(*)
            real(amrex_real), intent(inout) :: wdot(*)
        end subroutine

        subroutine cknu(kdim,nuki) bind(c,name='CKNU')
            integer, intent(in   ) :: kdim
            integer, intent(inout) :: nuki(*)
        end subroutine

        subroutine ckinu(i,nspec,ki,nu) bind(c,name='CKINU')
            integer, intent(in   ) :: i
            integer, intent(inout) :: nspec
            integer, intent(inout) :: ki(*)
            integer, intent(inout) :: nu(*)
        end subroutine

        subroutine ckncf(mdim,ncf) bind(c,name='CKNCF')
            integer, intent(in   ) :: mdim
            integer, intent(inout) :: ncf(*)
        end subroutine

        subroutine GET_REACTION_MAP(rmap) bind(c,name='GET_REACTION_MAP')
            integer, intent(inout) :: rmap(*)
        end subroutine

        subroutine get_imw(imw) bind(c,name='get_imw')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: imw(*)
        end subroutine

    end interface

contains

    subroutine chemistry_init()

      integer :: i, ic, ii
      integer :: names(NUM_SPECIES*L_spec_name)

      call ckinit()

      if (NUM_SPECIES > 1) then
          call cksyme(names, L_elem_name) 

          ic = 1
          do i = 1, NUM_ELEMENTS
             do ii=1, L_elem_name
                elem_names(i)(ii:ii) = char(names(ic))
                ic = ic + 1
             end do
          end do

          call cksyms(names, L_spec_name) 

          ic = 1
          do i = 1, NUM_SPECIES
             do ii=1, L_spec_name
                spec_names(i)(ii:ii) = char(names(ic))
                ic = ic+1
             end do
          end do

          call ckwt(molecular_weight)
      else
          elem_names(1) = "X"
          spec_names(1) = "X"
          molecular_weight(1) = mwt_scalar
      end if

      inv_mwt = 1.d0 / molecular_weight

      call ckrp(Ru, Ruc, Patm)

      chemistry_initialized = .true.

    end subroutine chemistry_init


    subroutine chemistry_close()
      call ckfinalize()
      chemistry_initialized = .false.
    end subroutine chemistry_close


    function get_species_index(name) result (iname)
      character(len=*), intent(in) :: name
      integer :: iname
      integer :: i
      iname = -1
      do i = 1, NUM_SPECIES
         if (trim(spec_names(i)) .eq. trim(name)) then
            iname = i
            exit
         end if
      end do
    end function get_species_index


    subroutine network_init

      use extern_probin_module, only: numaux, auxnamesin 
      use amrex_error_module

      implicit none

      integer :: iaux,nnames

      if (.not. chemistry_initialized)  call chemistry_init() 

      ! Get auxiliary variables (if any) 
      naux = numaux
      allocate(aux_names(naux))
      nnames = count(transfer(auxnamesin, 'a', len(auxnamesin)) == ",") +1
      if (naux .gt. 0) then
        if (naux .ne. nnames) then
          call amrex_error('simple ::actual_network_init wrong number of aux variable names')  
        end if
        read(auxnamesin,*) aux_names
        do iaux = 1,naux 
          aux_names(iaux) = trim(adjustl(aux_names(iaux)))
        end do
      end if

      ! Check to make sure, and if not, throw an error
      if ( NUM_SPECIES .lt. 1 ) then
              call bl_error("Network cannot have a nonpositive number of species.")
      endif

      if ( NUM_REACTIONS .lt. 0 ) then
              call bl_error("Network cannot have a negative number of reactions.")
      endif

      if ( naux .lt. 0 ) then
              call bl_error("Network cannot have a negative number of auxiliary variables.")
      endif

      network_initialized = .true.
    end subroutine network_init


    subroutine network_close
      call chemistry_close()

      deallocate(aux_names)

      network_initialized = .false.
    end subroutine network_close


    !function network_species_index(name) result(r)

    !  character(len=*) :: name
    !  integer :: r, n

    !  r = -1
    !  do n = 1, NUM_SPECIES
    !     if (name == spec_names(n)) then
    !        r = n
    !        return
    !     endif
    !  enddo

    !end function network_species_index

end module fuego_chemistry
