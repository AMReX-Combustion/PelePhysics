module fuego_chemistry

    interface
        subroutine egtransetLENIMC(LENIMC) bind(c,name='egtransetLENIMC')
            integer, intent(inout) :: LENIMC
        end subroutine
    !end interface

    !interface
        subroutine egtransetLENRMC(LENRMC) bind(c,name='egtransetLENRMC')
            integer, intent(inout) :: LENRMC
        end subroutine
    !end interface

    !interface
        subroutine egtransetNO(NO) bind(c,name='egtransetNO')
            integer, intent(inout) :: NO
        end subroutine
    !end interface

    !interface
        subroutine egtransetKK(KK) bind(c,name='egtransetKK')
            integer, intent(inout) :: KK
        end subroutine
    !end interface

    !interface
        subroutine egtransetNLITE(NLITE) bind(c,name='egtransetNLITE')
            integer, intent(inout) :: NLITE
        end subroutine
    !end interface

    !interface
        subroutine egtransetPATM(PATM) bind(c,name='egtransetPATM')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: PATM
        end subroutine
    !end interface

    !interface
        subroutine egtransetWT(WT) bind(c,name='egtransetWT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: WT(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetEPS(EPS) bind(c,name='egtransetEPS')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: EPS(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetSIG(SIG) bind(c,name='egtransetSIG')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: SIG(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetDIP(DIP) bind(c,name='egtransetDIP')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: DIP(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetPOL(POL) bind(c,name='egtransetPOL')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) ::POL(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetZROT(ZROT) bind(c,name='egtransetZROT')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: ZROT(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetNLIN(NLIN) bind(c,name='egtransetNLIN')
            use amrex_fort_module, only : amrex_real
            integer, intent(inout) :: NLIN(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetCOFETA(COFETA) bind(c,name='egtransetCOFETA')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFETA(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetCOFLAM(COFLAM) bind(c,name='egtransetCOFLAM')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFLAM(*)
        end subroutine
    !end interface

    !interface
        subroutine egtransetCOFD(COFD) bind(c,name='egtransetCOFD')
            use amrex_fort_module, only : amrex_real
            real(amrex_real), intent(inout) :: COFD(*)
        end subroutine
    !end interface

    !interface
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

    end interface

end module fuego_chemistry
