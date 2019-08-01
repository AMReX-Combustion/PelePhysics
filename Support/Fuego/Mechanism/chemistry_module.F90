module chemistry_module

  use amrex_fort_module, only : amrex_real
  use fuego_chemistry

#include "mechanism.h"

  implicit none

  integer, parameter :: nelements   = NUM_ELEMENTS  ! number of elements
  integer, parameter :: nspecies    = NUM_SPECIES   ! number of species
  integer, parameter :: nreactions  = NUM_REACTIONS ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), save :: elem_names(nelements)

  integer, parameter :: L_spec_name = 16 ! Each species name has at most 8 characters
  character*(L_spec_name), save :: spec_names(nspecies)

  real(amrex_real), save :: molecular_weight(nspecies), inv_mwt(nspecies)

  real(amrex_real), save :: Ru, Ruc, Patm, rwrk
  integer, save          :: iwrk

contains

  subroutine chemistry_init()
    integer :: nfit, i, ic, ii
    real(amrex_real) :: T0
    integer :: names(nspecies*L_spec_name)

    call ckinit()

    call cksyme(names, L_elem_name) 

    ic = 1
    do i = 1, nelements
       do ii=1, L_elem_name
          elem_names(i)(ii:ii) = char(names(ic))
          ic = ic + 1
       end do
    end do

    call cksyms(names, L_spec_name) 

    ic = 1
    do i = 1, nspecies
       do ii=1, L_spec_name
          spec_names(i)(ii:ii) = char(names(ic))
          ic = ic+1
       end do
    end do

    call ckwt(molecular_weight)
    inv_mwt = 1.d0 / molecular_weight

    call ckrp(Ru, Ruc, Patm)

    chemistry_initialized = .true.

  end subroutine chemistry_init


  subroutine chemistry_close()
    call ckfinalize()
  end subroutine chemistry_close


  function get_species_index(name) result (iname)
    character(len=*), intent(in) :: name
    integer :: iname
    integer :: i
    iname = -1
    do i = 1, nspecies
       if (trim(spec_names(i)) .eq. trim(name)) then
          iname = i
          exit
       end if
    end do
  end function get_species_index


end module chemistry_module
