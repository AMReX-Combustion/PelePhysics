module fuel_properties

     use meth_params_module
     use amrex_fort_module, only : amrex_real

!--------------------------------------------------------------------------------
! Fuel properties
!--------------------------------------------------------------------------------

      integer :: nspec_f
      integer, parameter :: max_nspec_f = 1
      real(amrex_real), dimension(max_nspec_f) :: fuel_density
      real(amrex_real), dimension(max_nspec_f) :: fuel_crit_temp
      real(amrex_real), dimension(max_nspec_f) :: fuel_latent
      real(amrex_real), dimension(max_nspec_f) :: fuel_mass_frac
      real(amrex_real), dimension(max_nspec_f) :: fuel_boil_temp
      real(amrex_real), dimension(max_nspec_f) :: fuel_cp
      real(amrex_real), dimension(max_nspec_f) :: fuel_molwt
      integer(amrex_real), dimension(max_nspec_f) :: fuel_indx

end module fuel_properties

