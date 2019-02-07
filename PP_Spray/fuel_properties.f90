module fuel_properties

     use meth_params_module

!--------------------------------------------------------------------------------
! Fuel properties
!--------------------------------------------------------------------------------

      integer :: nspec_f
      integer, parameter :: max_nspec_f = 1
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_density
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_crit_temp
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_latent
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_mass_frac
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_boil_temp
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_cp
      real(kind=dp_t), dimension(max_nspec_f) :: fuel_molwt
      integer(kind=dp_t), dimension(max_nspec_f) :: fuel_indx

end module fuel_properties

