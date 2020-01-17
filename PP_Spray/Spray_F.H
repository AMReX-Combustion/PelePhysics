#ifndef _Spray_F_H_
#define _Spray_F_H_
#include <AMReX_BLFort.H>

extern "C"
{
  void calc_vdot(const int* fill, const amrex::Real* vel_d, const amrex::Real* diameter_d,
		 const amrex::Real* rho_d, const amrex::Real* T_d, const amrex::Real* vel_g,
		 const amrex::Real* rho_g, const amrex::Real* temp_g, const amrex::Real* pressure_g,
		 const amrex::Real* molwt_g, const amrex::Real* vel_dot);

  void import_fuel_properties(const int* nfuel_species, amrex::Real* fuel_mass_frac,
			      amrex::Real* fuel_crit_temp, amrex::Real* fuel_latent,
			      amrex::Real* fuel_boil_temp, amrex::Real* fuel_cp,
			      amrex::Real* fuel_molwt, const int* fuel_indx);

  void import_control_parameters(const int* heat_transfer, const int* mass_transfer, const int* mom_transfer);

  void get_num_fuel_spec(int* n_fuel_spec);

  void calc_thermal_B(amrex::Real* Nu, const amrex::Real* Sh, const amrex::Real* Pr, const amrex::Real* Sc,
		      const amrex::Real* cp_fuel, const amrex::Real* cp_gas, const amrex::Real* B_m,
		      amrex::Real* B_T);

}

#endif
