#include "Transport.H"

TransParm *trans_parm_g;

void
transport_init()
{
   TransParm trans_parm;

   // Default
   trans_parm.Prandtl_number = 0.7;
   trans_parm.viscosity_mu_ref = 17.16;
   trans_parm.viscosity_T_ref = 273.15;
   trans_parm.viscosity_S = 110.4;
   trans_parm.const_bulk_viscosity = 0.0;
   trans_parm.const_diffusivity = 1.0;

   // User-specified
   amrex::ParmParse pp("transport");
   pp.query("Prandtl_number",       trans_parm.Prandtl_number);
   pp.query("viscosity_mu_ref",     trans_parm.viscosity_mu_ref);
   pp.query("viscosity_T_ref",      trans_parm.viscosity_T_ref);
   pp.query("viscosity_S",          trans_parm.viscosity_S);
   pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
   pp.query("const_diffusivity",    trans_parm.const_diffusivity);

   /* GPU */
   trans_parm_g = (TransParm *) amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_GPU
   amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
   std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
}

void
transport_close()
{
  amrex::The_Arena()->free(trans_parm_g);
}
