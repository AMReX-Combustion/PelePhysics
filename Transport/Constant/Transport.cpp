#include "Transport.H"

TransParm *trans_parm_g;

void
transport_init()
{
   TransParm trans_parm;

   // Default
   trans_parm.const_viscosity = 0.0;
   trans_parm.const_bulk_viscosity = 0.0;
   trans_parm.const_conductivity = 0.0;
   trans_parm.const_diffusivity = 0.0;

   // User-specified
   amrex::ParmParse pp("transport");
   pp.query("const_viscosity",      trans_parm.const_viscosity);
   pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
   pp.query("const_conductivity",   trans_parm.const_conductivity);
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
