#include <Transport.H>

TransParm *trans_parm_g;

void transport_init() 
{
    TransParm trans_parm;

    /* CPU */
    trans_parm.trans_wt   = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.trans_iwt  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.trans_eps  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);    
    trans_parm.trans_sig  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.trans_dip  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.trans_pol  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.trans_zrot = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);

    trans_parm.trans_fitmu   = (amrex::Real*) amrex::The_Arena()->alloc(
                              sizeof(amrex::Real) * trans_parm.array_size * trans_parm.fit_length);
    trans_parm.trans_fitlam  = (amrex::Real*) amrex::The_Arena()->alloc(
                              sizeof(amrex::Real) * trans_parm.array_size * trans_parm.fit_length);
    trans_parm.trans_fitdbin = (amrex::Real*) amrex::The_Arena()->alloc(
                               sizeof(amrex::Real) * trans_parm.array_size * trans_parm.array_size * trans_parm.fit_length);
    trans_parm.trans_nlin    = (int*) amrex::The_Arena()->alloc(sizeof(int) * trans_parm.array_size);

    egtransetWT(trans_parm.trans_wt);
    egtransetEPS(trans_parm.trans_eps);
    egtransetSIG(trans_parm.trans_sig);
    egtransetDIP(trans_parm.trans_dip);
    egtransetPOL(trans_parm.trans_pol);
    egtransetZROT(trans_parm.trans_zrot);
    egtransetNLIN(trans_parm.trans_nlin);
    egtransetCOFETA(trans_parm.trans_fitmu);
    egtransetCOFLAM(trans_parm.trans_fitlam);
    egtransetCOFD(trans_parm.trans_fitdbin);

    for (int i = 0; i < trans_parm.array_size; ++i) {
        trans_parm.trans_iwt[i] = 1. / trans_parm.trans_wt[i];
    }

    /* GPU */
    trans_parm_g = (TransParm *) amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_CUDA
    amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
    std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
}

void transport_close()
{
  amrex::The_Device_Arena()->free(trans_parm_g);
}
