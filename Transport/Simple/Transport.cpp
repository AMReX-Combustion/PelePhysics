#include <Transport.H>

TransParm *trans_parm_g;
TransParm trans_parm;

void transport_init() 
{
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

#ifdef PELEPHYSICS_NONIDEAL_EOS
    // Nonideal transport coefficients computed using Chung's method:
    // Chung, T.H., Ajlan, M., Lee, L.L. and Starling, K.E., 1988.
    // Generalized multiparameter correlation for nonpolar and polar
    // fluid transport properties. Industrial & engineering chemistry
    // research, 27(4), pp.671-679.

    trans_parm.Afac    = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * 10 * 4);
    trans_parm.Bfac    = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * 7 * 4);
    trans_parm.omega   = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.Kappai  = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size);
    trans_parm.Sigmaij = (amrex::Real*) amrex::The_Arena()->alloc(sizeof(amrex::Real) * trans_parm.array_size
								  * trans_parm.array_size);

    // Initialize coefficients of model
    {
      amrex::Real tmp1[trans_parm.array_size], tmp2[trans_parm.array_size], tmp3[trans_parm.array_size];
      GET_CRITPARAMS(tmp1, tmp2, tmp3, trans_parm.omega);
    }
    
    trans_parm.Afac[0]  = 6.32402;    trans_parm.Afac[1]  = 50.41190;    trans_parm.Afac[2]  = -51.68010;   trans_parm.Afac[3]  = 1189.020;
    trans_parm.Afac[4]  = 0.12102e-2; trans_parm.Afac[5]  = -0.11536e-2; trans_parm.Afac[6]  = -0.62571e-2; trans_parm.Afac[7]  = 0.37283e-1;
    trans_parm.Afac[8]  = 5.28346;    trans_parm.Afac[9]  = 254.209;     trans_parm.Afac[10] = -168.481;    trans_parm.Afac[11] = 3898.27;
    trans_parm.Afac[12] = 6.62263;    trans_parm.Afac[13] = 38.09570;    trans_parm.Afac[14] = -8.46414;    trans_parm.Afac[15] = 31.41780;
    trans_parm.Afac[16] = 19.74540;   trans_parm.Afac[17] = 7.63034;     trans_parm.Afac[18] = -14.35440;   trans_parm.Afac[19] = 31.52670;
    trans_parm.Afac[20] = -1.89992;   trans_parm.Afac[21] = -12.5367;    trans_parm.Afac[22] = 4.98529;     trans_parm.Afac[23] = -18.15070;
    trans_parm.Afac[24] = 24.27450;   trans_parm.Afac[25] = 3.44945;     trans_parm.Afac[26] = -11.29130;   trans_parm.Afac[27] = 69.34660;
    trans_parm.Afac[28] = 0.79716;    trans_parm.Afac[29] = 1.11764;     trans_parm.Afac[30] = 0.12348e-1;  trans_parm.Afac[31] = -4.11661;
    trans_parm.Afac[32] = -0.23816;   trans_parm.Afac[33] = 0.67695e-1;  trans_parm.Afac[34] = -0.81630;    trans_parm.Afac[35] = 4.02528;
    trans_parm.Afac[36] = 0.68629e-1; trans_parm.Afac[37] = 0.34793;     trans_parm.Afac[38] = 0.59256;     trans_parm.Afac[39] = -0.72663;
    
    trans_parm.Bfac[0]  = 2.41657;   trans_parm.Bfac[1]  = 0.74824;   trans_parm.Bfac[2]  = -0.91858;  trans_parm.Bfac[3]  = 121.721;
    trans_parm.Bfac[4]  = -0.50924;  trans_parm.Bfac[5]  = -1.50936;  trans_parm.Bfac[6]  = -49.99120; trans_parm.Bfac[7]  = 69.9834;
    trans_parm.Bfac[8]  = 6.61069;   trans_parm.Bfac[9]  = 5.62073;   trans_parm.Bfac[10] = 64.75990;  trans_parm.Bfac[11] = 27.0389;
    trans_parm.Bfac[12] = 14.54250;  trans_parm.Bfac[13] = -8.91387;  trans_parm.Bfac[14] = -5.63794;  trans_parm.Bfac[15] = 74.3435;
    trans_parm.Bfac[16] = 0.79274;   trans_parm.Bfac[17] = 0.82019;   trans_parm.Bfac[18] = -0.69369;  trans_parm.Bfac[19] = 6.31734;
    trans_parm.Bfac[20] = -5.86340;  trans_parm.Bfac[21] = 12.80050;  trans_parm.Bfac[22] = 9.58926;   trans_parm.Bfac[23] = -65.5292;
    trans_parm.Bfac[24] = 81.17100;  trans_parm.Bfac[25] = 114.15800; trans_parm.Bfac[26] = -60.84100; trans_parm.Bfac[27] = 466.7750;

    for (int i = 0; i < trans_parm.array_size; ++i) {
      for (int j = 0; j < trans_parm.array_size; ++j) {
	trans_parm.Sigmaij[i*trans_parm.array_size + j] = 0.5*(trans_parm.trans_sig[i]
							       + trans_parm.trans_sig[j]) * 1e-8; // converted to cm
      }
    }
    
    // Initialize Kappa, which has nonzero values only for specific polar species
    for (int i = 0; i < trans_parm.array_size; ++i) {
      trans_parm.Kappai[i] = 0.0;
    }
    {
      amrex::Vector<std::string> spec_names_kappa;
      spec_names_kappa.resize(1);
      spec_names_kappa[0] = "Null";
      CKSYMS_STR(spec_names_kappa);
      for (int i = 0; i < trans_parm.array_size; i++) {
	if (spec_names_kappa[i] == "H2O") {
	  trans_parm.Kappai[i] = 0.075908;
	  //trans_parm.Kappai[i] = 0.076;
	} else if (spec_names_kappa[i] == "CH3OH") {
	  trans_parm.Kappai[i] = 0.215175;
	} else if (spec_names_kappa[i] == "CH3CH2OH") {
	  trans_parm.Kappai[i] = 0.174823;
	}
      }
    }
#endif
    
    /* GPU */
    trans_parm_g = (TransParm *) amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
    std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
}

void transport_close()
{
    amrex::The_Device_Arena()->free(trans_parm_g);

    amrex::The_Arena()->free(trans_parm.trans_wt);
    amrex::The_Arena()->free(trans_parm.trans_iwt);
    amrex::The_Arena()->free(trans_parm.trans_eps);
    amrex::The_Arena()->free(trans_parm.trans_sig);
    amrex::The_Arena()->free(trans_parm.trans_dip);
    amrex::The_Arena()->free(trans_parm.trans_pol);
    amrex::The_Arena()->free(trans_parm.trans_zrot);
    amrex::The_Arena()->free(trans_parm.trans_fitmu);
    amrex::The_Arena()->free(trans_parm.trans_fitlam);
    amrex::The_Arena()->free(trans_parm.trans_fitdbin);
    amrex::The_Arena()->free(trans_parm.trans_nlin);
#ifdef PELEPHYSICS_NONIDEAL_EOS
    amrex::The_Arena()->free(trans_parm.Afac);
    amrex::The_Arena()->free(trans_parm.Bfac);
    amrex::The_Arena()->free(trans_parm.Kappai);
    amrex::The_Arena()->free(trans_parm.Sigmaij);
    amrex::The_Arena()->free(trans_parm.omega);
#endif
}
