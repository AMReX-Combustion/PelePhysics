#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include "TransportTypes.H"
#include "EOS.H"

namespace pele {
namespace physics {
namespace transport {

struct ConstTransport;
struct SimpleTransport;
struct SutherlandTransport;

template <typename EOSType, typename TransportType>
struct TransParm
{
};

template <typename EOSType>
struct TransParm<EOSType, ConstTransport>
{
  amrex::Real const_viscosity{0.0};
  amrex::Real const_bulk_viscosity{0.0};
  amrex::Real const_diffusivity{0.0};
  amrex::Real const_conductivity{0.0};
  TransParm() {}
  ~TransParm() {}
  void allocate()
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", const_viscosity);
    pp.query("const_bulk_viscosity", const_bulk_viscosity);
    pp.query("const_conductivity", const_conductivity);
    pp.query("const_diffusivity", const_diffusivity);
  }
  void deallocate() {}
};

template <typename EOSType>
struct TransParm<EOSType, SutherlandTransport>
{
  amrex::Real Prandtl_number{0.7};
  amrex::Real viscosity_mu_ref{17.16};
  amrex::Real viscosity_T_ref{273.15};
  amrex::Real viscosity_S{110.4};
  amrex::Real const_bulk_viscosity{0.0};
  amrex::Real const_diffusivity{1.0};
  TransParm() {}
  ~TransParm() {}
  void allocate()
  {
    amrex::ParmParse pp("transport");
    pp.query("Prandtl_number", Prandtl_number);
    pp.query("viscosity_mu_ref", viscosity_mu_ref);
    pp.query("viscosity_T_ref", viscosity_T_ref);
    pp.query("viscosity_S", viscosity_S);
    pp.query("const_bulk_viscosity", const_bulk_viscosity);
    pp.query("const_diffusivity", const_diffusivity);
  }
  void deallocate() {}
};

template <typename EOSType>
struct TransParm<EOSType, SimpleTransport>
{
  amrex::Real* trans_wt = nullptr;
  amrex::Real* trans_iwt = nullptr;
  amrex::Real* trans_eps = nullptr;
  amrex::Real* trans_sig = nullptr;
  amrex::Real* trans_dip = nullptr;
  amrex::Real* trans_pol = nullptr;
  amrex::Real* trans_zrot = nullptr;
  amrex::Real* trans_fitmu = nullptr;
  amrex::Real* trans_fitlam = nullptr;
  amrex::Real* trans_fitdbin = nullptr;
  int* trans_nlin = nullptr;
  bool m_allocated{false};
  TransParm() {}
  ~TransParm() {}
  void allocate()
  {
    if (!m_allocated) {
      trans_wt = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_iwt = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_eps = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_sig = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_dip = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_pol = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_zrot = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);

      trans_fitmu = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
      trans_fitlam = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
      trans_fitdbin = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES * NUM_FIT);
      trans_nlin = (int*)amrex::The_Arena()->alloc(sizeof(int) * NUM_SPECIES);

      egtransetWT(trans_wt);
      egtransetEPS(trans_eps);
      egtransetSIG(trans_sig);
      egtransetDIP(trans_dip);
      egtransetPOL(trans_pol);
      egtransetZROT(trans_zrot);
      egtransetNLIN(trans_nlin);
      egtransetCOFETA(trans_fitmu);
      egtransetCOFLAM(trans_fitlam);
      egtransetCOFD(trans_fitdbin);

      for (int i = 0; i < NUM_SPECIES; ++i) {
        trans_iwt[i] = 1. / trans_wt[i];
      }
      m_allocated = true;
    }
  }

  void deallocate()
  {
    if (m_allocated) {
      amrex::The_Arena()->free(trans_wt);
      amrex::The_Arena()->free(trans_iwt);
      amrex::The_Arena()->free(trans_eps);
      amrex::The_Arena()->free(trans_sig);
      amrex::The_Arena()->free(trans_dip);
      amrex::The_Arena()->free(trans_pol);
      amrex::The_Arena()->free(trans_zrot);
      amrex::The_Arena()->free(trans_fitmu);
      amrex::The_Arena()->free(trans_fitlam);
      amrex::The_Arena()->free(trans_fitdbin);
      amrex::The_Arena()->free(trans_nlin);
    }
  }
};

template <>
struct TransParm<eos::SRK, SimpleTransport>
{
  amrex::Real* trans_wt = nullptr;
  amrex::Real* trans_iwt = nullptr;
  amrex::Real* trans_eps = nullptr;
  amrex::Real* trans_sig = nullptr;
  amrex::Real* trans_dip = nullptr;
  amrex::Real* trans_pol = nullptr;
  amrex::Real* trans_zrot = nullptr;
  amrex::Real* trans_fitmu = nullptr;
  amrex::Real* trans_fitlam = nullptr;
  amrex::Real* trans_fitdbin = nullptr;
  int* trans_nlin = nullptr;
  amrex::Real* Afac = nullptr;
  amrex::Real* Bfac = nullptr;
  amrex::Real* Sigmaij = nullptr;
  amrex::Real* sqrtT2ij = nullptr;
  amrex::Real* sqrtEpsilonij = nullptr;
  amrex::Real* sqrtMWij = nullptr;
  amrex::Real* sqrtKappaij = nullptr;
  amrex::Real* Upsilonijk = nullptr;
  amrex::Real* Kappai = nullptr;
  amrex::Real* omega = nullptr;
  bool m_allocated{false};
  TransParm() {}
  ~TransParm() {}
  void allocate()
  {
    if (!m_allocated) {
      trans_wt = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_iwt = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_eps = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_sig = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_dip = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_pol = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      trans_zrot = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);

      trans_fitmu = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
      trans_fitlam = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
      trans_fitdbin = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES * NUM_FIT);
      trans_nlin = (int*)amrex::The_Arena()->alloc(sizeof(int) * NUM_SPECIES);

      egtransetWT(trans_wt);
      egtransetEPS(trans_eps);
      egtransetSIG(trans_sig);
      egtransetDIP(trans_dip);
      egtransetPOL(trans_pol);
      egtransetZROT(trans_zrot);
      egtransetNLIN(trans_nlin);
      egtransetCOFETA(trans_fitmu);
      egtransetCOFLAM(trans_fitlam);
      egtransetCOFD(trans_fitdbin);

      for (int i = 0; i < NUM_SPECIES; ++i) {
        trans_iwt[i] = 1. / trans_wt[i];
      }

      // Nonideal transport coefficients computed using Chung's method:
      // Chung, T.H., Ajlan, M., Lee, L.L. and Starling, K.E., 1988.
      // Generalized multiparameter correlation for nonpolar and polar
      // fluid transport properties. Industrial & engineering chemistry
      // research, 27(4), pp.671-679.

      Afac =
        (amrex::Real*)amrex::The_Arena()->alloc(sizeof(amrex::Real) * 10 * 4);
      Bfac =
        (amrex::Real*)amrex::The_Arena()->alloc(sizeof(amrex::Real) * 7 * 4);
      omega = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      Kappai = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES);
      Sigmaij = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);
      sqrtT2ij = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);
      sqrtEpsilonij = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);
      sqrtMWij = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);
      sqrtKappaij = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);
      Upsilonijk = (amrex::Real*)amrex::The_Arena()->alloc(
        sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES * NUM_SPECIES);

      // Initialize coefficients of model
      {
        amrex::Real tmp1[NUM_SPECIES], tmp2[NUM_SPECIES], tmp3[NUM_SPECIES];
        GET_CRITPARAMS(tmp1, tmp2, tmp3, omega);
      }

      Afac[0] = 6.32402;
      Afac[1] = 50.41190;
      Afac[2] = -51.68010;
      Afac[3] = 1189.020;
      Afac[4] = 0.12102e-2;
      Afac[5] = -0.11536e-2;
      Afac[6] = -0.62571e-2;
      Afac[7] = 0.37283e-1;
      Afac[8] = 5.28346;
      Afac[9] = 254.209;
      Afac[10] = -168.481;
      Afac[11] = 3898.27;
      Afac[12] = 6.62263;
      Afac[13] = 38.09570;
      Afac[14] = -8.46414;
      Afac[15] = 31.41780;
      Afac[16] = 19.74540;
      Afac[17] = 7.63034;
      Afac[18] = -14.35440;
      Afac[19] = 31.52670;
      Afac[20] = -1.89992;
      Afac[21] = -12.5367;
      Afac[22] = 4.98529;
      Afac[23] = -18.15070;
      Afac[24] = 24.27450;
      Afac[25] = 3.44945;
      Afac[26] = -11.29130;
      Afac[27] = 69.34660;
      Afac[28] = 0.79716;
      Afac[29] = 1.11764;
      Afac[30] = 0.12348e-1;
      Afac[31] = -4.11661;
      Afac[32] = -0.23816;
      Afac[33] = 0.67695e-1;
      Afac[34] = -0.81630;
      Afac[35] = 4.02528;
      Afac[36] = 0.68629e-1;
      Afac[37] = 0.34793;
      Afac[38] = 0.59256;
      Afac[39] = -0.72663;

      Bfac[0] = 2.41657;
      Bfac[1] = 0.74824;
      Bfac[2] = -0.91858;
      Bfac[3] = 121.721;
      Bfac[4] = -0.50924;
      Bfac[5] = -1.50936;
      Bfac[6] = -49.99120;
      Bfac[7] = 69.9834;
      Bfac[8] = 6.61069;
      Bfac[9] = 5.62073;
      Bfac[10] = 64.75990;
      Bfac[11] = 27.0389;
      Bfac[12] = 14.54250;
      Bfac[13] = -8.91387;
      Bfac[14] = -5.63794;
      Bfac[15] = 74.3435;
      Bfac[16] = 0.79274;
      Bfac[17] = 0.82019;
      Bfac[18] = -0.69369;
      Bfac[19] = 6.31734;
      Bfac[20] = -5.86340;
      Bfac[21] = 12.80050;
      Bfac[22] = 9.58926;
      Bfac[23] = -65.5292;
      Bfac[24] = 81.17100;
      Bfac[25] = 114.15800;
      Bfac[26] = -60.84100;
      Bfac[27] = 466.7750;

      for (int i = 0; i < NUM_SPECIES; ++i) {
        for (int j = 0; j < NUM_SPECIES; ++j) {
          Sigmaij[i * NUM_SPECIES + j] =
            0.5 * (trans_sig[i] + trans_sig[j]) * 1e-8; // converted to cm
        }
      }

      // Initialize Kappa, which has nonzero values only for specific polar
      // species
      for (int i = 0; i < NUM_SPECIES; ++i) {
        Kappai[i] = 0.0;
      }
      {
        amrex::Vector<std::string> spec_names_kappa;
        spec_names_kappa.resize(1);
        spec_names_kappa[0] = "Null";
        CKSYMS_STR(spec_names_kappa);
        for (int i = 0; i < NUM_SPECIES; i++) {
          if (spec_names_kappa[i] == "H2O") {
            Kappai[i] = 0.075908;
            // Kappai[i] = 0.076;
          } else if (spec_names_kappa[i] == "CH3OH") {
            Kappai[i] = 0.215175;
          } else if (spec_names_kappa[i] == "CH3CH2OH") {
            Kappai[i] = 0.174823;
          }
        }
      }
      for (int i = 0; i < NUM_SPECIES; ++i) {
        for (int j = 0; j < NUM_SPECIES; ++j) {
	  const int idx = i * NUM_SPECIES + j;
          sqrtT2ij[idx] = std::sqrt(trans_sig[i] + trans_sig[j]);
	  sqrtEpsilonij[idx] = std::sqrt(trans_eps[i] * trans_eps[j]);
	  sqrtMWij[idx] = std::sqrt(2.0 / (trans_iwt[i] + trans_iwt[j]));
	  sqrtKappaij[idx] = std::sqrt(Kappai[i] * Kappai[j]);
        }
      }

      for (int i = 0; i < NUM_SPECIES; ++i) {
        for (int j = 0; j < NUM_SPECIES; ++j) {
	  if (i != j) {
	    const int idx_ij = i + NUM_SPECIES * j;
	    const amrex::Real S_ij = Sigmaij[idx_ij];
	    const amrex::Real S_ij_inv = 1.0 / S_ij;
	    for (int k = 0; k < NUM_SPECIES; ++k) {
	      const amrex::Real S_ik = Sigmaij[i + NUM_SPECIES * k];
	      const amrex::Real S_jk = Sigmaij[j + NUM_SPECIES * k];
	      Upsilonijk[idx_ij*NUM_SPECIES + k] = 
	        trans_iwt[k] *
		(8.0 * (S_ik * S_ik * S_ik + S_jk * S_jk * S_jk) -
		 6.0 * (S_ik * S_ik + S_jk * S_jk) * S_ij -
		 3.0 *
                 ((S_ik * S_ik - S_jk * S_jk) * (S_ik * S_ik - S_jk * S_jk)) *
                 S_ij_inv +
		 S_ij * S_ij * S_ij);
	    }
	  }
	}
      }
      
      m_allocated = true;
    }
  }

  void deallocate()
  {
    if (m_allocated) {
      amrex::The_Arena()->free(trans_wt);
      amrex::The_Arena()->free(trans_iwt);
      amrex::The_Arena()->free(trans_eps);
      amrex::The_Arena()->free(trans_sig);
      amrex::The_Arena()->free(trans_dip);
      amrex::The_Arena()->free(trans_pol);
      amrex::The_Arena()->free(trans_zrot);
      amrex::The_Arena()->free(trans_fitmu);
      amrex::The_Arena()->free(trans_fitlam);
      amrex::The_Arena()->free(trans_fitdbin);
      amrex::The_Arena()->free(trans_nlin);

      amrex::The_Arena()->free(Afac);
      amrex::The_Arena()->free(Bfac);
      amrex::The_Arena()->free(Kappai);
      amrex::The_Arena()->free(Sigmaij);
      amrex::The_Arena()->free(sqrtT2ij);
      amrex::The_Arena()->free(sqrtEpsilonij);
      amrex::The_Arena()->free(sqrtMWij);
      amrex::The_Arena()->free(sqrtKappaij);
      amrex::The_Arena()->free(Upsilonijk);
      amrex::The_Arena()->free(omega);
    }
  }
};

template <typename TransportType>
class TransportParams
{
public:
  TransportParams() {}

  ~TransportParams() {}

  void allocate()
  {
    m_h_trans_parm.allocate();
    if (!m_device_allocated) {
      m_d_trans_parm =
        (TransParm<EosType, TransportType>*)amrex::The_Device_Arena()->alloc(
          sizeof(m_h_trans_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void deallocate()
  {
    m_h_trans_parm.deallocate();
    if (m_device_allocated)
      amrex::The_Device_Arena()->free(m_d_trans_parm);
  }

  void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {

      // would have liked to use this but apparently only works for trivial
      // types
      // amrex::Gpu::copy(
      //   amrex::Gpu::hostToDevice, &m_h_trans_parm, &m_h_trans_parm + 1,
      //   m_d_trans_parm);
#ifdef AMREX_USE_GPU
      amrex::Gpu::htod_memcpy(
        m_d_trans_parm, &m_h_trans_parm, sizeof(m_h_trans_parm));
#else
      std::memcpy(m_d_trans_parm, &m_h_trans_parm, sizeof(m_h_trans_parm));
#endif
    }
  }

  TransParm<EosType, TransportType>& host_trans_parm()
  {
    return m_h_trans_parm;
  }
  const TransParm<EosType, TransportType>* device_trans_parm()
  {
    return m_d_trans_parm;
  }

private:
  TransParm<EosType, TransportType> m_h_trans_parm;
  TransParm<EosType, TransportType>* m_d_trans_parm;
  bool m_device_allocated{false};
};
} // namespace transport
} // namespace physics
} // namespace pele
#endif
