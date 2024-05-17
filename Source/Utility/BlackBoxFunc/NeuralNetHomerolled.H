#ifndef _NEURAL_NET_HR_H_
#define _NEURAL_NET_HR_H_

#include "NeuralNetModelDef.H"
#include "ManifoldFunc.H"

namespace pele {
namespace physics {

struct NNFuncData : ManFuncData
{
  amrex::Real* nnrdata;          // Real data representing the NNModel
  int* nnidata;                  // int data representing the NNModel
  char nn_filename[len_str + 1]; // Path to neural network file
};

class NNFuncParams
  : public PeleParamsGeneric<
      ManFuncData>::Register<NNFuncParams, PeleParams<NNFuncData, ManFuncData>>
{
public:
  static std::string identifier() { return "NeuralNet"; }
};

template <>
struct InitParm<NNFuncData, ManFuncData>
{
  static void host_initialize(PeleParams<NNFuncData, ManFuncData>* parm_in)
  {
    parm_in->m_h_parm.manmodel = ManifoldModel::NEURAL_NET;

    std::string pp_pref =
      parm_in->m_host_only_parm.parm_parse_prefix.empty()
        ? "network"
        : parm_in->m_host_only_parm.parm_parse_prefix + ".network";
    amrex::ParmParse pp(pp_pref);

    int verbose = 2;
    pp.query("v", verbose);

    std::string nn_filename;
    pp.get("filename", nn_filename);
    if (verbose > 0) {
      amrex::Print() << "Loading neural network data from file: " << nn_filename
                     << std::endl;
    }
    std::string ext = amrex::trim(nn_filename).substr(nn_filename.length() - 4);
    if (ext != ".pnn") {
      amrex::Abort("Home-rolled neural net file should have .pnn extension.");
    }
    std::string info_filename;
    pp.get("info_filename", info_filename);
    if (verbose > 0) {
      amrex::Print() << "  with network metadata from file: " << info_filename
                     << std::endl;
    }

    bool cmlm_net = false;
    pp.query("cmlm_net", cmlm_net);

    NNModel nnmodel(nn_filename, cmlm_net);
    pack_model(nnmodel, parm_in->m_h_parm);
    read_metadata(info_filename, parm_in->m_h_parm);

    strncpy(
      parm_in->m_h_parm.nn_filename, nn_filename.c_str(),
      parm_in->m_h_parm.len_str * sizeof(char));
    parm_in->m_h_parm.nn_filename[parm_in->m_h_parm.len_str] = '\0';

    if (verbose > 1) {
      print(parm_in->m_h_parm);
      nnmodel.print();
    }
  }

  static void host_deallocate(PeleParams<NNFuncData, ManFuncData>* parm_in)
  {
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.nnrdata);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.nnidata);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.varnames);
    amrex::The_Pinned_Arena()->free(parm_in->m_h_parm.dimnames);
  }

  static void pack_model(NNModel& nnmodel, NNFuncData& m_h_nnf_data)
  {
    // Get sizes of buffers
    int nreals, nints;
    nnmodel.buffer_sizes_for_packing(nreals, nints);

    // Allocate memory for buffers so that they can be copied efficiently to the
    // GPU
    m_h_nnf_data.nnrdata = static_cast<amrex::Real*>(
      amrex::The_Pinned_Arena()->alloc(nreals * sizeof(amrex::Real)));
    m_h_nnf_data.nnidata =
      static_cast<int*>(amrex::The_Pinned_Arena()->alloc(nints * sizeof(int)));

    // Pack model
    nnmodel.pack(m_h_nnf_data.nnrdata, m_h_nnf_data.nnidata);
  }

  static void
  read_metadata(std::string& info_filename, NNFuncData& m_h_nnf_data)
  {
    std::ifstream fi(info_filename, std::ios::in);

    constexpr int READING_VAR = 0;
    constexpr int READING_VAL = 1;
    constexpr int READING_COMMENT = 2;
    constexpr int LINE_BREAK = 3;

    int state = READING_VAR;

    std::unordered_map<std::string, std::vector<std::string>> umap;
    std::string cur_var;
    std::string str;

    int str_size = m_h_nnf_data.len_str * sizeof(char);
    int str_len = m_h_nnf_data.len_str;

    if (!fi.is_open()) {
      amrex::Abort(
        "Unable to open supplementary info file for neural net manifold EOS.");
    }

    fi.seekg(0, std::ios::end);
    std::streampos length = fi.tellg();
    fi.seekg(0, std::ios::beg);

    std::vector<char> buffer(length);
    fi.read(buffer.data(), length);
    std::istringstream is(std::string(buffer.data(), length));

    while (!is.eof()) {
      if (is.peek() == '\n') {
        if (state == LINE_BREAK) {
          state = READING_VAL;
        } else {
          state = READING_VAR;
        }
      }

      is >> str;

      if (str == "#") {
        state = READING_COMMENT;
        continue;
      }

      if (state == READING_COMMENT) {
        continue;
      }

      if (str == "=") {
        state = READING_VAL;
        continue;
      }

      if (str == "\\") {
        state = LINE_BREAK;
        continue;
      }

      if (state == READING_VAR) {
        for (int i = 0; i < str.length(); i++) {
          str[i] = std::tolower(str[i]);
        }
        cur_var = str;
      }

      if (state == READING_VAL) {
        umap[cur_var].push_back(str);
      }
    }

    AMREX_ALWAYS_ASSERT(umap["model_name"].size() == 1);
    AMREX_ALWAYS_ASSERT(umap["ndim"].size() == 1);
    AMREX_ALWAYS_ASSERT(umap["nvar"].size() == 1);
    AMREX_ALWAYS_ASSERT(umap["nmanpar"].size() == 1);

    strncpy(m_h_nnf_data.model_name, umap["model_name"][0].c_str(), str_size);
    m_h_nnf_data.model_name[str_len] = '\0';

    m_h_nnf_data.Ndim = stoi(umap["ndim"][0]);
    m_h_nnf_data.Nvar = stoi(umap["nvar"][0]);
    m_h_nnf_data.Nmanpar = stoi(umap["nmanpar"][0]);

    // Hard coded limit on size of network (required for derivs function)
    AMREX_ALWAYS_ASSERT(m_h_nnf_data.Ndim < MAXD_NETWORK);

    AMREX_ALWAYS_ASSERT(umap["dimnames"].size() == m_h_nnf_data.Ndim);
    AMREX_ALWAYS_ASSERT(umap["varnames"].size() == m_h_nnf_data.Nvar);

    m_h_nnf_data.dimnames = static_cast<char*>(
      amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim * str_size));
    std::vector<std::string>& dimnames = umap["dimnames"];
    for (int i = 0; i < m_h_nnf_data.Ndim; i++) {
      // Pad string with whitespace as is done for the tabular data
      std::string sd =
        dimnames[i] + std::string(str_len - dimnames[i].length(), ' ');
      strncpy(&m_h_nnf_data.dimnames[i * str_len], sd.c_str(), str_size);
    }

    m_h_nnf_data.varnames = static_cast<char*>(
      amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Nvar * str_size));
    std::vector<std::string>& varnames = umap["varnames"];
    for (int i = 0; i < m_h_nnf_data.Nvar; i++) {
      // Pad string with whitespace as is done for the tabular data
      std::string sv =
        varnames[i] + std::string(str_len - varnames[i].length(), ' ');
      strncpy(&m_h_nnf_data.varnames[i * str_len], sv.c_str(), str_size);
    }
  }

  static void print(NNFuncData& m_h_nnf_data)
  {
    amrex::Print() << std::endl;
    amrex::Print() << "NEURAL NETWORK MODEL" << std::endl;
    amrex::Print() << "Neural net filename: " << m_h_nnf_data.nn_filename
                   << std::endl;
    amrex::Print() << "Model name: " << m_h_nnf_data.model_name << std::endl;
    amrex::Print() << "Ndim: " << m_h_nnf_data.Ndim << std::endl;
    amrex::Print() << "Nvar: " << m_h_nnf_data.Nvar << std::endl;
    amrex::Print() << "Nmanpar: " << m_h_nnf_data.Nmanpar << std::endl;
    amrex::Print() << std::endl;

    amrex::Print() << "Dimnames: Index | Variable Name" << std::endl;
    for (int i = 0; i < m_h_nnf_data.Ndim; i++) {
      std::string dimname(
        &m_h_nnf_data.dimnames[i * m_h_nnf_data.len_str], m_h_nnf_data.len_str);
      amrex::Print() << i << " | " << amrex::trim(dimname) << std::endl;
    }
    amrex::Print() << std::endl;

    amrex::Print() << "Varnames: Index | Variable Name" << std::endl;
    for (int i = 0; i < m_h_nnf_data.Nvar; i++) {
      std::string varname(
        &m_h_nnf_data.varnames[i * m_h_nnf_data.len_str], m_h_nnf_data.len_str);
      amrex::Print() << i << " | " << amrex::trim(varname) << std::endl;
    }
    amrex::Print() << std::endl;
  }
};

class NNFunc : public ManifoldFunc
{
public:
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc(const NNFuncData* nnf_data_in)
    : nnf_data{nnf_data_in}, nnmodel(nnf_data_in->nnrdata, nnf_data_in->nnidata)
  {
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void
  get_value(const int ivar, const amrex::Real indata[], amrex::Real& out)
  {
    out = nnmodel(indata)[ivar];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_values(
    const int nvar,
    const int ivar[],
    const amrex::Real indata[],
    amrex::Real out[])
  {
    const amrex::Real* outdata = nnmodel(indata);
    for (int i = 0; i < nvar; i++) {
      out[i] = (ivar[i] >= 0) ? outdata[ivar[i]] : 0.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_all_values(const amrex::Real indata[], amrex::Real out[])
  {
    const amrex::Real* outdata = nnmodel(indata);
    for (int i = 0; i < nnf_data->Nvar; i++) {
      out[i] = outdata[i];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void
  get_derivs(const int ivar, const amrex::Real indata[], amrex::Real derivs[])
  {
    // Use first-order finite difference approach for now
    constexpr amrex::Real eps = 1e-5;
    const amrex::Real outval0 = nnmodel(indata)[ivar];
    amrex::Real diff;

    // TODO: Better way to do this than hardcoding the size?
    amrex::Real indata_copy[MAXD_NETWORK];
    for (int i = 0; i < nnf_data->Ndim; i++) {
      indata_copy[i] = indata[i];
    }

    for (int i = 0; i < nnf_data->Ndim; i++) {
      diff = (indata[i] != 0.0) ? eps * indata[i] : eps;
      indata_copy[i] = indata[i] + diff;
      derivs[i] = (nnmodel(indata_copy)[ivar] - outval0) / diff;
      indata_copy[i] = indata[i];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual ManifoldModel model() { return nnf_data->manmodel; }

private:
  // Variables
  const NNFuncData* nnf_data;
  NNModel nnmodel;
  amrex::Real* indata_copy;

}; // class NNFunc

} // namespace physics
} // namespace pele

#endif
