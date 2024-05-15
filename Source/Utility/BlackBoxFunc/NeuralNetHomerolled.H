#ifndef _NEURAL_NET_HR_H_
#define _NEURAL_NET_HR_H_

#include "NeuralNetModelDef.H"

namespace pele {
namespace physics {
  
  struct NNFuncData : ManFuncData
  {
    amrex::Real* nnrdata;          // Real data representing the NNModel
    int* nnidata;                  // int data representing the NNModel
    char nn_filename[len_str + 1]; // Path to neural network file
  };
  
class NNFuncParams : public ManFuncParams
{
public:
  NNFuncParams() {}

  ~NNFuncParams() {}

  virtual void initialize()
  {
    m_h_nnf_data.manmodel = ManifoldModel::NEURAL_NET;

    amrex::ParmParse pp("manifold");

    std::string nn_filename;
    pp.get("filename", nn_filename);
    std::string ext = amrex::trim(nn_filename).substr(nn_filename.length() - 4);
    if (ext != ".pnn") {
      amrex::Abort("Home-rolled neural net file should have .pnn extension.");
    }

    bool cmlm_net = false;
    pp.query("cmlm_net", cmlm_net);

    NNModel nnmodel(nn_filename, cmlm_net);

    std::string info_filename;
    pp.get("info_filename", info_filename);

    strncpy(
      m_h_nnf_data.nn_filename, nn_filename.c_str(),
      m_h_nnf_data.len_str * sizeof(char));
    m_h_nnf_data.nn_filename[m_h_nnf_data.len_str] = '\0';

    pack_model(nnmodel);

    read_metadata(info_filename);

    pp.query("v", m_verbose);
    if (m_verbose >= 1) {
      print();
      nnmodel.print();
    }

    allocate();
  }

  void pack_model(NNModel& nnmodel)
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

  void read_metadata(std::string& info_filename)
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

  void print()
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

  virtual void allocate()
  {
    if (!m_device_allocated) {
      m_d_nnf_data = static_cast<NNFuncData*>(
        amrex::The_Device_Arena()->alloc(sizeof(m_h_nnf_data)));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  virtual void deallocate()
  {
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.nnrdata);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.nnidata);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.varnames);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.dimnames);
    if (m_device_allocated) {
      amrex::The_Device_Arena()->free(m_d_nnf_data);
    }
  }

  virtual void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_nnf_data, &m_h_nnf_data + 1,
        m_d_nnf_data);
    }
  }
  
  virtual ManFuncData& host_manfunc_data() { return m_h_nnf_data; }

  virtual ManFuncData* device_manfunc_data() { return m_d_nnf_data; }

  NNFuncData& host_nnfunc_data() { return m_h_nnf_data; }

  NNFuncData* device_nnfunc_data() { return m_d_nnf_data; }

private:
  int m_verbose = 0;
  NNFuncData m_h_nnf_data;
  NNFuncData* m_d_nnf_data;
  bool m_device_allocated{false};
};

class NNFunc : public ManifoldFunc
{
public:
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc(NNFuncData* nnf_data_in)
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
  NNFuncData* nnf_data;
  NNModel nnmodel;
  amrex::Real* indata_copy;

}; // class NNFunc

} // namespace physics
} // namespace pele

#endif
