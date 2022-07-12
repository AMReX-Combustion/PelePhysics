#ifndef _NEURAL_NET_HR_H_
#define _NEURAL_NET_HR_H_

#include "NeuralNetModelDef.H"

namespace pele {
namespace physics {
  
class NNFuncParams: public ManFuncParams
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
    std::string ext = amrex::trim(nn_filename).substr(nn_filename.length()-4);
    if(ext != ".pnn")
    {
      amrex::Abort("Home-rolled neural net file should have .pnn extension.");
    }
    
    bool cmlm_net = false;
    pp.query("cmlm_net", cmlm_net);

    NNModel nnmodel(nn_filename, cmlm_net);

    std::string info_filename;
    pp.get("info_filename", info_filename);

    strncpy(m_h_nnf_data.nn_filename, nn_filename.c_str(), m_h_nnf_data.len_str*sizeof(char));
    m_h_nnf_data.nn_filename[m_h_nnf_data.len_str] = '\0';

    pack_model(nnmodel);
    
    read_metadata(info_filename);

    pp.query("v", m_verbose);
    if(m_verbose >= 1)
    {
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
    
    // Allocate memory for buffers so that they can be copied efficiently to the GPU
    m_h_nnf_data.nnrdata = static_cast<amrex::Real*>(
        amrex::The_Pinned_Arena()->alloc(nreals*sizeof(amrex::Real)));
    m_h_nnf_data.nnidata = static_cast<int*>(
        amrex::The_Pinned_Arena()->alloc(nints*sizeof(int)));
        
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

    std::unordered_map<std::string, std::vector<std::string> > umap;
    std::string cur_var;
    std::string str;

    int str_size = m_h_nnf_data.len_str*sizeof(char);
    int str_len = m_h_nnf_data.len_str;

    if(!fi.is_open())
    {
      amrex::Abort("Unable to open supplementary info file for neural net manifold EOS.");
    }

    fi.seekg(0, std::ios::end);
    std::streampos length = fi.tellg();
    fi.seekg(0, std::ios::beg);

    std::vector<char> buffer(length);
    fi.read(buffer.data(), length);
    std::istringstream is(std::string(buffer.data(), length));

    while(!is.eof())
    {
      if(is.peek() == '\n')
      {
        if(state == LINE_BREAK)
        {
            state = READING_VAL;
        }
        else
        {
            state = READING_VAR;
        }
      }

      is >> str;

      if(str == "#")
      {
        state = READING_COMMENT;
        continue;
      }

      if(state == READING_COMMENT)
      {
        continue;
      }

      if(str == "=")
      {
        state = READING_VAL;
        continue;
      }

      if(str == "\\")
      {
        state = LINE_BREAK;
        continue;
      }

      if(state == READING_VAR)
      {
        for(int i = 0; i < str.length(); i++)
        {
          str[i] = std::tolower(str[i]);
        }
        cur_var = str;
      }

      if(state == READING_VAL)
      {
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

    AMREX_ALWAYS_ASSERT(umap["dimnames"].size() == m_h_nnf_data.Ndim);
    AMREX_ALWAYS_ASSERT(umap["varnames"].size() == m_h_nnf_data.Nvar);

    m_h_nnf_data.dimnames = static_cast<char*>(amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim*str_size));
    std::vector<std::string>& dimnames = umap["dimnames"];
    for(int i = 0; i < m_h_nnf_data.Ndim; i++)
    {
      // Pad string with whitespace as is done for the tabular data
      std::string sd = dimnames[i] + std::string(str_len - dimnames[i].length(), ' ');
      strncpy(&m_h_nnf_data.dimnames[i*str_len], sd.c_str(), str_size);
    }

    m_h_nnf_data.varnames = static_cast<char*>(amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Nvar*str_size));
    std::vector<std::string>& varnames = umap["varnames"];
    for(int i = 0; i < m_h_nnf_data.Nvar; i++)
    {
      // Pad string with whitespace as is done for the tabular data
      std::string sv = varnames[i] + std::string(str_len - varnames[i].length(), ' ');
      strncpy(&m_h_nnf_data.varnames[i*str_len], sv.c_str(), str_size);
    }

    m_h_nnf_data.Ncomb = umap["def_" + amrex::trim(dimnames[0])].size();
    m_h_nnf_data.comb_coeff = static_cast<amrex::Real*>(
        amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim*m_h_nnf_data.Ncomb*sizeof(amrex::Real)));
    m_h_nnf_data.comb_idx = static_cast<int*>(
        amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim*m_h_nnf_data.Ncomb*sizeof(int)));
    m_h_nnf_data.comb_src_idx = static_cast<int*>(
        amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim*m_h_nnf_data.Ncomb*sizeof(int)));
    for(int i = 0; i < m_h_nnf_data.Ndim; i++)
    {
      std::string dn = "def_" + dimnames[i];
      std::vector<std::string>& defn = umap[dn];
      AMREX_ALWAYS_ASSERT(defn.size() == m_h_nnf_data.Ncomb);
      for(int j = 0; j < defn.size(); j++)
      {
        auto start = 0u;
        auto end = defn[j].find("*");
        std::string coeff = defn[j].substr(start, end);
        std::string var = defn[j].substr(end+1, std::string::npos);
        std::string src_var;
        if(var.rfind("Y-") == 0)
        {
          src_var = "SRC_" + var.substr(2u, std::string::npos);
        }
        else
        {
          src_var = "SRC_" + var;
        }
        int idx = get_var_index(var.c_str(), &m_h_nnf_data, false);
        int src_idx = get_var_index(src_var.c_str(), &m_h_nnf_data, false);
        m_h_nnf_data.comb_coeff[i*m_h_nnf_data.Ncomb + j] = (amrex::Real)std::stod(coeff);
        m_h_nnf_data.comb_idx[i*m_h_nnf_data.Ncomb + j] = idx;
        m_h_nnf_data.comb_src_idx[i*m_h_nnf_data.Ncomb + j] = src_idx;
      }
    }

    m_h_nnf_data.man_bias = static_cast<amrex::Real*>(
        amrex::The_Pinned_Arena()->alloc(m_h_nnf_data.Ndim*sizeof(amrex::Real)));
    std::vector<std::string>& man_bias = umap["manibiases"];
    AMREX_ALWAYS_ASSERT(man_bias.size() == m_h_nnf_data.Ndim);
    for(int i = 0; i < m_h_nnf_data.Ndim; i++)
    {
      m_h_nnf_data.man_bias[i] = (amrex::Real)std::stod(man_bias[i]);
    }
  }

  void print()
  {
    amrex::Print() << std::endl;
    amrex::Print() << "NEURAL NETWORK MODEL" << std::endl;
    amrex::Print() << "Neural net filename: " << m_h_nnf_data.nn_filename << std::endl;
    amrex::Print() << "Model name: " << m_h_nnf_data.model_name << std::endl;
    amrex::Print() << "Ndim: " << m_h_nnf_data.Ndim << std::endl;
    amrex::Print() << "Nvar: " << m_h_nnf_data.Nvar << std::endl;
    amrex::Print() << "Nmanpar: " << m_h_nnf_data.Nmanpar << std::endl;
    amrex::Print() << std::endl;

    amrex::Print() << "Dimnames: Index | Variable Name" << std::endl;
    for(int i = 0; i < m_h_nnf_data.Ndim; i++)
    {
      std::string dimname(&m_h_nnf_data.dimnames[i*m_h_nnf_data.len_str], m_h_nnf_data.len_str);
      amrex::Print() << i << " | " << amrex::trim(dimname) << std::endl;
    }
    amrex::Print() << std::endl;

    amrex::Print() << "Varnames: Index | Variable Name" << std::endl;
    for(int i = 0; i < m_h_nnf_data.Nvar; i++)
    {
      std::string varname(&m_h_nnf_data.varnames[i*m_h_nnf_data.len_str], m_h_nnf_data.len_str);
      amrex::Print() << i << " | " << amrex::trim(varname) << std::endl;
    }
    amrex::Print() << std::endl;

    amrex::Print() << "Ncomb: " << m_h_nnf_data.Ncomb << std::endl;

    for(int i = 0; i < m_h_nnf_data.Ndim; i++)
    {
      std::string dimname(&m_h_nnf_data.dimnames[i*m_h_nnf_data.len_str], m_h_nnf_data.len_str);
      amrex::Print() << "Defn. " + amrex::trim(dimname) + ": Variable | Source Term | Coeff." << std::endl;
      for(int j = 0; j < m_h_nnf_data.Ncomb; j++)
      {
        const int& idx = m_h_nnf_data.comb_idx[i*m_h_nnf_data.Ncomb + j];
        const int& src_idx = m_h_nnf_data.comb_src_idx[i*m_h_nnf_data.Ncomb + j];
        const amrex::Real& coeff = m_h_nnf_data.comb_coeff[i*m_h_nnf_data.Ncomb + j];
        std::string varname(&m_h_nnf_data.varnames[idx*m_h_nnf_data.len_str], m_h_nnf_data.len_str);
        std::string src_varname;
        if(src_idx >= 0)
        {
          src_varname = std::string(&m_h_nnf_data.varnames[src_idx*m_h_nnf_data.len_str],
              m_h_nnf_data.len_str);
        }
        else
        {
          src_varname = "None";
        }
        amrex::Print() << amrex::trim(varname) << " | " << amrex::trim(src_varname) << " | " << coeff << std::endl;
      }

      amrex::Print() << "Bias: " << m_h_nnf_data.man_bias[i] << std::endl;
      amrex::Print();

      amrex::Print() << std::endl;
    }
  }

  virtual void allocate()
  {
    if(!m_device_allocated)
    {
      m_d_nnf_data = static_cast<NNFuncData*>(amrex::The_Device_Arena()->alloc(sizeof(m_h_nnf_data)));
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
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_src_idx);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_idx);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_coeff);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.man_bias);

    if(m_device_allocated)
    {
      amrex::The_Device_Arena()->free(m_d_nnf_data);
    }
  }

  virtual void sync_to_device()
  {
    if(!m_device_allocated)
    {
      amrex::Abort("Device params not allocated yet");
    }
    else
    {
      amrex::Gpu::copy(amrex::Gpu::hostToDevice, &m_h_nnf_data, &m_h_nnf_data + 1,
                       m_d_nnf_data);
    }
  }

  struct NNFuncData: ManFuncParams::ManFuncData
  {
    amrex::Real* nnrdata; // Real data representing the NNModel
    int* nnidata; // int data representing the NNModel
    char nn_filename[len_str+1]; // Path to neural network file
    int Ncomb; // Number of variables that the manifold params are linear combinations of
    amrex::Real* comb_coeff; // Coefficients in linear combination for each manifold param
    int* comb_idx; // Index of variable corresponding to each comb_coeff
    int* comb_src_idx; // Index of source term for variable corresponding to each comb_coeff
    amrex::Real* man_bias; // Offsets to add when calculating manifold parameters
  };

  virtual ManFuncData& host_manfunc_data()
  {
    return m_h_nnf_data;
  }

  virtual ManFuncData* device_manfunc_data()
  {
    return m_d_nnf_data;
  }

  NNFuncData& host_nnfunc_data()
  {
    return m_h_nnf_data;
  }

  NNFuncData* device_nnfunc_data()
  {
    return m_d_nnf_data;
  }

private:

  int m_verbose = 0;
  NNFuncData m_h_nnf_data;
  NNFuncData* m_d_nnf_data;
  bool m_device_allocated{false};

};

class NNFunc: public ManifoldFunc
{
public:

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc() {}

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNFunc(NNFuncParams::NNFuncData* nnf_data_in)
      : nnf_data{nnf_data_in}, nnmodel(nnf_data_in->nnrdata, nnf_data_in->nnidata)
  {

  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_value(const int ivar, const amrex::Real indata[], amrex::Real& out)
  {
    // amrex::Print() << indata[0] << " " << indata[1] << std::endl;
    out = nnmodel(indata)[ivar];
    // amrex::Print() << out << std::endl;
    // amrex::Print() << std::endl;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_values(const int nvar, const int ivar[], const amrex::Real indata[],
      amrex::Real out[])
  {
    const amrex::Real* outdata = nnmodel(indata);
    for(int i = 0; i < nvar; i++)
    {
      out[i] = (ivar[i] > 0) ? outdata[ivar[i]] : 0.0;  
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_derivs(const int ivar, const amrex::Real indata[], amrex::Real derivs[])
  {
    // Use first-order finite difference approach for now
    constexpr amrex::Real eps = 1e-5;
    const amrex::Real outval0 = nnmodel(indata)[ivar];
    amrex::Real diff;
    
    amrex::Real indata_copy[nnf_data->Ndim];
    for(int i = 0; i < nnf_data->Ndim; i++)
    {
      indata_copy[i] = indata[i];
    }
    
    for(int i = 0; i < nnf_data->Ndim; i++)
    {
      diff = (indata[i] != 0.0) ? eps*indata[i] : eps;
      indata_copy[i] = indata[i] + diff;
      derivs[i] = (nnmodel(indata_copy)[ivar] - outval0) / diff;
      indata_copy[i] = indata[i];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void calculate_Wdot(const int paramidx, const amrex::Real indata[], amrex::Real& out)
  {
    const amrex::Real* outdata = nnmodel(indata);
    
    amrex::Real rr = 0.0;

    for(int j = 0; j < nnf_data->Ncomb; j++)
    {
      const int& idx = nnf_data->comb_src_idx[paramidx*nnf_data->Ncomb + j];
      const amrex::Real& coeff = nnf_data->comb_coeff[paramidx*nnf_data->Ncomb + j];
      const amrex::Real src_term = (idx < 0) ? 0.0 : outdata[idx];
      rr += coeff*src_term;
    }

    out = rr;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual ManifoldModel model()
  {
    return nnf_data->manmodel;
  }

private:

  // Variables
  NNFuncParams::NNFuncData* nnf_data;
  NNModel nnmodel;

}; // class NNFunc
  
} // namespace physics
} // namespace pele

#endif
