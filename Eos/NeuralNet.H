#ifndef _NEURAL_NET_H_
#define _NEURAL_NET_H_

#include <torch/script.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include <vector>
#include "ManifoldFunc.H"

namespace pele {
namespace physics {
  
namespace TensorParams
{
  // Whether we are on a GPU or not
#ifdef AMREX_USE_GPUs
  const torch::Device device = torch::kGPU;
#else
  const torch::Device device = torch::kCPU;
#endif

  // Set data type
#ifdef BL_USE_FLOAT
  const torch::ScalarType dtype = torch::kFloat32;
#else
  const torch::ScalarType dtype = torch::kFloat64;
#endif

  // Tensor creation options
  const torch::TensorOptions tensoropt = torch::TensorOptions().dtype(dtype).device(device);
}

class NNFuncParams: public ManFuncParams
{
public:
  
  NNFuncParams()
  {
#ifdef AMREX_USE_GPU
    amrex::Abort("Neural network derived manifolds do not support GPUs at this time.");
#endif
  }
  
  ~NNFuncParams() {}

  virtual void initialize()
  {
    m_h_nnf_data.manmodel = ManifoldModel::NEURAL_NET;
    
    amrex::ParmParse pp("manifold");
    
    std::string nn_filename;
    pp.get("filename", nn_filename);
    
    torch::jit::script::Module nnmodel;
    try
    {
      // Deserialize the ScriptModule from a file using torch::jit::load().
      nnmodel = torch::jit::load(amrex::trim(nn_filename));
      // Convert to proper data type
      nnmodel.to(TensorParams::dtype);
      // Move to GPU if needed
// #ifdef AMREX_USE_GPU
//       nnmodel.to(device);
// #endif
    }
    catch(const c10::Error& e)
    {
      amrex::Error("Unable to load neural network model for manifold EOS.");
    }
    
    // The nnmodel should be freed up automatically once it goes out of scope
    // The code below clears the Torch GPU cache
    // c10::cuda::CUDACachingAllocator::empty_cache();
    
    std::string info_filename;
    pp.get("info_filename", info_filename);
    
    strncpy(m_h_nnf_data.nn_filename, nn_filename.c_str(), m_h_nnf_data.len_str*sizeof(char));
    m_h_nnf_data.nn_filename[m_h_nnf_data.len_str] = '\0';
    
    m_h_nnf_data.nnmodel = nnmodel;
    
    read_metadata(info_filename);
    
    pp.query("v", m_verbose);
    if(m_verbose >= 1) print();
    
    allocate();
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
    std::istringstream is(buffer.data(), length);
    
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
        umap[cur_var] = std::vector<std::string>();
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
        std::string src_var = "SRC_" + var;
        int idx = get_var_index(var.c_str(), &m_h_nnf_data);
        int src_idx = get_var_index(src_var.c_str(), &m_h_nnf_data);
        m_h_nnf_data.comb_coeff[i*m_h_nnf_data.Ncomb + j] = (amrex::Real)std::stod(coeff);
        m_h_nnf_data.comb_idx[i*m_h_nnf_data.Ncomb + j] = idx;
        m_h_nnf_data.comb_src_idx[i*m_h_nnf_data.Ncomb + j] = src_idx;
      }
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
    
    amrex::Print() << "Ncomb: " << std::endl;
    
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
        std::string src_varname(&m_h_nnf_data.varnames[src_idx*m_h_nnf_data.len_str], m_h_nnf_data.len_str);
        amrex::Print() << amrex::trim(varname) << " | " << amrex::trim(src_varname) << " | " << coeff << std::endl;
      }
      amrex::Print() << std::endl;
    }
  }

  virtual void allocate()
  {
    // if(!m_device_allocated)
    // {
    //   m_d_nnf_data = static_cast<NNFuncData*>(amrex::The_Device_Arena()->alloc(sizeof(m_h_nnf_data)));
    //   m_device_allocated = true;
    //   sync_to_device();
    // }
  }

  virtual void deallocate()
  {
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.varnames);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.dimnames);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_src_idx);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_idx);
    amrex::The_Pinned_Arena()->free(m_h_nnf_data.comb_coeff);
    
    // if(m_device_allocated)
    // {
    //   amrex::The_Device_Arena()->free(m_d_nnf_data);
    // }
  }

  virtual void sync_to_device()
  {
    // if(!m_device_allocated)
    // {
    //   amrex::Abort("Device params not allocated yet");
    // }
    // else
    // {
    //   amrex::Gpu::copy(amrex::Gpu::hostToDevice, &m_h_nnf_data, &m_h_nnf_data + 1,
    //                    m_d_nnf_data);
    // }
  }
  
  struct NNFuncData: ManFuncParams::ManFuncData
  {
    char nn_filename[len_str+1]; // Path to neural network file
    torch::jit::script::Module nnmodel;
    int Ncomb; // Number of variables that the manifold params are linear combinations of
    amrex::Real* comb_coeff; // Coefficients in linear combination for each manifold param
    int* comb_idx; // Index of variable corresponding to each comb_coeff
    int* comb_src_idx; // Index of source term for variable corresponding to each comb_coeff
  };
  
  virtual ManFuncData& host_manfunc_data()
  {
    return m_h_nnf_data;
  }
  
  virtual ManFuncData* device_manfunc_data()
  {
    return &m_h_nnf_data;
    // return m_d_nnf_data;
  }

  NNFuncData& host_nnfunc_data()
  {
    return m_h_nnf_data;
  }

  NNFuncData* device_nnfunc_data()
  {
    return &m_h_nnf_data;
    // return m_d_nnf_data;
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
  {
    nnf_data = nnf_data_in;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_value(const int ivar, const amrex::Real indata[], amrex::Real& out)
  {
    // Ostensibly cannot use from_blob because the indata array is const
    torch::Tensor intensor = torch::empty({1, nnf_data->Ndim}, TensorParams::tensoropt);
    std::memcpy(intensor.data_ptr(), indata, sizeof(amrex::Real) * intensor.numel());
    get_value(ivar, intensor, out);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_derivs(const int ivar, const amrex::Real indata[], amrex::Real derivs[])
  {
    // Ostensibly cannot use from_blob because the indata array is const
    torch::Tensor intensor = torch::empty({1, nnf_data->Ndim}, TensorParams::tensoropt.requires_grad(true));
    std::memcpy(intensor.data_ptr(), indata, sizeof(amrex::Real) * intensor.numel());
    get_derivs(ivar, intensor, derivs);
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void calculate_Wdot(const int paramidx, const amrex::Real indata[], amrex::Real& out)
  {
    // Ostensibly cannot use from_blob because the indata array is const
    torch::Tensor intensor = torch::empty({1, nnf_data->Ndim}, TensorParams::tensoropt);
    std::memcpy(intensor.data_ptr(), indata, sizeof(amrex::Real) * intensor.numel());
    calculate_Wdot(paramidx, intensor, out);
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

  // Private helper functions
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  torch::Tensor eval_model(const torch::Tensor& indata)
  {
    return nnf_data->nnmodel.forward({indata}).toTensor()[0];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  torch::Tensor eval_model(const int ivar, const torch::Tensor& indata)
  {
    return eval_model(indata)[ivar];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_value(const int ivar, const torch::Tensor& indata, amrex::Real& out)
  {
    out = eval_model(ivar, indata).item<amrex::Real>();
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_derivs(const int ivar, torch::Tensor& indata, amrex::Real derivs[])
  {
    torch::Tensor out_tensor = eval_model(ivar, indata);
    out_tensor.backward();

    auto grad = indata.grad();
    auto acc = grad.accessor<amrex::Real,1>();

    for(int i = 0; i < acc.size(0); i++)
    {
      derivs[i] = acc[i];
    }
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void calculate_Wdot(const int paramidx, const torch::Tensor& indata, amrex::Real& out)
  {
    torch::Tensor predvars = eval_model(indata);
    auto acc = predvars.accessor<amrex::Real,1>();
    
    amrex::Real rr = 0.0;
    
    for(int j = 0; j < nnf_data->Ncomb; j++)
    {
      const int& idx = nnf_data->comb_src_idx[paramidx*nnf_data->Ncomb + j];
      const amrex::Real& coeff = nnf_data->comb_coeff[paramidx*nnf_data->Ncomb + j];
      rr += coeff*acc[idx];
    }
    
    out = rr;
  }

}; // class NNFunc
}  // namespace physics
}  // namespace pele

#endif
