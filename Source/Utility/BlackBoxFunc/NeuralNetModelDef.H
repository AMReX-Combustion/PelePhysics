#ifndef _NEURAL_NET_MD_H_
#define _NEURAL_NET_MD_H_

#define MAXD_NETWORK 100
#define num_layers 7

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "NeuralNetLayerDef.H"

namespace pele {
namespace physics {

AMREX_FORCE_INLINE
void
read_fmt(
  const std::string& fmt,
  std::vector<char>& types,
  std::vector<std::vector<int>>& sizes)
{
  constexpr int READING_DTYPE = 0;
  constexpr int READING_PAREN = 1;
  constexpr int READING_SIZE = 2;
  constexpr int READING_DELIM = 3;
  constexpr int TRAP_STATE = 4;

  int state = READING_DTYPE;

  std::istringstream is(fmt);
  int char_idx = -1;

  while (is.peek() != EOF && is.peek() != '\0') {
    switch (state) {
    case READING_DTYPE:
      char c0;
      is >> c0;
      types.push_back(c0);
      state = READING_PAREN;
      break;
    case READING_PAREN:
      char c1;
      is >> c1;
      if (c1 != '(') {
        state = TRAP_STATE;
      } else {
        state = READING_SIZE;
      }
      break;
    case READING_SIZE:
      int i;
      is >> i;
      if (sizes.size() < types.size()) {
        sizes.push_back(std::vector<int>());
      }
      sizes[sizes.size() - 1].push_back(i);
      state = READING_DELIM;
      break;
    case READING_DELIM:
      char c3;
      is >> c3;
      if (c3 == ',') {
        state = READING_SIZE;
      } else if (c3 == ')') {
        state = READING_DTYPE;
      } else {
        state = TRAP_STATE;
      }
      break;
    default:
      std::string msg = "Error parsing neural net layer format. At character '";
      msg.append(std::string(1, fmt[char_idx]));
      msg.append("' and index " + std::to_string(char_idx) + ".");
      amrex::Abort(msg);
    }

    char_idx++;
  }
}

class NNModel
{
public:
  AMREX_GPU_HOST_DEVICE
  NNModel() {}

  NNModel(std::string& fname, bool cmlm_net = false) : cmlm_net{cmlm_net}
  {
    std::ifstream file(fname, std::ios::binary | std::ios::in);
    std::vector<NNLayer*> layer_vec;

    if (!file.is_open()) {
      throw std::runtime_error("Unable to open input file " + fname + " .");
    }

    file.seekg(0, std::ios::beg);

    int ssize, fsize, isize;
    file.read(reinterpret_cast<char*>(&ssize), sizeof(int));
    file.read(reinterpret_cast<char*>(&fsize), sizeof(int));
    file.read(reinterpret_cast<char*>(&isize), sizeof(int));
    file.read(reinterpret_cast<char*>(&inpsize), sizeof(int));

    while (file.peek() != EOF) {
      char name_arr[ssize];
      char fmt_arr[ssize];
      file.read(&name_arr[0], ssize * sizeof(char));
      file.read(&fmt_arr[0], ssize * sizeof(char));
      std::string name =
        amrex::trim(std::string(&name_arr[0], ssize), std::string(1, '\0'));
      std::string fmt =
        amrex::trim(std::string(&fmt_arr[0], ssize), std::string(1, '\0'));

      std::vector<char> types;
      std::vector<std::vector<int>> sizes;
      read_fmt(fmt, types, sizes);

      if (name == "Linear") {
        layer_vec.push_back(
          new LinearNNLayer(types, sizes, file, fsize, isize));
      } else if (name == "LeakyReLU") {
        layer_vec.push_back(
          new LeakyReluNNLayer(types, sizes, file, fsize, isize));
      } else if (name == "BatchNorm1d") {
        layer_vec.push_back(
          new BatchNorm1dNNLayer(types, sizes, file, fsize, isize));
      }
    }

    // num_layers = layer_vec.size();

    layers = new NNLayer*[num_layers];
    for (int i = 0; i < num_layers; i++) {
      layers[i] = layer_vec[i];
    }

    lsizes = new int[num_layers];
    lsizes[0] = (layers[0]->nout == -1) ? inpsize : layers[0]->nout;
    for (int i = 1; i < num_layers; i++) {
      lsizes[i] = (layers[i]->nout == -1) ? lsizes[i - 1] : layers[i]->nout;
    }

    int maxsize = 0;
    for (int i = 0; i < num_layers; i++) {
      layers[i]->nout = lsizes[i];
      maxsize = (lsizes[i] > maxsize) ? lsizes[i] : maxsize;
    }

    outvec0 = new amrex::Real[maxsize];
    outvec1 = new amrex::Real[maxsize];

    file.close();

    AMREX_ALWAYS_ASSERT(num_layers > 0);
    if (cmlm_net) {
      AMREX_ALWAYS_ASSERT(num_layers == 7);
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNModel(amrex::Real* nnrd_in, int* nnid_in)
  {
    cmlm_net = static_cast<bool>(nnid_in[0]);

    if (cmlm_net) {
      // num_layers = 7;
      layers = new NNLayer*[num_layers];
      lsizes = new int[num_layers];

      amrex::Real* nnrd = nnrd_in;
      int* nnid = nnid_in + 1; // Already read cmlm_net flag

      inpsize = nnid[0];
      nnid++;

      layers[0] = new LinearNNLayer(nnrd, nnid);
      nnrd += layers[0]->nfparams();
      nnid += layers[0]->niparams();

      layers[1] = new LeakyReluNNLayer(nnrd, nnid);
      nnrd += layers[1]->nfparams();
      nnid += layers[1]->niparams();

      layers[2] = new BatchNorm1dNNLayer(nnrd, nnid);
      nnrd += layers[2]->nfparams();
      nnid += layers[2]->niparams();

      layers[3] = new LinearNNLayer(nnrd, nnid);
      nnrd += layers[3]->nfparams();
      nnid += layers[3]->niparams();

      layers[4] = new LeakyReluNNLayer(nnrd, nnid);
      nnrd += layers[4]->nfparams();
      nnid += layers[4]->niparams();

      layers[5] = new BatchNorm1dNNLayer(nnrd, nnid);
      nnrd += layers[5]->nfparams();
      nnid += layers[5]->niparams();

      layers[6] = new LinearNNLayer(nnrd, nnid);

      lsizes[0] = (layers[0]->nout == -1) ? inpsize : layers[0]->nout;
      for (int i = 1; i < num_layers; i++) {
        lsizes[i] = (layers[i]->nout == -1) ? lsizes[i - 1] : layers[i]->nout;
      }

      int maxsize = 0;
      for (int i = 0; i < num_layers; i++) {
        layers[i]->nout = lsizes[i];
        maxsize = (lsizes[i] > maxsize) ? lsizes[i] : maxsize;
      }

      AMREX_ASSERT(maxsize < MAXD_NETWORK);
      outvec0 = new amrex::Real[MAXD_NETWORK];
      outvec1 = new amrex::Real[MAXD_NETWORK];
    } else {
      amrex::Abort("Only CMLM nets supported for now without owning data");
    }
  }

  AMREX_GPU_HOST_DEVICE
  ~NNModel()
  {
    for (int i = 0; i < num_layers; i++) {
      delete layers[i];
    }
    delete[] layers;
    delete[] lsizes;
    delete[] outvec0;
    delete[] outvec1;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  const amrex::Real* operator()(const amrex::Real* inputs)
  {
    amrex::Real* invec = outvec0;
    amrex::Real* outvec = outvec1;
    amrex::Real* temp;
    int nin = inpsize;

    (*layers[0])(inputs, outvec, nin);
    nin = lsizes[0];

    for (int i = 1; i < num_layers; i++) {
      // Swap input and output buffers
      temp = invec;
      invec = outvec;
      outvec = temp;
      // Evaluate layer
      (*layers[i])(invec, outvec, nin);
      nin = lsizes[i];
    }

    return outvec;
  }

  void print()
  {
    amrex::Print() << "Number of layers: " << num_layers << std::endl;
    amrex::Print() << "Input size: " << inpsize << std::endl;
    amrex::Print() << std::endl;
    amrex::Print() << "Outgoing Sizes:" << std::endl;
    for (int i = 0; i < num_layers; i++) {
      amrex::Print() << layers[i]->nout << std::endl;
    }
    amrex::Print() << std::endl;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  int nin() { return inpsize; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  int nout() { return lsizes[num_layers - 1]; }

  void buffer_sizes_for_packing(int& nreals_in, int& nints_in)
  {
    // Do not need to store Real parameters for net, just the layers
    int nreals = 0;
    // We want to store the cmlm_net flag and the input size
    int nints = 2;

    if (cmlm_net) {
      for (int i = 0; i < num_layers; i++) {
        nreals += layers[i]->nfparams();
        nints += layers[i]->niparams();
      }
    }

    nreals_in = nreals;
    nints_in = nints;
  }

  void pack(amrex::Real* nnrd, int* nnid)
  {
    int nreals = 0;
    int nints = 0;

    if (cmlm_net) {
      // Store cmlm_net flag and input size
      nnid[0] = static_cast<int>(cmlm_net);
      nints++;
      nnid[1] = inpsize;
      nints++;

      // Pack the layers into the arrays
      for (int i = 0; i < num_layers; i++) {
        layers[i]->pack(nnrd + nreals, nnid + nints);
        nreals += layers[i]->nfparams();
        nints += layers[i]->niparams();
      }
    }
  }

private:
  NNLayer** layers;
  int* lsizes;
  int inpsize;
  // static constexpr int num_layers = 7;
  amrex::Real* outvec0;
  amrex::Real* outvec1;
  bool cmlm_net;
};

} // namespace physics
} // namespace pele

#endif
