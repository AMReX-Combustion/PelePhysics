#include <pmf_data.H>
#include <pmf.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

using namespace amrex;

namespace PMF
{
   amrex::Vector<std::string> pmf_names;
}

void
PMF::close()
{
    if ( pmf_data_g != nullptr) {
        The_Arena()->free(pmf_data_g);
    }
}

static
std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

static
bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"')
      count++;
  }
  if ((count % 2) == 0)
    return true;
  else
    return false;
}

void
PMF::read_pmf(const std::string& myfile, bool do_average)
{

  PmfData pmf_data;

  std::string firstline, secondline, remaininglines;
  int pos1, pos2;
  int variable_count, line_count;

  std::ifstream infile(myfile);
  if (!infile.is_open()) {
    amrex::Abort("Unable to open pmf input file " + myfile);
  }
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  PMF::pmf_names.resize(variable_count);
  pos1 = 0;
  pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    PMF::pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file \n";

  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << PMF::pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  pmf_data.pmf_N = line_count;
  pmf_data.pmf_M = variable_count - 1;
  pmf_data.pmf_X = (amrex::Real *) The_Pinned_Arena()->alloc(pmf_data.pmf_N * sizeof(amrex::Real));
  pmf_data.pmf_Y = (amrex::Real *) The_Pinned_Arena()->alloc(pmf_data.pmf_N * pmf_data.pmf_M * sizeof(amrex::Real));

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < pmf_data.pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> pmf_data.pmf_X[i];
    for (int j = 0; j < pmf_data.pmf_M; j++) {
      sinput >> pmf_data.pmf_Y[j * pmf_data.pmf_N + i];
    }
  }

  pmf_data_g = (PmfData *) The_Arena()->alloc(sizeof(pmf_data));
#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(pmf_data_g,&pmf_data,sizeof(pmf_data));
#else
  pmf_data_g->pmf_N = line_count;
  pmf_data_g->pmf_M = variable_count - 1;
  pmf_data_g->pmf_do_average = do_average;
  std::memcpy(&pmf_data_g->pmf_X,&pmf_data.pmf_X, sizeof(pmf_data.pmf_X));
  std::memcpy(&pmf_data_g->pmf_Y,&pmf_data.pmf_Y, sizeof(pmf_data.pmf_Y));
#endif
}
