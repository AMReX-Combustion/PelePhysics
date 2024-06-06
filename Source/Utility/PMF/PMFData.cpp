#include <PMFData.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

static std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

static bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

namespace pele::physics {

void
InitParm<PMF::DataContainer>::read_pmf(
  const std::string& fname,
  const int a_doAverage,
  const int /*a_verbose*/,
  PMF::DataContainer& h_pmf_data)
{
  std::ifstream infile(fname);
  if (!infile.is_open()) {
    amrex::Abort("Unable to open pmf input file " + fname);
  }
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::string firstline;
  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::string secondline;
  std::getline(iss, secondline);
  unsigned long pos1 = 0;
  unsigned long pos2 = 0;
  int variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  pos1 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  // Check variable names
  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  if (variable_count != (NUM_SPECIES + 4)) {
    amrex::Abort("PMF file must have NUM_SPECIES+4 variables");
  }

  int line_count = 0;
  std::string remaininglines;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  h_pmf_data.m_nPoint = line_count;
  h_pmf_data.m_nVar = variable_count - 1;
  const int sizeYvec = line_count * (variable_count - 1);
  h_pmf_data.m_doAverage = a_doAverage;
  h_pmf_data.pmf_X = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
    line_count * sizeof(amrex::Real));
  h_pmf_data.pmf_Y = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
    sizeYvec * sizeof(amrex::Real));

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < h_pmf_data.m_nPoint; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> h_pmf_data.pmf_X[i];
    for (int j = 0; j < h_pmf_data.m_nVar; j++) {
      sinput >> h_pmf_data.pmf_Y[j * h_pmf_data.m_nPoint + i];
    }
  }
}
} // namespace pele::physics
