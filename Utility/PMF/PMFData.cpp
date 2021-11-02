#include <PMFData.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

using namespace amrex;

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

namespace pele {
namespace physics {
namespace PMF {

void
PmfData::read_pmf(const std::string& fname, int a_doAverage, int /*a_verbose*/)
{
  std::string firstline, secondline, remaininglines;
  unsigned long pos1, pos2;
  int variable_count, line_count;

  std::ifstream infile(fname);
  if (!infile.is_open()) {
    amrex::Abort("Unable to open pmf input file " + fname);
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

  pmf_names.resize(variable_count);
  pos1 = 0;
  pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  // Check variable names

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << PMF::pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  m_data_h.m_nPoint = line_count;
  m_data_h.m_nVar = variable_count - 1;
  m_data_h.m_doAverage = a_doAverage;
  m_data_h.pmf_X = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
    line_count * sizeof(amrex::Real));
  m_data_h.pmf_Y = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
    line_count * (variable_count - 1) * sizeof(amrex::Real));
  m_host_allocated = true;

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < m_data_h.m_nPoint; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> m_data_h.pmf_X[i];
    for (unsigned int j = 0; j < m_data_h.m_nVar; j++) {
      sinput >> m_data_h.pmf_Y[j * m_data_h.m_nPoint + i];
    }
  }
}
} // namespace PMF
} // namespace physics
} // namespace pele
