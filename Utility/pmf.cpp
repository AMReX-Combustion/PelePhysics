#include <pmf.H>

namespace PMF
{
  AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N;
  AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M;
  AMREX_GPU_DEVICE_MANAGED bool pmf_do_average;

  amrex::Gpu::ManagedVector<amrex::Real> pmf_X;
  amrex::Gpu::ManagedVector<amrex::Real> pmf_Y;

  amrex::Vector<std::string> pmf_names;
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
PMF::read_pmf(const std::string& myfile)
{
  std::string firstline, secondline, remaininglines;
  int pos1, pos2;
  int variable_count, line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline))
    amrex::Abort("PMF file variable quotes unbalanced");
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

  PMF::pmf_N = line_count;
  PMF::pmf_M = variable_count - 1;
  PMF::pmf_X.resize(PMF::pmf_N);
  PMF::pmf_Y.resize(PMF::pmf_N * PMF::pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < PMF::pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PMF::pmf_X[i];
    for (int j = 0; j < PMF::pmf_M; j++) {
      sinput >> PMF::pmf_Y[j * PMF::pmf_N + i];
    }
  }
}

