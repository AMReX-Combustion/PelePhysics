#include "EOS.H"

namespace EOS {

void init() 
{
  CKINIT();
}

void close()
{
  CKFINALIZE();
}

} // namespace EOS
