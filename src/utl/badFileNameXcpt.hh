#include "utl/xcpt.hh"

namespace utl
{

  class badFileNameXcpt :
    public xcpt
  {
  public:
    badFileNameXcpt()
      :
      xcpt("Error: No filename specified.  Proper usage is \"-f FILENAME\".")
    {}
  };
}
