#include "utl/xcpt.hh"

namespace utl
{

  class modelAlreadyLoaded :
    public xcpt
  {
  public:
    modelAlreadyLoaded()
      :
      xcpt("Error.  This class already has a loaded model.")
    {}
  };
}
