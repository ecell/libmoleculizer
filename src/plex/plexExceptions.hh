#ifndef PLEXEXCEPTIONS_HPP
#define PLEXEXCEPTIONS_HPP

#include "utl/xcpt.hh"

namespace plx
{
  class NonConstructableComplexOutputStateXcpt : public utl::xcpt
  {
  public:
    NonConstructableComplexOutputStateXcpt()
    :
      xcpt( "Error: Complex output state could not be constructed by the plexUnit.")
    {
    }
  };
}


#endif
