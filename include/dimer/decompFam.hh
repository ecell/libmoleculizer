#ifndef DECOMPFAM_H
#define DECOMPFAM_H

#include "mzr/reactionFamily.hh"
#include "dimer/decompRxnGen.hh"

namespace dimer
{
  class decompFam :
    public mzr::reactionFamily
  {
    decompRxnGen decompGen;

  public:
    decompFam(mzr::mzrUnit& rMzrUnit,
	      plx::plexUnit& rPlexUnit,
	      decomposeExtrapolator* pExtrap) :
      decompGen(this,
		rMzrUnit,
		rPlexUnit,
		pExtrap)
    {}

    decompRxnGen*
    getRxnGen(void)
    {
      return &decompGen;
    }
  };
}

#endif
