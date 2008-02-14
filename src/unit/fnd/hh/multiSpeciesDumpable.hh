/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef FND_MULTISPECIESDUMPABLE_H
#define FND_MULTISPECIESDUMPABLE_H

#include "fnd/dumpable.hh"
#include "fnd/sensitive.hh"
#include "fnd/newSpeciesStimulus.hh"

namespace fnd
{
  // This doeesn't include any output (state dump) functionality.  Decided to
  // keep xml parsing and output out of template components like this.

  template<class speciesT,
	   class dumpArgT>
  class multiSpeciesDumpable :
    public dumpable<dumpArgT>,
    public fnd::sensitive<fnd::newSpeciesStimulus<speciesT> >
  {
  protected:
    std::vector<const speciesT*> dumpedSpecies;

  public:
    multiSpeciesDumpable(const std::string& rName) :
      dumpable<dumpArgT>(rName)
    {}

    virtual
    ~multiSpeciesDumpable(void)
    {}

    virtual int
    getTotalPop(const dumpArgT& rDumpArg) const;

    virtual void
    doDump(const dumpArgT& rDumpArg) const
    {
      rDumpArg.getOstream() << getTotalPop(rDumpArg);
    }

    virtual
    void
    respond(const fnd::newSpeciesStimulus<speciesT>& rStimulus)
    {
      dumpedSpecies.push_back(rStimulus.getSpecies());
    }
  };

  template<class speciesT>
  class addDumpedSpeciesPop :
    public std::unary_function<const speciesT*, void>
  {
    int& rTotalPop;
    
  public:
    addDumpedSpeciesPop(int& refTotalPop) :
      rTotalPop(refTotalPop)
    {}
    
    void
    operator()(const speciesT* pSpecies) const
    {
      rTotalPop += pSpecies->getPop();
    }
  };

  // The automatically defined version doesn't use the dumpArg
  // in the process of getting the total population of the species.
  template<class speciesT,
	   class dumpArgT>
  int
  multiSpeciesDumpable<speciesT,
		       dumpArgT>::
  getTotalPop(const dumpArgT& rDumpArg) const
  {
    int totalPop = 0;
    std::for_each(dumpedSpecies.begin(),
		  dumpedSpecies.end(),
		  addDumpedSpeciesPop<speciesT>(totalPop));
    return totalPop;
  }
}

#endif // FND_MULTISPECIESDUMPABLE_H
