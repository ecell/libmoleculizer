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

#include <ctime>
#include "mzr/species.hh"
#include "mzr/reaction.hh"
#include "mzr/dumpable.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/moleculizer.hh"
#include "mzr/mzrEltName.hh"

namespace mzr
{
  void
  speciesDumpable::insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
    throw(std::exception)
  {
    xmlpp::Element* pTaggedSpeciesStreamRefElt
      = pParent->add_child(eltName::taggedSpeciesStreamRef);

    pTaggedSpeciesStreamRefElt
      ->set_attribute(eltName::taggedSpeciesStreamRef_nameAttr,
		      getName());
  }
  
  void
  singleSpeciesDumpable::doDump(std::ostream& rOs) const
  {
    rOs << rSpecies.getPop();
  }

  void
  singleSpeciesDumpable::insertDumpedSpeciesTags(xmlpp::Element* pParentElt)
    const throw(std::exception)
  {
    xmlpp::Element* pSpeciesRefElt
      = pParentElt->add_child(eltName::taggedSpeciesRef);

    pSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
				  rSpecies.getTag());
  }

  // Function class to sum populations of dumped species.
  class multiSpeciesDumpable::accumPop :
    public std::unary_function<species*, void>
  {
    int& rPop;

  public:
    accumPop(int& rPopulation) :
      rPop(rPopulation)
    {}

    void
    operator()(const species* pSpecies) const
    {
      rPop += pSpecies->getPop();
    }
  };

  void
  multiSpeciesDumpable::doDump(std::ostream& rOs) const
  {
    int pop = 0;
    std::for_each(dumpedSpecies.begin(),
		  dumpedSpecies.end(),
		  accumPop(pop));
    rOs << pop;
  }

  // Function class to insert tag of one of the dumped species.
  class multiSpeciesDumpable::insertTag :
    public std::unary_function<species*, void>
  {
    xmlpp::Element* pParentElt;
    
  public:
    insertTag(xmlpp::Element* pParentElement) :
      pParentElt(pParentElement)
    {}
    
    void
    operator()(const species* pSpecies) const
      throw(std::exception)
    {
      xmlpp::Element* pTaggedSpeciesRefElt
	= pParentElt->add_child(eltName::taggedSpeciesRef);

      pTaggedSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
					  pSpecies->getTag());
    }
  };

  void
  multiSpeciesDumpable::
  insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
    throw(std::exception)
  {
    std::for_each(dumpedSpecies.begin(),
		  dumpedSpecies.end(),
		  insertTag(pParentElt));
  }
  

  void
  simTimeDumpable::doDump(std::ostream& rOs) const
  {
    rOs << rMolzer.eventQ.getSimTime();
  }

  void
  clockDumpable::doDump(std::ostream& rOs) const
  {
    rOs << (((double) clock()) / ((double) CLOCKS_PER_SEC));
  }

  void
  secondsDumpable::doDump(std::ostream& rOs) const
  {
    time_t secondsNow = time(0);
    rOs << secondsNow - rMzrUnit.startSeconds;
  }

  void
  reactEventCountDumpable::doDump(std::ostream& rOs) const
  {
    rOs << reaction::reactionEventCount;
  }

  void
  reactCountDumpable::doDump(std::ostream& rOs) const
  {
    rOs << reaction::reactionCount;
  }

  void
  activationCountDumpable::doDump(std::ostream& rOs) const
  {
    rOs << reaction::activationCount;
  }

  void
  volumeDumpable::doDump(std::ostream& rOs) const
  {
    rOs << rMzrUnit.getMolarFactor().getVolume();
  }

  void
  speciesCountDumpable::doDump(std::ostream& rOs) const
  {
    rOs << species::speciesCount;
  }
}
