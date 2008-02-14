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

#include "utl/platform.hh"
#include "mzr/mzrReaction.hh"
#include "mzr/mzrUnit.hh"
#include "mzr/mzrEltName.hh"
#include "mzr/dupSpeciesNameXcpt.hh"
#include "stoch/stochUnit.hh"
#include "stoch/unkStochSpeciesXcpt.hh"

namespace stoch
{
  bool
  stochUnit::
  addStochSpecies(const std::string& rSpeciesName,
		  stochSpecies* pStochSpecies)
  {
   
    bool result = rMzrUnit.addSpecies(rSpeciesName,
				      pStochSpecies);

    // stochSpecies are memory-managed via stochSpeciesByName.
    return result && stochSpeciesByName.addEntry(rSpeciesName,
						 pStochSpecies);
  }

  void
  stochUnit::
  mustAddStochSpecies(const std::string& rSpeciesName,
		      stochSpecies* pStochSpecies,
		      xmlpp::Node* pRequestingNode)
    throw(utl::xcpt)
  {
    if(! addStochSpecies(rSpeciesName,
			 pStochSpecies))
      throw mzr::dupSpeciesNameXcpt(rSpeciesName,
				    pRequestingNode);
  }
  
  stochSpecies*
  stochUnit::
  findStochSpecies(const std::string& rSpeciesName)
  {
    return stochSpeciesByName.findEntry(rSpeciesName);
  }

  stochSpecies*
  stochUnit::
  mustGetStochSpecies(xmlpp::Node* pRequestingNode,
		      const std::string& rSpeciesName)
    throw(unkStochSpeciesXcpt)
  {
    stochSpecies* pSpecies = findStochSpecies(rSpeciesName);
    if(0 == pSpecies) throw unkStochSpeciesXcpt(rSpeciesName,
						pRequestingNode);
    return pSpecies;
  }

  void
  stochUnit::
  addNoSubstrateArrow(mzr::mzrReaction* pReaction)
  {
    noSubstrateArrows.push_back(pReaction);
  }

  class insertStochSpeciesElt :
    public std::unary_function<utl::autoCatalog<stochSpecies>::value_type, void>
  {
    xmlpp::Element* pTaggedSpeciesElt;
    double molFact;
    
  public:
    insertStochSpeciesElt(xmlpp::Element* pTaggedSpeciesElement,
			  double molarFactor) :
      pTaggedSpeciesElt(pTaggedSpeciesElement),
      molFact(molarFactor)
    {}

    void
    operator()(const argument_type& rCatalogEntry) const throw(std::exception)
    {
      const stochSpecies* pStochSpecies = rCatalogEntry.second;

      pStochSpecies->insertElt(pTaggedSpeciesElt,
			       molFact);
    }
  };

  void
  stochUnit::
  insertStateElts(xmlpp::Element* pRootElt) throw(std::exception)
  {
    xmlpp::Element* pModelElt
      = utl::dom::mustGetUniqueChild(pRootElt,
				     mzr::eltName::model);

    // Ensure that the tagged-species element is there.
    xmlpp::Element* pTaggedSpeciesElt
      = utl::dom::mustGetUniqueChild(pModelElt,
				     mzr::eltName::taggedSpecies);

    // Insert stoch-species nodes.
    std::for_each(stochSpeciesByName.begin(),
		  stochSpeciesByName.end(),
		  insertStochSpeciesElt(pTaggedSpeciesElt,
					rMzrUnit.getMolarFactor().getFactor()));

    // Leaving out the "no-substrate reactions" for now.  Don't know
    // why I put them in this unit to begin with anyway...
  }
}
