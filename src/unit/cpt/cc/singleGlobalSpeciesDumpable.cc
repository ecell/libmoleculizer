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

#include "cpt/cptEltName.hh"
#include "cpt/singleGlobalSpeciesDumpable.hh"
#include "cpt/globalSpeciesDumpableAux.hh"

namespace cpt
{
  void
  singleGlobalSpeciesDumpable::
  doDump(const globalDumpArg& rDumpArg) const
  {
    std::ostream& rOstream = rDumpArg.getOstream();

    if(rDumpArg.dumpTotal)
      {
	// Emit the population of the single species totaled over
	// the dumped compartments.
	int totalPop = 0;

	std::for_each(rDumpArg.dumpCompartments.begin(),
		      rDumpArg.dumpCompartments.end(),
		      addDumpedCompartmentPop(getSpecies(),
					      totalPop));

	rOstream << totalPop;
      }
    else
      {
	// Emit the populations of the single species as a tab-separated
	// list of species in each compartment.
	std::vector<int>::const_iterator 
	  iCptNdx = rDumpArg.dumpCompartments.begin();

	if(rDumpArg.dumpCompartments.end() != iCptNdx)
	  {
	    int compartmentPop
	      = getSpecies()->getCompartmentSpecies(*iCptNdx)->getPop();
	      
	    rOstream << compartmentPop;

	    while(rDumpArg.dumpCompartments.end() != ++iCptNdx)
	      {
		int compartmentPop
		  = getSpecies()->getCompartmentSpecies(*iCptNdx)->getPop();
	      
		rOstream << "\t"
			 << compartmentPop;
	      }
	  }
      }
  }

  void
  singleGlobalSpeciesDumpable::
  dumpHeader(const globalDumpArg& rDumpArg) const
  {
    std::ostream& rOstream = rDumpArg.getOstream();
    const compartmentGraph& rGraph = getSpecies()->getCompartmentGraph();

    if(rDumpArg.dumpTotal)
      {
	rOstream << getName();
      }
    else
      {
	std::vector<int>::const_iterator
	  iCptNdx = rDumpArg.dumpCompartments.begin();

	if(rDumpArg.dumpCompartments.end() != iCptNdx)
	  {
	    std::string compartmentName
	      = rGraph.compartments[*iCptNdx]->getName();

	    rOstream << getName()
		     << ":"
		     << compartmentName;

	    while(rDumpArg.dumpCompartments.end() != ++iCptNdx)
	      {
		std::string compartmentName
		  = rGraph.compartments[*iCptNdx]->getName();

		rOstream << "\t"
			 << getName()
			 << ":"
			 << compartmentName;
	      }
	  }
      }
  }

  void
  singleGlobalSpeciesDumpable::
  insertTaggedSpeciesStreamRef(xmlpp::Element* pParent) const
    throw(std::exception)
  {
    xmlpp::Element* pTaggedSpeciesStreamRefElt
      = pParent->add_child(eltName::taggedSpeciesStreamRef);

    pTaggedSpeciesStreamRefElt
      ->set_attribute(eltName::taggedSpeciesStreamRef_nameAttr,
		      getName());
  }

  // For now, this inserts just the tag of the global species.
  // Presumably, the tags of compartment species will not appear
  // in state dump.
  void
  singleGlobalSpeciesDumpable::
  insertDumpedSpeciesTags(xmlpp::Element* pParentElt) const
    throw(std::exception)
  {
    xmlpp::Element* pSpeciesRefElt
      = pParentElt->add_child(eltName::taggedSpeciesRef);

    pSpeciesRefElt->set_attribute(eltName::taggedSpeciesRef_tagAttr,
				  getSpecies()->getTag());
  }
}
