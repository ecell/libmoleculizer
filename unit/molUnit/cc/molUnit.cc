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

#include "mol/molUnit.hh"

namespace bnd
{
  class insertModElt :
    public std::unary_function<mzr::autoCatalog<const modification>::value_type, void>
  {
    xmlpp::Element* pModificationsElt;
  public:
    insertModElt(xmlpp::Element* pModificationsElement) :
      pModificationsElt(pModificationsElement)
    {}

    void
    operator()(const argument_type& rEntry) const throw(std::exception)
    {
      rEntry.second->insertElt(pModificationsElt);
    }
  };

  class insertMolElt :
    public std::unary_function<mzr::autoCatalog<mol>::value_type, void>
  {
    xmlpp::Element* pMolsElt;
  public:
    insertMolElt(xmlpp::Element* pMolsElement) :
      pMolsElt(pMolsElement)
    {}

    void
    operator()(const argument_type& rEntry) const throw(std::exception)
    {
      rEntry.second->insertElt(pMolsElt);
    }
  };

  void
  molUnit::insertStateElts(xmlpp::Element* pRootElt)
    throw(std::exception)
  {
    // Get the model element
    xmlpp::Element* pModelElt
      = domUtils::mustGetUniqueChild(pRootElt,
				     mzr::eltName::model);

    // Get the units-states element.
    xmlpp::Element* pUnitsStatesElt
      = domUtils::mustGetUniqueChild(pModelElt,
				     mzr::eltName::unitsStates);
    
    // Insert the modifications element.
    xmlpp::Element* pModificationsElt
      = pUnitsStatesElt->add_child(eltName::modifications);

    // Have each modification insert itself into the modifications section.
    std::for_each(knownMods.begin(),
		  knownMods.end(),
		  insertModElt(pModificationsElt));

    // Insert the mols element.
    xmlpp::Element* pMolsElt
      = pUnitsStatesElt->add_child(eltName::mols);

    // Have each mol insert itself into the mols section.
    // At this time, there is actually only one kind of mol.
    std::for_each(molsByName.begin(),
		  molsByName.end(),
		  insertMolElt(pMolsElt));
  }

  mol*
  molUnit::mustFindMol(xmlpp::Node* pRequestingNode,
		       const std::string& rMolName)
    throw(unknownMolXcpt)
  {
    mol* pMol = findMol(rMolName);
    if(0 == pMol) throw unknownMolXcpt(pRequestingNode,
				       rMolName);
    return pMol;
  }

  modMol*
  molUnit::mustFindModMol(xmlpp::Node* pRequestingNode,
			  const std::string& rModMolName)
    throw(unknownMolXcpt, badModMolCastXcpt)
  {
    return mustBeModMolPtr(pRequestingNode,
			   mustFindMol(pRequestingNode,
				       rModMolName));
  }

  smallMol*
  molUnit::mustFindSmallMol(xmlpp::Node* pRequestingNode,
			    const std::string& rSmallMolName)
    throw(unknownMolXcpt, badSmallMolCastXcpt)
  {
    return mustBeSmallMolPtr(pRequestingNode,
			     mustFindMol(pRequestingNode,
					 rSmallMolName));
  }
}
