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

#ifndef MOLUNIT_H
#define MOLUNIT_H

#include "domUtils/domUtils.hh"
#include "mzr/util.hh"
#include "mzr/mzrUnit.hh"
#include "mol/modification.hh"
#include "mol/mol.hh"
// For unknownModXcpt.
#include "mol/molDomParse.hh"
#include "mol/molEltName.hh"

namespace bnd
{
  class molUnit : public mzr::unit
  {
    class getValueMod :
      public std::unary_function
    <std::pair<std::string, std::string>,
      std::pair<std::string, const modification*> >
    {
      molUnit& rMolUnit;
      
    public:
      getValueMod(molUnit& refMolUnit) :
	rMolUnit(refMolUnit)
      {}
      
      std::pair<std::string, const modification*>
      operator()(const std::pair<std::string, std::string>& rEntry) const
      {
	return std::make_pair(rEntry.first,
			      rMolUnit.mustGetMod(rEntry.second));
      }
    };
    
    mzr::autoCatalog<const modification> knownMods;

    // Sites an mols are purely user objects: they are created only by
    // the user, and they are destroyed when the userData goes away.
    mzr::autoCatalog<mol> molsByName;

  public:

    bool
    addMol(const std::string& rMolName,
	   mol* pMol)
    {
      return molsByName.addEntry(rMolName,
				 pMol);
    }

    mol*
    findMol(const std::string& rMolName)
    {
      return molsByName.findEntry(rMolName);
    }

    mol*
    mustFindMol(xmlpp::Node* pRequestingNode,
		const std::string& rMolName)
      throw(unknownMolXcpt);

    modMol*
    mustFindModMol(xmlpp::Node* pRequestingNode,
		   const std::string& rModMolName)
      throw(unknownMolXcpt, badModMolCastXcpt);

    bool
    addMod(const modification* pModification)
    {
      return knownMods.addEntry(pModification->name,
				pModification);
    }

    void
    mustAddParsedMod(xmlpp::Node* pModNode) throw(std::exception)
    {
      modification* pMod = new modification(pModNode);

      if(! addMod(pMod))
	throw(duplicateModNameXcpt(pModNode,
				   pMod->name));
    }

    // Unification of modifications by name.  This "intern" style
    // of operation avoids the possibility of overwriting an existing
    // modification in the catalog.
//     bool
//     internMod(const std::string& rModName,
// 	      double molWeightDelta,
// 	      const modification*& rModificationPointer)
//     {
//       rModificationPointer = knownMods.findEntry(rModName);

//       bool modWasKnown = (0 != rModificationPointer);
//       if(! modWasKnown)
// 	{
// 	  rModificationPointer = new modification(rModName,
// 						  molWeightDelta);
// 	  knownMods.addEntry(rModificationPointer->name,
// 			     rModificationPointer);
// 	}
//       return modWasKnown;
//     }

    // Same as above; discards pointer to found or created modifier.
//     bool
//     internMod(const std::string& rModName,
// 	      double molWeightDelta)
//     {
//       // Pointer to discard.
//       const modification* pMod = 0;

//       return internMod(rModName,
// 		       molWeightDelta,
// 		       pMod);
//     }

    // Gets modification with the given name, or returns null.
    const modification*
    getMod(const std::string& rModName)
    {
      return knownMods.findEntry(rModName);
    }

    // Get modification with the given name, or throws unknownModXcpt.
    // This is used in getModMap below, where there is no DOM parsing
    // context.
    const modification*
    mustGetMod(const std::string& rModificationName) throw(unknownModXcpt)
    {
      const modification* pMod = getMod(rModificationName);
      if(! pMod) throw unknownModXcpt(rModificationName);
      return pMod;
    }

    // Do I really want tons of xmlpp-dependent routines of this
    // kind to assist with parsing?
    const modification*
    mustGetMod(xmlpp::Node* pNodeWithName,
	       const std::string& rModificationName)
      throw(unknownModXcpt)
    {
      const modification* pMod = getMod(rModificationName);
      if(0 == pMod) throw unknownModXcpt(pNodeWithName,
					 rModificationName);
      return pMod;
    }

    // Converts a map from modification site names to modification names
    // into a map from modification site names to modification pointers.
    void
    getModMap(const std::map<std::string, std::string>& rMapToModNames,
	      std::map<std::string, const modification*>& rMapToMods)
    {
      transform(rMapToModNames.begin(),
		rMapToModNames.end(),
		inserter(rMapToMods,
			 rMapToMods.begin()),
		getValueMod(*this));
    }

    mzr::mzrUnit& rMzrUnit;
  
    molUnit(mzr::moleculizer& rMoleculizer,
	    mzr::mzrUnit& refMzrUnit) :
      mzr::unit("mol",
		rMoleculizer),
      rMzrUnit(refMzrUnit)
    {
      // Register the model elements that this unit is
      // responsible for.
      inputCap.addModelContentName(eltName::modifications);
      inputCap.addModelContentName(eltName::mols);

      // This unit isn't responsible for any reaction generators,
      // explicit species, species streams, or events.
    }

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement) throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt) throw(std::exception);
  };
}

#endif // MOLUNIT_H
