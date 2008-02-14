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

#ifndef CML_CMLUNIT_H
#define CML_CMLUNIT_H

#include "utl/autoVector.hh"
#include "utl/autoCatalog.hh"
#include "utl/dom.hh"
#include "cpt/cptUnit.hh"
#include "cpx/modification.hh"
#include "cml/cptSmallMol.hh"
#include "cml/cptModMol.hh"

namespace cml
{
  class cmlUnit : 
    public cpt::unit
  {
    utl::autoCatalog<const cpx::modification> knownMods;

    utl::autoCatalog<cptMol> molsByName;

  public:

    cmlUnit(cpt::cptApp& rCptApp);

    // Attempts to add the given mol to the catalog.  Fails, returning false,
    // if there is a mol with the same name already in the catalog.
    bool
    addMol(cptMol* pMol)
    {
      return molsByName.addEntry(pMol->getName(),
				 pMol);
    }

    // Attempts to add the given mol to the catalog.  Fails, throwing an
    // exception, if there is a mol with the same name already in the catalog.
    // Exception message incorporates Xpath of requesting node, if given.
    void
    mustAddMol(cptMol* pMol,
	       const xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // Looks up mol with the given name.  Returns 0 if there is no mol
    // with the given name.
    cptMol*
    findMol(const std::string& rMolName) const
    {
      return molsByName.findEntry(rMolName);
    }

    // Looks up mol with the given name.  Throws an exception if there
    // is no mol with the given name.  Exception message includes Xpath
    // of requesting node, if given.
    cptMol*
    mustFindMol(const std::string& rMolName,
		const xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    // Looks up mol with the given name.  Throws an exception if there
    // is no mol with the given name or if the mol with that name is not
    // a mzrModMol.  Exception message includes Xpath of requesting node
    // if given.
    cptModMol*
    mustFindModMol(const std::string& rModMolName,
		   const xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    // Looks up mol with the given name.  Throws an exception if there
    // is no mol with the given name or if the mol with that name is not
    // a smallMol.  Exception message includes Xpath of requesting node
    // if given.
    cptSmallMol*
    mustFindSmallMol(const std::string& rSmallMolName,
		     const xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    // Adds a modification to the catalog.  If there is a modification
    // having the same name as the given modification, fails and returns
    // false.
    bool
    addMod(const cpx::modification* pModification)
    {
      return knownMods.addEntry(pModification->getName(),
				pModification);
    }

    // Adds a modification to the catalog.  If there is a modification
    // having the same name as the given modification, fails and thows an
    // exception.  If a node is given, its Xpath is included in the
    // exception message.
    void
    mustAddMod(const cpx::modification* pModification,
	       const xmlpp::Node* pRequestingNode = 0)
      throw(utl::xcpt);

    // Gets modification with the given name, or fails by returning null.
    const cpx::modification*
    getMod(const std::string& rModName) const
    {
      return knownMods.findEntry(rModName);
    }

    // Get modification with the given name, or throws an exception.
    // If a node pointer is provided, then the Xpath of the node is
    // included in the exception message.
    const cpx::modification*
    mustGetMod(const std::string& rModificationName,
	       xmlpp::Node* pRequestingNode = 0) const
      throw(utl::xcpt);

    // Converts a map from modification site names to modification names
    // into a map from modification site names to modification pointers.
    //
    // Is this ever used?
    void
    getModMap(const std::map<std::string, std::string>& rMapToModNames,
	      std::map<std::string, const cpx::modification*>& rMapToMods) const
      throw(utl::xcpt);

    void
    parseDomInput(xmlpp::Element* pRootElement,
		  xmlpp::Element* pModelElement,
		  xmlpp::Element* pStreamsElement,
		  xmlpp::Element* pEventsElement)
      throw(std::exception);

    void
    insertStateElts(xmlpp::Element* pRootElt)
      throw(std::exception);
  };
}

#endif // CML_CMLUNIT_H
