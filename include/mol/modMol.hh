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

#ifndef MODMOL_H
#define MODMOL_H

#include "mol/molState.hh"
#include "mol/mol.hh"
#include "mol/modMixin.hh"
#include "plex/prm.hh"

namespace bnd
{
  class modMolState :
    public molState, public modStateMixin
  {
  public:
    modMolState(double molecularWeight,
		const modification* pModification,
		int count) :
      molState(molecularWeight),
      modStateMixin(pModification, count)
    {}

    modMolState(double molecularWeight,
		const modStateMixin& rModStateMixin) :
      molState(molecularWeight),
      modStateMixin(rModStateMixin)
    {}

    bool operator< (const modMolState& rRight) const
    {
      const molState& rThisMolState = *this;
      const molState& rRightMolState = rRight;

      if(rThisMolState < rRightMolState) return true;
      if(rRightMolState < rThisMolState) return false;

      const modStateMixin& rThisModState = *this;
      const modStateMixin& rRightModState = rRight;

      return rThisModState < rRightModState;
    }

    virtual
    ~modMolState(void)
    {}

    virtual double
    getMolWeight(void) const
    {
      return baseWeight + totalWeightDelta();
    }
  };

  class modMol :
    public alloMol<modMolState>, public modMolMixin
  {
  public:

    // Use modification::getModMap to convert a
    // map<string, string> into a map<string, const modification*>
    // as a preliminary to using this function.
    modMol(const std::string& rName,
	   const std::vector<bindingSite>& rSites,
	   double molecularWeight,
	   const std::map<std::string, const modification*>& rModMap) :
      alloMol<modMolState>(rName,
			   rSites),
      // Generates the map from modification site name to index.
      modMolMixin(rModMap)
    {
      // Uses the map from modification site name to index to
      // index the const modification*'s.
      //
      // Make this an accessor call, setDefaultState(...)?
      pDefaultState = internState(modMolState(molecularWeight,
					      indexModMap(rModMap)));
    }

    // Still using the default state to get the molecular weight.
    // 
    // Use modification::getModMap to convert a
    // map<string, string> into a map<string, const modification*>
    // as a preliminary to using this function.
    const modMolState*
    internModMap(const std::map<std::string, const modification*>& rModMap)
    {
      const modMolState& rDflt = *getDefaultState();

      // This seems like an insane way to get the base molecular weight.
      const molState& rBaseState = rDflt;
      double baseWeight = rBaseState.getMolWeight();
    
      return internState(modMolState(baseWeight,
				     substituteModMap(rModMap,
						      rDflt)));
    }

    // Returns the index of the modification site with the given name,
    // or throws a parsing exception.
    int
    mustGetModSiteNdx(xmlpp::Node* pRequestingNode,
		      const std::string& rModSiteName) const
      throw(unknownModSiteXcpt);

    // Output to state dump.
    xmlpp::Element*
    insertInstanceState(xmlpp::Element* pInstanceStatesElt,
			int molInstanceNdx,
			molParam param) const;

    // Used to generate instance names for mod-mols in complexes in
    // state dump.
    std::string
    genInstanceName(int molInstanceNdx) const;

    // Output to state dump.
    xmlpp::Element*
    insertElt(xmlpp::Element* pMolsElt) const throw(std::exception);
  };

  modMol*
  mustBeModMol(xmlpp::Node* pRequestingNode,
	       mol* pMol)
    throw(badModMolCastXcpt);

  const modMol*
  mustBeModMol(xmlpp::Node* pRequestingNode,
	       const mol* pMol)
    throw(badModMolCastXcpt);
}

#endif // MODMOL_H
