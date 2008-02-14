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

#include "clx/cptRecognizer.hh"
#include "clx/cptPlexFamily.hh"
#include "clx/clxUnit.hh"

namespace clx
{
  cptPlexFamily*
  cptRecognizer::
  makePlexFamily(const cptPlex& rPlex) const
  {
    return new cptPlexFamily(rPlex,
			     rClxUnit.bindingFeatures,
			     rClxUnit.omniPlexFamilies,
			     rCptUnit,
			     rGraph);
  }
  
  class insertFamilySpecies :
    public std::unary_function<std::map<int, cptPlexFamily*>::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesElt;
  public:
    insertFamilySpecies(xmlpp::Element* pExplicitSpeciesElement) :
      pExplicitSpeciesElt(pExplicitSpeciesElement)
    {}

    void
    operator()(const argument_type& rHasherEntry) const
      throw(std::exception)
    {
      const cptPlexFamily* pFamily = rHasherEntry.second;
      pFamily->insertSpecies(pExplicitSpeciesElt);
    }
  };

  void
  cptRecognizer::
  insertSpecies(xmlpp::Element* pExplicitSpeciesElt) const
    throw(std::exception)
  {
    std::for_each(plexHasher.begin(),
		  plexHasher.end(),
		  insertFamilySpecies(pExplicitSpeciesElt));
  }
}
