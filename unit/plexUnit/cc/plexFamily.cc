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

#include <algorithm>
#include "sampleDist/sampleDist.hh"
#include "mol/molState.hh"
#include "mol/mol.hh"
#include "plex/plexFamily.hh"
#include "plex/prm.hh"
#include "plex/plexQuery.hh"

namespace plx
{
  // Constructs minimal plexFamily for a given plex.
  plexFamily::plexFamily(const plex& rParadigm,
			 plexUnit& refPlexUnit) :
    paradigm(rParadigm),
    rPlexUnit(refPlexUnit)
  {}

  class zeroUpdateSpecies : public
  std::unary_function<plexFamily::value_type, void>
  {
    std::set<mzr::reaction*>& rAffectedReactions;
  public:
    zeroUpdateSpecies(std::set<mzr::reaction*>& rAffected) :
      rAffectedReactions(rAffected)
    {}
    
    void
    operator()(const argument_type& rParamSpeciesPair) const
    {
      plexSpecies* pPlexSpecies = rParamSpeciesPair.second;
      pPlexSpecies->update(rAffectedReactions,
			   0);
    }
  };

  plexSpecies*
  plexFamily::makeMember(const std::vector<bnd::molParam>& rMolParams)
  {
    plexSpecies* pNewSpecies = new plexSpecies(*this,
					       allostery(rMolParams));

    std::cerr << "plexFamily::makeMember: plexFamily "
	      << this
	      << " finished making plexSpecies "
	      << pNewSpecies
	      << std::endl;
    
    return pNewSpecies;
  }

  namespace
  {
    class omniHasParent :
      public std::unary_function<omniPlex*, bool>
    {
      xmlpp::Node* pParent;
    public:
      omniHasParent(xmlpp::Node* pParentNode) :
	pParent(pParentNode)
      {}

      bool
      operator()(const omniPlex* pOmniPlex) const
      {
	return pParent == (pOmniPlex->getParentNode());
      }
    };
  }
  omniPlex*
  plexFamily::
  mustFindOmniForNode(xmlpp::Node* pParentNode) const
    throw(mzr::mzrXcpt)
  {
    // Search through the omniPlexes in this family,
    // looking for the one that was found by the plexUnit
    // under the given parent node.
    mzr::autoVector<omniPlex>::const_iterator iOmni
      = std::find_if(omniPlexes.begin(),
		     omniPlexes.end(),
		     omniHasParent(pParentNode));

    // Throw exception if none was found.
    if(omniPlexes.end() == iOmni)
      throw noOmniForNodeXcpt(pParentNode);

    // Return the first (and should be only) one that was found.
    return *iOmni;
  }

  class insertPlexSpecies :
    public std::unary_function<plexFamily::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesElt;
    double molFact;
    
  public:
    insertPlexSpecies(xmlpp::Element* pExplicitSpeciesElement,
		      double molarFactor) :
      pExplicitSpeciesElt(pExplicitSpeciesElement),
      molFact(molarFactor)
    {}

    void
    operator()(const argument_type& rPlexFamilyEntry) const
      throw(std::exception)
    {
      plexSpecies* pSpecies = rPlexFamilyEntry.second;
      pSpecies->insertElt(pExplicitSpeciesElt,
			  molFact);
    }
  };

  void
  plexFamily::insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
			    double molarFactor) const
    throw(std::exception)
  {
    std::for_each(begin(),
		  end(),
		  insertPlexSpecies(pExplicitSpeciesElt,
				    molarFactor));
  }

//   // Predicate used in getOmniImbeddings below.
//   class omniFeatureInFamily :
//     public std::unary_function<mzr::featureMap<plexSpecies, subPlexSpec>::value_type, bool>
//   {
//     const plexFamily* pFamily;
//   public:
//     omniFeatureInFamily(const plexFamily* pOmniFamily) :
//       pFamily(pOmniFamily)
//     {
//     }

//     bool
//     operator()(const argument_type& rFeatureMapEntry) const
//     {
//       const subPlexSpec& rSubPlexSpec = rFeatureMapEntry.first;
//       return pFamily == rSubPlexSpec.getFamilyPtr();
//     }
//   };

  class selectImbeddingsOfOmni :
    public std::unary_function<mzr::featureMap<plexSpecies, subPlexSpec>::value_type,
	  void>
  {
    const plexFamily* pOmniFamily;
    std::vector<subPlexSpec>& rResult;

  public:
    selectImbeddingsOfOmni(const plexFamily* pSearchOmniFamily,
			   std::vector<subPlexSpec>& rSelectedImbeddings) :
      pOmniFamily(pSearchOmniFamily),
      rResult(rSelectedImbeddings)
    {}

    void
    operator()(const argument_type& rFeatureMapEntry) const
    {
      const subPlexSpec& rSubPlexSpec = rFeatureMapEntry.first;

      // If the omniplex feature map entry points to the omniplex family
      // we're hunting for, add the imbedding to the result.
      const omniPlex* pOmni = rSubPlexSpec.getOmni();
      
      if(pOmniFamily == pOmni->getFamily())
	rResult.push_back(rSubPlexSpec);
    }
  };

  // This may be obsolete, or due for upgrade.  It's used in nucExRxnGen,
  // and we should have the omniPlex, rather than its associated family.
  std::vector<subPlexSpec>
  plexFamily::getOmniImbeddings(const plexFamily* pOmniFamily) const
  {
    std::vector<subPlexSpec> result;

    std::for_each(omniFeatures.begin(),
		  omniFeatures.end(),
		  selectImbeddingsOfOmni(pOmniFamily,
					 result));
    return result;
  }
}
