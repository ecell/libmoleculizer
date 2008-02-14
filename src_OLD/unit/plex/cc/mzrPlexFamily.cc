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

#include "plex/mzrOmniPlex.hh"
#include "plex/mzrPlexFamily.hh"

namespace plx
{
  mzrPlexFamily::
  mzrPlexFamily(const mzrPlex& rParadigm,
		cpx::knownBindings<bnd::mzrMol, fnd::feature<cpx::cxBinding<mzrPlexSpecies, mzrPlexFamily> > >& refKnownBindings,
		std::set<mzrPlexFamily*>& refOmniplexFamilies) :
    cpx::plexFamily<bnd::mzrMol,
		    mzrPlex,
		    mzrPlexSpecies,
		    mzrPlexFamily,
		    mzrOmniPlex>(rParadigm,
				 refKnownBindings,
				 refOmniplexFamilies)
  {}

  mzrPlexSpecies*
  mzrPlexFamily::
  constructSpecies(const cpx::siteToShapeMap& rSiteParams,
		   const std::vector<cpx::molParam>& rMolParams)
  {
    return new mzrPlexSpecies(*this,
			      rSiteParams,
			      rMolParams);
  }

  class insertMzrPlexSpecies :
    public std::unary_function<mzrPlexFamily::value_type, void>
  {
    xmlpp::Element* pExplicitSpeciesElt;
    double molFact;
    
  public:
    insertMzrPlexSpecies(xmlpp::Element* pExplicitSpeciesElement,
			 double molarFactor) :
      pExplicitSpeciesElt(pExplicitSpeciesElement),
      molFact(molarFactor)
    {}

    void
    operator()(const argument_type& rPlexFamilyEntry) const
      throw(std::exception)
    {
      mzrPlexSpecies* pSpecies = rPlexFamilyEntry.second;
      pSpecies->insertElt(pExplicitSpeciesElt,
			  molFact);
    }
  };

  void
  mzrPlexFamily::insertSpecies(xmlpp::Element* pExplicitSpeciesElt,
			       double molarFactor) const
    throw(std::exception)
  {
    std::for_each(begin(),
		  end(),
		  insertMzrPlexSpecies(pExplicitSpeciesElt,
				       molarFactor));
  }
}
