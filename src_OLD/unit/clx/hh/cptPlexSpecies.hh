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

#ifndef CLX_CPTPLEXSPECIES_H
#define CLX_CPTPLEXSPECIES_H

#include "utl/dom.hh"
#include "cpx/plexSpcsMixin.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/queryGlobalSpeciesDumpable.hh"
#include "cpt/multiGlobalSpeciesDumpable.hh"
#include "clx/cptPlexSpecies.hh"

namespace clx
{
  class cptPlexFamily;

  class cptPlexSpecies : 
    public cpt::globalSpecies,
    public cpx::plexSpeciesMixin<cptPlexFamily>,
    public fnd::onceNotifier
  {
  public:
    typedef cpt::queryGlobalSpeciesDumpable<cptPlexSpecies> queryDumpableType;

    typedef cpt::multiGlobalSpeciesDumpable<cptPlexSpecies> msDumpableType;

    cptPlexSpecies(cptPlexFamily& rContainingFamily,
		   const cpx::siteToShapeMap& rSiteParams,
		   const std::vector<cpx::molParam>& rMolParams,
		   const cpt::compartmentGraph& rCompartmentGraph,
		   const std::vector<double>& rDiffusionRates);

    ~cptPlexSpecies(void)
    {}

    // This fulfils the pure virtual function in the ancestor class
    // fnd::massive.
    double
    getWeight(void) const;

    void
    notify(int notifyDepth);

    void
    respond(const fnd::newSpeciesStimulus<cpt::compartmentSpecies>& rStim);

    std::string
    getTag(void) const;

    // This overrides basicSpecies::getName(), which just returns a tag.
    std::string
    getName(void) const;

    xmlpp::Element*
    insertElt(xmlpp::Element* pExplicitSpeciesElt) const
      throw(std::exception);
  };
}

#endif // CLX_CPTPLEXSPECIES_H
