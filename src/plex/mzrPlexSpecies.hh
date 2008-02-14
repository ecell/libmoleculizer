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

#ifndef PLEXSPECIES_H
#define PLEXSPECIES_H

/*! \file plexSpecies.hh
  \ingroup plexStructGroup
  \brief Defines plexSpecies, a species of protein complex. */

#include "utl/dom.hh"
#include "fnd/basicDumpable.hh"
#include "cpx/plexSpcsMixin.hh"
#include "mzr/mzrSpeciesDumpable.hh"

namespace plx
{
  class mzrPlexFamily;

  class mzrPlexSpecies : 
    public mzr::mzrSpecies,
    public cpx::plexSpeciesMixin<mzrPlexFamily>
  {
  public:

    typedef mzr::querySpeciesDumpable<mzrPlexSpecies> queryDumpableType;

    typedef mzr::multiSpeciesDumpable<mzrPlexSpecies> msDumpableType;
    
    mzrPlexSpecies(mzrPlexFamily& rContainingFamily,
		   const cpx::siteToShapeMap& rSiteParams,
		   const std::vector<cpx::molParam>& rMolParams) :
      cpx::plexSpeciesMixin<mzrPlexFamily>(rContainingFamily,
					   rSiteParams,
					   rMolParams)
    {}

    ~mzrPlexSpecies(void)
    {}

    // This fulfils the pure virtual function in the ancestor class
    // fnd::massive.
    double
    getWeight(void) const;

    // This fulfills the pure virtual function in the ancestor class
    // cpx::notifier.
    void
    notify(int generateDepth);

    // This overrides basicSpecies::getName(), which just returns a tag.
    std::string
    getName(void) const;

    xmlpp::Element*
    insertElt(xmlpp::Element* pExplicitSpeciesElt,
	      double molarFactor) const
      throw(std::exception);
  };
}

#endif
