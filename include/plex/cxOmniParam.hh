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

#ifndef CXOMNIPARAM_H
#define CXOMNIPARAM_H

#include "mzr/feature.hh"
#include "plex/plexFamily.hh"
#include "plex/subPlexSpec.hh"

namespace plx
{
  class cxOmni : public mzr::featureContext<plexSpecies, subPlexSpec>
  {
  public:
    cxOmni(const mzr::featureContext<plexSpecies, subPlexSpec>& rfc) :
      mzr::featureContext<plexSpecies, subPlexSpec>(rfc)
    {}

    // Retrieves the embedding of the omniplex into the complex where
    // it was recognized.
    const plexIsoPair&
    getEmbedding(void) const
    {
      return getSpec().getInjection();
    }

    // Uses the embedding to translate mol indices in the omniplex
    // to the corresponding mol indices in the complex where the
    // omniplex was recognized.
    plexMolSpec
    translateMolSpec(plexMolSpec specInOmni) const
    {
      return getEmbedding().forward.molMap[specInOmni];
    }

    // Uses the embedding to translate binding indices in the omniplex
    // to the corresponding binding indices in the complex where the
    // omniplex was recognized.
    plexBindingSpec
    translateBindingSpec(plexBindingSpec specInOmni) const
    {
      return getEmbedding().forward.bindingMap[specInOmni];
    }
  
    // Uses the embedding to translate a site spec in the omniplex into
    // a site spec in the complex where the omniplex was recognized.
    plexSiteSpec
    translateSiteSpec(plexSiteSpec specInOmni) const
    {
      return plexSiteSpec(translateMolSpec(specInOmni.molNdx()),
			  specInOmni.siteNdx());
    }

    // Get the plex family in which the omniplex was recognized.
    plexFamily&
    getPlexFamily(void) const
    {
      // Same as the template function.
      return getSpecies()->getFamily();
    }

    // Returns the full plexParam of the particular complex where the omni
    // occurs.  This is (frequently) used to specify the precise species
    // of complex that is a substrate of a reaction, usually to decrement
    // it, because the reaction has transformed it into something(s) else.
    const plexParam&
    getPlexParam(void) const
    {
      // Same as the template function.
      return getSpecies()->getParam();
    }

    // Returns the population of the particular species of complex where
    // the omni was found.
    int
    getPop(void) const
    {
      return getSpecies()->getPop();
    }

    // Used in many propensity calculations.
    double
    getPlexWeight(void) const
    {
      return getSpecies()->getWeight();
    }

    // Used in sensitizeToSubstrates.
    void
    addSensitiveReaction(mzr::reaction* pReaction) const
    {
      getSpecies()->addSensitiveReaction(pReaction);
    }

    // Returns just the molParams of the particular complex where the
    // omni occurs.  These are usually used to "cook up" the vector of
    // molParams of a reaction product.  The remaining parameters of the
    // reaction product are usually created by the product complex's
    // allostery function.
    const std::vector<bnd::molParam>&
    getMolParams(void) const
    {
      return getPlexParam().molParams;
    }

    // Note that this does NOT translate the molSpec.
    bnd::molParam
    getMolParam(plx::plexMolSpec molSpec) const
    {
      const std::vector<bnd::molParam>& rMolParams = getMolParams();
      return rMolParams[molSpec];
    }

    // Get the vector of mols of the plex.  Note that the indexing in
    // this vector is that given by the plex, not by the omniplex.  Use
    // translateMolNdx above to overcome this.
    //
    // Not building the translation routines into any of these, since
    // it's frequently a case of translate once, use many times.
    const std::vector<bnd::mol*>&
    getMols(void) const
    {
      return getPlexFamily().getParadigm().mols;
    }
  };
}

#endif
