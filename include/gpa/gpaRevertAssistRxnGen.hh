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

#ifndef GPAREVERTASSISTRXNGEN_H
#define GPAREVERTASSISTRXNGEN_H

#include "mzr/reactionFamily.hh"
#include "mol/gpaMol.hh"
#include "plex/plexFeature.hh"
#include "plex/plexSpec.hh"
#include "plex/cxOmniParam.hh"

namespace gpa
{
  template<class gpaMolClass>
  class gpaRevertAssistRxnGen :
    public omniInContext::rxnGen
  {
    reactionFamily* pFamily;

    gpaMolClass* pGpaMol;
    plexMolSpec gpaMolSpec;

    // A unary reaction.
    double revertRate;

  public:

    gpaRevertAssistRxnGen(gpaMolClass* pTheGpaMol,
			  plexMolSpec theGpaMolSpec,
			  double reactionRate,
			  reactionFamily* pReactionFamily) :
      pFamily(pReactionFamily),
      pGpaMol(pTheGpaMol),
      gpaMolSpec(theGpaMolSpec),
      revertRate(reactionRate)
    {}

    void
    makeReactions(const omniInContext& rContext) const;
  };

  template<class gpaMolClass> void gpaRevertAssistRxnGen<gpaMolClass>::
  makeReactions(const omniInContext& rContext) const
  {
    // Wrap the new species (containing the enabling omniplex) in accessors.
    cxOmni cxGpaComplex(rContext);

    // Get the state of gpa1.
    int gpaNdx = cxGpaComplex.translateMolSpec(gpaMolSpec);
    molParam gpaParam = cxGpaComplex.getMolParam(gpaNdx);
    const typename gpaMolClass::molStateClass* pState
      = dynamic_cast<const typename gpaMolClass::molStateClass*>(gpaParam);
    if(! pState)
      {
	fatal("gpaRevertAssistRxnGne::makeReactions: "
	      "mol is not a GPA mol.");
      }

    // Generate a reaction only if the mol is GTP-bound.
    if(pState->getGpaForm() == gtpForm)
      {
	// Create the reaction.
	reaction* pReaction = new reaction();
	pFamily->addReaction(pReaction);

	// Add the complex containing the gpa protein as a substrate.
	pReaction->addSubstrate(cxGpaComplex.getSpecies(),
				1);
	pReaction->sensitizeToSubstrates();

	// Set the reaction rate.  Since this is essentially a unary
	// reaction, the given rate is the real rate.
	pReaction->setRate(revertRate);

	// Make the new (gdp) state of the gpaMol.
	typename gpaMolClass::molStateClass nuState(*pState);
	nuState.setGpaForm(gdpForm);
	molParam nuParam = pGpaMol->internState(nuState);

	// Make the product species generator.
	std::vector<molParam> molParams(cxGpaComplex.getMolParams());
	molParams[gpaNdx] = nuParam;
	plexSpeciesGen* pGen
	  = new plexSpeciesGen(molParams,
			       &(cxGpaComplex.getPlexFamily()));
	pReaction->addProductGenerator(pGen,
				       1);
      }
  }
}

#endif // GPAREVERTASSISTRXNGEN_H
