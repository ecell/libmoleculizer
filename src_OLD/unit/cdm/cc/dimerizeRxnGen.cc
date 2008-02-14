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

#include "cpt/globalReaction.hh"
#include "fnd/pchem.hh"
#include "cpt/cptApp.hh"
#include "cml/cptMol.hh"
#include "cpx/cxSite.hh"
#include "clx/clxUnit.hh"
#include "cdm/dimerizeRxnGen.hh"
#include "cdm/cdmUnit.hh"

namespace cdm
{
  //! Auxiliary class for joining plexes in a dimerization.
  class offsetBinding :
    public std::unary_function<cpx::binding, cpx::binding>
  {
    int offset;
  public:
    offsetBinding(int molOffset) :
      offset(molOffset)
    {}

    cpx::siteSpec
    offsetSite(const cpx::siteSpec& rSpec)
    {
      return cpx::siteSpec(rSpec.molNdx() + offset,
			   rSpec.siteNdx());
    }

    cpx::binding
    operator()(const cpx::binding& rBinding)
    {
      return cpx::binding(offsetSite(rBinding.leftSite()),
			  offsetSite(rBinding.rightSite()));
    }
  };

  // The generateDepth argument here came directly from the featureStimulus
  // that promped reaction the reaction generator; it has not yet
  // been decremented.
  void
  dimerizeRxnGenPair::
  makeBinaryReactions
  (const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rLeftContext,
   const cpx::cxSite<clx::cptPlexSpecies, clx::cptPlexFamily>& rRightContext,
   int generateDepth) const
  {
    if(rCptUnit.getGenerateOk())
      {
	// Construct the (only) new reaction and install it in the family
	// of dimerization reactions.
	cpt::globalReaction* pReaction
	  = new cpt::globalReaction(rCptUnit.getCompartmentGraph());
	pFamily->addEntry(pReaction);

	// Add the reactants.
	pReaction->addReactant(rLeftContext.getSpecies(),
				1);
	pReaction->addReactant(rRightContext.getSpecies(),
				1);

	// Extrapolate the rate for the new reaction.
	pReaction->setRate(pExtrap->getRate(rLeftContext,
					    rRightContext));

	// Now start cooking up the product species.
	std::vector<cpx::molParam> resultMolParams;
	clx::cptPlexFamily* pResultFamily;
    
	// We join two distinct virtual molecules together.
	// Begin constructing the "joined" plex by copying the
	// paradigm of the left isomorphism class.
	const clx::cptPlex& rLeftParadigm
	  = rLeftContext.getPlexFamily().getParadigm();
	clx::cptPlex joined(rLeftParadigm);

	// Insert the rightIsoClass paradigm's mols at the end.
	const clx::cptPlex& rRightParadigm
	  = rRightContext.getPlexFamily().getParadigm();
	joined.mols.insert(joined.mols.end(),
			   rRightParadigm.mols.begin(),
			   rRightParadigm.mols.end());

	// Insert the rest of the "joined" plex's sites and
	// bindings, offsetting the mol indices.
	//
	// First, make the function that will do the offsetting.
	int leftMolCount = rLeftParadigm.mols.size();
	offsetBinding offset(leftMolCount);

	// Now do the offsetting.
	transform
	  (rRightParadigm.bindings.begin(),
	   rRightParadigm.bindings.end(),
	   std::back_insert_iterator<std::vector<cpx::binding> >(joined.bindings),
	   offset);
  
	// Make a binding that binds the sites corresponding to the two
	// original free binding sites.
	cpx::binding
	  joint(rLeftContext.getSiteSpec(),
		offset.offsetSite(rRightContext.getSiteSpec()));

	// Add the new binding to the joined plex.
	joined.bindings.push_back(joint);

	// Intern the joined plex.
	cpx::plexIso joinedToParadigm;
	pResultFamily = rClxUnit.recognize(joined,
					    joinedToParadigm);

	// Reorder the joinedPlex's parameters using recognition
	// permutation.
	for(int paradigmNdx = 0;
	    paradigmNdx < (int) joined.mols.size();
	    paradigmNdx++)
	  {
	    // Construct the joinedPlex's (mol) paramters by
	    // virtually appending the left and right plexes' (mol)
	    // parameters.
	    //
	    // This needs to be cleaned up somehow.
	    int joinedNdx = joinedToParadigm.backward.molMap[paradigmNdx];
	    if(joinedNdx < leftMolCount)
	      {
		resultMolParams.push_back
		  (rLeftContext.getMolParams()[joinedNdx]);
	      }
	    else
	      {
		resultMolParams.push_back(rRightContext.getMolParams()
					  [joinedNdx - leftMolCount]);
	      }
	  }

	// Construct result species.
	clx::cptPlexSpecies* pResult
	  = pResultFamily->getMember(resultMolParams);

	// Add it as a product of multiplicity 1.
	pReaction->addProduct(pResult,
			      1);

	// Continue reaction network generation at one lower depth.
	int notifyDepth = generateDepth - 1;
	if(0 < notifyDepth) pResult->ensureNotified(notifyDepth);

	// Finalize the reaction, creating and "hooking up" the compartment
	// reactions with the compartment species.
	pReaction->finalizeCompartments(rCptApp.getPropensities());
      }
  }
}
