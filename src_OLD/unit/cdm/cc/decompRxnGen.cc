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

#include "cpt/cptApp.hh"
#include "cpt/globalReaction.hh"
#include "cpx/cxBinding.hh"
#include "cml/cptBndSite.hh"
#include "clx/clxUnit.hh"
#include "cdm/decompRxnGen.hh"
#include "cdm/cdmUnit.hh"

namespace cdm
{
  void
  decompRxnGen::
  respond(const fnd::featureStimulus<cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily> >& rStimulus)
  {
    if(rCptUnit.getGenerateOk())
      {
	// Get the depth to notify product species for further reaction
	// generation.
	int notificationDepth = rStimulus.getNotificationDepth() - 1;

	const cpx::cxBinding<clx::cptPlexSpecies, clx::cptPlexFamily>&
	  rNewContext = rStimulus.getContext();

	// Create the new reaction, and install it in the family of
	// decomposition reactions for memory management.
	cpt::globalReaction* pReaction
	  = new cpt::globalReaction(rCptUnit.getCompartmentGraph());
	pFamily->addEntry(pReaction);

	// Add the substrate and sensitize to it.
	pReaction->addReactant(rNewContext.getSpecies(),
			       1);

	// Extrapolate the rate of the reaction.
	pReaction->setRate(pExtrap->getRate(rNewContext));
    
	// Ingredients for the mandatory result species.
	std::vector<cpx::molParam> mndtryResultParams;
	clx::cptPlexFamily* pMndtryResultFamily;

	const clx::cptPlex& rWholePlex
	  = rNewContext.getPlexFamily().getParadigm();

	// The index of the binding that is decomposing.
	cpx::bindingSpec breakingBindingNdx
	  = rNewContext.getBindingSpec();
  
	// Make a copy of the chosen plexSpecies's paradigm, but omitting
	// the binding that is to be broken.  At the same time, we'll get
	// the molNdx's of the ends of the binding that is to be broken.
	clx::cptPlex brokenPlex;
	brokenPlex.mols = rWholePlex.mols;

	// Need the two mols on the ends of the specified binding as
	// "seeds" for computing the (1 or 2) connected components of the
	// brokenPlex.
	int leftMolNdx = -1;
	int rightMolNdx = -1;
	// Copy over all except the specified binding.
	for(int bindingNdx = 0;
	    bindingNdx < (int) rWholePlex.bindings.size();
	    bindingNdx++)
	  {
	    const cpx::binding& rBinding
	      = rWholePlex.bindings[bindingNdx];

	    if(bindingNdx == breakingBindingNdx)
	      {
		leftMolNdx = rBinding.leftSite().molNdx();
		rightMolNdx = rBinding.rightSite().molNdx();
	      }
	    else
	      {
		brokenPlex.bindings.push_back(rBinding);
	      }
	  }

	// Get the connected component of the left (mandatory) half of the
	// broken plex.  Note that this way of finding connected components
	// is quite general; actually more than we need, and most of the
	// work (reversing the edge->vertex map) gets repeated here.
	clx::cptPlex leftComponent;
	cpx::plexIso leftIso(brokenPlex.mols.size(),
			     brokenPlex.bindings.size());
	brokenPlex.makeTrackedComponent(leftMolNdx,
					leftComponent,
					leftIso);

	// Find the species of the left component, along with an isomorphism
	// to the species's paradigm.
	cpx::plexIso toLeftParadigm;
	pMndtryResultFamily = rClxUnit.recognize(leftComponent,
						  toLeftParadigm);

	// Construct the parameters for the left component.
	// 
	// First, construct the vector of molParams by permuting the
	// molParams from the original plex.
	for(int leftParaMolNdx = 0;
	    leftParaMolNdx < (int) leftComponent.mols.size();
	    leftParaMolNdx++)
	  {
	    // Bring along the old mol params, using the two tracking maps.
	    int molNdx = toLeftParadigm.backward.molMap[leftParaMolNdx];
	    int preImgNdx = leftIso.backward.molMap[molNdx];
	    mndtryResultParams.push_back
	      (rNewContext.getMolParams()[preImgNdx]);
	  }

	// Construct the mandatory result species.
	clx::cptPlexSpecies* pMndtryResult
	  = pMndtryResultFamily->getMember(mndtryResultParams);

	// Add it as a product of multiplicity one.
	pReaction->addProduct(pMndtryResult,
			      1);

	// Do reaction network generation on this product
	// species at one lower depth.
	if(0 <= notificationDepth)
	  pMndtryResult->ensureNotified(notificationDepth);

	// Now start looking at the other result component, if any.
	bool resultIsConnected
	  = (leftComponent.mols.size() == brokenPlex.mols.size());
	if(! resultIsConnected)
	  {
	    // Ingredients for the optional result species, which only is
	    // needed if the molecule decomposes into two parts.
	    std::vector<cpx::molParam> optResultParams;
	    clx::cptPlexFamily* pOptResultFamily;
  
	    // Extract the other connected component.  Note that we really
	    // might want to do both of these at the same time.
	    clx::cptPlex rightComponent;
	    cpx::plexIso rightIso(brokenPlex.mols.size(),
				  brokenPlex.bindings.size());
	    brokenPlex.makeTrackedComponent(rightMolNdx,
					    rightComponent,
					    rightIso);

	    // Intern the right connected component.
	    cpx::plexIso toRightParadigm;
	    pOptResultFamily = rClxUnit.recognize(rightComponent,
						   toRightParadigm);

	    // Construct the parameters for the right component.
	    // 
	    // First, construct the vector of molParams by permuting the
	    // molParams from the original plex.
	    for(int rightParaMolNdx = 0;
		rightParaMolNdx < (int) rightComponent.mols.size();
		rightParaMolNdx++)
	      {
		// Bring along the old mol params, using the two tracking maps.
		int molNdx = toRightParadigm.backward.molMap[rightParaMolNdx];
		int preImgNdx = rightIso.backward.molMap[molNdx];
	  
		optResultParams.push_back
		  (rNewContext.getMolParams()[preImgNdx]);
	      }

	    // Construct the optional result species.
	    clx::cptPlexSpecies* pOptResult
	      = pOptResultFamily->getMember(optResultParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pOptResult,
				  1);

	    // Do reaction network generation on this product
	    // species at one lower depth.
	    if(0 <= notificationDepth)
	      pOptResult->ensureNotified(notificationDepth);

	    // Finalize the reaction, creating and "hooking up" the
	    // compartment reactions with the compartment species.
	    pReaction->finalizeCompartments(rCptApp.getPropensities());
	  }
      }
  }
}
