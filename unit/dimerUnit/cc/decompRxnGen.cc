#include "mzr/moleculizer.hh"
#include "mzr/reaction.hh"
#include "plex/cxBindingParam.hh"
#include "plex/plexUnit.hh"
#include "dimer/decompRxnGen.hh"
#include "dimer/dimerUnit.hh"

namespace dimer
{
  void
  decompRxnGen::makeReactions(const plx::bindingInContext& rNewContext) const
  {
    if(rMzrUnit.generateOk)
      {
	
	// Wrap the new bindingInContext with accessor functions.
	plx::cxBinding cxSubstrate(rNewContext);

	// Create the new reaction, and install it in the family of
	// decomposition reactions.
	mzr::reaction* pReaction = new mzr::reaction();
	pFamily->addReaction(pReaction,
			     rMzrUnit);

	// Add the substrate and sensitize to it.
	pReaction->addSubstrate(cxSubstrate.getSpecies(),
				1);
	pReaction->sensitizeToSubstrates(rMzrUnit);

	// Extrapolate the rate of the reaction.
	pReaction->setRate(pExtrap->getRate(rNewContext));
    
	// Set the rate of reaction.  Since this is a unary reaction,
	// no adjustment for molecular weight is done.
	//     pReaction->setRate(cxSubstrate.getBindingParam());

	// Now we start cooking up the product species, which from the old
	// routine "doReactionFirstTime".

	// Ingredients for the mandatory result species.
	std::vector<bnd::molParam> mndtryResultParams;
	plx::plexFamily* pMndtryResultFamily;

	const plx::plex& rWholePlex = cxSubstrate.getPlexFamily().getParadigm();
	plx::plexBindingSpec breakingBindingNdx = cxSubstrate.getBindingSpec();
  
	// Make a copy of the chosen plexSpecies's paradigm, but omitting
	// the binding that is to be broken.  At the same time, we'll get
	// the molNdx's of the ends of the binding that is to be broken.
	plx::plex brokenPlex;
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
	    const plx::plexBinding& rBinding
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
	plx::plex leftComponent;
	plx::plexIsoPair leftIso(brokenPlex.mols.size(),
				 brokenPlex.bindings.size());
	brokenPlex.makeTrackedComponent(leftMolNdx,
					leftComponent,
					leftIso);

	// Find the species of the left component, along with an isomorphism
	// to the species's paradigm.
	plx::plexIsoPair toLeftParadigm;
	pMndtryResultFamily = rPlexUnit.recognize(leftComponent,
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
	      (cxSubstrate.getMolParams()[preImgNdx]);
	  }

	// Construct the mandatory result species.
	plx::plexSpecies* pMndtryResult
	  = pMndtryResultFamily->getMember(mndtryResultParams);

	// Add it as a product of multiplicity one.
	pReaction->addProduct(pMndtryResult,
			      1);

	// Now start looking at the other result component, if any.
	bool resultIsConnected
	  = (leftComponent.mols.size() == brokenPlex.mols.size());
	if(! resultIsConnected)
	  {
	    // Ingredients for the optional result species, which only is
	    // needed if the molecule decomposes into two parts.
	    std::vector<bnd::molParam> optResultParams;
	    plx::plexFamily* pOptResultFamily;
  
	    // Extract the other connected component.  Note that we really
	    // might want to do both of these at the same time.
	    plx::plex rightComponent;
	    plx::plexIsoPair rightIso(brokenPlex.mols.size(),
				      brokenPlex.bindings.size());
	    brokenPlex.makeTrackedComponent(rightMolNdx,
					    rightComponent,
					    rightIso);

	    // Intern the right connected component.
	    plx::plexIsoPair toRightParadigm;
	    pOptResultFamily = rPlexUnit.recognize(rightComponent,
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
		  (cxSubstrate.getMolParams()[preImgNdx]);
	      }

	    // Construct the optional result species.
	    plx::plexSpecies* pOptResult
	      = pOptResultFamily->getMember(optResultParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pOptResult,
				  1);
	  }
      }
  }
}
