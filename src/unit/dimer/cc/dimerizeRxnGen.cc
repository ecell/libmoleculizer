#include "mzr/mzrReaction.hh"
#include "fnd/pchem.hh"
#include "mzr/moleculizer.hh"
#include "mol/mzrMol.hh"
#include "cpx/cxSite.hh"
#include "plex/plexUnit.hh"
#include "dimer/dimerizeRxnGen.hh"
#include "dimer/dimerUnit.hh"

namespace dimer
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

  void
  dimerizeRxnGenPair::
  makeBinaryReactions
  (const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rLeftContext,
   const cpx::cxSite<plx::mzrPlexSpecies, plx::mzrPlexFamily>& rRightContext,
   int generateDepth) const
  {
    if(rMzrUnit.getGenerateOk())
      {
	// Construct the (only) new reaction and install it in the family
	// of dimerization reactions.
	mzr::mzrReaction* pReaction
	  = new mzr::mzrReaction(rMzrUnit.globalVars.begin(),
				 rMzrUnit.globalVars.end());
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
	plx::mzrPlexFamily* pResultFamily;
    
	// We join two distinct virtual molecules together.
	// Begin constructing the "joined" plex by copying the
	// paradigm of the left isomorphism class.
	const plx::mzrPlex& rLeftParadigm
	  = rLeftContext.getPlexFamily().getParadigm();
	plx::mzrPlex joined(rLeftParadigm);

	// Insert the rightIsoClass paradigm's mols at the end.
	const plx::mzrPlex& rRightParadigm
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
	pResultFamily = rPlexUnit.recognize(joined,
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
	plx::mzrPlexSpecies* pResult
	  = pResultFamily->getMember(resultMolParams);

	// Add it as a product of multiplicity 1.
	pReaction->addProduct(pResult,
			      1);

	// Continue reaction generation at one lower depth.
	int notificationDepth = generateDepth - 1;
	if(0 <= notificationDepth)
	  {
	    pResult->ensureNotified(notificationDepth);
	  }
      }
  }
}
