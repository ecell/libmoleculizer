#include "mzr/reaction.hh"
#include "mzr/pchem.hh"
#include "mzr/moleculizer.hh"
#include "mol/mol.hh"
#include "plex/cxSiteParam.hh"
#include "plex/plexUnit.hh"
#include "dimer/dimerizeRxnGen.hh"
#include "dimer/dimerUnit.hh"

namespace dimer
{
  //! Auxiliary class for joining plexes in a dimerization.
  class offsetBinding :
    public std::unary_function<plx::plexBinding, plx::plexBinding>
  {
    int offset;
  public:
    offsetBinding(int molOffset) :
      offset(molOffset)
    {}

    plx::plexSiteSpec offsetSite(const plx::plexSiteSpec& rSpec)
    {
      return plx::plexSiteSpec(rSpec.molNdx() + offset,
			  rSpec.siteNdx());
    }

    plx::plexBinding operator()(const plx::plexBinding& rBinding)
    {
      return plx::plexBinding(offsetSite(rBinding.leftSite()),
			 offsetSite(rBinding.rightSite()));
    }
  };

  void
  dimerizeRxnGenPair::
  makeBinaryReactions(const plx::siteInContext& rLeftContext,
		      const plx::siteInContext& rRightContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the siteInContext's with accessor functions.
	plx::cxSite cxLeft(rLeftContext);
	plx::cxSite cxRight(rRightContext);
    
	// Construct the (only) new reaction and install it in the family
	// of dimerization reactions.
	mzr::reaction* pReaction = new mzr::reaction();
	pFamily->addReaction(pReaction,
			     rMzrUnit);

	// Add the substrates.
	pReaction->addSubstrate(cxLeft.getSpecies(),
				1);
	pReaction->addSubstrate(cxRight.getSpecies(),
				1);
	pReaction->sensitizeToSubstrates(rMzrUnit);

	// Extrapolate the rate for the new reaction.
	pReaction->setRate(pExtrap->getRate(cxLeft,
					    cxRight));

	// Now start cooking up the product species.
	std::vector<bnd::molParam> resultMolParams;
	plx::plexFamily* pResultFamily;
    
	// We join two distinct virtual molecules together.
	// Begin constructing the "joined" plex by copying the
	// paradigm of the left isomorphism class.
	const plx::plex& rLeftParadigm
	  = cxLeft.getPlexFamily().getParadigm();
	plx::plex joined(rLeftParadigm);

	// Insert the rightIsoClass paradigm's mols at the end.
	const plx::plex& rRightParadigm
	  = cxRight.getPlexFamily().getParadigm();
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
	   std::back_insert_iterator<std::vector<plx::plexBinding> >(joined.bindings),
	   offset);
  
	// Make a binding that binds the sites corresponding to the two
	// original free binding sites.
	plx::plexBinding
	  joint(cxLeft.getSiteSpec(),
		offset.offsetSite(cxRight.getSiteSpec()));

	// Add the new binding to the joined plex.
	joined.bindings.push_back(joint);

	// Intern the joined plex.
	plx::plexIsoPair joinedToParadigm;
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
		  (cxLeft.getMolParams()[joinedNdx]);
	      }
	    else
	      {
		resultMolParams.push_back(cxRight.getMolParams()
					  [joinedNdx - leftMolCount]);
	      }
	  }

	// Construct result species.
	plx::plexSpecies* pResult
	  = pResultFamily->getMember(resultMolParams);

	// Add it as a product of multiplicity 1.
	pReaction->addProduct(pResult,
			      1);
      }
  }
}
