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

#include "mzr/pchem.hh"
#include "mol/molUnit.hh"
#include "plex/cxOmniParam.hh"
#include "gpa/gpaExchangeRxnGen.hh"

namespace gpa
{
  // One alternative that really might be more correct would be to generate
  // four reactions for every occurrence of the g-protein: binding and
  // unbinding reactions for both GTP and GDP.  The rates would be different
  // depending on whether the g-protein was in the "enabling" omni or not.
  //
  // This would require a significantly different kind of reaction generator,
  // as it would have to be "sensitive" to the mol, but decide whether the mol
  // was in the enabling complex at reaction generation time.  When the
  // reaction generator is notified of a new species of complex containing the
  // g-protein, it would have to see if the enabling subcomplex is present by
  // looking at the plexFamily's subPlex features.  If any of the subPlex's
  // featured by the plexFamily are the enabling subcomplex, then it only
  // remains to verify that the "reporting" target mol is really the image of
  // the target mol under the one of the subPlex's inclusions (supposing for
  // the hell of it that there might be more than one.)  If so, then a
  // binding, unbinding, or no reaction is generated depending on the status
  // of the "nucleotide-binding" modification site.
  //
  // On the other hand, if we fail for any reason to get to this point,
  // then a binding, unbinding, or no reaction is generated, again depending on
  // the status of the "nucleotide-binding" modification site, but with the
  // other set of rates.

  // Simple nucleotide-exchange reaction generator, see above for a possibly
  // better implementation.
  void
  gpaExchangeRxnGen::makeReactions(const plx::omniInContext& rContext) const
  {
    if(rMzrUnit.generateOk)
      {
	// Wrap the new species of receptor complex in accessors.
	plx::cxOmni cxReceptor(rContext);

	// Get the state of the gpa mol.
	//
	// First, translate its index in the "enabling subcomplex" into
	// the corresponding index in the species containing the subcomplex.
	int gpaMolSpec = cxReceptor.translateMolSpec(gpaSpec);
	bnd::molParam pState = cxReceptor.getMolParam(gpaMolSpec);

	// Get the modification at the index where we expect to find
	// GDP, GTP, or nothing.  We only construct an exchange reaction
	// when we find GDP.
	const bnd::modMolState& rModMolState = pMol->externState(pState);
	if(rModMolState[gtpSiteNdx] == pGdpBound)
	  {
	    // Create the new reaction and install it in the reaction family
	    // for memory management.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add the receptor complex containing the gpa protein as a substrate
	    // of the reaction.
	    pReaction->addSubstrate(cxReceptor.getSpecies(),
				    1);
	    // Add GTP as a substrate.
	    pReaction->addSubstrate(pGTP,
				    1);

	    // Sensitize to all substrates.
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // Set the reaction rate.
	    pReaction->setRate(pExtrap->getRate(rContext));

	    // Add GDP as a product.
	    pReaction->addProduct(pGDP,
				  1);

	    // Make the new state of the mol that is being converted from gdp
	    // form to gtp form.
	    bnd::modMolState nuState(rModMolState);
	    nuState[gtpSiteNdx] = pGtpBound;
	    bnd::molParam nuParam = pMol->internState(nuState);

	    // Make the molParams for the result plexSpecies.
	    std::vector<bnd::molParam> molParams(cxReceptor.getMolParams());
	    molParams[gpaMolSpec] = nuParam;

	    // Make the result plexSpecies.
	    plx::plexSpecies* pResult
	      = cxReceptor.getPlexFamily().getMember(molParams);

	    // Add it as a product of multiplicity 1.
	    pReaction->addProduct(pResult,
				  1);
	  }
      }
  }
}
