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

#include "mzr/mzrUnit.hh"
#include "mol/modMol.hh"
#include "plex/cxMolParam.hh"
#include "ftr/uniMolGen.hh"

namespace ftr
{
  class exchangeMods :
    public std::unary_function<molModExchange, void>
  {
    bnd::modMolState& rState;

  public:
    exchangeMods(bnd::modMolState& rTargetState) :
      rState(rTargetState)
    {}

    void
    operator()(const molModExchange& rExchange) const
    {
      rState[rExchange.modSiteNdx] = rExchange.pReplacementMod;
    }
  };
  
  void
  uniMolRxnGen::
  makeReactions(const plx::molInContext& rContext) const
  {
    // See if reaction generation is enabled.
    if(rMzrUnit.generateOk)
      {
	// Wrap the enabling mol in accessors.
	plx::cxMol cxEnabling(rContext);

	// Check if the mol's state passes the tests for reaction generation.
	bnd::andMolQueries& rMolQueries = *pMolQueries;
	bnd::molParam enablingParam = cxEnabling.getMolParam();
	if(rMolQueries(enablingParam))
	  {
	    // Construct the reaction and intern it for memory management.
	    mzr::reaction* pReaction = new mzr::reaction();
	    pFamily->addReaction(pReaction,
				 rMzrUnit);

	    // Add the triggering complex as a reactant of multiplicity 1.
	    pReaction->addSubstrate(cxEnabling.getSpecies(),
				    1);

	    // If an auxiliary reactant was specified, add it with
	    // multiplicity 1.
	    if(pAdditionalReactant)
	      {
		pReaction->addSubstrate(pAdditionalReactant,
					1);
	      }

	    // Sensitize the reaction to its reactants (i.e. add the reaction
	    // to the sensitivity list of each reaction, and to the
	    // sensitivity list of a few other state variables, such as
	    // volume.)
	    pReaction->sensitizeToSubstrates(rMzrUnit);

	    // If an auxiliary product was specified, add it with multiplicity
	    // 1.
	    if(pAdditionalProduct)
	      {
		pReaction->addProduct(pAdditionalProduct,
				      1);
	      }

	    // Set the rate of the reaction.
	    pReaction->setRate(pExtrapolator->getRate(rContext));

	    // Construct the primary product species.
	    // 
	    // Start by constructing the state of the enabling mol
	    // after modification exchange.
	    bnd::modMolState enablingState
	      = pEnablingMol->externState(enablingParam);
	    std::for_each(molModExchanges.begin(),
			  molModExchanges.end(),
			  exchangeMods(enablingState));

	    // Put the post-modification state of the mol into the product
	    // species's molParams.
	    std::vector<bnd::molParam> productMolParams 
	      = cxEnabling.getMolParams();
	    productMolParams[cxEnabling.getMolSpec()]
	      = pEnablingMol->internState(enablingState);

	    // Determine the primary product species.
	    plx::plexFamily& rPlexFamily
	      = cxEnabling.getPlexFamily();
	    plx::plexSpecies* pProductSpecies
	      = rPlexFamily.makeMember(productMolParams);

	    // Add the primary product species to the reaction, with
	    // multiplicity 1.
	    pReaction->addProduct(pProductSpecies,
				  1);
	  }
      }
  }
}
