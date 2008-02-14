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

#include "cft/uniMolGen.hh"

namespace cft
{
  class exchangeMods :
    public std::unary_function<molModExchange, void>
  {
    cpx::modMolState& rState;

  public:
    exchangeMods(cpx::modMolState& rTargetState) :
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
  respond(const fnd::featureStimulus<cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily> >& rStimulus)
  {
    // See if reaction generation is enabled.
    if(rCptUnit.getGenerateOk())
      {
	const cpx::cxMol<clx::cptPlexSpecies, clx::cptPlexFamily>& rNewContext
	  = rStimulus.getContext();
	
	// Check if the mol's state passes the tests for reaction generation.
	cpx::andMolStateQueries& rMolQueries = *pMolQueries;
	cpx::molParam enablingParam = rNewContext.getMolParam();
	if(rMolQueries(enablingParam))
	  {
	    // Construct the reaction and intern it for memory management.
	    cpt::globalReaction* pReaction
	      = new cpt::globalReaction(rCptUnit.getCompartmentGraph());
	    pFamily->addEntry(pReaction);

	    // Add the triggering complex as a reactant of multiplicity 1.
	    pReaction->addReactant(rNewContext.getSpecies(),
				    1);

	    // If an auxiliary reactant was specified, add it with
	    // multiplicity 1.
	    if(pAdditionalReactant)
	      {
		pReaction->addReactant(pAdditionalReactant,
					1);
	      }

	    // If an auxiliary product was specified, add it with multiplicity
	    // 1.
	    if(pAdditionalProduct)
	      {
		pReaction->addProduct(pAdditionalProduct,
				      1);
	      }

	    // Set the rate of the reaction.
	    pReaction->setRate(pExtrapolator->getRate(rNewContext));

	    // Construct the primary product species.
	    // 
	    // Start by constructing the state of the enabling mol
	    // after modification exchange.
	    cpx::modMolState enablingState
	      = pEnablingMol->externState(enablingParam);
	    std::for_each(molModExchanges.begin(),
			  molModExchanges.end(),
			  exchangeMods(enablingState));

	    // Put the post-modification state of the mol into the product
	    // species's molParams.
	    std::vector<cpx::molParam> productMolParams 
	      = rNewContext.getMolParams();
	    productMolParams[rNewContext.getMolSpec()]
	      = pEnablingMol->internState(enablingState);

	    // Determine the primary product species.
	    clx::cptPlexFamily& rPlexFamily
	      = rNewContext.getPlexFamily();
	    clx::cptPlexSpecies* pProductSpecies
	      = rPlexFamily.makeMember(productMolParams);

	    // Add the primary product species to the reaction, with
	    // multiplicity 1.
	    pReaction->addProduct(pProductSpecies,
				  1);

	    // Continue reaction network generation at depth one lower.
	    int notificationDepth = rStimulus.getNotificationDepth() - 1;
	    if(0 <= notificationDepth)
	      pProductSpecies->ensureNotified(notificationDepth);
	    
	    // Finalize the reaction, creating and "hooking up" the
	    // compartment reactions withe the compartment species.
	    pReaction->finalizeCompartments(rCptApp.getPropensities());
	  }
      }
  }
}
