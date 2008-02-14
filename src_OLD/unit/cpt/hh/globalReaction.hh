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

#ifndef GLOBALREACTION_H
#define GLOBALREACTION_H

#include "fnd/basicReaction.hh"
#include "cpt/propensityDistro.hh"
#include "cpt/globalSpecies.hh"
#include "cpt/cptReaction.hh"

namespace cpt
{
  class globalReaction :
    public fnd::basicReaction<globalSpecies>
  {
    // This arguably should not be here.  If it is here, it should arguably
    // be static, since it should be the same for all globalReactions.
    const compartmentGraph& rGraph;

    // Not clear that globalReaction aught to own all of its compartment
    // reactions, as this autoVector implies.
    utl::autoVector<cptReaction> cptReactionVector;

    static int 
    globalReactionCount;

  public:
    globalReaction(const compartmentGraph& rCompartmentGraph,
		   double rate = 0) :
      fnd::basicReaction<globalSpecies>(rate),
      rGraph(rCompartmentGraph)
    {
      ++globalReactionCount;
    }

    const compartmentGraph&
    getCompartmentGraph(void) const
    {
      return rGraph;
    }

    static
    int
    getGlobalReactionCount(void)
    {
      return globalReactionCount;
    }

    // After all the reactants and products have been added, rate set, etc.
    // call this method to generate the actual reactions.
    //
    // Note the similarity between this and the "after part" of the
    // globalReaction constructor.
    //
    // The only elegant way I can think of to avoid this is to make this class
    // into something like "parserPlex", and make this step something like
    // extraction of the final globalReaction from the "parserGlobalReaction,"
    // whose purpose is to facilitate stepwise construction, as during
    // parsing.  It's hard for me to see a use for this class once simulation
    // is up and running.
    //
    // An alternative location for much of what goes on in this method
    // is compartmentReaction's constructor.
    void
    finalizeCompartments(propensityDistro& rPropensityDistro);

    // Override basic_reaction<globalSpecies>::addReactant to sensitize
    // to the added reactant.
    void
    addReactant(globalSpecies* pSpecies,
		int multiplicity);

    // Output generation for state dump.
    xmlpp::Element*
    insertElt(xmlpp::Element* pParentElt) const 
      throw(std::exception);
  };
}

#endif // GLOBALREACTION_H
