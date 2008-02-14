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

#include "odi/odieApp.hh"
#include "odi/odieParse.hh"

namespace odie
{
  // Multiplies derivative of a reaction amount by the concntration
  // of a substrate raised to the power of the multiplicity of the
  // substrate.
  //
  // The derivative of a reaction amount is a monomial in the species
  // concentrations, so this is just adding the multiplicity of the substrate
  // to the exponent of the substrate's concentration.
  class timesSubstrateRate :
    public std::unary_function<const std::pair<int, int>&, void>
  {
    rk4util::monomial& rReactionDeriv;
  public:
    timesSubstrateRate(rk4util::monomial& rReactionDerivative) :
      rReactionDeriv(rReactionDerivative)
    {}

    void
    operator()(const std::pair<int, int>& rSubstrateMultPair) const
    {
      rReactionDeriv.addToExponent(rSubstrateMultPair.first,
				   rSubstrateMultPair.second);
    }
  };

  class subtractSubstrateMult : public
  std::unary_function<std::map<int, int>::value_type, void>
  {
    std::map<int, int>& rDeltas;
  public:
    subtractSubstrateMult(std::map<int, int>& rDeltaMap) :
      rDeltas(rDeltaMap)
    {
    }

    void
    operator()(const argument_type& rSpeciesMultPair) const
    {
      // Substrate species are decremented by their multiplicity
      // as a substrate.
      int delta = - rSpeciesMultPair.second;

      // Attempt to insert, in case this species hasn't yet appeared
      // in the deltas map.
      std::pair<std::map<int, int>::iterator, bool> insertResult
	= rDeltas.insert(std::make_pair(rSpeciesMultPair.first,
					delta));

      // If the species has already appeared in the deltas map, then
      // update its entry.
      if(! insertResult.second)
	{
	  insertResult.first->second += delta;
	}
    }
  };

  class insertReactionTerm :
    public std::unary_function<std::map<int, int>::value_type, void>
  {
    const rk4util::monomial& rReactDeriv;
    double rate;
    rk4util::polymap<double>& rSpeciesDerivs;
    
  public:
    insertReactionTerm(const rk4util::monomial& rReactionDeriv,
		       double reactionRate,
		       rk4util::polymap<double>& rSpeciesDerivatives) :
      rReactDeriv(rReactionDeriv),
      rate(reactionRate),
      rSpeciesDerivs(rSpeciesDerivatives)
    {}

    void
    operator()(const argument_type& rSpeciesDeltaPair) const
    {
      // Multiply the monomial associated to this reaction
      // by the reaction rate and the delta; this is the term
      // that this reaction contributes to the derivative of the
      // affected species.
      rk4util::polynomial<double> term;
      term.setCoeff(rReactDeriv,
		    rate * ((double) rSpeciesDeltaPair.second));

      // Add the term to the affected species's derivative.
      term.addTo(rSpeciesDerivs[rSpeciesDeltaPair.first]);
    }
  };

  class insertReactionDerivatives :
    public std::unary_function<parserReaction&, void>
  {
    rk4util::polymap<double>& rDerivatives;
  public:
    insertReactionDerivatives(rk4util::polymap<double>& rSpeciesDerivs) :
      rDerivatives(rSpeciesDerivs)
    {}

    void
    operator()(const parserReaction& rParserReaction) const
    {
      // Make the monomial associated to the substrates of the reaction.
      // This gives the rate at which the reaction will happen at any time.
      rk4util::monomial reactionDeriv;
      std::for_each(rParserReaction.substrates.begin(),
	       rParserReaction.substrates.end(),
	       timesSubstrateRate(reactionDeriv));

      // Run through both the substrates and the products to accumulate
      // the deltas.
      std::map<int, int> deltas;
      deltas.insert(rParserReaction.products.begin(),
		    rParserReaction.products.end());
      std::for_each(rParserReaction.substrates.begin(),
		    rParserReaction.substrates.end(),
		    subtractSubstrateMult(deltas));

      // Insert a multiple of the above monomial to the derivative
      // of each of the species affected by this reaction.
      std::for_each(deltas.begin(),
		    deltas.end(),
		    insertReactionTerm(reactionDeriv,
				       rParserReaction.rate,
				       rDerivatives));
    }
  };

  void
  odieApp::makeDerivatives(std::vector<parserReaction>& rReactions)
  {
    for_each(rReactions.begin(),
	     rReactions.end(),
	     insertReactionDerivatives(derivative));
  }
}
