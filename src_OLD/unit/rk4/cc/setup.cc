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

#include "rk4/rk4tau.hh"
#include "rk4/tauApp.hh"
#include "rk4/tauParse.hh"

namespace rk4tau
{
  // Perhaps a general generic function to convert a scalar type to a
  // polynomial over that scalar type would be better.
  static polynomial<int>
  int2intPoly(int anInt)
  {
    polynomial<int> result;
    monomial one;
    result.setCoeff(one,
		    anInt);
    return result;
  }

  // Performs function similar to that of int2intPoly above, but
  // converts integer population to concentration using "molar factor."
  class pop2concPoly :
    public std::unary_function<int, polynomial<double> >
  {
    double molarConst;

  public:
    pop2concPoly(double molarConstant) :
      molarConst(molarConstant)
    {}

    polynomial<double>
    operator()(int population) const
    {
      polynomial<double> result;
      monomial one;
      result.setCoeff(one,
		      ((double) population) / molarConst);
      return result;
    }
  };

  template<class scalarType>
  class insertDelta :
    public std::unary_function<const std::pair<int, int>&, void>
  {
    polymap<scalarType>& rConverter;
    monomial linMonom;
  
  public:
    insertDelta(polymap<scalarType>& rConverterMap,
		int reactionIndex) :
      rConverter(rConverterMap)
    {
      linMonom.setExponent(reactionIndex,
			   1);
    }

    void
    operator()(const std::pair<int, int>& rDelta) const
    {
      // An alternative would be to construct a polynomial and add
      // polynomials.
      scalarType newCoeff
	= rConverter[rDelta.first].getCoeff(linMonom)
	+ ((scalarType) rDelta.second);

      rConverter[rDelta.first].setCoeff(linMonom,
					newCoeff);
    }
  };

  // Multiplies derivative of a reactionAmount by the concntration
  class timesSubstrateRate :
    public std::unary_function<const std::pair<int, int>&, void>
  {
    const polymap<double>& rConverter;
    polynomial<double>& rDeriv;
  public:
    timesSubstrateRate(const polymap<double>& rRateConverter,
		       polynomial<double>& rDerivative) :
      rConverter(rRateConverter),
      rDeriv(rDerivative)
    {}

    void
    operator()(const std::pair<int, int>& rSubstrateMultPair) const
    {
      rConverter[rSubstrateMultPair.first].
	powerMultiply(rSubstrateMultPair.second,
		      rDeriv);
    }
  };

  class makeReactionDeriv :
    public std::unary_function<parserReaction&, polynomial<double> >
  {
    // The derivative is basically a product of these polynomials.
    // Each one gives the molar concentration of a reactant
    // in terms of the "molar concentration" of each reaction,
    // the reaction amount.
    const polymap<double>& rConverter;
  public:
    makeReactionDeriv(const polymap<double>& rRateConverter) :
      rConverter(rRateConverter)
    {}

    polynomial<double>
    operator()(const parserReaction& rParserReaction) const
    {
      // Start out with by making the derivative a constant,
      // equal to the reaction rate constant.
      polynomial<double> deriv;
      deriv.setCoeff(monomial(),
		     rParserReaction.rate);

      // Now multiply the derivative by powers of the linear polynomials
      // that give the current molar concentrations.  The power is given by
      // the multiplicity of the population in the substrates list of the
      // reaction.
      for_each(rParserReaction.substrates.begin(),
	       rParserReaction.substrates.end(),
	       timesSubstrateRate(rConverter,
				  deriv));
      return deriv;
    }
  };

  void
  tauApp::makeDerivatives(std::vector<parserReaction>& rReactions)
  {
    // Construct the linear polynomial that converts reaction concentrations
    // to species concentrations.
    //
    // Note that this polynomial is related to the integral linear
    // polynomial that converts reaction counts to population counts,
    // but they are not scalar multiples:  The constant terms are related
    // by the molar factor, but the coefficients of the linear terms are
    // the same (integers).
    polymap<double> rateConverterPoly;
    makeRateConverter(rateConverterPoly,
		      rReactions);

    transform(rReactions.begin(),
	      rReactions.end(),
	      back_inserter(derivatives),
	      makeReactionDeriv(rateConverterPoly));
  }

  void
  tauApp::makeRateConverter(polymap<double>& rConverterPoly,
			    std::vector<parserReaction>& rReactions)
  {
    // The constant terms are given by the inital concentrations,
    // and we have the initial populations.  Hence, we need to scale
    // while we extend scalars.
    transform(populations.begin(),
	      populations.end(),
	      back_inserter(rConverterPoly),
	      pop2concPoly(molarConst));

    // Traverse the reaction deltas, inserting linear terms into the
    // integer polynomial.  This has to be done by index.
    int reactionNdx = rReactions.size();
    while(0 < reactionNdx--)
      {
	for_each(rReactions[reactionNdx].deltas.begin(),
		 rReactions[reactionNdx].deltas.end(),
		 insertDelta<double>(rConverterPoly,
				     reactionNdx));
      }
  }

  void
  tauApp::makeCountConverter(std::vector<parserReaction>& rReactions)
  {
    // The constant terms are given by the initial populations.
    //
    // This is the first time that the number of species becomes available,
    // so the converterPoly polynomials are pushed onto the vector.
    transform(populations.begin(),
	      populations.end(),
	      back_inserter(popsForCounts),
	      int2intPoly);

    // Traverse the reaction deltas, inserting linear terms into the
    // integer polynomial.  This has to be done by index.
    int reactionNdx = rReactions.size();
    while(0 < reactionNdx--)
      {
	for_each(rReactions[reactionNdx].deltas.begin(),
		 rReactions[reactionNdx].deltas.end(),
		 insertDelta<int>(popsForCounts,
				  reactionNdx));
      }
  }
}

