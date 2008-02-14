/////////////////////////////////////////////////////////////////////////////
// rk4tau - an accelerated stochastic simulator.
// Copyright (C) 2002  Walter Lawrence (Larry) Lok.
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

#ifndef MONOMIAL_H
#define MONOMIAL_H

#include <vector>
#include <algorithm>
#include <map>

namespace rk4util
{
  // Forward definition of polynomial; partial derivatives of
  // monomials are polynomials.
  template<class scalarType> class polynomial;
  
  // Monomials represented as map from int (variable index) to int
  // (exponent).
  class monomial :
    public std::map<int, int>
  {
    template<class scalarType>
    class evalVarb :
      public std::unary_function<std::pair<int, int>, void>
    {
      const std::vector<scalarType>& rValues;
      scalarType& rAccum;
    public:
      evalVarb(const std::vector<scalarType>& rVarValues,
	       scalarType& rAccumulatedValue) :
	rValues(rVarValues),
	rAccum(rAccumulatedValue)
      {}

      void
      operator()(const std::pair<int, int>& rMonomialEntry) const
      {
	// Note that this will fail if there are variables that are out
	// whose indices are out of range of the vector of values.
	int exponent = rMonomialEntry.second;
	while(0 < exponent--) rAccum = rAccum * rValues[rMonomialEntry.first];
      }
    };

    // Auxiliary function for multiplying one monomial by another.
    class multiplyByVar :
      public std::unary_function<const std::pair<int, int>&, void>
    {
      monomial& rTarget;
    public:
      multiplyByVar(monomial& rTargetMonomial) :
	rTarget(rTargetMonomial)
      {}

      void
      operator()(const std::pair<int, int>& rVar) const
      {
	rTarget.setExponent(rVar.first,
			    rTarget.getExponent(rVar.first) + rVar.second);
      }
    };

  public:

    int
    getExponent(int varNdx) const
    {
      const_iterator iEntry = find(varNdx);
      return iEntry == end()
	? 0
	: iEntry->second;
    }

    void
    setExponent(int varNdx,
		int exponent)
    {
      std::pair<iterator, bool> insertResult
	= insert(std::make_pair(varNdx, exponent));

      if(! insertResult.second)
	insertResult.first->second = exponent;
    }

    void
    addToExponent(int varNdx,
		  int exponent)
    {
      std::pair<iterator, bool> insertResult
	= insert(std::make_pair(varNdx, exponent));

      if(! insertResult.second)
	insertResult.first->second += exponent;
    }

    template<class scalarType>
    scalarType
    evaluate(const std::vector<scalarType>& rVarValues) const
    {
      scalarType result = (scalarType) 1;
      for_each(begin(),
	       end(),
	       evalVarb<scalarType>(rVarValues,
				    result));
      return result;
    }

    // Destructively multiplies the target monomial by this
    // monomial.
    //
    // It's actually a little trickier than you'd think to multiply two
    // monomials and put the result into a third monomial; it's actually
    // easier to copy the other vactor into the target, then multiply the
    // target this way.
    //
    // Note that, again when push comes to shove, addition of polynomials
    // does something just like this.
    void
    multiplyBy(monomial& rTargetMonomial) const
    {
      std::for_each(begin(),
		    end(),
		    multiplyByVar(rTargetMonomial));
    }

    // For completeness.
    monomial
    multiply(const monomial& rOtherFactor) const
    {
      monomial result(rOtherFactor);

      multiplyBy(result);

      return result;
    }

    void
    multiply(const monomial& rOtherFactor,
	     monomial& rTargetMonomial) const
    {
      rTargetMonomial = rOtherFactor;
      multiplyBy(rTargetMonomial);
    }

    // Assumes that the target polynomial is zero to begin with.
    //
    // Note that dx/dx = 1x^0, not just dx/dx = 1; variables with
    // zero exponent are still included in the map.
    //
    // This is defined in polynomial.hh for now.
    template<class scalarType>
    void
    partialDerivative(int varNdx,
		      polynomial<scalarType>& rTarget) const;
  };
}

#endif // MONOMIAL_H
