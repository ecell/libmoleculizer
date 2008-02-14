/////////////////////////////////////////////////////////////////////////////
// rk4 - a 4th order adaptive Runge Kutta solver for polynomial ODE's.
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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <map>
#include <vector>
#include <functional>
#include "rk4util/monomial.hh"

namespace rk4util
{
  template<class scalarType> class polymap;
  
  template<class scalarType>
  class polynomial :
    public std::map<monomial, scalarType>
  {
  public:

    typedef
    typename std::map<monomial, scalarType>::iterator
    iterator;

    typedef
    typename std::map<monomial, scalarType>::const_iterator
    const_iterator;

    typedef
    typename std::map<monomial, scalarType>::value_type
    term_type;

  private:

    // Evaluates a term of this polynomial, given the values of the variables.
    template<class extensionType>
    class evalTerm :
      public std::unary_function<term_type, void>
    {
      const std::vector<extensionType>& rValues;
      extensionType& rResult;
    public:
      evalTerm(const std::vector<extensionType>& rVarValues,
	       extensionType& rAccumResult) :
	rValues(rVarValues),
	rResult(rAccumResult)
      {}

      void
      operator()(const term_type& rTerm) const
      {
	rResult += rTerm.second * rTerm.first.evaluate(rValues);
      }
    };

    // Adds a term (as a polynomial) to a polynomial with extended scalars.
    template<class extensionType>
    class addTermToTarget :
      public std::unary_function<term_type, void>
    {
      polynomial<extensionType>& rTarget;
    public:
      addTermToTarget(polynomial<extensionType>& rTargetPoly) :
	rTarget(rTargetPoly)
      {}

      void
      operator()(const term_type& rTerm) const
      {
	rTarget.setCoeff(rTerm.first,
			 rTarget.getCoeff(rTerm.first) + rTerm.second);
      }
    };

    // For multiplying a term of this polynomial by a term of a polynomial
    // with extended scalars.
    template<class extensionType>
    class multiplyTermByTerm :
      public std::unary_function<typename polynomial<extensionType>::term_type, void>
    {
      typedef
      typename polynomial<extensionType>::term_type
      extensionTerm;
      
      const term_type& rLeft;
      polynomial<extensionType>& rTarget;
    public:
      multiplyTermByTerm(const term_type& rLeftTerm,
			 polynomial<extensionType>& rTargetPoly) :
	rLeft(rLeftTerm),
	rTarget(rTargetPoly)
      {}

      void
      operator()(const extensionTerm& rRight) const
      {
	monomial productMonomial(rRight.first);
	rLeft.first.multiplyBy(productMonomial);

	extensionType coeff = rTarget.getCoeff(productMonomial)
	  + (rLeft.second * rRight.second);

	rTarget.setCoeff(productMonomial,
			 coeff);
      }
    };

    // This is used in taking partial derivatives of polynomials below.
    template<class extensionType>
    class multiplyTermByScalar :
      public std::unary_function<
      std::pair<monomial, scalarType>,
      std::pair<monomial, extensionType> >
    {

      extensionType scalar;
    public:
      multiplyTermByScalar(extensionType theScalar) :
	scalar(theScalar)
      {}

      std::pair<monomial, extensionType>
      operator()(const std::pair<monomial, scalarType>& rTerm) const
      {
	return std::make_pair(rTerm.first,
			      scalar * (extensionType) rTerm.second);
      }
    };

    // For multiplying a term of this polynomial by a polynomial with extended
    // scalars.
    template<class extensionType>
    class multiplyByTerm :
      public std::unary_function<term_type&, void>
    {
      const polynomial<extensionType>& rOther;
      polynomial<extensionType>& rTarget;
    public:
      multiplyByTerm(const polynomial<extensionType>& rOtherFactor,
		     polynomial<extensionType>& rTargetPolynomial) :
	rOther(rOtherFactor),
	rTarget(rTargetPolynomial)
      {}

      void
      operator()(const term_type& rTerm) const
      {
	for_each(rOther.begin(),
		 rOther.end(),
		 multiplyTermByTerm<extensionType>(rTerm,
						   rTarget));
      }
    };

    // Extend scalars on one term of this polynomial.
    template<class extensionType>
    class extendTermScalars :
      public std::unary_function<term_type, void>
    {
      polynomial<extensionType>& rTarget;
    
    public:
      extendTermScalars(polynomial<extensionType>& rTargetPoly) :
	rTarget(rTargetPoly)
      {}

      void
      operator()(const term_type& rTerm) const
      {
	rTarget.insert(std::pair<monomial, extensionType>(rTerm.first,
							  rTerm.second));
      }
    };

  public:

    scalarType
    getCoeff(const monomial& rMonomial) const
    {
      const_iterator iEntry = this->find(rMonomial);

      return iEntry == this->end()
	? (scalarType) 0
	: iEntry->second;
    }

    void
    setCoeff(const monomial& rMonomial,
	     scalarType coefficient)
    {
      std::pair<iterator, bool> insertResult
	= insert(std::make_pair(rMonomial, coefficient));

      if(! insertResult.second)
	{
	  insertResult.first->second = coefficient;
	}
    }

    void
    addToCoeff(const monomial& rMonomial,
	       scalarType coefficient)
    {
      std::pair<iterator, bool> insertResult
	= insert(std::make_pair(rMonomial, coefficient));

      // Note that this requires scalarType to implement operator+=.
      if(! insertResult.second)
	insertResult.first->second += coefficient;
    }

    // Note that we are still missing some forms of the following functions:
    // those that apply this polynomial to a vector (or whatever) of with
    // MORE RESTRICTED scalars, returning a vector (or whatever) with
    // scalars of type scalarType.

    template<class extensionType>
    extensionType
    evaluate(const std::vector<extensionType>& rVarValues) const
    {
      extensionType result(0);
      
      for_each(this->begin(),
	       this->end(),
	       evalTerm<extensionType>(rVarValues,
				       result));
      return result;
    }

    // Destructively add this polynomial to another polynomial.
    template<class extensionType>
    void
    addTo(polynomial<extensionType>& rTargetPoly) const
    {
      for_each(this->begin(),
	       this->end(),
	       addTermToTarget<extensionType>(rTargetPoly));
    }

    template<class extensionType>
    polynomial<extensionType>
    add(const polynomial<extensionType>& rOtherPoly) const
    {
      polynomial<extensionType> result(rOtherPoly);
      addTo(result);
      return result;
    }

    // This assumes that the target polynomial is zero to begin with.
    template<class extensionType>
    void
    multiply(const polynomial<extensionType>& rOtherFactor,
	     polynomial<extensionType>& rTargetPoly) const
    {
      for_each(this->begin(),
	       this->end(),
	       multiplyByTerm<extensionType>(rOtherFactor,
					     rTargetPoly));
    }

    template<class extensionType>
    polynomial<extensionType>
    multiply(const polynomial<extensionType>& rOtherFactor) const
    {
      polynomial<extensionType> result;
      multiply(rOtherFactor,
	       result);
      return result;
    }

    template<class extensionType>
    void
    multiplyByScalar(extensionType theScalar,
		     polynomial<extensionType>& rTarget) const
    {
      multiplyTermByScalar<extensionType> mm(theScalar);
      
      std::transform(this->begin(),
		     this->end(),
		     std::inserter(rTarget, rTarget.begin()),
		     mm);
    }

    template<class extensionType>
    polynomial<extensionType>
    multiplyByScalar(extensionType theScalar)
    {
      polynomial<extensionType> result;
      multiplyByScalar(theScalar,
		       result);
      return result;
    }

    // Destructively multiplies rTargetPoly by a power, assumed to be
    // >=0 of course, of this polynomial.
    //
    // Note that this operation is on the order of the exponent in
    // time. In this situation, we expect the exponent almost always to
    // be 1, so simplicity is favored over worst-case speed.
    template<class extensionType>
    void
    powerMultiply(int thePower,
		  polynomial<extensionType>& rTargetPoly) const
    {
      while(0 < thePower--) rTargetPoly = multiply(rTargetPoly);
    }

    // This is to support conversion from polynomial<int> to
    // polynomial<double>.
    template<class extensionType>
    void
    extendScalars(polynomial<extensionType>& rTargetPoly) const
    {
      for_each(this->begin(),
	       this->end(),
	       extendTermScalars<extensionType>(rTargetPoly));
    }

    template<class extensionType>
    polynomial<extensionType>
    extendScalars(void) const
    {
      polynomial<extensionType> result;
      extendScalars(result);
      return result;
    }

    template<class extensionType>
    void
    partialDerivative(int varNdx,
		      polynomial<extensionType>& rTarget) const;

    template<class extensionType>
    polynomial<extensionType>
    partialDerivative(int varNdx)
    {
      polynomial<extensionType> result;
      partialDerivative(varNdx,
			result);
      return result;
    }

    // Defined in polymap.hh, so that a special inline header isn't needed
    // for polynomials.
    template<class extensionType>
    void
    gradient(int varCount,
	     polymap<extensionType>& rTarget) const;

//     template<class extensionType>
//     std::vector<polynomial<extensionType> >
//     gradient(varCount)
//     {
//       std::vector<polynomial<extensionType> > result(varCount);
//       gradient(varCount,
// 	       result);
//       return result;
//     }
  };

  // Putting this here avoids having to have a special inline header for
  // monomials.
  template<class scalarType>
  void
  monomial::partialDerivative(int varNdx,
			      polynomial<scalarType>& rTarget) const
  {
    int varExponent = getExponent(varNdx);
    if(0 < varExponent)
      {
	monomial derivedMonomial(*this);
	derivedMonomial.setExponent(varNdx, varExponent - 1);

	scalarType coefficient = (scalarType) varExponent;

	rTarget.setCoeff(derivedMonomial,
			 coefficient);
      }
  }

  template<class scalarType, class extensionType>
  class partialOfPolyTerm :
    public std::unary_function<std::pair<const monomial, scalarType>, void>
  {
    int varNdx;
    polynomial<extensionType>& rTarget;

  public:
    partialOfPolyTerm(int variableIndex,
		      polynomial<extensionType>& rTargetPoly) :
      varNdx(variableIndex),
      rTarget(rTargetPoly)
    {}

    void
    operator()(const std::pair<const monomial, scalarType>& rEntry) const
    {
      polynomial<scalarType> monomDeriv;
      rEntry.first.partialDerivative(varNdx,
				     monomDeriv);

      polynomial<extensionType> termDeriv;
      monomDeriv.multiplyByScalar((extensionType) rEntry.second,
				  termDeriv);

      termDeriv.addTo(rTarget);
    }
  };

  template<class scalarType>
  template<class extensionType>
  void
  polynomial<scalarType>
  ::partialDerivative(int varNdx,
		      polynomial<extensionType>& rTarget) const
  {
    std::for_each(this->begin(),
		  this->end(),
		  partialOfPolyTerm<scalarType, extensionType>(varNdx,
							       rTarget));
  }
}

#endif // POLYNOMIAL_H
