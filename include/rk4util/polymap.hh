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

#ifndef POLYMAP_H
#define POLYMAP_H

#include "rk4util/polynomial.hh"

namespace rk4util
{
  // A polynomial mapping from several variables to several variables.
  template<class scalarType>
  class polymap :
    public std::vector<polynomial<scalarType> >
  {
    // Auxiliary class used to evaluate a polymap.
    template<class extensionType>
    class evaluatePoly :
      public std::unary_function<polynomial<scalarType>, extensionType>
    {
      const std::vector<extensionType>& rArgument;
    public:
      evaluatePoly(const std::vector<extensionType>& rArgumentVector) :
	rArgument(rArgumentVector)
      {}

      extensionType
      operator()(const polynomial<scalarType>& rPoly) const
      {
	return rPoly.evaluate(rArgument);
      }
    };

    // Auxiliar class used to extend scalars of a polymap.
    template<class extensionType>
    class extendPolyScalars : public
    std::unary_function<polynomial<scalarType>, polynomial<extensionType> >
    {
    public:
      polynomial<extensionType>
      operator()(const polynomial<scalarType>& rSourcePoly)
      {
	polynomial<extensionType> result;
	rSourcePoly.extendScalars(result);
	return result;
      }
    };

  public:

    // Default constructor has no range variables; polynomials
    // must be inserted.
    polymap()
    {}

    // For constructing a polynomial map with a specific
    // number of variables in the range.
    polymap(int varCount) :
      std::vector<polynomial<scalarType> >(varCount)
    {}

    // Note that this assumes that the target vector has the same length
    // as the number of polynomials.  It does not push the coordinates.
    template<class extensionType>
    void
    evaluate(const std::vector<extensionType>& rVariableValues,
	     std::vector<extensionType>& rTargetVector) const
    {
      transform(this->begin(),
		this->end(),
		rTargetVector.begin(),
		evaluatePoly<extensionType>(rVariableValues));
    }

    template<class extensionType>
    std::vector<extensionType>
    evaluate(const std::vector<extensionType>& rVariableValues) const
    {
      // The result has as many coordinates as there are polynomials
      // in this vector.
      std::vector<extensionType> result(this->size());
      
      evaluate(rVariableValues,
	       result);

      return result;
    }

    // This assumes that the target vector of polynomials is empty, so
    // that each coordinate polynomial is pushed onto the end of the vector.
    template<class extensionType>
    void
    extendScalars(polymap<extensionType>& rConverted) const
    {
      // There seems to be a g++ bug that prevents taking a pointer
      // to a template member function, so using
      //
      // mem_fun_ref(&polynomial<scalarType>::extendScalars<extensionType>)
      //
      // signals a parse error.
      std::transform(this->begin(),
		     this->end(),
		     std::back_inserter(rConverted),
		     extendPolyScalars<extensionType>());
    }

    template<class extensionType>
    polymap<extensionType>
    extendScalars(void) const
    {
      polymap<extensionType> result;
      extendScalars(result);
      return result;
    }

    template<class extensionType>
    std::vector<polymap<extensionType> >
    jacobian(int varCount);
  };

  // This is here, rather than in polynomial.hh so that I don't have
  // to have a special inline header for polynomials.
  //
  // Note that the target polymap has to have the correct size.
  template<class scalarType>
  template<class extensionType>
  void
  polynomial<scalarType>::gradient(int varCount,
				   polymap<extensionType>& rTarget) const
  {
    for(int varNdx = 0;
	varNdx < varCount;
	varNdx++)
      {
	partialDerivative(varNdx,
			  rTarget[varNdx]);
      }
  }

  template<class scalarType, class extensionType>
  class takeGradient :
    public std::unary_function<polynomial<scalarType>, polymap<extensionType> >
  {
    int varCount;
    
  public:
    takeGradient(int variableCount) :
      varCount(variableCount)
    {}

    polymap<extensionType>
    operator()(const polynomial<scalarType>& rPoly) const
    {
      // There should be a combined function to do this all at once.
      polymap<extensionType> result(varCount);
      rPoly.gradient(varCount,
			    result);
      return result;
    }
  };

  template<class scalarType>
  template<class extensionType>
  std::vector<polymap<extensionType> >
  polymap<scalarType>::jacobian(int varCount)
  {
    // The Jacobian should have as many rows as the polymap does.
    // Each row of the Jacobian is the gradient of the corresponding
    // row of the polymap.
    std::vector<polymap<extensionType> >
      result(this->size(),
	     polymap<extensionType>(varCount));

    std::transform(this->begin(),
		   this->end(),
		   result.begin(),
		   takeGradient<scalarType, extensionType>(varCount));

    return result;
  }
}


#endif // POLYMAP_H
