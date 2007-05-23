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

#ifndef LINALG_H
#define LINALG_H

#include<vector>
#include<functional>
#include<algorithm>
#include<numeric>
#include <vector>

namespace rk4util
{
  template<class scalarType, class extensionType>
  void
  extendScalars(const std::vector<scalarType>& rVector,
		std::vector<extensionType>& rTargetVector)
  {
    std::copy(rVector.begin(),
	      rVector.end(),
	      rTargetVector.begin());
  }

  template<class scalarType, class extensionType>
  std::vector<extensionType>
  extendScalars(const std::vector<scalarType>& rVector)
  {
    std::vector<extensionType> result(rVector.size());

    extendScalars(rVector,
		  result);

    return result;
  }

  template<class scalarType, class extensionType>
  void
  scalarMult(const std::vector<scalarType>& rVector,
	     extensionType extendedScalar,
	     std::vector<extensionType>& rTargetVector)
  {
    std::transform(rVector.begin(),
		   rVector.end(),
		   rTargetVector.begin(),
		   std::bind2nd(std::multiplies<extensionType>(),
				extendedScalar));
  }

  template<class scalarType, class extensionType>
  std::vector<extensionType>
  scalarMult(const std::vector<scalarType>& rVector,
	     extensionType extendedScalar)
  {
    std::vector<extensionType> result(rVector.size());

    scalarMult(rVector,
	       extendedScalar,
	       result);

    return result;
  }

  template<class scalarType, class extensionType> 
  void
  scalarMult(scalarType scalar,
	     const std::vector<extensionType>& rVector,
	     std::vector<extensionType>& rTargetVector)
  {
    std::transform(rVector.begin(),
		   rVector.end(),
		   rTargetVector.begin(),
		   std::bind2nd(std::multiplies<extensionType>(),
				scalar));
  }

  template<class scalarType, class extensionType>
  std::vector<extensionType>
  scalarMult(scalarType aScalar,
	     const std::vector<extensionType>& rVector)
  {
    std::vector<extensionType> result(rVector.size());

    scalarMult(aScalar,
	       rVector,
	       result);

    return result;
  }
		     
  template<class scalarType, class extensionType>
  void
  vectorAdd(const std::vector<scalarType>& rLeftVector,
	    const std::vector<extensionType>& rRightVector,
	    std::vector<extensionType>& rResultVector)
  {
    std::transform(rLeftVector.begin(),
		   rLeftVector.end(),
		   rRightVector.begin(),
		   rResultVector.begin(),
		   std::plus<extensionType>());
  }

  template<class scalarType, class extensionType>
  std::vector<extensionType>
  vectorAdd(const std::vector<scalarType>& rLeftVector,
	    const std::vector<extensionType>& rRightVector)
  {
    std::vector<extensionType> result(rLeftVector.size());

    vectorAdd(rLeftVector,
	      rRightVector,
	      result);

    return result;
  }

  // Due to an apparent bug in g++ resolution of overloaded operators,
  // something funny about templates and overloaded operators, or possibly
  // something funny about namespace std, letting std::inner_product
  // try to find the (possibly overloaded) operators + and * seems not
  // to work.
  //
  // For example, if I define multiplication (operator*)of a complex and an
  // integer in the global namespace, the invocation of std::inner_product
  // doesn't seem to find the def of operator*.  On the other hand if it
  // define the operator* in the std namespace, the template does seem to find
  // the def.  This is apparently erroneous, since if it is visible in std, it
  // should be visible in the global namespace.
  template<class scalarType, class extensionType>
  class dotProductOp :
    public std::binary_function<std::vector<scalarType>,
				std::vector<extensionType>,
				extensionType>
  {
    class doPlus :
      std::binary_function<extensionType, extensionType, extensionType>
    {
    public:
      extensionType
      operator()(const extensionType& rLeft,
		 const extensionType& rRight) const
      {
	return rLeft + rRight;
      }
    };

    class doTimes :
      public std::binary_function<scalarType, extensionType, extensionType>
    {
    public:
      extensionType
      operator()(const scalarType& rLeft,
		 const extensionType& rRight) const
      {
	return rLeft * rRight;
      }
    };
  
  public:
    extensionType
    operator()(const std::vector<scalarType>& rLeftVector,
	       const std::vector<extensionType>& rRightVector) const
    {
      return std::inner_product(rLeftVector.begin(),
				rLeftVector.end(),
				rRightVector.begin(),
				extensionType(0),
				doPlus(),
				doTimes());
    }
  };

  // This quite poor, but rarely used in this code, and I have better, but
  // more cumbersome, matrix and tensor types elsewhere.
  template<class scalarType>
  class linTrans : public std::vector<std::vector<scalarType> >
  {
    class extractColumn :
      public std::unary_function<std::vector<scalarType>, scalarType>
    {
      int colNdx;
    public:
      extractColumn(int columnNdx) :
	colNdx(columnNdx)
      {}

      scalarType
      operator()(const std::vector<scalarType>& rVector) const
      {
	return rVector[colNdx];
      }
    };

    class makeConstantRow :
      public std::unary_function<scalarType, std::vector<scalarType> >
    {
      int colCount;
	
    public:
      makeConstantRow(int columnCount) :
	colCount(columnCount)
      {}
	
      std::vector<scalarType>
      operator()(scalarType aScalar) const
      {
	return std::vector<scalarType>(aScalar, colCount);
      }
    };
    
  public:
    // Untrapped error if rowCnt <= 0 or colCnt <= 0.
    linTrans(int rowCnt,
	     int colCnt,
	     scalarType aScalar = 0) :
      std::vector<std::vector<scalarType> >
    (rowCnt,
     std::vector<scalarType>(colCnt, aScalar))
    {}

    linTrans(const std::vector<scalarType>& rRow,
	     int rowCnt) :
      std::vector<std::vector<scalarType> >(rRow, rowCount)
    {}

    linTrans(int colCnt,
	     const std::vector<scalarType>& rCol) :
      std::vector<std::vector<scalarType> >(rCol.size())
    {
      std::transform(rCol.begin(),
		     rCol.end(),
		     this->begin(),
		     makeConstantRow(colCnt));
    }

    int
    rowCount(void) const
    {
      return this->size();
    }

    int
    colCount(void) const
    {
      return (*this)[0].size();
    }

    const std::vector<scalarType>&
    row(int rowNdx) const
    {
      return (*this)[rowNdx];
    }

    void
    column(int columnNdx,
	   std::vector<scalarType>& rTargetVector) const
    {
      std::transform(this->begin(),
		     this->end(),
		     rTargetVector.begin(),
		     extractColumn(columnNdx));
    }

    std::vector<scalarType>
    column(int columnNdx) const
    {
      std::vector<scalarType> result(rowCount());

      column(columnNdx,
	     result);

      return result;
    }

    template<class extensionType>
    void
    apply(const std::vector<extensionType>& rVector,
	  std::vector<extensionType>& rTargetVector) const
    {
      std::transform(this->begin(),
		     this->end(),
		     rTargetVector.begin(),
		     std::bind2nd(dotProductOp<scalarType, extensionType>(),
				  rVector));
    }

    template<class extensionType>
    std::vector<extensionType>
    apply(const std::vector<extensionType>& rVector) const
    {
      std::vector<extensionType> result(rowCount());

      apply(rVector,
	    result);

      return result;
    }
  };
}

#endif // LINALG_H

