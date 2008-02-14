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

#ifndef UTL_GSL_H
#define UTL_GSL_H

#include <vector>

// Unfortunately, gsl_vector and gsl_matrix are typedefs, so that
// it's not possible to use a forward declaration here.  Since the whole
// headers have to be included, all of this could as well be inlined.
//
// Leaving it as it would be if things were better and we didn't have
// to include the whole headers here.
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

namespace utl
{
  namespace gsl
  {
    class autoGslVector
    {
      gsl_vector *pVector;

    public:
      autoGslVector(int coordCount,
		    bool initZero = false);

      autoGslVector(const autoGslVector& rOriginal);
      
      ~autoGslVector(void);

      int
      size(void) const;

      // Note that this does not resize *this; the vectors must already be
      // of the same size.  This is lame and probably useless.
      autoGslVector&
      operator=(const autoGslVector& rSource);

      gsl_vector*
      getGslVector(void)
      {
	return pVector;
      }

      const gsl_vector*
      getGslVector(void) const
      {
	return pVector;
      }

      double&
      at(int ndx);

      const double&
      at(int ndx) const;

      void
      add(const autoGslVector& rOtherVector,
	  autoGslVector& rSum) const;

      autoGslVector
      add(const autoGslVector& rOtherVector) const
      {
	autoGslVector result(size());

	add(rOtherVector,
	    result);

	return result;
      }

      void
      scale(double scalar,
	    autoGslVector& rScaled) const;

      autoGslVector
      scale(double scalar) const
      {
	autoGslVector result(size());

	scale(scalar,
	      result);

	return result;
      }

      double
      boxNorm(void) const;
    };

    class autoGslMatrix
    {
      gsl_matrix* pMatrix;

    public:
      autoGslMatrix(int rowCount,
		    int colCount,
		    bool initZero = false);

      autoGslMatrix(const autoGslMatrix& rOriginal);

      ~autoGslMatrix(void);

      autoGslMatrix&
      operator=(const autoGslMatrix& rOriginal);

      int
      getRowCount(void) const;

      int
      getColCount(void) const;

      gsl_matrix*
      getGslMatrix(void)
      {
	return pMatrix;
      }

      const gsl_matrix*
      getGslMatrix(void) const
      {
	return pMatrix;
      }

      void
      setIdentity(void);

      void
      getColumnVector(int colNdx,
		      autoGslVector& rTarget) const;

      autoGslVector
      getColumnVector(int colNdx) const;

      void
      getRowVector(int rowNdx,
		   autoGslVector& rTarget) const;

      autoGslVector
      getRowVector(int rowNdx) const;

      double&
      at(int rowNdx,
	 int colNdx);
    
      const double&
      at(int rowNdx,
	 int colNdx) const;
    
      double
      boxNorm(void) const;

      void
      apply(const autoGslVector& rVector,
	    autoGslVector& rResult) const;

      autoGslVector
      apply(const autoGslVector& rVector) const
      {
	autoGslVector result(getRowCount(),
			     true);
	apply(rVector,
	      result);

	return result;
      }

      void
      plus(const autoGslMatrix& rOtherTerm,
	   autoGslMatrix& rResult) const;

      autoGslMatrix
      plus(const autoGslMatrix& rOtherTerm)
      {
	autoGslMatrix result(getRowCount(),
			     getColCount(),
			     true);
	plus(rOtherTerm,
	     result);

	return result;
      }

      void
      scale(double scalar,
	    autoGslMatrix& rResult) const;

      autoGslMatrix
      scale(double scalar) const
      {
	autoGslMatrix result(getRowCount(),
			     getColCount());
	scale(scalar,
	      result);

	return result;
      }

      void
      times(const autoGslMatrix& rOtherFactor,
	    autoGslMatrix& rResult) const;

      autoGslMatrix
      times(const autoGslMatrix& rOtherFactor) const
    
      {
	autoGslMatrix result(getRowCount(),
			     rOtherFactor.getColCount());
	times(rOtherFactor,
	      result);

	return result;
      }

      void
      taylorExp(double t,
		double epsilon,
		autoGslMatrix& rResultMatrix) const;

      autoGslMatrix
      taylorExp(double t,
		double epsilon) const
      {
	autoGslMatrix result(getRowCount(),
			     getColCount());
	taylorExp(t,
		  epsilon,
		  result);

	return result;
      }

      void
      exp(double t,
	  double epsilon,
	  autoGslMatrix& rResultMatrix) const;

      autoGslMatrix
      exp(double t,
	  double epsilon) const
      {
	autoGslMatrix result(getRowCount(),
			     getColCount());
	exp(t,
	    epsilon,
	    result);

	return result;
      }
    };

    class autoGslRng
    {
      gsl_rng* pRng;

    public:
      autoGslRng(void);

      ~autoGslRng(void);

      void
      seed(unsigned long int seed);

      gsl_rng*
      getGslRng(void)
      {
	return pRng;
      }
    };

    // For the time being, going along with my earlier uniform sampler's
    // interface: returns uniform on (0, 1], as if one were about to take a
    // logarithm.
    class uniformSampler
    {
      autoGslRng& rRng;
    public:
      uniformSampler(autoGslRng& rAutoGslRng) :
	rRng(rAutoGslRng)
      {}

      double
      sample(void);
    };

    class exponentialSampler
    {
      autoGslRng& rRng;
    public:
      exponentialSampler(autoGslRng& rAutoGslRng) :
	rRng(rAutoGslRng)
      {}

      double
      sample(double rate);
    };

    class poissonSampler
    {
      autoGslRng& rRng;
    public:
      poissonSampler(autoGslRng& rAutoGslRng) :
	rRng(rAutoGslRng)
      {}

      unsigned int
      sample(double rate);
    };

    // This approach avoids using variable-length arrays for p and s in the
    // code for "sample," which would be okay in g++, but...
    class multinomialSampler
    {
      autoGslRng& rRng;

      int altCount;
      double* p;
      unsigned int* s;

    public:
      multinomialSampler(int alternativeCount,
			 autoGslRng& rAutoGslRng) :
	rRng(rAutoGslRng),
	altCount(alternativeCount),
	p(new double[alternativeCount]),
	s(new unsigned int[alternativeCount])
      {}

      ~multinomialSampler(void)
      {
	delete [] p;
	delete [] s;
      }

      // rProbabilities should have size altCount, given in the
      // constructor above, as should rSample.  This restriction is unchecked.
      void
      sample(const autoGslVector& rProbabilities,
	     int sampleSize,
	     std::vector<int>& rSample) const;
    };
  }
}

#endif // UTL_GSL_H
