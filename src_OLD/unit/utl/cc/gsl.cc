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

#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include "utl/gsl.hh"

namespace utl
{
  namespace gsl
  {
    autoGslVector::
    autoGslVector(int coordCount,
		  bool initZero)
    {
      if(initZero)
	{
	  pVector = gsl_vector_alloc(coordCount);
	}
      else
	{
	  pVector = gsl_vector_calloc(coordCount);
	}
    }

    autoGslVector::
    autoGslVector(const autoGslVector& rOriginal)
    {
      pVector = gsl_vector_alloc(rOriginal.size());
      gsl_vector_memcpy(pVector, rOriginal.getGslVector());
    }

    autoGslVector::
    ~autoGslVector(void)
    {
      gsl_vector_free(pVector);
    }

    int
    autoGslVector::
    size(void) const
    {
      return pVector->size;
    }

    autoGslVector&
    autoGslVector::
    operator=(const autoGslVector& rSource)
    {
      gsl_vector_memcpy(pVector,
			rSource.getGslVector());
      return *this;
    }

    double&
    autoGslVector::
    at(int ndx)
    {
      return *(gsl_vector_ptr(pVector,
			      ndx));
    }

    const double&
    autoGslVector::
    at(int ndx) const
    {
      return *(gsl_vector_const_ptr(pVector,
				    ndx));
    }

    void
    autoGslVector::
    add(const autoGslVector& rOtherVector,
	autoGslVector& rSum) const
    {
      // For some reason, gsl gives destination first, then source.
      gsl_vector_memcpy(rSum.getGslVector(),
			rOtherVector.getGslVector());

      // Again, destination first, source second.
      gsl_vector_add(rSum.getGslVector(),
		     pVector);
    }

    void
    autoGslVector::
    scale(double scalar,
	  autoGslVector& rScaled) const
    {
      gsl_vector_memcpy(rScaled.getGslVector(),
			getGslVector());

      gsl_vector_scale(rScaled.getGslVector(),
		       scalar);
    }

    double
    autoGslVector::
    boxNorm(void) const
    {
      
      double maxCoord = 0.0;
      double minCoord = 0.0;

      gsl_vector_minmax(pVector, &minCoord, &maxCoord);

      return ((0.0 < minCoord) || (0.0 < maxCoord))
	? maxCoord
	: - minCoord;
    }

    autoGslMatrix::
    autoGslMatrix(int rowCount,
		  int colCount,
		  bool initZero)
    {
      if(initZero)
	{
	  pMatrix = gsl_matrix_calloc(rowCount,
				      colCount);
	}
      else
	{
	  pMatrix = gsl_matrix_alloc(rowCount,
				     colCount);
	}
    }

    autoGslMatrix::
    autoGslMatrix(const autoGslMatrix& rOriginal)
    {
      pMatrix = gsl_matrix_alloc(rOriginal.getRowCount(),
				 rOriginal.getColCount());
      gsl_matrix_memcpy(pMatrix,
			rOriginal.getGslMatrix());
    }

    autoGslMatrix::
    ~autoGslMatrix(void)
    {
      gsl_matrix_free(pMatrix);
    }

    autoGslMatrix&
    autoGslMatrix::
    operator=(const autoGslMatrix& rOriginal)
    {
      gsl_matrix_memcpy(pMatrix,
			rOriginal.getGslMatrix());
      return *this;
    }

    int
    autoGslMatrix::
    getRowCount(void) const
    {
      return pMatrix->size1;
    }

    int
    autoGslMatrix::
    getColCount(void) const
    {
      return pMatrix->size2;
    }

    void
    autoGslMatrix::
    setIdentity(void)
    {
      gsl_matrix_set_identity(pMatrix);
    }

    void
    autoGslMatrix::
    getColumnVector(int colNdx,
		    autoGslVector& rTarget) const
    {
      gsl_vector_const_view colView
	= gsl_matrix_const_column(pMatrix,
				  colNdx);

      gsl_vector_memcpy(rTarget.getGslVector(),
			&colView.vector);
    }

    autoGslVector
    autoGslMatrix::
    getColumnVector(int colNdx) const
    {
      autoGslVector result(getRowCount());
      getColumnVector(colNdx,
		      result);
      return result;
    }

    void
    autoGslMatrix::
    getRowVector(int rowNdx,
		 autoGslVector& rTarget) const
    {
      gsl_vector_const_view rowView
	= gsl_matrix_const_row(pMatrix,
			       rowNdx);

      gsl_vector_memcpy(rTarget.getGslVector(),
			&rowView.vector);
    }

    autoGslVector
    autoGslMatrix::
    getRowVector(int rowNdx) const
    {
      autoGslVector result(getColCount());
      getRowVector(rowNdx,
		   result);
      return result;
    }

    double&
    autoGslMatrix::
    at(int rowNdx,
       int colNdx)
    {
      return *(gsl_matrix_ptr(pMatrix,
			      rowNdx,
			      colNdx));
    }

    const double&
    autoGslMatrix::
    at(int rowNdx,
       int colNdx) const
    {
      return *(gsl_matrix_ptr(pMatrix,
			      rowNdx,
			      colNdx));
    }

    double
    autoGslMatrix::
    boxNorm(void) const
    {
      
      double maxCoord = 0.0;
      double minCoord = 0.0;

      gsl_matrix_minmax(getGslMatrix(),
			&minCoord,
			&maxCoord);

      return ((0.0 < minCoord) || (0.0 < maxCoord))
	? maxCoord
	: - minCoord;
    }

    void
    autoGslMatrix::
    apply(const autoGslVector& rVector,
	  autoGslVector& rResult) const
    {
      gsl_blas_dgemv(CblasNoTrans,
		     1.0,
		     pMatrix,
		     rVector.getGslVector(),
		     1.0,
		     rResult.getGslVector());
    }

    void
    autoGslMatrix::
    plus(const autoGslMatrix& rOtherTerm,
	 autoGslMatrix& rResult) const
    {
      rResult = rOtherTerm;
      gsl_matrix_add (rResult.getGslMatrix(),
		      getGslMatrix());
    }

    void
    autoGslMatrix::
    scale(double scalar,
	  autoGslMatrix& rResult) const
    {
      rResult = *this;
      gsl_matrix_scale(rResult.getGslMatrix(),
		       scalar) ;
    }

    void
    autoGslMatrix::
    times(const autoGslMatrix& rOtherFactor,
	  autoGslMatrix& rResult) const
    {
      gsl_blas_dgemm(CblasNoTrans,
		     CblasNoTrans,
		     1.0,
		     getGslMatrix(),
		     rOtherFactor.getGslMatrix(),
		     0.0,
		     rResult.getGslMatrix());
    }

    void
    autoGslMatrix::
    taylorExp(double t,
	      double epsilon,
	      autoGslMatrix& rResultMatrix) const
    {
      rResultMatrix.setIdentity();
      autoGslMatrix termMatrix(rResultMatrix);
      autoGslMatrix scaledTermMatrix(getRowCount(),
				     getColCount());
      int n = 0;
      double termMatrixNorm = termMatrix.boxNorm();
      double resultMatrixNorm = rResultMatrix.boxNorm();
      while(epsilon <  (termMatrixNorm / resultMatrixNorm))
	{
	  termMatrix.scale(t / ((double) ++n),
			   scaledTermMatrix);
	  times(scaledTermMatrix,
		termMatrix);
	  rResultMatrix = rResultMatrix.plus(termMatrix);

	  termMatrixNorm = termMatrix.boxNorm();
	  resultMatrixNorm = rResultMatrix.boxNorm();
	}
    }

    void
    autoGslMatrix::
    exp(double t,
	double epsilon,
	autoGslMatrix& rResultMatrix) const
    {
      // Get a power of 2 large enough to reduce our matrix
      // to a small size.
      int base2exponent = 0;
      frexp(boxNorm(),
	    &base2exponent);
      if(base2exponent < 1) base2exponent = 1;

      double scalar = ldexp(1.0, -base2exponent);

      // Compute the exponential of the smaller matrix using Taylor series.
      autoGslMatrix smallerMatrix = scale(scalar);
      rResultMatrix = smallerMatrix.taylorExp(t,
					      epsilon);

      // Square the exponential of the smaller matrix repeatedly
      // to get the exponential of the original matrix.
      while(0 < base2exponent--)
	{
	  rResultMatrix = rResultMatrix.times(rResultMatrix);
	}
    }

    autoGslRng::
    autoGslRng(void) :
      pRng(gsl_rng_alloc(gsl_rng_taus))
    {}

    autoGslRng::
    ~autoGslRng(void)
    {
      gsl_rng_free(pRng);
    }

    void
    autoGslRng::
    seed(unsigned long int seed)
    {
      gsl_rng_set(pRng, seed);
    }

    double
    uniformSampler::
    sample(void)
    {
      return 1.0 - gsl_ran_flat(rRng.getGslRng(),
				0.0, 
				1.0);
    }

    double
    exponentialSampler::
    sample(double rate)
    {
      return gsl_ran_exponential(rRng.getGslRng(), 
				 rate);
    }

    unsigned int
    poissonSampler::
    sample(double rate)
    {
      return gsl_ran_poisson(rRng.getGslRng(),
			     rate);
    }

    void
    multinomialSampler::
    sample(const autoGslVector& rProbabilities,
	   int sampleSize,
	   std::vector<int>& rSample) const
    {
      int altNdx = altCount;
      while(0 < altNdx--)
	{
	  p[altNdx] = rProbabilities.at(altNdx);
	}

      gsl_ran_multinomial(rRng.getGslRng(),
			  altCount,
			  sampleSize,
			  p,
			  s);
      std::copy(s,
		s + altCount,
		rSample.begin());
    }
  }
}
