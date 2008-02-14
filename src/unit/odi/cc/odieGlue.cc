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

#include "odi/odieGlue.hh"

namespace odie
{
  // Global function always used by gsl ode solver.  The real function
  // is in the params, of course.  We're assuming the ode is autonomous,
  // so that we don't use t here.
  int
  evalDeriv(double t,
	    const double arg[],
	    double returnValue[],
	    void* pParams)
  {
    const odieParams* pOdieParams = (const odieParams*) pParams;
    const rk4util::polymap<double>& rDerivative = pOdieParams->rDeriv;
    
    // Construct an STL vector from the input array.
    const std::vector<double> argVector(arg,
					arg + rDerivative.size());

    // Apply the derivative polymap, which returns its value in a vector.
    std::vector<double> returnVector = rDerivative.evaluate(argVector);

    // Copy the return value's coordinates into the output array.
    std::copy(returnVector.begin(),
	      returnVector.end(),
	      returnValue);

    return GSL_SUCCESS;
  }

  int
  evalJacobian(double t,
	       const double arg[],
	       double* pJac,
	       double dfdt[],
	       void* pParams)
  {
    const odieParams* pOdieParams = (const odieParams*) pParams;
    const std::vector<rk4util::polymap<double> >& rJacobian
      = pOdieParams->rJac;

    int dimension = rJacobian.size();

    // Construct an STL vector from the input array.
    const std::vector<double> argVector(arg, arg + dimension);

    // Evaluate the Jacobian and copy components into the target array.
    for(int i = 0;
	i < dimension;
	++i)
      {
	// Evaluate the ith gradient map of the argument.
	std::vector<double> ithGradient
	  = rJacobian[i].evaluate(argVector);

	// Copy the result into the output array.
	std::copy(ithGradient.begin(),
		  ithGradient.end(),
		  pJac + (i * dimension));
      }

    // Zero the partial derivative wrt t.
    double* pComponent = dfdt + dimension;
    while(dfdt < pComponent--) *pComponent = 0.0;

    return GSL_SUCCESS;
  }
}
