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

#ifndef ODIEGLUE_H
#define ODIEGLUE_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Since I won't be using rk4 technology any more, perhaps I need
// to reorganize this?  Still want to leave tau leaping stuff in?
#include "rk4util/polymap.hh"

namespace odie
{
  class odieParams
  {
  public:
    const rk4util::polymap<double>& rDeriv;
    const std::vector<rk4util::polymap<double> >& rJac;

    odieParams(const rk4util::polymap<double>& rDerivative,
	       const std::vector<rk4util::polymap<double> >& rJacobian) :
      rDeriv(rDerivative),
      rJac(rJacobian)
    {}
  };
  
  int
  evalDeriv(double t,
	    const double arg[],
	    double returnValue[],
	    void* pParams);

  int
  evalJacobian(double t,
	       const double arg[],
	       double* pJac,
	       double dfdt[],
	       void* pParams);
}

#endif // ODIEGLUE_H
