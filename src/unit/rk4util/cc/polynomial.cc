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

#include <algorithm>
#include "rk4util/polynomial.hh"

// This is really just a test to force instantiation of the two classes
// of polynomials that are important to tau, and it does the extension
// of scalars that is used in tau.
void
convertIntPoly(const polynomial<int>& rIntPoly,
	       polynomial<double>& rTargetDoublePoly)
{
  rIntPoly.extendScalars(rTargetDoublePoly);
}

// Again, just testing.
void
extendIntPolyVector(const std::vector<polynomial<int> >& rSourcePolyVector,
		    std::vector<polynomial<double> >& rTargetPolyVector)
{
  extendVectorScalars(rSourcePolyVector,
		      rTargetPolyVector);
}



