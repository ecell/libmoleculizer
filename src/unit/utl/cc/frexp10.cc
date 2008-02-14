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
#include "utl/frexp10.hh"

namespace utl
{
  double
  frexp10(double num,
	  int& rExponent)
  {
    // Get the proper decomposition for base 2.
    int exponent2 = 0;
    double fraction2 = frexp(num, &exponent2);

    // Convert the exponent to an exponent base 10, leaving a remainder.
    double exponent10 = 0.0;
    double remainder10 = modf(exponent2 * M_LN2 / M_LN10,
			      &exponent10);

    // Floor shouldn't do anything here, I think.  If it DOES do anything,
    // then the correction below should take care of it.
    double exponent10Floor = floor(exponent10);
    rExponent = (int) exponent10Floor;

    return fraction2 * exp((remainder10 - exponent10Floor + exponent10)
			   * M_LN10);
  }
}
