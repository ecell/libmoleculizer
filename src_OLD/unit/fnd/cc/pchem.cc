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
#include "fnd/pchem.hh"

namespace fnd
{
  // Reaction kinetics formulae that are independent of molecular
  // weight.  The formulae for bindingInvariant and bindingRate are
  // actually the basis for all the rest.

  // Rob says that this has something the effective mass of the center
  // of mass of a 2-body system if you figure the kinetic energy of
  // the system as a whole, or something like that.  Don't know.
  double
  reducedMass(double leftMass,
	      double rightMass)
  {
    return std::sqrt((leftMass * rightMass)
		     /(leftMass + rightMass));
  }

  double
  bindingInvariant(double bindingRate,
		   double leftMass,
		   double rightMass)
  {
    return bindingRate * reducedMass(leftMass,
				     rightMass);
  }

  double
  bindingRate(double bindingInvariant,
	      double leftMass,
	      double rightMass)
  {
    return bindingInvariant / reducedMass(leftMass,
					  rightMass);
  }

  // Implied formulae for dissociation constants.
  double
  dissocInvariant(double dissocConst,
		  double leftMass,
		  double rightMass)
  {
    return dissocConst / reducedMass(leftMass,
				     rightMass);
  }

  double
  dissocConst(double dissocInvariant,
	      double leftMass,
	      double rightMass)
  {
    return dissocInvariant * reducedMass(leftMass,
					 rightMass);
  }

  // Implied formulae for Michaelis constants.
  double
  michaelisInvariant(double michaelisConst,
		     double leftMass,
		     double rightMass)
  {
    return michaelisConst / reducedMass(leftMass,
					rightMass);
  }

  double
  michaelisConst(double michaelisInvariant,
		 double leftMass,
		 double rightMass)
  {
    return michaelisInvariant * reducedMass(leftMass,
					    rightMass);
  }
}
