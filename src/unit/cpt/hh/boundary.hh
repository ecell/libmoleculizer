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

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "cpt/compartment.hh"

namespace cpt
{
  class boundary
  {
    std::pair<compartment*, compartment*> neighbors;
    double area;
    double thickness;

  public:
    boundary(compartment* pFirst,
	     compartment* pSecond,
	     double boundaryArea,
	     double boundaryThickness) :
      neighbors(pFirst,
		pSecond),
      area(boundaryArea),
      thickness(boundaryThickness)
    {}

    // For now, I will assume that each edge in the graph is stored
    // just once as an ordered pair (boundary), and that searches
    // will look at both the first and second compartments as appropriate.
    compartment*
    getFirstCpt(void) const
    {
      return neighbors.first;
    }

    compartment*
    getSecondCpt(void) const
    {
      return neighbors.second;
    }

    double
    getArea(void) const
    {
      return area;
    }

    double
    getThickness(void) const
    {
      return thickness;
    }

    void
    insertState(xmlpp::Element* pBoundariesElt) const
      throw(std::exception);
  };
}

#endif // BOUNDARY_H
