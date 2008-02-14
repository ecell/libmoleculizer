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

#ifndef CPT_GROWEVENT_H
#define CPT_GROWEVENT_H

#include "cpt/cptEvent.hh"

namespace cpt
{
  class cptUnit;
  class cptApp;

  // I guess that, for now, this will grow (or shrink) all compartments
  // by the same factor.
  //
  // In effect, this issue has already arisen in the calculation of reaction
  // propensity.  There, I settled it by (for now?) confining each reaction to
  // a single compartment, which the reaction knows.  This defeats (for now)
  // the possibility of inter-compartment reactions.

  class growEvent :
    public cptEvent
  {
    double factor;
    double period;

  public:
    growEvent(double growthFactor,
	      double schedulingPeriod);

    fnd::eventResult
    happen(cptApp& rCptApp)
      throw(std::exception);
  };
}

#endif // CPT_GROWEVENT_H
