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

#include "mol/molState.hh"
#include "plex/prm.hh"
#include "plex/plex.hh"
#include "mzr/moleculizer.hh"

namespace plx
{
  // This constructor doesn't construct the default parameter for a
  // given plex; that's done in the constructor for plexFamily.  This
  // just makes sure that the two vectors, molParams and bindingParams
  // have the right length, so that it's okay to index into them when
  // the plexFamily is constructed.
  plexParam::plexParam(const plex& rParadigm) :
    molParams(rParadigm.mols.size())
  {}

  class molWeightAccum : public std::unary_function<bnd::molParam, void>
  {
    double& rTotal;
  public:
    molWeightAccum(double& refTotal) :
      rTotal(refTotal)
    {}
    void operator()(const bnd::molParam prm) const
    {
      rTotal += prm->getMolWeight();
    }
  };
  double
  plexParam::getWeight(void) const
  {
    double mw = 0.0;
    for_each(molParams.begin(),
	     molParams.end(),
	     molWeightAccum(mw));
    return mw;
  }
}
